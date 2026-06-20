/**
 * @file sdft.h
 * @brief Sliding Discrete Fourier Transform (SDFT) in fixed-point.
 *
 * The classical block FFT computes the spectrum of N samples in O(N log N)
 * time but introduces latency spikes: you must wait for a full block before
 * producing output, and the worst-case compute time is non-deterministic
 * relative to the sample rate. This breaks real-time control loops.
 *
 * The Sliding DFT updates the spectrum incrementally on every new sample in
 * O(N) time per tick. After a one-block warm-up, every sample produces a
 * fresh spectrum with flat, deterministic latency. The MPTC paper reports
 * 2.6 ms/tick vs 5.1 ms worst-case for block FFT.
 *
 * Recurrence (paper eq. 5):
 *
 *   X_k[new] = (X_k[old] - x_old + x_new) * exp(j * 2*pi*k/N)
 *
 * In fixed-point, we store Re and Im separately as Q31, with the twiddle
 * factor stored as a Q15 cos/sin pair. The complex multiply uses Q31
 * intermediate to preserve precision.
 *
 * @par Numerical stability
 * The naive SDFT recurrence is marginally stable — round-off errors can
 * accumulate. We use the "modulated SDFT" variant (Jacobsen & Lyons 2003)
 * which guarantees bounded error growth over arbitrary run lengths.
 */
#ifndef MPTC_SDFT_H
#define MPTC_SDFT_H

#include <stdint.h>
#include <stdlib.h>
#include "mptc_config.h"
#include "mptc_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief SDFT state, Q31 fixed-point.
 *
 * Use one instance per signal stream. Not thread-safe; guard with mutex if
 * shared.
 */
typedef struct {
    uint16_t n_bins;         /**< Number of frequency bins (FFT size N).       */
    uint16_t idx;            /**< Circular index into the delay line.          */
    q31_t   *re;             /**< Real part of each bin, size n_bins.          */
    q31_t   *im;             /**< Imaginary part of each bin, size n_bins.      */
    q15_t   *twiddle_cos;    /**< cos(2*pi*k/N) for each k, Q15.                */
    q15_t   *twiddle_sin;    /**< sin(2*pi*k/N) for each k, Q15.                */
    q15_t   *delay;          /**< Circular delay line, size n_bins.             */
    q15_t   *cos_n;          /**< Per-step cos(2*pi/N), Q15 (modulated SDFT).   */
    q15_t   *sin_n;          /**< Per-step sin(2*pi/N), Q15 (modulated SDFT).   */
    q31_t    cf_re;          /**< Modulation carrier, real.                     */
    q31_t    cf_im;          /**< Modulation carrier, imag.                     */
    uint32_t sample_rate;    /**< For metadata only; does not affect compute.   */
} sdft_q31_t;

/**
 * @brief Initialize an SDFT instance.
 *
 * Allocates internal buffers via a single malloc(). For bare-metal systems
 * without malloc, you may pre-allocate the buffers and assign the pointers
 * manually (see ::sdft_q31_attach).
 *
 * @param s           State to initialize.
 * @param n_bins      Number of frequency bins (FFT size N). Must be >= 4.
 * @param sample_rate Sample rate in Hz. Used only for metadata.
 * @return 0 on success, -1 on allocation failure or invalid argument.
 */
int  sdft_q31_init(sdft_q31_t *s, uint16_t n_bins, uint32_t sample_rate);

/**
 * @brief Attach pre-allocated buffers (no malloc).
 *
 * For bare-metal systems. The caller must provide buffers of the right size:
 *   - re, im:        n_bins * sizeof(q31_t) each
 *   - twiddle_cos/sin: n_bins * sizeof(q15_t) each
 *   - delay:         n_bins * sizeof(q15_t)
 *   - cos_n, sin_n:  n_bins * sizeof(q15_t) each
 *
 * After attaching, call ::sdft_q31_reset to populate the twiddle tables.
 */
void sdft_q31_attach(sdft_q31_t *s,
                     q31_t *re, q31_t *im,
                     q15_t *twiddle_cos, q15_t *twiddle_sin,
                     q15_t *delay, q15_t *cos_n, q15_t *sin_n,
                     uint16_t n_bins, uint32_t sample_rate);

/**
 * @brief Reset state (zero all bins and delay line, rebuild twiddles).
 */
void sdft_q31_reset(sdft_q31_t *s);

/**
 * @brief Free buffers allocated by ::sdft_q31_init.
 *
 * Safe to call on an attached (non-owning) instance — it's a no-op.
 */
void sdft_q31_free(sdft_q31_t *s);

/**
 * @brief Process one new sample, updating the entire spectrum.
 *
 * O(N) per call. Deterministic latency independent of input data.
 *
 * @param s     SDFT state.
 * @param sample New sample in Q15.
 */
void sdft_q31_update(sdft_q31_t *s, q15_t sample);

/**
 * @brief Read the real part of bin k.
 * @note Value is in Q31. Magnitude is roughly N * (input amplitude)/2.
 */
MPTC_INLINE q31_t sdft_q31_real(const sdft_q31_t *s, uint16_t k) {
    MPTC_ASSERT(k < s->n_bins);
    return s->re[k];
}

/** @brief Read the imaginary part of bin k. */
MPTC_INLINE q31_t sdft_q31_imag(const sdft_q31_t *s, uint16_t k) {
    MPTC_ASSERT(k < s->n_bins);
    return s->im[k];
}

/**
 * @brief Read |X_k| (magnitude) using branchless integer sqrt.
 *
 * We use a 4-iteration Newton sqrt for Q31, which gives < 1 % error.
 * The result is in Q31; convert to Q15 by >> 16 if needed.
 */
q31_t sdft_q31_mag(const sdft_q31_t *s, uint16_t k);

/**
 * @brief Read |X_k|^2 (power) — cheaper than ::sdft_q31_mag, no sqrt.
 *
 * Result is in Q31 (squared magnitude, scaled). Useful for threshold-based
 * detectors and for argmax across bins.
 */
MPTC_INLINE q31_t sdft_q31_power(const sdft_q31_t *s, uint16_t k) {
    MPTC_ASSERT(k < s->n_bins);
    int64_t r = s->re[k];
    int64_t i = s->im[k];
    int64_t p = r*r + i*i;
    /* p is in Q62. Shift to Q31. */
    return (q31_t)(p >> 31);
}

/**
 * @brief Find the bin with the largest magnitude (argmax).
 *
 * O(N) branchless scan. Useful for pitch detection / dominant frequency.
 */
uint16_t sdft_q31_argmax(const sdft_q31_t *s);

/**
 * @brief Estimate the dominant frequency in Hz.
 *
 * @return sample_rate * argmax / n_bins.
 */
MPTC_INLINE float sdft_q31_dominant_hz(const sdft_q31_t *s) {
    uint16_t k = sdft_q31_argmax(s);
    return (float)s->sample_rate * (float)k / (float)s->n_bins;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MPTC_SDFT_H */
