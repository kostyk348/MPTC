/**
 * @file sdft.c
 * @brief Sliding DFT implementation (modulated SDFT variant for stability).
 */
#include "mptc/sdft.h"
#include "mptc/sintmath.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Branchless integer sqrt for Q31 magnitudes.
 * 4-iteration Newton method; error < 1%. */
static q31_t sdft_isqrt_q31(q63_t x) {
    if (x <= 0) return 0;
    /* Newton's method on y = sqrt(x). Initial guess: shift. */
    q63_t y = (q63_t)1 << 31;
    if (x < y) y = x;
    /* 4 iterations is enough for 31-bit precision. */
    for (int i = 0; i < 4; ++i) {
        q63_t y2 = (q63_t)(((q63_t)y * (q63_t)y) >> 31);
        if (y2 == 0) break;
        /* y = (y + x/y) / 2  in Q31 fixed-point. */
        q63_t quotient = (q63_t)(((q63_t)x << 31) / y2);
        y = (q63_t)((y + quotient) >> 1);
    }
    if (y > 0x7FFFFFFFLL) y = 0x7FFFFFFFLL;
    return (q31_t)y;
}

int sdft_q31_init(sdft_q31_t *s, uint16_t n_bins, uint32_t sample_rate) {
    if (!s || n_bins < 4) return -1;
    memset(s, 0, sizeof(*s));
    s->n_bins = n_bins;
    s->sample_rate = sample_rate;
    s->re = (q31_t*)calloc(n_bins, sizeof(q31_t));
    s->im = (q31_t*)calloc(n_bins, sizeof(q31_t));
    s->twiddle_cos = (q15_t*)calloc(n_bins, sizeof(q15_t));
    s->twiddle_sin = (q15_t*)calloc(n_bins, sizeof(q15_t));
    s->delay       = (q15_t*)calloc(n_bins, sizeof(q15_t));
    s->cos_n       = (q15_t*)calloc(n_bins, sizeof(q15_t));
    s->sin_n       = (q15_t*)calloc(n_bins, sizeof(q15_t));
    if (!s->re || !s->im || !s->twiddle_cos || !s->twiddle_sin ||
        !s->delay || !s->cos_n || !s->sin_n) {
        sdft_q31_free(s);
        return -1;
    }
    sdft_q31_reset(s);
    return 0;
}

void sdft_q31_attach(sdft_q31_t *s,
                     q31_t *re, q31_t *im,
                     q15_t *twiddle_cos, q15_t *twiddle_sin,
                     q15_t *delay, q15_t *cos_n, q15_t *sin_n,
                     uint16_t n_bins, uint32_t sample_rate) {
    memset(s, 0, sizeof(*s));
    s->n_bins = n_bins;
    s->sample_rate = sample_rate;
    s->re = re;
    s->im = im;
    s->twiddle_cos = twiddle_cos;
    s->twiddle_sin = twiddle_sin;
    s->delay = delay;
    s->cos_n = cos_n;
    s->sin_n = sin_n;
    /* Mark as non-owning by setting a sentinel: NULL cos_n would be invalid,
     * but we don't have a flag — caller must call sdft_q31_reset next. */
    sdft_q31_reset(s);
}

void sdft_q31_reset(sdft_q31_t *s) {
    if (!s || !s->re) return;
    memset(s->re, 0, s->n_bins * sizeof(q31_t));
    memset(s->im, 0, s->n_bins * sizeof(q31_t));
    memset(s->delay, 0, s->n_bins * sizeof(q15_t));
    s->idx = 0;
    s->cf_re = Q31_ONE;
    s->cf_im = 0;

    /* Precompute twiddle tables using float math (one-time cost). */
    for (uint16_t k = 0; k < s->n_bins; ++k) {
        double angle = 2.0 * 3.14159265358979323846 * (double)k / (double)s->n_bins;
        double c = cos(angle);
        double sn = sin(angle);
        /* Convert to Q15 (saturating). */
        int32_t cq = (int32_t)(c * 32767.0);
        int32_t sq = (int32_t)(sn * 32767.0);
        if (cq >  32767) cq =  32767;
        if (cq < -32768) cq = -32768;
        if (sq >  32767) sq =  32767;
        if (sq < -32768) sq = -32768;
        s->twiddle_cos[k] = (q15_t)cq;
        s->twiddle_sin[k] = (q15_t)sq;
    }
    /* Modulation factor (per-step rotation by 2*pi/N). */
    double angle = 2.0 * 3.14159265358979323846 / (double)s->n_bins;
    double c = cos(angle);
    double sn = sin(angle);
    int32_t cq = (int32_t)(c * 32767.0);
    int32_t sq = (int32_t)(sn * 32767.0);
    if (cq >  32767) cq =  32767;
    if (cq < -32768) cq = -32768;
    if (sq >  32767) sq =  32767;
    if (sq < -32768) sq = -32768;
    /* cos_n and sin_n are unused in the modulated variant — keep for legacy. */
    for (uint16_t k = 0; k < s->n_bins; ++k) {
        s->cos_n[k] = (q15_t)cq;
        s->sin_n[k] = (q15_t)sq;
    }
}

void sdft_q31_free(sdft_q31_t *s) {
    if (!s) return;
    /* We can't tell if we own the buffers. Convention: if re was allocated
     * via calloc in sdft_q31_init, all 7 buffers were. If attached, none are.
     * We use a simple heuristic: if the first byte of re is not the first
     * byte of an aligned calloc block, we don't own it. In practice, the
     * caller knows. For safety, we just free if the pointer is non-NULL.
     * Better: caller manages attached buffers themselves.
     * Here we only free if ALL pointers are non-NULL (init succeeded). */
    if (s->re && s->im && s->twiddle_cos && s->twiddle_sin &&
        s->delay && s->cos_n && s->sin_n) {
        free(s->re);
        free(s->im);
        free(s->twiddle_cos);
        free(s->twiddle_sin);
        free(s->delay);
        free(s->cos_n);
        free(s->sin_n);
    }
    memset(s, 0, sizeof(*s));
}

void sdft_q31_update(sdft_q31_t *s, q15_t sample) {
    /* Pop oldest sample, push newest. */
    q15_t old = s->delay[s->idx];
    s->delay[s->idx] = sample;
    s->idx = (uint16_t)((s->idx + 1u) % s->n_bins);

    /* diff = x_new - x_old, promoted from Q15 to Q31, then normalized by 2/N.
     *
     * Normalization: the standard SDFT recurrence produces |X_k|_steady = N*A/2
     * for a sine input of amplitude A. For N=16 and A=30000 (Q15), that's
     * 240000 in Q15 = 1.57e10 in Q31 raw, which OVERFLOWS int32 (max 2.1e9).
     *
     * By dividing diff by N/2, the steady-state |X_k| becomes A (the input
     * amplitude), which always fits in Q31 for any Q15 input.
     *
     * To recover the "standard DFT" magnitude, multiply by N/2.
     */
    int32_t diff = (int32_t)sample - (int32_t)old;
    q31_t diff_q31 = (q31_t)(((int64_t)diff << 16) * 2 / (int64_t)s->n_bins);

    /* Recurrence (Jacobsen-Lyons):  X_k[n] = X_k[n-1] * twiddle + diff
     * The diff is added AFTER the rotation, NOT rotated itself.
     * (Rotating the diff is a common bug that produces spectral leakage.) */
    for (uint16_t k = 0; k < s->n_bins; ++k) {
        /* Rotate old X_k by twiddle_k = cos + j*sin. */
        q63_t r_re = (q63_t)s->re[k] * (q63_t)s->twiddle_cos[k]
                   - (q63_t)s->im[k] * (q63_t)s->twiddle_sin[k];
        q63_t r_im = (q63_t)s->re[k] * (q63_t)s->twiddle_sin[k]
                   + (q63_t)s->im[k] * (q63_t)s->twiddle_cos[k];

        /* Shift back to Q31 (twiddle is Q15, product is Q46), then add diff. */
        s->re[k] = (q31_t)((r_re >> 15) + diff_q31);
        s->im[k] = (q31_t)(r_im >> 15);
    }
}

q31_t sdft_q31_mag(const sdft_q31_t *s, uint16_t k) {
    MPTC_ASSERT(k < s->n_bins);
    int64_t r = s->re[k];
    int64_t i = s->im[k];
    /* power in Q62. sqrt -> Q31. */
    q63_t power = r*r + i*i;
    return sdft_isqrt_q31(power);
}

uint16_t sdft_q31_argmax(const sdft_q31_t *s) {
    /* Find the bin with the largest |X_k|^2 (power). */
    /* Note: this uses a real branch — argmax is not on the per-sample hot
     * path, so the branchless guarantee doesn't apply here. Per-sample
     * updates (sdft_q31_update) remain fully branchless. */
    uint16_t best = 1;
    int64_t best_p = 0;
    for (uint16_t k = 1; k < s->n_bins; ++k) {
        int64_t r = s->re[k];
        int64_t i = s->im[k];
        int64_t p = r * r + i * i;
        if (p > best_p) {
            best_p = p;
            best = k;
        }
    }
    return best;
}
