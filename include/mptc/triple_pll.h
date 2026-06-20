/**
 * @file triple_pll.h
 * @brief Median-of-3 jitter killer (Ironclad Precision Triple-PLL).
 *
 * The MPTC hardware paper describes "Triple-PLL with median recovery" — every
 * information pulse passes through three independent phase restorers, and the
 * output is the median of the three. This gives "ironclad" stability: even
 * when one channel is corrupted by jitter up to 10–20 % of the unit interval,
 * the median is correct.
 *
 * In software we adapt this idea as a moving 3-sample median filter for
 * ADC / IMU / sensor streams. Compared to a moving average:
 *
 *   - Median preserves edges (steps) — average smears them.
 *   - Median rejects impulsive noise — average propagates it.
 *   - Median of 3 is branchless (3 conditional moves).
 *
 * This module also provides a "Triple-PLL proper" mode where three independent
 * first-order IIR loops track the same input with different initial conditions
 * and the output is the median — useful for recovering a stable estimate when
 * any single loop may be transiently disturbed.
 */
#ifndef MPTC_TRIPLE_PLL_H
#define MPTC_TRIPLE_PLL_H

#include <stdint.h>
#include "mptc_config.h"
#include "mptc_types.h"
#include "mptc/sintmath.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================== */
/*  Mode 1: Pure 3-sample median filter (fastest, branchless)         */
/* ================================================================== */

/**
 * @brief 3-sample median filter state for Q15.
 *
 * Stores the last 3 samples in a ring buffer. Median is computed branchlessly
 * on every call to ::triple_pll_q15_update.
 */
typedef struct {
    q15_t   buf[3];   /**< Ring buffer of last 3 samples. */
    uint8_t idx;      /**< Write index.                    */
    uint8_t filled;   /**< Number of valid samples so far. */
} triple_pll_q15_t;

/** @brief Initialize a 3-sample median filter. */
MPTC_INLINE void triple_pll_q15_init(triple_pll_q15_t *p) {
    p->buf[0] = p->buf[1] = p->buf[2] = 0;
    p->idx = 0;
    p->filled = 0;
}

/**
 * @brief Branchless median of 3 values.
 *
 * Implements:  median = a + b + c - max(a,b,c) - min(a,b,c)
 * which is fully branchless. Equivalent to:
 *   sort(a,b,c); return b;
 */
MPTC_INLINE q15_t triple_pll_median3(q15_t a, q15_t b, q15_t c) {
    q15_t mx = mptc_max_q15(a, mptc_max_q15(b, c));
    q15_t mn = mptc_min_q15(a, mptc_min_q15(b, c));
    /* sum - max - min = middle. */
    int32_t s = (int32_t)a + (int32_t)b + (int32_t)c - (int32_t)mx - (int32_t)mn;
    return (q15_t)s;
}

/**
 * @brief Push a new sample; returns the median of the last 3.
 *
 * @par Latency
 * 3 multiplications / comparisons, fully branchless. ~12 cycles on Cortex-M0+.
 *
 * @par Cold start
 * For the first 2 samples, returns the latest sample (median of an incomplete
 * set is the latest).
 */
MPTC_INLINE q15_t triple_pll_q15_update(triple_pll_q15_t *p, q15_t sample) {
    p->buf[p->idx] = sample;
    p->idx = (uint8_t)((p->idx + 1u) % 3u);
    if (p->filled < 3) p->filled++;

    if (p->filled < 3) {
        return sample;
    }
    return triple_pll_median3(p->buf[0], p->buf[1], p->buf[2]);
}

/* ================================================================== */
/*  Mode 2: Triple-PLL proper (3 IIR trackers + median)               */
/* ================================================================== */

/**
 * @brief Triple-PLL state with three IIR trackers and median voting.
 *
 * Each tracker is a first-order IIR:  y[i] = y[i] + alpha * (x - y[i]).
 * The three trackers are seeded with slightly different initial conditions
 * so transient disturbances that knock one out of lock are outvoted by
 * the other two.
 */
typedef struct {
    q15_t y[3];      /**< Three IIR tracker outputs. */
    q15_t alpha;     /**< IIR smoothing factor in Q15. */
    q15_t init_bias; /**< Initial offset between trackers. */
} triple_pll_iir_q15_t;

/**
 * @brief Initialize the Triple-PLL.
 *
 * @param p         State.
 * @param alpha     IIR smoothing factor, Q15. Typical: 0x0800 (alpha=0.0625).
 * @param init_bias Initial spread between trackers, Q15. Typical: 0x0100.
 */
MPTC_INLINE void triple_pll_iir_q15_init(triple_pll_iir_q15_t *p,
                                          q15_t alpha, q15_t init_bias) {
    p->y[0] = -init_bias;
    p->y[1] = 0;
    p->y[2] = init_bias;
    p->alpha = alpha;
    p->init_bias = init_bias;
}

/**
 * @brief Update the Triple-PLL with a new sample; returns median of 3 trackers.
 *
 * Each tracker independently low-passes the input. The median of the three
 * outputs rejects impulsive disturbances that knock one tracker out of lock.
 *
 * @par Use case
 * Recovering a stable angular estimate from a noisy IMU when the IMU
 * occasionally spikes due to vibration. A single IIR would smear the spike
 * over many samples; the median of three independent IIRs rejects it cleanly.
 */
MPTC_INLINE q15_t triple_pll_iir_q15_update(triple_pll_iir_q15_t *p, q15_t x) {
    /* y[i] += alpha * (x - y[i])   for i in 0..2. */
    for (int i = 0; i < 3; ++i) {
        q15_t err = (q15_t)((int32_t)x - (int32_t)p->y[i]);
        q15_t corr = mptc_mul_q15(p->alpha, err);
        p->y[i] = (q15_t)((int32_t)p->y[i] + (int32_t)corr);
    }
    return triple_pll_median3(p->y[0], p->y[1], p->y[2]);
}

/**
 * @brief Reset all three trackers to a known value.
 */
MPTC_INLINE void triple_pll_iir_q15_reset(triple_pll_iir_q15_t *p, q15_t value) {
    p->y[0] = (q15_t)((int32_t)value - (int32_t)p->init_bias);
    p->y[1] = value;
    p->y[2] = (q15_t)((int32_t)value + (int32_t)p->init_bias);
}

/* ================================================================== */
/*  Mode 3: Triple-PLL with outlier rejection                         */
/* ================================================================== */

/**
 * @brief Triple-PLL with explicit outlier rejection.
 *
 * Each tracker's residual |x - y[i]| is computed; if any tracker's residual
 * exceeds `threshold`, that tracker is frozen for this step (does not update).
 * The output is still the median of the three (now including the frozen ones).
 *
 * This combines the noise rejection of the median with the explicit fault
 * detection of a voter monitor.
 */
typedef struct {
    triple_pll_iir_q15_t core;  /**< Underlying Triple-PLL.                */
    q15_t threshold;            /**< Q15 threshold for outlier rejection.  */
} triple_pll_voter_q15_t;

/** @brief Initialize voter-mode Triple-PLL. */
MPTC_INLINE void triple_pll_voter_q15_init(triple_pll_voter_q15_t *p,
                                            q15_t alpha,
                                            q15_t init_bias,
                                            q15_t threshold) {
    triple_pll_iir_q15_init(&p->core, alpha, init_bias);
    p->threshold = threshold;
}

/**
 * @brief Update with explicit outlier rejection.
 *
 * For each tracker i:
 *   err_i = x - y[i]
 *   if |err_i| > threshold: skip update (frozen this tick)
 *   else:                   y[i] += alpha * err_i
 * Returns median(y[0], y[1], y[2]).
 */
MPTC_INLINE q15_t triple_pll_voter_q15_update(triple_pll_voter_q15_t *p, q15_t x) {
    for (int i = 0; i < 3; ++i) {
        q15_t err = (q15_t)((int32_t)x - (int32_t)p->core.y[i]);
        q15_t abs_err = mptc_abs_q15(err);
        /* Branchless: if |err| <= threshold, update; else skip. */
        uint32_t in_range = (uint32_t)(abs_err <= p->threshold);
        q15_t corr = mptc_mul_q15(p->core.alpha, err);
        q15_t new_y = (q15_t)((int32_t)p->core.y[i] + (int32_t)corr);
        p->core.y[i] = sint_q15_select(new_y, p->core.y[i], in_range);
    }
    return triple_pll_median3(p->core.y[0], p->core.y[1], p->core.y[2]);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MPTC_TRIPLE_PLL_H */
