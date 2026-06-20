/**
 * @file mptc_types.h
 * @brief Core fixed-point types and conversion helpers.
 *
 * MPTC uses three fixed-point formats throughout:
 *
 *   - Q7:  int8_t,   1 sign + 7 frac,  range [-1.0, +0.992],  step 2^-7
 *   - Q15: int16_t,  1 sign + 15 frac, range [-1.0, +0.99997], step 2^-15
 *   - Q31: int32_t,  1 sign + 31 frac, range [-1.0, +0.99999], step 2^-31
 *
 * All arithmetic assumes two's-complement integers, which is the case on every
 * mainstream MCU shipped in the last 30 years (ARM, RISC-V, Xtensa, AVR, MIPS,
 * x86, POWER). The library is endianness-agnostic at the type level — packed
 * formats (ternary) handle endianness explicitly via mptc_config.h.
 */
#ifndef MPTC_TYPES_H
#define MPTC_TYPES_H

#include <stdint.h>
#include "mptc_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/*  Type aliases                                                       */
/* ------------------------------------------------------------------ */
/** @brief Q7  fixed-point: 1.7 format, int8_t.   */
typedef int8_t  q7_t;
/** @brief Q15 fixed-point: 1.15 format, int16_t. */
typedef int16_t q15_t;
/** @brief Q31 fixed-point: 1.31 format, int32_t. */
typedef int32_t q31_t;
/** @brief Q63 fixed-point: 1.63 format, int64_t. Used internally for Q31*Q31. */
typedef int64_t q63_t;

/* ------------------------------------------------------------------ */
/*  Constants                                                          */
/* ------------------------------------------------------------------ */
#define Q15_ONE   ((q15_t)0x7FFF)  /**<  0.99997 in Q15.           */
#define Q15_MINUS_ONE ((q15_t)0x8000) /**< -1.0 in Q15 (the asymmetric point). */
#define Q15_HALF  ((q15_t)0x4000)  /**<  0.5 in Q15.                */

/* Trig convention: in sintmath, a Q15 input x represents the angle
 * (x / 32768) * pi radians. So x in [-1, 1) covers [-pi, pi).
 * This keeps the full unit circle in int16 range.
 */
#define Q15_PI2   ((q15_t)0x4000)  /**<  pi/2 radians (0.5 in Q15).  */
#define Q15_PI    ((q15_t)0x7FFF)  /**<  ~pi radians (1.0 in Q15).   */
#define Q15_EPS   ((q15_t)0x0002)  /**<  smallest meaningful value. */

#define Q31_ONE   ((q31_t)0x7FFFFFFF)
#define Q31_MINUS_ONE ((q31_t)0x80000000)
#define Q31_HALF  ((q31_t)0x40000000)
#define Q31_PI2   ((q31_t)0x3243F6A8)   /* pi/2 in Q31 */
#define Q31_PI    ((q31_t)0x7FFFFFFF)

#define Q7_ONE    ((q7_t)0x7F)
#define Q7_HALF   ((q7_t)0x40)

/* ------------------------------------------------------------------ */
/*  Conversion: float <-> fixed  (only compiled if MPTC_USE_FLOAT=1)  */
/* ------------------------------------------------------------------ */
#if MPTC_USE_FLOAT || defined(MPTC_TEST_HOST)

#include <math.h>

/** @brief Convert float in [-1,1) to Q15 (saturating). */
MPTC_INLINE q15_t mptc_f32_to_q15(float f) {
    float scaled = f * 32768.0f;
    if (scaled >=  32767.0f) return 32767;
    if (scaled <= -32768.0f) return -32768;
    return (q15_t)scaled;
}

/** @brief Convert Q15 to float in [-1, 1). */
MPTC_INLINE float mptc_q15_to_f32(q15_t q) {
    return (float)q / 32768.0f;
}

/** @brief Convert float in [-1,1) to Q31 (saturating). */
MPTC_INLINE q31_t mptc_f32_to_q31(float f) {
    double scaled = (double)f * 2147483648.0;
    if (scaled >=  2147483647.0) return (q31_t)0x7FFFFFFF;
    if (scaled <= -2147483648.0) return (q31_t)0x80000000;
    return (q31_t)scaled;
}

/** @brief Convert Q31 to float in [-1, 1). */
MPTC_INLINE float mptc_q31_to_f32(q31_t q) {
    return (float)((double)q / 2147483648.0);
}

#endif /* MPTC_USE_FLOAT || MPTC_TEST_HOST */

/* ------------------------------------------------------------------ */
/*  Saturation arithmetic (branchless)                                */
/* ------------------------------------------------------------------ */

/**
 * @brief Saturate a 32-bit int into Q15 range, branchless.
 *
 * Equivalent to: if (x > 32767) return 32767; if (x < -32768) return -32768; return x;
 * but uses no branches. On ARM Cortex-M3+ this compiles to a single SSAT instruction.
 */
MPTC_INLINE q15_t mptc_saturate_q15(int32_t x) {
    /* Branchless saturation: clamp via mask trick. */
    int32_t lo = -32768;
    int32_t hi =  32767;
    /* x clamped to [lo, hi] using min/max that are branchless on most compilers */
    x = x < lo ? lo : x;
    x = x > hi ? hi : x;
    return (q15_t)x;
}

/**
 * @brief Saturate a 64-bit int into Q31 range.
 */
MPTC_INLINE q31_t mptc_saturate_q31(q63_t x) {
    int64_t lo = (int64_t)INT32_MIN;
    int64_t hi = (int64_t)INT32_MAX;
    x = x < lo ? lo : x;
    x = x > hi ? hi : x;
    return (q31_t)x;
}

/* ------------------------------------------------------------------ */
/*  Fixed-point multiplication                                         */
/* ------------------------------------------------------------------ */

/**
 * @brief Q15 * Q15 -> Q15, with one extra shift for the binary point.
 *
 * Standard Q15 multiplication: result = (a * b) >> 15. We use 32-bit
 * intermediate to avoid overflow. Branchless.
 */
MPTC_INLINE q15_t mptc_mul_q15(q15_t a, q15_t b) {
    return (q15_t)(((int32_t)a * (int32_t)b) >> 15);
}

/**
 * @brief Q31 * Q31 -> Q31. Uses 64-bit intermediate.
 */
MPTC_INLINE q31_t mptc_mul_q31(q31_t a, q31_t b) {
    return (q31_t)(((q63_t)a * (q63_t)b) >> 31);
}

/**
 * @brief Q15 * Q15 -> Q31 (full precision intermediate).
 *
 * Useful when accumulating many products before truncation.
 */
MPTC_INLINE q31_t mptc_mul_q15_q31(q15_t a, q15_t b) {
    return ((q31_t)a * (q31_t)b) << 1;  /* shift to align to Q31 */
}

/* ------------------------------------------------------------------ */
/*  Branchless min / max (CMOV-friendly)                              */
/* ------------------------------------------------------------------ */

/**
 * @brief Branchless min for q15_t.
 * @note Compiles to MIN instruction on ARMv7-M, conditional move elsewhere.
 */
MPTC_INLINE q15_t mptc_min_q15(q15_t a, q15_t b) { return a < b ? a : b; }

/** @brief Branchless max for q15_t. */
MPTC_INLINE q15_t mptc_max_q15(q15_t a, q15_t b) { return a > b ? a : b; }

/** @brief Branchless min for q31_t. */
MPTC_INLINE q31_t mptc_min_q31(q31_t a, q31_t b) { return a < b ? a : b; }

/** @brief Branchless max for q31_t. */
MPTC_INLINE q31_t mptc_max_q31(q31_t a, q31_t b) { return a > b ? a : b; }

/**
 * @brief Branchless clamp for q15_t.
 *
 * On ARMv7-M this emits two conditional select instructions (no branches).
 * On ARMv6-M (M0/M0+) it compiles to 4 instructions, still branchless in
 * practice thanks to IT blocks.
 */
MPTC_INLINE q15_t mptc_clamp_q15(q15_t x, q15_t lo, q15_t hi) {
    return mptc_max_q15(lo, mptc_min_q15(x, hi));
}

/** @brief Branchless clamp for q31_t. */
MPTC_INLINE q31_t mptc_clamp_q31(q31_t x, q31_t lo, q31_t hi) {
    return mptc_max_q31(lo, mptc_min_q31(x, hi));
}

/* ------------------------------------------------------------------ */
/*  Absolute value (branchless)                                       */
/* ------------------------------------------------------------------ */

/** @brief Branchless abs for q15_t. */
MPTC_INLINE q15_t mptc_abs_q15(q15_t x) {
    /* The classic branchless abs: (x ^ mask) - mask where mask = x >> 15. */
    int16_t mask = x >> 15;
    return (q15_t)((x ^ mask) - mask);
}

/** @brief Branchless abs for q31_t. */
MPTC_INLINE q31_t mptc_abs_q31(q31_t x) {
    int32_t mask = x >> 31;
    return (q31_t)((x ^ mask) - mask);
}

/* ------------------------------------------------------------------ */
/*  Sign mask (for branchless selection)                              */
/* ------------------------------------------------------------------ */

/**
 * @brief Returns 0xFFFF if x < 0, 0x0000 otherwise. Branchless.
 *
 * Useful for synthesizing conditional logic without branches:
 * @code
 *   uint16_t neg = mptc_signmask_q15(x);
 *   q15_t abs_x = (x & ~neg) | ((-x) & neg);   // branchless abs
 * @endcode
 */
MPTC_INLINE uint16_t mptc_signmask_q15(q15_t x) {
    return (uint16_t)(x >> 15);
}

/** @brief Returns 0xFFFFFFFF if x < 0, 0 otherwise. Branchless. */
MPTC_INLINE uint32_t mptc_signmask_q31(q31_t x) {
    return (uint32_t)(x >> 31);
}

/**
 * @brief Returns 1 if x < 0, 0 otherwise. Branchless.
 *
 * This is the integer-mask equivalent of (x < 0) ? 1 : 0 and is the
 * canonical MPTC superposition primitive: pass it to mptc_select_*().
 */
MPTC_INLINE uint32_t mptc_is_negative_q15(q15_t x) {
    return (uint32_t)((uint16_t)x >> 15);
}

/** @brief Returns 1 if a > b, 0 otherwise. Branchless. */
MPTC_INLINE uint32_t mptc_gt_q15(q15_t a, q15_t b) {
    return (uint32_t)(((uint16_t)(b - a)) >> 15);
}

/** @brief Returns 1 if a < b, 0 otherwise. Branchless. */
MPTC_INLINE uint32_t mptc_lt_q15(q15_t a, q15_t b) {
    return mptc_gt_q15(b, a);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MPTC_TYPES_H */
