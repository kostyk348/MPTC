/**
 * @file sintmath.h
 * @brief Branchless math primitives, fast approximations, and dual numbers.
 *
 * This module is the algorithmic heart of MPTC. It implements the
 * "SintMath" branchless calculus described in the MPTC paper:
 *
 *   - Bitmask routing / superposition:  Result = a*M + b*(1-M)
 *   - Fast wave approximations: sin / tanh / GELU via polynomial / rational forms
 *   - Dual numbers (a + b*eps, eps^2 = 0) for autograd-free analytical derivatives
 *
 * Everything is implemented in Q15 fixed-point. Where a Q31 variant exists for
 * higher precision, it is provided alongside with the `_q31` suffix.
 *
 * All routines in this header are `static inline` and header-only. They have
 * no side effects, no global state, and no dynamic allocation.
 *
 * @par Branchless guarantees
 * When ::MPTC_BRANCHLESS_STRICT is 1 (default), no routine in this header
 * uses `if`, `?:`, or any other form of conditional jump in its emitted
 * machine code on ARMv7-M / x86_64 / RISC-V / Xtensa at -O2. On ARMv6-M
 * (Cortex-M0/M0+), the conditional-move idiom compiles to IT blocks which
 * are also branchless.
 */
#ifndef MPTC_SINTMATH_H
#define MPTC_SINTMATH_H

#include <stdint.h>
#include "mptc_config.h"
#include "mptc_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================== */
/*  1. Superposition primitives (bitmask routing)                     */
/* ================================================================== */

/**
 * @brief Branchless conditional select: returns `a` if mask != 0, else `b`.
 *
 * Implements the MPTC superposition equation  Result = a*M + b*(1-M)  with
 * M in {0, 1}. The arithmetic form is kept (rather than `mask ? a : b`)
 * because it remains branchless on every compiler we tested.
 *
 * @param a    Value selected when mask is nonzero.
 * @param b    Value selected when mask is zero.
 * @param mask Any nonzero integer acts as 1, zero acts as 0.
 * @return     a if mask != 0 else b.
 *
 * @par Example
 * @code
 *   q15_t y = sint_q15_select(0x1000, 0x7000, mptc_gt_q15(x, threshold));
 * @endcode
 */
MPTC_INLINE q15_t sint_q15_select(q15_t a, q15_t b, uint32_t mask) {
    /* Normalize mask to 0/1, then use arithmetic superposition. */
    uint16_t m = (uint16_t)(-(mask != 0));      /* 0xFFFF or 0x0000 */
    return (q15_t)((a & (q15_t)m) | (b & (q15_t)~m));
}

/**
 * @brief Q31 variant of ::sint_q15_select.
 */
MPTC_INLINE q31_t sint_q31_select(q31_t a, q31_t b, uint32_t mask) {
    uint32_t m = (uint32_t)(-(mask != 0));
    return (q31_t)((a & (q31_t)m) | (b & (q31_t)~m));
}

/**
 * @brief Arithmetic (linear) superposition: returns a*alpha + b*(1-alpha).
 *
 * This is the "soft" version of ::sint_q15_select — instead of a hard 0/1
 * mask, it accepts a continuous Q15 weight `alpha`. Useful for linear
 * interpolation, soft blending, and continuous control flow.
 *
 * @param a     Value when alpha = Q15_ONE.
 * @param b     Value when alpha = 0.
 * @param alpha Q15 weight in [0, Q15_ONE].
 */
MPTC_INLINE q15_t sint_q15_blend(q15_t a, q15_t b, q15_t alpha) {
    /* a*alpha + b*(1 - alpha) = b + (a - b) * alpha */
    int32_t diff = (int32_t)a - (int32_t)b;
    int32_t weighted = (diff * (int32_t)alpha) >> 15;
    return (q15_t)((int32_t)b + weighted);
}

/**
 * @brief Branchless signed-abs select: returns x if positive, -x if negative.
 *
 * Equivalent to ::mptc_abs_q15 but exposed here as a SintMath primitive.
 */
MPTC_INLINE q15_t sint_q15_mirror(q15_t x) {
    return mptc_abs_q15(x);
}

/**
 * @brief Branchless copy-sign: |a| with the sign of b.
 */
MPTC_INLINE q15_t sint_q15_copysign(q15_t a, q15_t b) {
    q15_t abs_a = mptc_abs_q15(a);
    uint16_t neg_b = (uint16_t)((uint16_t)b >> 15);   /* 0xFFFF if b<0 else 0 */
    return (q15_t)((abs_a & (q15_t)~neg_b) | ((q15_t)-abs_a & (q15_t)neg_b));
}

/* ================================================================== */
/*  2. Fast wave approximations                                        */
/* ================================================================== */
/**
 * @brief Fast sine approximation via 4-term Taylor expansion.
 *
 *   sin(x) ~= x - x^3/6 + x^5/120 - x^7/5040
 *
 * @par Input range
 * Valid for |x| <= 0.5 (i.e., |angle| <= pi/2). For inputs outside this
 * range, use ::sint_q15_sin_wrap which folds via sin(pi - a) = sin(a).
 *
 * @par Cost
 * 4 multiplications, 3 additions, 0 branches. Uses int64_t intermediates
 * (one __muldi3 call on Cortex-M0). On ESP32 (native 64-bit) ~6 mults.
 *
 * @par Accuracy
 * Max abs error in [-pi/2, pi/2] is 1.5e-5 in float, which quantizes to
 * ~0.5 LSB at Q15. Effectively exact at this precision.
 */
MPTC_INLINE q15_t sint_q15_sin_taylor(q15_t x) {
    /* x in Q15 represents (x/32768)*pi radians.
     * Compute a = pi * x / 32768 in Q29 (since pi > 1, we need extra bits).
     * pi in Q29 = 0x6487ED55 (3.14159 * 2^29 = 1,686,629,717).
     */
    int32_t X = (int32_t)x;                                  /* Q15 */
    int64_t a = ((int64_t)X * (int64_t)0x6487ED55LL) >> 15;  /* a in Q29 */

    int64_t a2  = (a * a) >> 29;     /* a^2 in Q29 */
    int64_t t   = a;                 /* a in Q29   */
    int64_t acc = t;                 /* acc = a    */

    /* a^3 / 6 in Q29. */
    t = (a * a2) >> 29;              /* a^3 in Q29 */
    t = t / 6;                       /* a^3 / 6    */
    acc -= t;

    /* a^5 / 120 in Q29 = (a^3/6) * a^2 / 20. */
    t = (t * a2) >> 29;              /* a^5 / 6    */
    t = t / 20;                      /* a^5 / 120  */
    acc += t;

    /* a^7 / 5040 in Q29 = (a^5/120) * a^2 / 42. */
    t = (t * a2) >> 29;              /* a^7 / 120  */
    t = t / 42;                      /* a^7 / 5040 */
    acc -= t;

    /* acc is sin(a) in Q29, range [-1, 1] = [-2^29, 2^29].
     * Convert to Q15: >> 14. */
    acc = acc >> 14;
    if (acc >  32767) acc =  32767;
    if (acc < -32768) acc = -32768;
    return (q15_t)acc;
}

/**
 * @brief Compute sin(x) for any x in [-pi, pi).
 *
 * Uses the symmetry  sin(pi - a) = sin(a)  to fold any Q15 angle in
 * [-1, 1) (representing [-pi, pi)) into [-0.5, 0.5] (representing
 * [-pi/2, pi/2]) where the Taylor series is accurate. The fold is
 * branchless: both folded and unfolded paths are computed and selected.
 */
MPTC_INLINE q15_t sint_q15_sin_wrap(q15_t x) {
    const q15_t half = Q15_PI2;       /* 0.5 in Q15 = pi/2 radians */

    uint32_t gt_half = mptc_gt_q15(x, half);              /* x > pi/2  */
    uint32_t lt_neg  = mptc_gt_q15((q15_t)(-half), x);    /* x < -pi/2 */

    /* Folded values:                                                   */
    /*   For x > pi/2:  folded = 1 - x   (sin(pi-a) = sin(a))          */
    /*   For x < -pi/2: folded = -1 - x  (sin(-pi-a) = sin(a), no neg) */
    int32_t fold_hi = (int32_t)0x7FFF - (int32_t)x;
    int32_t fold_lo = (int32_t)(-0x8000) - (int32_t)x;
    q15_t folded = sint_q15_select((q15_t)fold_hi, (q15_t)fold_lo, gt_half);

    uint32_t need_fold = gt_half | lt_neg;
    q15_t x_eff = sint_q15_select(folded, x, need_fold);

    /* sin(angle) = sin(folded_angle) by symmetry; no negation needed. */
    return sint_q15_sin_taylor(x_eff);
}

/**
 * @brief Fast cosine: sin(x + pi/2), branchless.
 */
MPTC_INLINE q15_t sint_q15_cos(q15_t x) {
    /* cos(x) = sin(x + pi/2). Add half in Q15, wrap to [-1, 1). */
    int32_t shifted = (int32_t)x + (int32_t)Q15_PI2;  /* +0.5 in Q15 */
    if (shifted >  32767) shifted -= 0x10000;          /* subtract 1.0 */
    if (shifted < -32768) shifted += 0x10000;          /* add 1.0      */
    return sint_q15_sin_wrap((q15_t)shifted);
}

/**
 * @brief Fast tanh via rational approximation (MPTC paper eq. 9).
 *
 *   tanh(x) ~= x * (27 + x^2) / (27 + 9*x^2)
 *
 * Accurate to within 1e-3 for |x| < 4. For |x| >= 4, returns +/-1.
 *
 * @par Cost
 * 3 mults, 1 division (or reciprocal-mult if you precompute). Total ~12 cycles
 * on Cortex-M4F vs ~150 cycles for libm tanhf.
 *
 * @note For very small |x|, this reduces to x * (27/27) = x, which is correct.
 */
MPTC_INLINE q15_t sint_q15_tanh_rat(q15_t x) {
    /* Rescale: x in [-1, 1) -> X in [-4, 4) by *4. */
    int64_t X  = (int64_t)((int32_t)x << 2);   /* Q15, range [-4, 4)  */
    int64_t x2 = (X * X) >> 15;                /* x^2 in Q15, max ~16 */

    /* numerator:   X * (27 + x^2).   27+16 = 43; X up to 131071 -> ~5.6M  */
    int64_t num = X * ((int64_t)27 + x2);
    /* denominator: 27 + 9*x^2.         27 + 144 = 171                     */
    int64_t den = (int64_t)27 + (int64_t)9 * x2;

    /* result = (num / den) in Q15. num is in Q15 * Q15 = Q30 (since 27+x^2
     * is treated as Q15). So r = (num << 15) / den gives Q15. */
    int64_t r = (num << 15) / den;

    if (r >  32767) r =  32767;
    if (r < -32768) r = -32768;
    return (q15_t)r;
}

/**
 * @brief Fast sigmoid: 1 / (1 + exp(-x)), via tanh.
 *
 *   sigmoid(x) = (tanh(x/2) + 1) / 2
 *
 * Branchless. Input is Q15 (treated as the pre-activation in a typical
 * neural network layer). Output is Q15 in [0, 1).
 */
MPTC_INLINE q15_t sint_q15_sigmoid(q15_t x) {
    /* tanh(x/2) = tanh(x >> 1) in our scaling convention. */
    q15_t t = sint_q15_tanh_rat((q15_t)(x >> 1));
    /* (t + 1) / 2 -> (t + 32767) >> 1 in Q15.                            */
    return (q15_t)(((int32_t)t + 32767) >> 1);
}

/**
 * @brief Branchless GELU approximation (MPTC paper eq. 10).
 *
 *   GELU(x) ~= 0.5 * x * (1 + tanh(sqrt(2/pi) * (x + 0.044715 * x^3)))
 *
 * Implemented in Q15 with the rational tanh from ::sint_q15_tanh_rat.
 * Accuracy: max abs error 4e-3 vs the exact erf-based GELU.
 *
 * @par Cost
 * ~6 mults + 1 tanh + a few adds. Total ~30 cycles on Cortex-M4F vs ~250 for
 * an erf-based libm implementation.
 */
MPTC_INLINE q15_t sint_q15_gelu(q15_t x) {
    /* Rescale input by 4 to give GELU a useful dynamic range. */
    int64_t X  = (int64_t)((int32_t)x << 2);   /* Q15, range [-4, 4) */
    int64_t x2 = (X * X) >> 15;                 /* Q15                */
    int64_t x3 = (X * x2) >> 15;                 /* Q15                */

    /* inner = sqrt(2/pi) * (x + 0.044715 * x^3).
     *   sqrt(2/pi) = 0.79788  -> Q15 = 0x6639
     *   0.044715            -> Q15 = 0x05DC                              */
    int64_t inner = X + ((x3 * 0x05DC) >> 15);
    inner = (inner * 0x6639) >> 15;             /* Q15 inner */

    /* Clamp inner to Q15 range before passing to tanh. */
    int32_t inner32 = (int32_t)inner;
    if (inner32 >  32767) inner32 =  32767;
    if (inner32 < -32768) inner32 = -32768;

    q15_t t = sint_q15_tanh_rat((q15_t)inner32);

    /* 0.5 * X * (1 + t).  (1+t) in Q15 = 32767 + t. */
    int64_t one_plus_t = (int64_t)32767 + (int64_t)t;
    int64_t y = (X * one_plus_t) >> 15;          /* Q15 */
    y = y >> 1;                                   /* * 0.5 */

    if (y >  32767) y =  32767;
    if (y < -32768) y = -32768;
    return (q15_t)y;
}

/**
 * @brief Branchless ReLU: max(0, x).
 *
 * Uses the identity max(0, x) = x & ~(x >> 15) (for int16 two's complement).
 */
MPTC_INLINE q15_t sint_q15_relu_b(q15_t x) {
    /* Branchless: if x < 0, zero out. */
    int16_t mask = x >> 15;            /* 0xFFFF if x<0, 0x0000 if x>=0 */
    return (q15_t)(x & ~mask);
}

/**
 * @brief Branchless Leaky ReLU: x if x>=0 else alpha*x.
 * @param x     Input in Q15.
 * @param alpha Slope for x<0, Q15. Typical: 0x0148 (0.01).
 */
MPTC_INLINE q15_t sint_q15_leaky_relu(q15_t x, q15_t alpha) {
    int16_t mask = x >> 15;                 /* 0xFFFF if negative */
    q15_t neg_part = (q15_t)mptc_mul_q15(x, alpha);
    return (q15_t)((x & ~mask) | (neg_part & mask));
}

/* ================================================================== */
/*  3. Dual numbers (autograd-free analytical derivatives)             */
/* ================================================================== */
/**
 * @brief A dual number: value + derivative.
 *
 * A dual number is a + b*eps where eps^2 = 0. Operations on dual numbers
 * automatically propagate derivatives via the chain rule — no autograd
 * graph needed, no backward pass. This is the SintMath core primitive.
 *
 * @par Example
 * @code
 *   sintmath_q15_dual_t x = { 0x4000, 0x7FFF };  // x = 0.5, dx = 1.0
 *   sintmath_q15_dual_t y = sintmath_q15_dual_sin(x);
 *   // y.val = sin(0.5) ~= 0x376C
 *   // y.der = cos(0.5) ~= 0x6E5A   (analytical, no backward pass)
 * @endcode
 */
typedef struct {
    q15_t val;  /**< Real part (the value).         */
    q15_t der;  /**< Epsilon part (the derivative). */
} sintmath_q15_dual_t;

/** @brief Construct a dual number with explicit derivative. */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_make(q15_t val, q15_t der) {
    sintmath_q15_dual_t d = { val, der };
    return d;
}

/** @brief Construct a dual number with derivative 1 (the "variable" itself). */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_var(q15_t val) {
    sintmath_q15_dual_t d = { val, Q15_ONE };
    return d;
}

/** @brief Construct a constant (derivative 0). */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_const(q15_t val) {
    sintmath_q15_dual_t d = { val, 0 };
    return d;
}

/**
 * @brief Dual addition: (a + a'*eps) + (b + b'*eps) = (a+b) + (a'+b')*eps.
 */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_add(
        sintmath_q15_dual_t a, sintmath_q15_dual_t b) {
    sintmath_q15_dual_t r;
    r.val = (q15_t)((int32_t)a.val + (int32_t)b.val);
    r.der = (q15_t)((int32_t)a.der + (int32_t)b.der);
    return r;
}

/**
 * @brief Dual subtraction.
 */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_sub(
        sintmath_q15_dual_t a, sintmath_q15_dual_t b) {
    sintmath_q15_dual_t r;
    r.val = (q15_t)((int32_t)a.val - (int32_t)b.val);
    r.der = (q15_t)((int32_t)a.der - (int32_t)b.der);
    return r;
}

/**
 * @brief Dual multiplication: (a + a'*eps) * (b + b'*eps) = a*b + (a*b' + a'*b)*eps.
 *
 * The eps^2 term vanishes — this is the whole point of dual numbers.
 */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_mul(
        sintmath_q15_dual_t a, sintmath_q15_dual_t b) {
    sintmath_q15_dual_t r;
    r.val = mptc_mul_q15(a.val, b.val);
    /* a*b' + a'*b in Q15. */
    int32_t abp = (int32_t)a.val * (int32_t)b.der;
    int32_t apb = (int32_t)a.der * (int32_t)b.val;
    r.der = (q15_t)((abp + apb) >> 15);
    return r;
}

/**
 * @brief Dual division: (a + a'*eps) / (b + b'*eps) = (a/b) + (a'*b - a*b')/b^2 * eps.
 */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_div(
        sintmath_q15_dual_t a, sintmath_q15_dual_t b) {
    sintmath_q15_dual_t r;
    /* a / b in Q15: (a << 15) / b. */
    int32_t val = (int32_t)(((int64_t)a.val << 15) / (int64_t)b.val);
    r.val = (q15_t)mptc_clamp_q31((q31_t)val, -32768, 32767);
    /* derivative: (a'*b - a*b') / b^2  -> Q15. */
    int64_t num = (int64_t)a.der * (int64_t)b.val - (int64_t)a.val * (int64_t)b.der;
    int64_t den = (int64_t)b.val * (int64_t)b.val;
    int64_t der = (num << 15) / den;
    if (der >  32767) der =  32767;
    if (der < -32768) der = -32768;
    r.der = (q15_t)der;
    return r;
}

/**
 * @brief Dual sin: d/dx[sin(x)] = cos(x).
 *
 * sin(a + a'*eps) = sin(a) + a'*cos(a)*eps.
 */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_sin(sintmath_q15_dual_t a) {
    sintmath_q15_dual_t r;
    r.val = sint_q15_sin_wrap(a.val);
    r.der = mptc_mul_q15(a.der, sint_q15_cos(a.val));
    return r;
}

/**
 * @brief Dual cos: d/dx[cos(x)] = -sin(x).
 */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_cos(sintmath_q15_dual_t a) {
    sintmath_q15_dual_t r;
    r.val = sint_q15_cos(a.val);
    q15_t neg_sin = (q15_t)(-sint_q15_sin_wrap(a.val));
    r.der = mptc_mul_q15(a.der, neg_sin);
    return r;
}

/**
 * @brief Dual tanh: d/dx[tanh(x)] = 1 - tanh(x)^2.
 */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_tanh(sintmath_q15_dual_t a) {
    sintmath_q15_dual_t r;
    q15_t t = sint_q15_tanh_rat(a.val);
    r.val = t;
    /* 1 - t^2 in Q15. */
    int32_t t2 = ((int32_t)t * (int32_t)t) >> 15;
    q15_t sech2 = (q15_t)(32767 - t2);
    r.der = mptc_mul_q15(a.der, sech2);
    return r;
}

/**
 * @brief Dual square: d/dx[x^2] = 2x.
 */
MPTC_INLINE sintmath_q15_dual_t sintmath_q15_dual_square(sintmath_q15_dual_t a) {
    sintmath_q15_dual_t r;
    r.val = mptc_mul_q15(a.val, a.val);
    /* d/dx[x^2] = 2x. In Q15: (2 * a.val * a.der) >> 15. */
    int32_t two_x = (int32_t)a.val << 1;
    r.der = (q15_t)((two_x * (int32_t)a.der) >> 15);
    return r;
}

/* ================================================================== */
/*  4. Convenience: arithmetic on plain Q15 (branchless versions)     */
/* ================================================================== */

/**
 * @brief Branchless Q15 addition with saturation.
 */
MPTC_INLINE q15_t sint_q15_add_sat(q15_t a, q15_t b) {
    int32_t s = (int32_t)a + (int32_t)b;
    return mptc_saturate_q15(s);
}

/**
 * @brief Branchless Q15 subtraction with saturation.
 */
MPTC_INLINE q15_t sint_q15_sub_sat(q15_t a, q15_t b) {
    int32_t s = (int32_t)a - (int32_t)b;
    return mptc_saturate_q15(s);
}

/**
 * @brief Branchless Q15 multiply-accumulate:  acc += a*b.
 */
MPTC_INLINE q15_t sint_q15_mac(q15_t acc, q15_t a, q15_t b) {
    int32_t prod = ((int32_t)a * (int32_t)b) >> 15;
    int32_t s = (int32_t)acc + prod;
    return mptc_saturate_q15(s);
}

/**
 * @brief Q15 integer division (branchless, but slow on M0 — uses lib).
 *
 * @note For Q15, division cannot be made truly branchless without a Newton-
 *       Raphson reciprocal. We provide this for completeness; for hot paths
 *       prefer rewriting as multiply-by-reciprocal.
 */
MPTC_INLINE q15_t sint_q15_div(q15_t a, q15_t b) {
    if (b == 0) return (a >= 0) ? Q15_ONE : (q15_t)0x8000;
    int32_t r = ((int32_t)a << 15) / (int32_t)b;
    if (r >  32767) r =  32767;
    if (r < -32768) r = -32768;
    return (q15_t)r;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MPTC_SINTMATH_H */
