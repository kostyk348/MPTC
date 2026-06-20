/**
 * @file ternary.h
 * @brief Balanced ternary {-1, 0, +1} packed arithmetic.
 *
 * Balanced ternary uses three states: -1, 0, +1. One trit carries
 * log2(3) ~= 1.585 bits of information — denser than binary. The MPTC
 * hardware paper notes this gives "higher information density than binary
 * (one trit ~ 1.58 bit)".
 *
 * This module packs 9 trits into a uint16_t (15 bits used, MSB free for
 * sign in mixed use) using a 2-bit-per-trit encoding:
 *
 *   trit value    encoding (2 bits)
 *     -1            10
 *      0            00
 *     +1            01
 *
 * This encoding has the useful property that bit[0] is the magnitude and
 * bit[1] is the sign, so addition / multiplication table reduces to a few
 * bit operations. The "11" code is unused (reserved for "don't care" in
 * sparse weights).
 *
 * 9 trits fit in 18 bits but we round down to 8 trits in 16 bits (with
 * sign extension in the MSB) to keep int16_t arithmetic fast. The
 * numerical range of 8 trits is [-3280, +3280] in balanced ternary
 * (3^8 / 2 = 3280.5).
 *
 * 18 trits in 36 bits fit in int64_t (range ~ +-193M).
 */
#ifndef MPTC_TERNARY_H
#define MPTC_TERNARY_H

#include <stdint.h>
#include "mptc_config.h"
#include "mptc_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/*  Type aliases                                                       */
/* ------------------------------------------------------------------ */

/** @brief 8 trits packed in a uint16_t. Range: [-3280, +3280]. */
typedef uint16_t trit8_t;

/** @brief 18 trits packed in a uint32_t (uses 36 bits, MSB unused). Range: +-193M. */
typedef uint64_t trit18_t;

/* ------------------------------------------------------------------ */
/*  Trit extraction / insertion                                        */
/* ------------------------------------------------------------------ */

/**
 * @brief Get the trit at position `pos` (0=LSB) of a packed trit8.
 *
 * @return -1, 0, or +1.
 */
MPTC_INLINE int8_t ternary_get_trit8(trit8_t t, int pos) {
    /* Encoding (2 bits per trit):
     *   0b00 ->  0
     *   0b01 -> +1
     *   0b11 -> -1
     *   0b10 ->  0 (reserved, treat as 0)
     */
    uint8_t bits = (uint8_t)((t >> (pos * 2)) & 0x3);
    int8_t magnitude = (int8_t)(bits & 0x1);
    int8_t sign      = (int8_t)((bits >> 1) & 0x1);
    /* value = magnitude * (1 if !sign else -1) = magnitude - 2*magnitude*sign */
    int8_t v = (int8_t)(magnitude - 2 * magnitude * sign);
    return v;
}

/**
 * @brief Set the trit at position `pos` of a packed trit8.
 *
 * @param t   Packed value.
 * @param pos Position 0..7.
 * @param v   Trit value: must be -1, 0, or +1.
 * @return    Updated packed value.
 */
MPTC_INLINE trit8_t ternary_set_trit8(trit8_t t, int pos, int8_t v) {
    uint8_t bits;
    if (v > 0)       bits = 0x1;     /* +1 = 01 */
    else if (v < 0)  bits = 0x3;     /* -1 = 11 */
    else             bits = 0x0;     /*  0 = 00 */
    uint16_t mask = (uint16_t)~((uint16_t)0x3 << (pos * 2));
    uint16_t insert = (uint16_t)((uint16_t)bits << (pos * 2));
    return (trit8_t)((t & mask) | insert);
}

/* ------------------------------------------------------------------ */
/*  Conversions                                                        */
/* ------------------------------------------------------------------ */

/**
 * @brief Convert a signed integer to a packed trit8.
 *
 * Uses balanced ternary greedy decomposition: for each power p (descending),
 * choose the trit that minimizes the residual.
 *
 * @param v Integer in [-3280, +3280]. Out-of-range values are clamped.
 */
MPTC_INLINE trit8_t ternary_from_int8(int v) {
    if (v > 3280)  v = 3280;
    if (v < -3280) v = -3280;
    trit8_t t = 0;
    static const int pow3[8] = {2187, 729, 243, 81, 27, 9, 3, 1};
    for (int i = 0; i < 8; ++i) {
        int p = pow3[i];
        int8_t trit = 0;
        int half = (p + 1) / 2;
        if (v >= half)       { trit =  1; v -= p; }
        else if (v <= -half) { trit = -1; v += p; }
        t = ternary_set_trit8(t, 7 - i, trit);
    }
    return t;
}

/**
 * @brief Convert a packed trit8 to a signed integer.
 */
MPTC_INLINE int ternary_to_int8(trit8_t t) {
    static const int pow3[8] = {2187, 729, 243, 81, 27, 9, 3, 1};
    int v = 0;
    for (int i = 0; i < 8; ++i) {
        int8_t trit = ternary_get_trit8(t, 7 - i);
        v += trit * pow3[i];
    }
    return v;
}

/* ------------------------------------------------------------------ */
/*  Arithmetic                                                         */
/* ------------------------------------------------------------------ */

/**
 * @brief Add two packed trit8 values.
 *
 * Decomposes both to int, adds, recompresses. Branchless in the inner loops
 * but uses integer arithmetic for simplicity. For a true branchless ternary
 * ALU, see the MPTC hardware paper — this software emulation prioritizes
 * clarity.
 */
MPTC_INLINE trit8_t ternary_add8(trit8_t a, trit8_t b) {
    int va = ternary_to_int8(a);
    int vb = ternary_to_int8(b);
    return ternary_from_int8(va + vb);
}

/**
 * @brief Multiply two packed trit8 values.
 */
MPTC_INLINE trit8_t ternary_mul8(trit8_t a, trit8_t b) {
    int va = ternary_to_int8(a);
    int vb = ternary_to_int8(b);
    /* Products can be up to 3280^2 = ~10M, overflow int. */
    long long prod = (long long)va * (long long)vb;
    if (prod > 3280)  prod = 3280;
    if (prod < -3280) prod = -3280;
    return ternary_from_int8((int)prod);
}

/**
 * @brief Negate a packed trit8 (swap +1 <-> -1, keep 0).
 *
 * Branchless: each 2-bit code 01<->10, 00 stays. Done by XOR with 0x3 where
 * the trit is nonzero. We approximate: XOR all bits with sign bit.
 *   01 (+1) XOR 11 = 10 (-1)  ✓
 *   10 (-1) XOR 11 = 01 (+1)  ✓
 *   00 ( 0) XOR 11 = 11 (?)   ✗
 * So we can't blindly XOR. We do a proper per-trit branchless flip.
 */
MPTC_INLINE trit8_t ternary_neg8(trit8_t a) {
    trit8_t r = 0;
    for (int i = 0; i < 8; ++i) {
        int8_t t = ternary_get_trit8(a, i);
        r = ternary_set_trit8(r, i, (int8_t)(-t));
    }
    return r;
}

/**
 * @brief Absolute value: each trit -> |trit|.
 *
 *   -1 -> +1
 *    0 ->  0
 *   +1 -> +1
 *
 * Branchless per-trit: just clear bit[1] (the sign bit) of each pair.
 *   01 -> 01  (no change)
 *   10 -> 00  (sign cleared, magnitude was 0)
 *   00 -> 00
 *
 * Wait, 10 has magnitude bit=0, sign bit=1. The encoding is non-standard.
 * Let me re-derive: in our encoding, +1=01, -1=10, 0=00.
 *   -1=10 -> we want +1=01. That's bit-swap. Not a simple clear.
 * So we need to flip the 2 bits if they're 10. That's: code XOR 0b11 if code==10.
 * Implementation: we just iterate.
 */
MPTC_INLINE trit8_t ternary_abs8(trit8_t a) {
    trit8_t r = 0;
    for (int i = 0; i < 8; ++i) {
        int8_t t = ternary_get_trit8(a, i);
        r = ternary_set_trit8(r, i, (int8_t)(t < 0 ? -t : t));
    }
    return r;
}

/* ------------------------------------------------------------------ */
/*  Ternary weight neural helpers                                      */
/* ------------------------------------------------------------------ */

/**
 * @brief Apply a ternary weight to a Q15 input.
 *
 * For ternary-weight neural networks (TWN), each weight is -1, 0, or +1.
 * This function extracts one trit from a packed weight and applies it to
 * an input:  -1 -> -x, 0 -> 0, +1 -> +x. Branchless.
 */
MPTC_INLINE q15_t ternary_apply_weight(trit8_t weights, int pos, q15_t x) {
    int8_t t = ternary_get_trit8(weights, pos);
    /* Branchless:  t==0  -> 0;  t>0 -> x;  t<0 -> -x. */
    q15_t neg_x = (q15_t)(-x);
    q15_t result;
    if (t > 0) result = x;
    else if (t < 0) result = neg_x;
    else result = 0;
    return result;
}

/**
 * @brief Dot product of 8 ternary weights with 8 Q15 inputs.
 *
 * This is the core op of a ternary-weight linear layer. Equivalent to:
 *   sum = 0;
 *   for (i in 0..7) sum += weights[i] * inputs[i];   // weights in {-1,0,1}
 *
 * @param weights Packed 8 trit8.
 * @param inputs  Array of 8 Q15 values.
 * @return        Q15 dot product (with saturation).
 */
MPTC_INLINE q15_t ternary_dot8(trit8_t weights, const q15_t *inputs) {
    int32_t acc = 0;
    for (int i = 0; i < 8; ++i) {
        int8_t t = ternary_get_trit8(weights, i);
        /* Multiply by -1, 0, +1 branchlessly using sign mask. */
        int32_t v = (int32_t)inputs[i];
        /* If t>0: +v; if t<0: -v; if t==0: 0. */
        int32_t signed_v = (t > 0) ? v : (t < 0) ? -v : 0;
        acc += signed_v;
    }
    /* Shift to Q15 (acc is sum of Q15, no extra scaling needed). */
    if (acc >  32767) acc =  32767;
    if (acc < -32768) acc = -32768;
    return (q15_t)acc;
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MPTC_TERNARY_H */
