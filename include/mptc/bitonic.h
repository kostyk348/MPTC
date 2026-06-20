/**
 * @file bitonic.h
 * @brief Branchless bitonic sorting networks for N = 4, 8, 16, 32, 64.
 *
 * A bitonic sorting network is a fixed sequence of compare-and-swap
 * operations that sorts any input of size N (power of 2). Because the
 * sequence is data-independent, the implementation has:
 *
 *   - No branches in the hot path
 *   - Constant execution time (data-independent)
 *   - No data-dependent memory accesses (side-channel safe)
 *
 * This is critical for:
 *   1. Real-time systems — deterministic latency.
 *   2. Cryptography — timing-attack resistance.
 *   3. SIMD / pipelined cores — no pipeline stalls.
 *
 * We also provide median-of-N helpers built on top of the sorting networks.
 *
 * @par References
 *   - K. E. Batcher, "Sorting networks and their applications", 1968.
 *   - The MPTC paper §2.2 "Phase-Ordered Sorting (Bitonic Wave)".
 */
#ifndef MPTC_BITONIC_H
#define MPTC_BITONIC_H

#include <stdint.h>
#include "mptc_config.h"
#include "mptc_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------ */
/*  Branchless compare-and-swap primitive                              */
/* ------------------------------------------------------------------ */

/**
 * @brief Branchless compare-and-swap: ensures a <= b after call.
 *
 * Implemented as:
 *   q15_t mn = min(a, b);
 *   q15_t mx = max(a, b);
 *   *pa = mn; *pb = mx;
 *
 * The min/max themselves use conditional moves (CMOV on x86, SEL on ARMv7-M,
 * IT blocks on ARMv6-M). On all targets this is branchless.
 */
MPTC_INLINE void bitonic_cas_q15(q15_t *pa, q15_t *pb) {
    q15_t a = *pa;
    q15_t b = *pb;
    *pa = mptc_min_q15(a, b);
    *pb = mptc_max_q15(a, b);
}

/**
 * @brief Branchless compare-and-swap, descending order.
 */
MPTC_INLINE void bitonic_cas_q15_desc(q15_t *pa, q15_t *pb) {
    q15_t a = *pa;
    q15_t b = *pb;
    *pa = mptc_max_q15(a, b);
    *pb = mptc_min_q15(a, b);
}

/* ------------------------------------------------------------------ */
/*  Sorting networks (power-of-two sizes)                              */
/* ------------------------------------------------------------------ */

/**
 * @brief Branchless sort of 4 Q15 values in ascending order.
 *
 * Uses Batcher's bitonic network for N=4 (3 stages, 5 CASes total).
 */
MPTC_INLINE void bitonic_q15_sort4(q15_t a[4]) {
    bitonic_cas_q15(&a[0], &a[1]);
    bitonic_cas_q15(&a[2], &a[3]);
    bitonic_cas_q15(&a[0], &a[2]);
    bitonic_cas_q15(&a[1], &a[3]);
    bitonic_cas_q15(&a[1], &a[2]);
}

/**
 * @brief Branchless sort of 8 Q15 values.
 *
 * Uses odd-even transposition sort (a valid sorting network with N stages
 * of N/2 CASes). Total: 8 stages * 4 CASes = 32 CASes. Slightly slower than
 * optimal bitonic (19 CASes) but trivially correct and still fully branchless.
 */
MPTC_INLINE void bitonic_q15_sort8(q15_t a[8]) {
    /* 8 stages of odd-even transposition. */
    for (int stage = 0; stage < 8; ++stage) {
        if (stage & 1) {
            /* Odd stage: pairs (1,2), (3,4), (5,6) */
            bitonic_cas_q15(&a[1], &a[2]);
            bitonic_cas_q15(&a[3], &a[4]);
            bitonic_cas_q15(&a[5], &a[6]);
        } else {
            /* Even stage: pairs (0,1), (2,3), (4,5), (6,7) */
            bitonic_cas_q15(&a[0], &a[1]);
            bitonic_cas_q15(&a[2], &a[3]);
            bitonic_cas_q15(&a[4], &a[5]);
            bitonic_cas_q15(&a[6], &a[7]);
        }
    }
}

/**
 * @brief Branchless sort of 16 Q15 values.
 *
 * Odd-even transposition sort: 16 stages of 8 CASes each = 128 CASes.
 */
MPTC_INLINE void bitonic_q15_sort16(q15_t a[16]) {
    for (int stage = 0; stage < 16; ++stage) {
        for (int i = (stage & 1); i + 1 < 16; i += 2) {
            bitonic_cas_q15(&a[i], &a[i + 1]);
        }
    }
}

/**
 * @brief Branchless sort of 32 Q15 values.
 *
 * Odd-even transposition sort: 32 stages of 16 CASes each = 512 CASes.
 */
MPTC_INLINE void bitonic_q15_sort32(q15_t a[32]) {
    for (int stage = 0; stage < 32; ++stage) {
        for (int i = (stage & 1); i + 1 < 32; i += 2) {
            bitonic_cas_q15(&a[i], &a[i + 1]);
        }
    }
}

/**
 * @brief Branchless sort of 64 Q15 values.
 *
 * Odd-even transposition sort: 64 stages of 32 CASes each = 2048 CASes.
 */
MPTC_INLINE void bitonic_q15_sort64(q15_t a[64]) {
    for (int stage = 0; stage < 64; ++stage) {
        for (int i = (stage & 1); i + 1 < 64; i += 2) {
            bitonic_cas_q15(&a[i], &a[i + 1]);
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Median filters (built on top of sorting networks)                 */
/* ------------------------------------------------------------------ */

/**
 * @brief 3x3 median filter for a 9-element window.
 *
 * Sorts the 9 values (using the 8-element network plus one insert) and
 * returns the middle (index 4) value. Used for image denoising.
 */
MPTC_INLINE q15_t bitonic_q15_median9(q15_t in[9]) {
    /* Sort first 8 with the network. */
    q15_t buf[8];
    for (int i = 0; i < 8; ++i) buf[i] = in[i];
    bitonic_q15_sort8(buf);
    /* Insert in[8] into the sorted 8 to get a sorted 9. */
    q15_t x = in[8];
    /* Branchless insertion: find position via 8 CASes. */
    /* We insert at the end and bubble down. */
    q15_t sorted9[9];
    for (int i = 0; i < 8; ++i) sorted9[i] = buf[i];
    sorted9[8] = x;
    for (int i = 8; i > 0; --i) {
        bitonic_cas_q15(&sorted9[i - 1], &sorted9[i]);
    }
    return sorted9[4];
}

/**
 * @brief 5x5 median filter (25 elements).
 *
 * Uses a hybrid: bucket-sort-friendly since Q15 values fit in int16.
 * Here we use the 32-element network with padding.
 */
MPTC_INLINE q15_t bitonic_q15_median25(q15_t in[25]) {
    q15_t buf[32];
    /* Fill 25 values, pad remaining 7 with INT16_MAX so they sort to the end. */
    for (int i = 0; i < 25; ++i) buf[i] = in[i];
    for (int i = 25; i < 32; ++i) buf[i] = 32767;
    bitonic_q15_sort32(buf);
    return buf[12];  /* median of 25 is index 12 (0-indexed). */
}

/**
 * @brief Top-k selection via sort. Returns the k-th largest element.
 *
 * Sorts all N and returns element[N - k]. N must be 4, 8, 16, 32, or 64.
 * For large N or repeated queries, a partial sort would be faster —
 * but for branchless guarantees, the full sort is preferable.
 *
 * @param arr  Array of N Q15 values (will be modified in place).
 * @param N    Must be 4, 8, 16, 32, or 64.
 * @param k    1-indexed: 1 = largest, N = smallest.
 * @return     The k-th largest value.
 */
MPTC_INLINE q15_t bitonic_q15_topk(q15_t *arr, int N, int k) {
    switch (N) {
        case 4:  bitonic_q15_sort4(arr);  break;
        case 8:  bitonic_q15_sort8(arr);  break;
        case 16: bitonic_q15_sort16(arr); break;
        case 32: bitonic_q15_sort32(arr); break;
        case 64: bitonic_q15_sort64(arr); break;
        default: return 0;
    }
    return arr[N - k];
}

/* ------------------------------------------------------------------ */
/*  Generic (non-power-of-two) sorting via network embedding           */
/* ------------------------------------------------------------------ */

/**
 * @brief Branchless sort of N Q15 values, N in [1, 64].
 *
 * Pads to next power of two with INT16_MAX, sorts, returns.
 * Useful when N is known only at runtime but bounded.
 */
MPTC_INLINE void bitonic_q15_sort_n(q15_t *arr, int N) {
    if (N <= 1) return;
    /* Pad to next power of two. */
    q15_t buf[64];
    for (int i = 0; i < N; ++i) buf[i] = arr[i];
    for (int i = N; i < 64; ++i) buf[i] = 32767;

    if      (N <= 4)  bitonic_q15_sort4(buf);
    else if (N <= 8)  bitonic_q15_sort8(buf);
    else if (N <= 16) bitonic_q15_sort16(buf);
    else if (N <= 32) bitonic_q15_sort32(buf);
    else              bitonic_q15_sort64(buf);

    for (int i = 0; i < N; ++i) arr[i] = buf[i];
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MPTC_BITONIC_H */
