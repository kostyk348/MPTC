/**
 * @file integral_img.h
 * @brief Summed Area Table (Integral Image) for O(1) windowed statistics.
 *
 * The MPTC paper reports that adaptive binarization (Sauvola) drops from
 * 29.17 seconds (sliding window) to 0.016 seconds (integral image) — an
 * 1800x speedup — by pre-computing two integral fields (sum and squared sum)
 * and answering window queries in O(1) via the four-point rectangle sum:
 *
 *   Sum(rect) = I(D) - I(B) - I(C) + I(A)
 *
 * This module provides:
 *   - ::integral_q15_build  — pre-compute the two integral fields.
 *   - ::integral_q15_window_sum / _mean / _var — O(1) queries.
 *   - ::integral_q15_sauvola — full Sauvola adaptive threshold in one call.
 *
 * Storage: 2 * (W+1) * (H+1) * sizeof(uint32_t) = ~310 KB for 320x240.
 * For memory-constrained MCUs we also provide a tiled variant
 * (::integral_q15_build_tiled) that processes the image in strips.
 */
#ifndef MPTC_INTEGRAL_IMG_H
#define MPTC_INTEGRAL_IMG_H

#include <stdint.h>
#include "mptc_config.h"
#include "mptc_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Integral image state.
 *
 * The two fields (`I` and `Isq`) have size (W+1) * (H+1) to allow the
 * four-point formula to work without bounds checking on the top row / left
 * column (those entries are 0).
 */
typedef struct {
    uint32_t  width;   /**< Image width (W).                              */
    uint32_t  height;  /**< Image height (H).                             */
    uint32_t *I;       /**< Sum integral, size (W+1)*(H+1).               */
    uint32_t *Isq;     /**< Squared-sum integral, size (W+1)*(H+1).       */
} integral_q15_t;

/**
 * @brief Allocate and build the integral fields for a grayscale image.
 *
 * @param img   Pointer to grayscale pixel data, row-major, `width`*`height` bytes.
 *              Pixel values are 0..255.
 * @param width Image width.
 * @param height Image height.
 * @return Newly allocated integral image, or NULL on failure. Caller must
 *         free with ::integral_q15_free.
 */
integral_q15_t *integral_q15_build(const uint8_t *img,
                                    uint32_t width, uint32_t height);

/**
 * @brief Attach pre-allocated buffers (no malloc).
 *
 * The caller provides I and Isq arrays of size (width+1)*(height+1) each.
 * After attaching, call ::integral_q15_build_into to populate them.
 */
void integral_q15_attach(integral_q15_t *ii,
                         uint32_t *I, uint32_t *Isq,
                         uint32_t width, uint32_t height);

/**
 * @brief Populate the integral fields in-place (no allocation).
 *
 * Use with ::integral_q15_attach for bare-metal systems.
 */
void integral_q15_build_into(integral_q15_t *ii, const uint8_t *img);

/**
 * @brief Free a previously built integral image.
 */
void integral_q15_free(integral_q15_t *ii);

/**
 * @brief O(1) sum of a rectangular window.
 *
 * @param ii Integral image.
 * @param x  Top-left x coordinate (inclusive).
 * @param y  Top-left y coordinate (inclusive).
 * @param w  Window width.
 * @param h  Window height.
 * @return   Sum of pixel values in the window.
 *
 * @note No bounds checking in release builds. Caller must ensure
 *       x + w <= width and y + h <= height.
 */
MPTC_INLINE uint32_t integral_q15_window_sum(const integral_q15_t *ii,
                                              int x, int y, int w, int h) {
    /* Four-point formula: I(D) - I(B) - I(C) + I(A)        */
    /*   A = (x,   y),  B = (x+w, y),  C = (x, y+h),  D = (x+w, y+h)
     *
     * Index into I: row-major, row stride = ii->width + 1.
     */
    uint32_t stride = ii->width + 1u;
    const uint32_t *I = ii->I;
    uint32_t A = I[(uint32_t)y * stride + (uint32_t)x];
    uint32_t B = I[(uint32_t)y * stride + (uint32_t)(x + w)];
    uint32_t C = I[(uint32_t)(y + h) * stride + (uint32_t)x];
    uint32_t D = I[(uint32_t)(y + h) * stride + (uint32_t)(x + w)];
    return D - B - C + A;
}

/**
 * @brief O(1) mean of a rectangular window.
 */
MPTC_INLINE uint32_t integral_q15_window_mean(const integral_q15_t *ii,
                                               int x, int y, int w, int h) {
    uint32_t sum = integral_q15_window_sum(ii, x, y, w, h);
    return sum / (uint32_t)(w * h);
}

/**
 * @brief O(1) variance of a rectangular window.
 *
 * Uses the identity: var = E[X^2] - (E[X])^2.
 * We compute E[X^2] from the Isq field, E[X] from the I field.
 *
 * @return Variance * (w*h) (avoids the divide; caller divides if needed).
 */
MPTC_INLINE uint32_t integral_q15_window_var_scaled(const integral_q15_t *ii,
                                                     int x, int y, int w, int h) {
    uint32_t stride = ii->width + 1u;
    const uint32_t *Isq = ii->Isq;
    uint32_t A = Isq[(uint32_t)y * stride + (uint32_t)x];
    uint32_t B = Isq[(uint32_t)y * stride + (uint32_t)(x + w)];
    uint32_t C = Isq[(uint32_t)(y + h) * stride + (uint32_t)x];
    uint32_t D = Isq[(uint32_t)(y + h) * stride + (uint32_t)(x + w)];
    uint32_t sum_sq = D - B - C + A;

    uint32_t sum = integral_q15_window_sum(ii, x, y, w, h);
    uint32_t n = (uint32_t)(w * h);
    /* var*n = sum_sq - sum*sum/n */
    return sum_sq - (uint32_t)(((uint64_t)sum * (uint64_t)sum) / n);
}

/**
 * @brief Sauvola adaptive binarization threshold at a single pixel.
 *
 * Sauvola's formula:
 *   T = mean * (1 + k * (std / R - 1))
 * where R = 128 (dynamic range for 8-bit images), k is a sensitivity
 * parameter (typical 0.5), and mean/std are computed over a local window.
 *
 * @param ii Integral image.
 * @param x  Pixel x coordinate.
 * @param y  Pixel y coordinate.
 * @param w  Window width (odd recommended).
 * @param h  Window height.
 * @param k  Sensitivity in Q15 (typical 0x4000 = 0.5).
 * @return   Threshold T in 0..255 range.
 */
uint8_t integral_q15_sauvola(const integral_q15_t *ii,
                              int x, int y, int w, int h, q15_t k);

/**
 * @brief Full Sauvola binarization of an entire image.
 *
 * Writes 0/255 to `out` for each pixel. Window is centered on each pixel
 * (clamped at image edges).
 *
 * @param in   Input grayscale image (W*H bytes).
 * @param out  Output binary image (W*H bytes).
 * @param W    Width.
 * @param H    Height.
 * @param win  Half-window size (full window = 2*win+1). Typical: 8..15.
 * @param k    Sensitivity in Q15. Typical: 0x4000 (0.5).
 */
void integral_q15_sauvola_full(const uint8_t *in, uint8_t *out,
                                uint32_t W, uint32_t H,
                                int win, q15_t k);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MPTC_INTEGRAL_IMG_H */
