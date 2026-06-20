/**
 * @file integral_img.c
 * @brief Summed Area Table (Integral Image) implementation.
 */
#include "mptc/integral_img.h"
#include "mptc/sintmath.h"
#include <stdlib.h>
#include <string.h>

integral_q15_t *integral_q15_build(const uint8_t *img,
                                    uint32_t width, uint32_t height) {
    if (!img || width == 0 || height == 0) return NULL;
    integral_q15_t *ii = (integral_q15_t*)malloc(sizeof(integral_q15_t));
    if (!ii) return NULL;
    ii->width = width;
    ii->height = height;
    size_t n = (size_t)(width + 1u) * (size_t)(height + 1u);
    ii->I   = (uint32_t*)calloc(n, sizeof(uint32_t));
    ii->Isq = (uint32_t*)calloc(n, sizeof(uint32_t));
    if (!ii->I || !ii->Isq) {
        free(ii->I); free(ii->Isq); free(ii);
        return NULL;
    }
    integral_q15_build_into(ii, img);
    return ii;
}

void integral_q15_attach(integral_q15_t *ii,
                         uint32_t *I, uint32_t *Isq,
                         uint32_t width, uint32_t height) {
    ii->width = width;
    ii->height = height;
    ii->I = I;
    ii->Isq = Isq;
}

void integral_q15_build_into(integral_q15_t *ii, const uint8_t *img) {
    uint32_t W = ii->width;
    uint32_t H = ii->height;
    uint32_t stride = W + 1u;
    uint32_t *I   = ii->I;
    uint32_t *Isq = ii->Isq;

    /* Zero the first row and column (already zeroed if calloc'd). */
    for (uint32_t x = 0; x <= W; ++x) { I[x] = 0; Isq[x] = 0; }

    for (uint32_t y = 1; y <= H; ++y) {
        I[y * stride]   = 0;
        Isq[y * stride] = 0;
        uint32_t row_sum   = 0;
        uint32_t row_sum_sq = 0;
        for (uint32_t x = 1; x <= W; ++x) {
            uint8_t p = img[(y - 1u) * W + (x - 1u)];
            uint32_t psq = (uint32_t)p * (uint32_t)p;
            row_sum    += p;
            row_sum_sq += psq;
            I[y * stride + x]   = I[(y - 1u) * stride + x] + row_sum;
            Isq[y * stride + x] = Isq[(y - 1u) * stride + x] + row_sum_sq;
        }
    }
}

void integral_q15_free(integral_q15_t *ii) {
    if (!ii) return;
    free(ii->I);
    free(ii->Isq);
    free(ii);
}

uint8_t integral_q15_sauvola(const integral_q15_t *ii,
                              int x, int y, int w, int h, q15_t k) {
    uint32_t sum = integral_q15_window_sum(ii, x, y, w, h);
    uint32_t n = (uint32_t)(w * h);
    uint32_t mean = sum / n;
    /* Variance (scaled by n to avoid the divide): var*n = sum_sq - sum^2/n. */
    uint32_t stride = ii->width + 1u;
    const uint32_t *Isq = ii->Isq;
    uint32_t A = Isq[(uint32_t)y * stride + (uint32_t)x];
    uint32_t B = Isq[(uint32_t)y * stride + (uint32_t)(x + w)];
    uint32_t C = Isq[(uint32_t)(y + h) * stride + (uint32_t)x];
    uint32_t D = Isq[(uint32_t)(y + h) * stride + (uint32_t)(x + w)];
    uint32_t sum_sq = D - B - C + A;
    /* var = sum_sq/n - mean^2 */
    int64_t var = (int64_t)sum_sq - (int64_t)sum * (int64_t)sum / (int64_t)n;
    if (var < 0) var = 0;
    /* std = sqrt(var).  Use integer sqrt. */
    int32_t std = 0;
    int32_t v = (int32_t)var;
    if (v > 0) {
        int32_t r = v;
        int32_t t = 0;
        /* Classic integer sqrt. */
        for (int32_t bit = 1 << 15; bit > 0; bit >>= 1) {
            int32_t nt = t + bit;
            if (nt * nt <= r) t = nt;
        }
        std = t;
    }
    /* Sauvola: T = mean * (1 + k * (std / R - 1)),  R = 128. */
    /* (std / 128 - 1) in Q15: (std << 7) - 32768 ... but we need it as Q15.   */
    int32_t ratio = (int32_t)std - 128;   /* std - R */
    int32_t term  = mptc_mul_q15(k, (q15_t)mptc_clamp_q15(ratio, -32768, 32767));
    int32_t factor = 32767 + term;        /* (1 + k*(std/R - 1)) in Q15 */
    int32_t T = (int32_t)mean * factor;   /* Q15 * 0..255 -> needs shift */
    T = T >> 15;
    if (T > 255) T = 255;
    if (T < 0)   T = 0;
    return (uint8_t)T;
}

void integral_q15_sauvola_full(const uint8_t *in, uint8_t *out,
                                uint32_t W, uint32_t H,
                                int win, q15_t k) {
    integral_q15_t ii_storage;
    integral_q15_t *ii = &ii_storage;
    size_t n = (size_t)(W + 1u) * (size_t)(H + 1u);
    /* Allocate on heap to avoid stack overflow for large images. */
    ii->I   = (uint32_t*)calloc(n, sizeof(uint32_t));
    ii->Isq = (uint32_t*)calloc(n, sizeof(uint32_t));
    ii->width = W;
    ii->height = H;
    if (!ii->I || !ii->Isq) {
        free(ii->I); free(ii->Isq);
        return;
    }
    integral_q15_build_into(ii, in);

    int wh = 2 * win + 1;  /* full window size */
    for (uint32_t y = 0; y < H; ++y) {
        int wy = (int)y - win;
        int hy = wh;
        if (wy < 0) { hy += wy; wy = 0; }
        if (wy + hy > (int)H) hy = (int)H - wy;
        for (uint32_t x = 0; x < W; ++x) {
            int wx = (int)x - win;
            int hw = wh;
            if (wx < 0) { hw += wx; wx = 0; }
            if (wx + hw > (int)W) hw = (int)W - wx;
            uint8_t T = integral_q15_sauvola(ii, wx, wy, hw, hy, k);
            uint8_t p = in[y * W + x];
            out[y * W + x] = (p > T) ? 255 : 0;
        }
    }
    free(ii->I);
    free(ii->Isq);
}
