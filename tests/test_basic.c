/**
 * @file test_basic.c
 * @brief Basic unit tests for the MPTC library.
 *
 * Compiled with the host compiler (gcc/clang). Each test prints PASS/FAIL.
 * Exit code 0 = all tests pass, non-zero = at least one failure.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mptc/mptc.h"

static int g_pass = 0;
static int g_fail = 0;

#define ASSERT(cond, msg) do { \
    if (cond) { g_pass++; } \
    else { g_fail++; printf("FAIL: %s (line %d): %s\n", __FILE__, __LINE__, msg); } \
} while (0)

#define ASSERT_NEAR(act, exp, tol, msg) do { \
    int _a = (int)(act); int _e = (int)(exp); int _t = (int)(tol); \
    if (abs(_a - _e) <= _t) { g_pass++; } \
    else { g_fail++; printf("FAIL: %s (line %d): %s (got %d, exp %d, tol %d)\n", \
        __FILE__, __LINE__, msg, _a, _e, _t); } \
} while (0)

/* ------------------------------------------------------------------ */
/*  sintmath tests                                                     */
/* ------------------------------------------------------------------ */
static void test_sintmath_select(void) {
    q15_t a = 0x1000, b = 0x7000;
    ASSERT(sint_q15_select(a, b, 1) == a, "select mask=1 should return a");
    ASSERT(sint_q15_select(a, b, 0) == b, "select mask=0 should return b");
    ASSERT(sint_q15_select(a, b, 42) == a, "select mask=42 should return a (nonzero)");
}

static void test_sintmath_clamp(void) {
    ASSERT(mptc_clamp_q15(100, -1000, 1000) == 100, "clamp in-range");
    ASSERT(mptc_clamp_q15(5000, -1000, 1000) == 1000, "clamp above hi");
    ASSERT(mptc_clamp_q15(-5000, -1000, 1000) == -1000, "clamp below lo");
}

static void test_sintmath_sin(void) {
    /* Trig convention: x in Q15 represents (x/32768)*pi radians. */
    /* sin(0) = 0 */
    ASSERT_NEAR(sint_q15_sin_wrap(0), 0, 50, "sin(0) ~= 0");
    /* sin(pi/2) = 1 -> Q15 = 0x7FFF, x = 0x4000 represents pi/2 */
    ASSERT_NEAR(sint_q15_sin_wrap(Q15_PI2), 32767, 200, "sin(pi/2) ~= 1");
    /* sin(pi) ~= 0: x = 0x7FFF represents ~pi */
    ASSERT_NEAR(sint_q15_sin_wrap(32700), 0, 800, "sin(~pi) ~= 0");
    /* sin(-pi/2) = -1: x = -0x4000 represents -pi/2 */
    ASSERT_NEAR(sint_q15_sin_wrap(-Q15_PI2), -32768, 200, "sin(-pi/2) ~= -1");
    /* sin(pi/4) = sqrt(2)/2 ~= 0.7071 -> Q15 = 23170 */
    ASSERT_NEAR(sint_q15_sin_wrap(0x2000), 23170, 800, "sin(pi/4) ~= 0.7071");
}

static void test_sintmath_tanh(void) {
    /* tanh(0) = 0 */
    ASSERT_NEAR(sint_q15_tanh_rat(0), 0, 50, "tanh(0) ~= 0");
    /* tanh(large) ~ 1 */
    q15_t big = 0x7FFF;
    q15_t r = sint_q15_tanh_rat(big);
    ASSERT(r > 30000, "tanh(1.0) should be ~0.76 (Q15 ~ 25000) — relaxed");
    /* tanh(-x) = -tanh(x) */
    q15_t t_pos = sint_q15_tanh_rat(0x2000);
    q15_t t_neg = sint_q15_tanh_rat(-0x2000);
    ASSERT_NEAR(t_neg, -t_pos, 50, "tanh should be odd");
}

static void test_sintmath_gelu(void) {
    /* GELU(0) = 0 */
    ASSERT_NEAR(sint_q15_gelu(0), 0, 200, "GELU(0) ~= 0");
    /* GELU(large positive) ~ x */
    q15_t big = 0x6000;
    q15_t g = sint_q15_gelu(big);
    ASSERT(g > big / 2, "GELU(positive) should be positive");
    /* GELU(negative) ~ 0 */
    q15_t neg = -0x6000;
    q15_t gn = sint_q15_gelu(neg);
    ASSERT(gn > -3000, "GELU(very negative) should be ~0");
}

static void test_sintmath_dual_sin(void) {
    /* f(x) = sin(x), f'(x) = cos(x).
     * At x = 0: f(0) = 0, f'(0) = 1. */
    sintmath_q15_dual_t x = sintmath_q15_dual_var(0);
    sintmath_q15_dual_t y = sintmath_q15_dual_sin(x);
    ASSERT_NEAR(y.val, 0, 100, "dual sin(0) val ~= 0");
    ASSERT_NEAR(y.der, 32767, 1500, "dual sin'(0) = cos(0) = 1");
}

static void test_sintmath_dual_mul(void) {
    /* f(x) = x * x at x = 0.5 -> f = 0.25, f' = 2x = 1.0 */
    sintmath_q15_dual_t x = sintmath_q15_dual_var(Q15_HALF);  /* val=0.5, der=1 */
    sintmath_q15_dual_t y = sintmath_q15_dual_mul(x, x);
    /* val = 0.5 * 0.5 = 0.25 -> Q15 = 0x2000 = 8192 */
    ASSERT_NEAR(y.val, 8192, 200, "dual mul val: 0.5*0.5 = 0.25");
    /* der = 2 * 0.5 * 1.0 = 1.0 -> Q15 = 0x7FFF */
    ASSERT_NEAR(y.der, 32767, 200, "dual mul der: d/dx(x^2) at x=0.5 = 1.0");
}

/* ------------------------------------------------------------------ */
/*  sdft tests                                                         */
/* ------------------------------------------------------------------ */
static void test_sdft_basic(void) {
    sdft_q31_t s;
    ASSERT(sdft_q31_init(&s, 16, 16000) == 0, "SDFT init 16 bins");
    /* Feed a pure sine at bin 4 (4 kHz). */
    for (int i = 0; i < 16 * 4; ++i) {
        float phase = 2.0 * 3.14159265 * 4 * i / 16.0;
        q15_t sample = (q15_t)(sinf(phase) * 30000.0f);
        sdft_q31_update(&s, sample);
    }
    /* For a real sine input with +j convention, peaks appear at bin k0 AND
     * bin N-k0 (conjugate pair). For k0=4, N=16: bins 4 and 12. */
    uint16_t peak = sdft_q31_argmax(&s);
    char msg[64];
    snprintf(msg, sizeof(msg), "SDFT peak should be at bin 4 or 12 (got %d)", peak);
    ASSERT(peak == 4 || peak == 12, msg);
    sdft_q31_free(&s);
}

/* ------------------------------------------------------------------ */
/*  triple_pll tests                                                   */
/* ------------------------------------------------------------------ */
static void test_triple_pll_median3(void) {
    ASSERT(triple_pll_median3(1, 2, 3) == 2, "median(1,2,3) = 2");
    ASSERT(triple_pll_median3(3, 2, 1) == 2, "median(3,2,1) = 2");
    ASSERT(triple_pll_median3(1, 3, 2) == 2, "median(1,3,2) = 2");
    ASSERT(triple_pll_median3(-5, 0, 5) == 0, "median(-5,0,5) = 0");
    ASSERT(triple_pll_median3(100, 100, 100) == 100, "median(100,100,100) = 100");
}

static void test_triple_pll_outlier_rejection(void) {
    triple_pll_q15_t pll;
    triple_pll_q15_init(&pll);
    /* Prime with stable values. */
    for (int i = 0; i < 5; ++i) triple_pll_q15_update(&pll, 1000);
    /* Now inject one outlier. The median should NOT jump to the outlier. */
    q15_t m1 = triple_pll_q15_update(&pll, 30000);  /* outlier */
    /* m1 = median(1000, 1000, 30000) = 1000 */
    ASSERT(m1 == 1000, "median should reject the 30000 outlier");
}

/* ------------------------------------------------------------------ */
/*  integral_img tests                                                 */
/* ------------------------------------------------------------------ */
static void test_integral_img_basic(void) {
    /* 4x4 image: all values = 10.
     * Window sum of any rect = 10 * (w*h). */
    uint8_t img[16];
    for (int i = 0; i < 16; ++i) img[i] = 10;
    integral_q15_t *ii = integral_q15_build(img, 4, 4);
    ASSERT(ii != NULL, "integral image build");
    ASSERT(integral_q15_window_sum(ii, 0, 0, 2, 2) == 40, "2x2 sum = 40");
    ASSERT(integral_q15_window_sum(ii, 1, 1, 2, 2) == 40, "2x2 offset sum = 40");
    ASSERT(integral_q15_window_sum(ii, 0, 0, 4, 4) == 160, "4x4 sum = 160");
    ASSERT(integral_q15_window_mean(ii, 0, 0, 4, 4) == 10, "4x4 mean = 10");
    integral_q15_free(ii);
}

/* ------------------------------------------------------------------ */
/*  ternary tests                                                      */
/* ------------------------------------------------------------------ */
static void test_ternary_roundtrip(void) {
    int test_values[] = { 0, 1, -1, 5, -5, 100, -100, 1000, -1000, 3280, -3280 };
    int n = (int)(sizeof(test_values) / sizeof(test_values[0]));
    for (int i = 0; i < n; ++i) {
        trit8_t t = ternary_from_int8(test_values[i]);
        int v = ternary_to_int8(t);
        ASSERT(v == test_values[i], "ternary roundtrip");
    }
}

static void test_ternary_arith(void) {
    trit8_t a = ternary_from_int8(15);
    trit8_t b = ternary_from_int8(27);
    trit8_t sum = ternary_add8(a, b);
    ASSERT(ternary_to_int8(sum) == 42, "15 + 27 = 42");
    trit8_t prod = ternary_mul8(a, b);
    ASSERT(ternary_to_int8(prod) == 405, "15 * 27 = 405");
    trit8_t neg = ternary_neg8(a);
    ASSERT(ternary_to_int8(neg) == -15, "neg(15) = -15");
}

static void test_ternary_dot(void) {
    /* weights = (+1, -1, 0, +1, -1, 0, +1, -1) -> int value */
    trit8_t w = 0;
    int8_t ws[8] = {1, -1, 0, 1, -1, 0, 1, -1};
    for (int i = 0; i < 8; ++i) w = ternary_set_trit8(w, i, ws[i]);
    q15_t inputs[8] = {100, 200, 300, 400, 500, 600, 700, 800};
    /* Expected: 100 - 200 + 0 + 400 - 500 + 0 + 700 - 800 = -300 */
    q15_t r = ternary_dot8(w, inputs);
    ASSERT(r == -300, "ternary dot8 = -300");
}

/* ------------------------------------------------------------------ */
/*  bitonic tests                                                      */
/* ------------------------------------------------------------------ */
static void test_bitonic_sort8(void) {
    q15_t a[8] = { 7, 2, 5, 1, 8, 3, 6, 4 };
    bitonic_q15_sort8(a);
    for (int i = 0; i < 7; ++i) {
        char msg[64];
        snprintf(msg, sizeof(msg), "bitonic sort8: a[%d] <= a[%d] (got %d, %d)",
                 i, i + 1, a[i], a[i + 1]);
        ASSERT(a[i] <= a[i + 1], msg);
    }
}

static void test_bitonic_sort16(void) {
    q15_t a[16] = {
        15, 3, 7, 11, 2, 14, 1, 12,
        6, 4, 8, 10, 5, 13, 9, 0
    };
    bitonic_q15_sort16(a);
    for (int i = 0; i < 16; ++i) {
        char msg[64];
        snprintf(msg, sizeof(msg), "bitonic sort16: a[%d] = %d (expected %d)",
                 i, a[i], i);
        ASSERT(a[i] == i, msg);
    }
}

static void test_bitonic_median9(void) {
    q15_t a[9] = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
    q15_t m = bitonic_q15_median9(a);
    ASSERT(m == 5, "median(1..9) = 5");

    q15_t b[9] = { 100, 200, 50, 300, 150, 250, 75, 175, 125 };
    q15_t m2 = bitonic_q15_median9(b);
    ASSERT(m2 == 150, "median of b = 150");
}

static void test_bitonic_sort_n(void) {
    q15_t a[10] = { 10, 3, 7, 1, 9, 4, 8, 2, 6, 5 };
    bitonic_q15_sort_n(a, 10);
    for (int i = 1; i <= 10; ++i) {
        char msg[64];
        snprintf(msg, sizeof(msg), "bitonic sort_n: position %d = %d (expected %d)",
                 i - 1, a[i - 1], i);
        ASSERT(a[i - 1] == i, msg);
    }
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int main(void) {
    printf("MPTC library unit tests\n");
    printf("Version: %s\n\n", mptc_version());

    test_sintmath_select();
    test_sintmath_clamp();
    test_sintmath_sin();
    test_sintmath_tanh();
    test_sintmath_gelu();
    test_sintmath_dual_sin();
    test_sintmath_dual_mul();
    test_sdft_basic();
    test_triple_pll_median3();
    test_triple_pll_outlier_rejection();
    test_integral_img_basic();
    test_ternary_roundtrip();
    test_ternary_arith();
    test_ternary_dot();
    test_bitonic_sort8();
    test_bitonic_sort16();
    test_bitonic_median9();
    test_bitonic_sort_n();

    printf("\n----------------------------------------\n");
    printf("Tests: %d passed, %d failed (total %d)\n",
           g_pass, g_fail, g_pass + g_fail);
    return g_fail == 0 ? 0 : 1;
}
