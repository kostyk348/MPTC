/**
 * @file main.cpp
 * @brief ESP32 demo for the MPTC library.
 *
 * Runs each module on a synthetic signal and prints timing / results to UART.
 * Connect at 115200 baud to view output.
 *
 * The demo is structured as a series of "showcases" — one per module —
 * so you can pick the ones relevant to your project and ignore the rest.
 */
#include <Arduino.h>
#include <math.h>
#include "mptc/mptc.h"

static void print_sep(void) {
    Serial.println("------------------------------------------------");
}

/* ------------------------------------------------------------------ */
/*  Showcase 1: sintmath — fast sin/tanh/GELU                         */
/* ------------------------------------------------------------------ */
static void showcase_sintmath(void) {
    Serial.println("\n[1] sintmath: branchless approximations");
    print_sep();

    /* Sweep sin over [-pi/2, pi/2]. */
    Serial.println("  x       sin(x)      tanh(x*4)   GELU(x*4)");
    for (int i = -8; i <= 8; ++i) {
        q15_t x = (q15_t)(i * 4096);  /* -32768..32768, mapped to [-pi/2, pi/2] */
        q15_t s = sint_q15_sin_wrap(x);
        q15_t t = sint_q15_tanh_rat(x);
        q15_t g = sint_q15_gelu(x);
        Serial.printf("  %+6d  %+6d  %+6d  %+6d\n", x, s, t, g);
    }

    /* Dual numbers: derivative of sin at x = pi/4. */
    sintmath_q15_dual_t x = sintmath_q15_dual_var(Q15_PI2 >> 1);
    sintmath_q15_dual_t y = sintmath_q15_dual_sin(x);
    Serial.printf("\n  dual sin(pi/4): val = %d, der = %d (cos(pi/4) = 23170)\n",
                  y.val, y.der);

    /* Benchmark: how many microseconds per sin call? */
    uint32_t t0 = micros();
    volatile q15_t sink = 0;
    for (int i = 0; i < 10000; ++i) {
        sink = sint_q15_sin_wrap((q15_t)(i & 0x7FFF));
    }
    uint32_t t1 = micros();
    Serial.printf("  10000 sint_q15_sin_wrap calls: %u us (%.2f us/call)\n",
                  t1 - t0, (float)(t1 - t0) / 10000.0f);
    (void)sink;
}

/* ------------------------------------------------------------------ */
/*  Showcase 2: sdft — Sliding DFT                                     */
/* ------------------------------------------------------------------ */
static void showcase_sdft(void) {
    Serial.println("\n[2] sdft: real-time spectral analysis");
    print_sep();

    /* 32-bin SDFT at 8 kHz sample rate -> 250 Hz/bin resolution. */
    sdft_q31_t s;
    if (sdft_q31_init(&s, 32, 8000) != 0) {
        Serial.println("  ERROR: sdft init failed");
        return;
    }

    /* Feed 64 samples of a 1 kHz sine (bin 4). */
    for (int i = 0; i < 64; ++i) {
        float phase = 2.0f * 3.14159265f * 1000.0f * i / 8000.0f;
        q15_t sample = (q15_t)(sinf(phase) * 30000.0f);
        sdft_q31_update(&s, sample);
    }

    /* Find peak bin. */
    uint16_t peak = sdft_q31_argmax(&s);
    float peak_hz = (float)s.sample_rate * (float)peak / (float)s.n_bins;
    Serial.printf("  32-bin SDFT at 8 kHz, fed 1 kHz sine:\n");
    Serial.printf("  Peak bin: %d (%.1f Hz)\n", peak, peak_hz);
    Serial.printf("  Expected: bin 4 (1000 Hz)\n");

    /* Print the spectrum bar chart. */
    Serial.println("  Spectrum (bins 1..15):");
    for (uint16_t k = 1; k <= 15; ++k) {
        q31_t mag = sdft_q31_mag(&s, k);
        int bars = (int)(mag >> 22);  /* scale to 0..32 range */
        if (bars > 32) bars = 32;
        Serial.printf("  bin %2d |", k);
        for (int b = 0; b < bars; ++b) Serial.print("#");
        Serial.println();
    }
    sdft_q31_free(&s);
}

/* ------------------------------------------------------------------ */
/*  Showcase 3: triple_pll — outlier rejection                         */
/* ------------------------------------------------------------------ */
static void showcase_triple_pll(void) {
    Serial.println("\n[3] triple_pll: outlier-rejecting filter");
    print_sep();

    triple_pll_q15_t pll;
    triple_pll_q15_init(&pll);

    /* Stable signal with one outlier. */
    q15_t samples[] = { 1000, 1000, 1000, 1000, 32000, 1000, 1000, 1000 };
    Serial.println("  Input  ->  Output (median of 3)");
    for (int i = 0; i < 8; ++i) {
        q15_t out = triple_pll_q15_update(&pll, samples[i]);
        Serial.printf("  %+6d  ->  %+6d %s\n",
                      samples[i], out,
                      (abs(samples[i] - out) > 1000) ? "(outlier rejected)" : "");
    }

    /* Triple-PLL IIR mode for noisy IMU-like signal. */
    Serial.println("\n  Triple-PLL IIR mode (alpha=0.1, noisy signal):");
    triple_pll_iir_q15_t iir;
    triple_pll_iir_q15_init(&iir, 0x0800, 0x0100);
    /* Simulate a slowly varying signal with random spikes. */
    int32_t base = 5000;
    for (int i = 0; i < 16; ++i) {
        /* Add a small drift + occasional spike. */
        base += (i % 5 == 0) ? 20000 : (i % 3 - 1) * 100;
        q15_t in = (q15_t)base;
        q15_t out = triple_pll_iir_q15_update(&iir, in);
        Serial.printf("  in: %+6d  out: %+6d\n", in, out);
    }
}

/* ------------------------------------------------------------------ */
/*  Showcase 4: integral_img — fast windowed statistics                */
/* ------------------------------------------------------------------ */
static void showcase_integral_img(void) {
    Serial.println("\n[4] integral_img: O(1) windowed statistics");
    print_sep();

    /* Build a 16x16 test image. */
    static uint8_t img[16 * 16];
    for (int y = 0; y < 16; ++y)
        for (int x = 0; x < 16; ++x)
            img[y * 16 + x] = (uint8_t)((x + y) * 8);

    integral_q15_t *ii = integral_q15_build(img, 16, 16);
    if (!ii) {
        Serial.println("  ERROR: integral_img build failed");
        return;
    }

    /* Query a few windows. */
    Serial.println("  Window (x,y,w,h)  ->  Sum  Mean  Var(scaled)");
    int queries[][4] = {
        {0, 0, 4, 4}, {4, 4, 8, 8}, {8, 8, 4, 4}, {0, 0, 16, 16}
    };
    for (int q = 0; q < 4; ++q) {
        int x = queries[q][0], y = queries[q][1];
        int w = queries[q][2], h = queries[q][3];
        uint32_t sum = integral_q15_window_sum(ii, x, y, w, h);
        uint32_t mean = integral_q15_window_mean(ii, x, y, w, h);
        uint32_t var = integral_q15_window_var_scaled(ii, x, y, w, h);
        Serial.printf("  (%d,%d,%d,%d)       -> %5u %5u %5u\n",
                      x, y, w, h, sum, mean, var);
    }

    /* Benchmark: time to binarize the 16x16 image. */
    static uint8_t out[16 * 16];
    uint32_t t0 = micros();
    integral_q15_sauvola_full(img, out, 16, 16, 2, 0x4000);
    uint32_t t1 = micros();
    Serial.printf("\n  Sauvola binarize 16x16: %u us\n", t1 - t0);

    integral_q15_free(ii);
}

/* ------------------------------------------------------------------ */
/*  Showcase 5: ternary — packed balanced ternary                      */
/* ------------------------------------------------------------------ */
static void showcase_ternary(void) {
    Serial.println("\n[5] ternary: packed balanced ternary arithmetic");
    print_sep();

    int test_values[] = { 0, 1, -1, 27, -27, 100, 500, -500 };
    Serial.println("  int -> trit8 -> int (roundtrip)");
    for (int i = 0; i < 8; ++i) {
        trit8_t t = ternary_from_int8(test_values[i]);
        int v = ternary_to_int8(t);
        Serial.printf("  %+5d -> 0x%04X -> %+5d %s\n",
                      test_values[i], t, v,
                      (v == test_values[i]) ? "OK" : "MISMATCH");
    }

    /* Dot product of 8 ternary weights with 8 inputs. */
    trit8_t w = 0;
    int8_t ws[8] = {1, -1, 0, 1, -1, 0, 1, -1};
    for (int i = 0; i < 8; ++i) w = ternary_set_trit8(w, i, ws[i]);
    q15_t inputs[8] = {100, 200, 300, 400, 500, 600, 700, 800};
    q15_t r = ternary_dot8(w, inputs);
    Serial.printf("\n  ternary_dot8([+1,-1,0,+1,-1,0,+1,-1], [100..800]) = %d\n", r);
    Serial.printf("  Expected: 100-200+0+400-500+0+700-800 = -300\n");
}

/* ------------------------------------------------------------------ */
/*  Showcase 6: bitonic — branchless sorting                           */
/* ------------------------------------------------------------------ */
static void showcase_bitonic(void) {
    Serial.println("\n[6] bitonic: branchless sorting networks");
    print_sep();

    /* Sort 16 values. */
    q15_t a[16] = {
        15, 3, 7, 11, 2, 14, 1, 12,
        6, 4, 8, 10, 5, 13, 9, 0
    };
    Serial.println("  Before: ");
    Serial.print("  ");
    for (int i = 0; i < 16; ++i) { Serial.printf("%d ", a[i]); }
    Serial.println();

    uint32_t t0 = micros();
    for (int i = 0; i < 1000; ++i) {
        /* Re-shuffle to prevent hoisting. */
        a[0] = (q15_t)(i & 0x7FFF);
        bitonic_q15_sort16(a);
    }
    uint32_t t1 = micros();
    Serial.println("  After: ");
    Serial.print("  ");
    for (int i = 0; i < 16; ++i) { Serial.printf("%d ", a[i]); }
    Serial.println();
    Serial.printf("\n  1000x bitonic_q15_sort16: %u us (%.2f us/call)\n",
                  t1 - t0, (float)(t1 - t0) / 1000.0f);

    /* Median of 9 (3x3 image filter). */
    q15_t win[9] = { 100, 200, 50, 300, 150, 250, 75, 175, 125 };
    q15_t m = bitonic_q15_median9(win);
    Serial.printf("\n  median9(100,200,50,300,150,250,75,175,125) = %d (expected 150)\n", m);
}

/* ------------------------------------------------------------------ */
/*  Setup & loop                                                       */
/* ------------------------------------------------------------------ */
void setup(void) {
    Serial.begin(115200);
    delay(500);

    Serial.println("\n========================================");
    Serial.println("MPTC library — ESP32 demo");
    Serial.printf("Library version: %s\n", mptc_version());
    Serial.printf("CPU freq: %u MHz\n", ESP.getCpuFreqMHz());
    Serial.println("========================================");

    showcase_sintmath();
    showcase_sdft();
    showcase_triple_pll();
    showcase_integral_img();
    showcase_ternary();
    showcase_bitonic();

    Serial.println("\n========================================");
    Serial.println("All showcases complete.");
    Serial.println("========================================");
}

void loop(void) {
    /* Nothing to do in the loop — demo runs once on boot. */
    delay(1000);
}
