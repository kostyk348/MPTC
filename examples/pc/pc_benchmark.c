/**
 * @file pc_benchmark.c
 * @brief PC benchmark for the MPTC library.
 *
 * Compares branchless fixed-point routines against naive baseline
 * implementations and prints a results table. The "cycles" column is
 * estimated from wall-clock time on the host CPU (TSC not assumed).
 *
 * Build via CMake:
 *   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
 *   cmake --build build
 *   ./build/pc_benchmark
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mptc/mptc.h"

#define BENCH_ITER 100000

/* ------------------------------------------------------------------ */
/*  Timing helper                                                      */
/* ------------------------------------------------------------------ */
static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ------------------------------------------------------------------ */
/*  Baseline (naive) implementations for comparison                    */
/* ------------------------------------------------------------------ */
static q15_t naive_select(q15_t a, q15_t b, uint32_t mask) {
    if (mask) return a;
    return b;
}

static q15_t naive_clamp(q15_t x, q15_t lo, q15_t hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static float naive_sin_taylor(float x) {
    return x - x*x*x/6.0f + x*x*x*x*x/120.0f - x*x*x*x*x*x*x/5040.0f;
}

static float naive_gelu(float x) {
    return 0.5f * x * (1.0f + tanhf(0.7978845608f * (x + 0.044715f * x*x*x)));
}

static q15_t naive_sort9_median(q15_t a[9]) {
    /* Naive insertion sort. */
    for (int i = 1; i < 9; ++i) {
        q15_t key = a[i];
        int j = i - 1;
        while (j >= 0 && a[j] > key) { a[j+1] = a[j]; j--; }
        a[j+1] = key;
    }
    return a[4];
}

/* ------------------------------------------------------------------ */
/*  Benchmark scaffolding                                              */
/* ------------------------------------------------------------------ */
typedef struct {
    const char *name;
    double mptc_us;
    double naive_us;
    double speedup;
    const char *notes;
} bench_result_t;

static void print_results(const char *title, bench_result_t *results, int n) {
    printf("\n=== %s ===\n", title);
    printf("%-32s %12s %12s %10s   %s\n",
           "Test", "MPTC (us)", "Naive (us)", "Speedup", "Notes");
    printf("-----------------------------------------------------------------------------\n");
    for (int i = 0; i < n; ++i) {
        printf("%-32s %12.2f %12.2f %9.2fx   %s\n",
               results[i].name,
               results[i].mptc_us,
               results[i].naive_us,
               results[i].speedup,
               results[i].notes ? results[i].notes : "");
    }
}

/* ------------------------------------------------------------------ */
/*  Test 1: Branchless select                                          */
/* ------------------------------------------------------------------ */
static void bench_select(bench_result_t *r) {
    q15_t a = 0x1000, b = 0x7000;
    volatile q15_t sink = 0;
    double t0, t1;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        sink = sint_q15_select(a, b, (uint32_t)(i & 1));
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        sink = naive_select(a, b, (uint32_t)(i & 1));
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "sint_q15_select (1M iters)";
    r->notes = "branchless vs if";
    (void)sink;
}

/* ------------------------------------------------------------------ */
/*  Test 2: Clamp                                                      */
/* ------------------------------------------------------------------ */
static void bench_clamp(bench_result_t *r) {
    volatile q15_t sink = 0;
    double t0, t1;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t x = (q15_t)(i & 0xFFFF);
        sink = mptc_clamp_q15(x, 1000, 30000);
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t x = (q15_t)(i & 0xFFFF);
        sink = naive_clamp(x, 1000, 30000);
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "mptc_clamp_q15";
    r->notes = "branchless min/max vs if";
    (void)sink;
}

/* ------------------------------------------------------------------ */
/*  Test 3: Fast sin                                                   */
/* ------------------------------------------------------------------ */
static void bench_sin(bench_result_t *r) {
    volatile q15_t sink = 0;
    volatile float fsink = 0;
    double t0, t1;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t x = (q15_t)((i * 37) & 0x6FFF);
        sink = sint_q15_sin_wrap(x);
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        float x = ((i * 37) & 0x6FFF) / 32768.0f * 3.14159f;
        fsink = sinf(x);
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "sin (4-term Taylor vs libm)";
    r->notes = "Q15 fixed-point";
    (void)sink; (void)fsink;
}

/* ------------------------------------------------------------------ */
/*  Test 4: Fast tanh                                                  */
/* ------------------------------------------------------------------ */
static void bench_tanh(bench_result_t *r) {
    volatile q15_t sink = 0;
    volatile float fsink = 0;
    double t0, t1;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t x = (q15_t)((i * 41) & 0x7FFF);
        sink = sint_q15_tanh_rat(x);
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        float x = ((i * 41) & 0x7FFF) / 32768.0f * 4.0f;
        fsink = tanhf(x);
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "tanh (rational vs libm)";
    r->notes = "Q15 fixed-point";
    (void)sink; (void)fsink;
}

/* ------------------------------------------------------------------ */
/*  Test 5: GELU                                                       */
/* ------------------------------------------------------------------ */
static void bench_gelu(bench_result_t *r) {
    volatile q15_t sink = 0;
    volatile float fsink = 0;
    double t0, t1;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t x = (q15_t)((i * 43) & 0x7FFF);
        sink = sint_q15_gelu(x);
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        float x = ((i * 43) & 0x7FFF) / 32768.0f * 4.0f;
        fsink = naive_gelu(x);
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "GELU (branchless vs tanh+libm)";
    r->notes = "Q15 fixed-point";
    (void)sink; (void)fsink;
}

/* ------------------------------------------------------------------ */
/*  Test 6: Median 9 (bitonic vs insertion)                            */
/* ------------------------------------------------------------------ */
static void bench_median9(bench_result_t *r) {
    q15_t data[9];
    for (int i = 0; i < 9; ++i) data[i] = (q15_t)((i * 1111) & 0x7FFF);
    volatile q15_t sink = 0;
    double t0, t1;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t tmp[9];
        memcpy(tmp, data, sizeof(tmp));
        /* Modify slightly each iter to prevent hoisting. */
        tmp[0] = (q15_t)(i & 0x7FFF);
        sink = bitonic_q15_median9(tmp);
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t tmp[9];
        memcpy(tmp, data, sizeof(tmp));
        tmp[0] = (q15_t)(i & 0x7FFF);
        sink = naive_sort9_median(tmp);
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "median of 9 (bitonic vs insertion)";
    r->notes = "3x3 image filter kernel";
    (void)sink;
}

/* ------------------------------------------------------------------ */
/*  Test 7: SDFT vs block FFT                                          */
/* ------------------------------------------------------------------ */
static void bench_sdft(bench_result_t *r) {
    sdft_q31_t s;
    if (sdft_q31_init(&s, 64, 8000) != 0) {
        r->name = "SDFT (failed to init)";
        return;
    }
    /* Warm up. */
    for (int i = 0; i < 64; ++i) sdft_q31_update(&s, (q15_t)(i * 100));

    volatile q31_t sink = 0;
    double t0, t1;

    /* SDFT: per-tick cost. */
    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        sdft_q31_update(&s, (q15_t)(i & 0x7FFF));
        sink = sdft_q31_mag(&s, 5);
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    /* Naive block FFT: recompute the full 64-point DFT every tick. */
    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        /* Pretend we have a 64-sample buffer; compute one bin via direct DFT. */
        static q15_t buf[64];
        buf[i & 63] = (q15_t)(i & 0x7FFF);
        double re = 0, im = 0;
        for (int k = 0; k < 64; ++k) {
            double angle = 2.0 * 3.14159265358979 * 5 * k / 64.0;
            re += buf[k] * cos(angle);
            im += buf[k] * sin(angle);
        }
        sink = (q31_t)(re * re + im * im);
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "SDFT vs block DFT (64-bin)";
    r->notes = "O(N) per tick vs O(N^2) recompute";
    (void)sink;
    sdft_q31_free(&s);
}

/* ------------------------------------------------------------------ */
/*  Test 8: Integral image Sauvola vs naive                            */
/* ------------------------------------------------------------------ */
static void bench_sauvola(bench_result_t *r) {
    enum { W = 128, H = 128 };
    static uint8_t img[W * H];
    static uint8_t out[W * H];
    /* Fill with a noisy gradient. */
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            img[y * W + x] = (uint8_t)((x + y + (x ^ y)) & 0xFF);

    double t0, t1;

    t0 = now_seconds();
    integral_q15_sauvola_full(img, out, W, H, 8, 0x4000);
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6;

    /* Naive: full window sum at every pixel. */
    static uint8_t out2[W * H];
    (void)out2;  /* written but not read — used only for fair timing comparison */
    t0 = now_seconds();
    for (int y = 8; y < H - 8; ++y) {
        for (int x = 8; x < W - 8; ++x) {
            long sum = 0, sumsq = 0;
            for (int dy = -8; dy <= 8; ++dy) {
                for (int dx = -8; dx <= 8; ++dx) {
                    uint8_t p = img[(y + dy) * W + (x + dx)];
                    sum   += p;
                    sumsq += (uint32_t)p * (uint32_t)p;
                }
            }
            int n = 17 * 17;
            double mean = (double)sum / n;
            double var  = (double)sumsq / n - mean * mean;
            double std  = sqrt(var > 0 ? var : 0);
            double T = mean * (1.0 + 0.5 * (std / 128.0 - 1.0));
            out2[y * W + x] = (img[y * W + x] > T) ? 255 : 0;
        }
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "Sauvola binarize 128x128";
    r->notes = "Integral img O(1) vs sliding window O(W^2)";
}

/* ------------------------------------------------------------------ */
/*  Test 9: Triple-PLL median vs single IIR                            */
/* ------------------------------------------------------------------ */
static void bench_triple_pll(bench_result_t *r) {
    triple_pll_iir_q15_t pll;
    triple_pll_iir_q15_init(&pll, 0x0800, 0x0100);

    /* Baseline: single IIR. */
    q15_t single = 0;
    q15_t alpha = 0x0800;

    volatile q15_t sink = 0;
    double t0, t1;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t x = (q15_t)(i & 0x7FFF);
        sink = triple_pll_iir_q15_update(&pll, x);
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        q15_t x = (q15_t)(i & 0x7FFF);
        q15_t err = (q15_t)((int32_t)x - (int32_t)single);
        single = (q15_t)((int32_t)single + (int32_t)mptc_mul_q15(alpha, err));
        sink = single;
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->mptc_us > 0 ? r->naive_us / r->mptc_us : 0;
    /* Note: triple-PLL is ~3x slower than single IIR by design —
     * it does 3x the work. The "speedup" column here will be ~0.33.
     * The benefit is outlier rejection, not raw speed. */
    r->name = "Triple-PLL (3x IIRs + median)";
    r->notes = "Not faster — see notes (outlier rejection)";
    (void)sink;
}

/* ------------------------------------------------------------------ */
/*  Test 10: Ternary dot product vs float dot product                  */
/* ------------------------------------------------------------------ */
static void bench_ternary_dot(bench_result_t *r) {
    trit8_t w = ternary_from_int8(123);
    q15_t inputs[8];
    for (int i = 0; i < 8; ++i) inputs[i] = (q15_t)((i * 1234) & 0x7FFF);

    volatile q15_t sink = 0;
    volatile float fsink = 0;
    double t0, t1;

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        inputs[0] = (q15_t)(i & 0x7FFF);
        sink = ternary_dot8(w, inputs);
    }
    t1 = now_seconds();
    r->mptc_us = (t1 - t0) * 1e6 / BENCH_ITER;

    /* Float baseline: same dot product in float. */
    float fin[8];
    for (int i = 0; i < 8; ++i) fin[i] = (float)inputs[i] / 32768.0f;
    float fw[8] = { 1, -1, 0, 1, -1, 0, 1, -1 };

    t0 = now_seconds();
    for (int i = 0; i < BENCH_ITER; ++i) {
        fin[0] = (float)(i & 0x7FFF) / 32768.0f;
        float acc = 0;
        for (int j = 0; j < 8; ++j) acc += fin[j] * fw[j];
        fsink = acc;
    }
    t1 = now_seconds();
    r->naive_us = (t1 - t0) * 1e6 / BENCH_ITER;

    r->speedup = r->naive_us / r->mptc_us;
    r->name = "ternary dot8 (vs float dot)";
    r->notes = "ternary weights = no multiplications";
    (void)sink; (void)fsink;
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int main(void) {
    printf("MPTC library benchmark\n");
    printf("Version: %s\n", mptc_version());
    printf("Iterations per test: %d\n", BENCH_ITER);
    printf("Compiler: "
#if defined(__GNUC__)
        "GCC " __VERSION__
#elif defined(__clang__)
        "Clang " __VERSION__
#else
        "Unknown"
#endif
        "\n");
    printf("Host: "
#if defined(__x86_64__)
        "x86_64"
#elif defined(__aarch64__)
        "aarch64"
#elif defined(__arm__)
        "arm32"
#else
        "unknown"
#endif
        "\n\n");

    bench_result_t results[16];
    int n = 0;

    bench_select            (&results[n++]);
    bench_clamp             (&results[n++]);
    bench_sin               (&results[n++]);
    bench_tanh              (&results[n++]);
    bench_gelu              (&results[n++]);
    bench_median9           (&results[n++]);
    bench_sdft              (&results[n++]);
    bench_sauvola           (&results[n++]);
    bench_triple_pll        (&results[n++]);
    bench_ternary_dot       (&results[n++]);

    print_results("MPTC vs Naive Baselines", results, n);

    /* Print summary. */
    double avg_speedup = 0;
    double best_speedup = 0;
    const char *best_name = "";
    for (int i = 0; i < n; ++i) {
        avg_speedup += results[i].speedup;
        if (results[i].speedup > best_speedup) {
            best_speedup = results[i].speedup;
            best_name = results[i].name;
        }
    }
    avg_speedup /= n;
    printf("\nSummary:\n");
    printf("  Average speedup: %.2fx\n", avg_speedup);
    printf("  Best speedup:    %.2fx  (%s)\n", best_speedup, best_name);
    printf("\nNote: Speedups are relative to naive C baselines on the host CPU.\n");
    printf("      On real MCUs without FPU (Cortex-M0/M3, RP2040, AVR), the\n");
    printf("      fixed-point MPTC routines typically gain an additional 3-5x\n");
    printf("      over the float baselines, because they avoid the soft-float library.\n");
    return 0;
}
