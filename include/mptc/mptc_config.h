/**
 * @file mptc_config.h
 * @brief Compile-time configuration for the MPTC library.
 *
 * Override any of these defaults by defining them on the compiler command line
 * before including any mptc/ header. For example:
 *
 *     cc -DMPTC_USE_FLOAT=1 -DMPTC_ASSERT=my_assert ...
 *
 * All MPTC headers include this file first; user code does not need to include
 * it explicitly.
 */
#ifndef MPTC_CONFIG_H
#define MPTC_CONFIG_H

/* ------------------------------------------------------------------ */
/*  Numeric back-ends                                                  */
/* ------------------------------------------------------------------ */
/**
 * @def MPTC_USE_FLOAT
 * @brief If 1, expose float32 overloads alongside Q15/Q31 fixed-point.
 *
 * Default is 0 (fixed-point only). Enable on Cortex-M4F/M7/ESP32 where the
 * hardware FPU makes float competitive with Q31 on energy.
 */
#ifndef MPTC_USE_FLOAT
#define MPTC_USE_FLOAT 0
#endif

/* ------------------------------------------------------------------ */
/*  Strictness                                                        */
/* ------------------------------------------------------------------ */
/**
 * @def MPTC_BRANCHLESS_STRICT
 * @brief If 1, hot-path routines are guaranteed branch-free.
 *
 * Disable only for debugging — when off, some routines may use early-return
 * for clarity. Default is 1 (production).
 */
#ifndef MPTC_BRANCHLESS_STRICT
#define MPTC_BRANCHLESS_STRICT 1
#endif

/* ------------------------------------------------------------------ */
/*  Implementation choices                                            */
/* ------------------------------------------------------------------ */
/**
 * @def MPTC_TABLE_SIN
 * @brief If 1, sint_q15_sin uses a 256-entry LUT instead of Taylor series.
 *
 * Trades 512 bytes of Flash for ~30 % speedup on cores with slow multipliers
 * (AVR, RP2040). Default is 0 (Taylor — branchless, no Flash cost).
 */
#ifndef MPTC_TABLE_SIN
#define MPTC_TABLE_SIN 0
#endif

/**
 * @def MPTC_TABLE_COS
 * @brief Same as MPTC_TABLE_SIN but for cos.
 */
#ifndef MPTC_TABLE_COS
#define MPTC_TABLE_COS 0
#endif

/* ------------------------------------------------------------------ */
/*  Endianness (auto-detect, override if needed)                      */
/* ------------------------------------------------------------------ */
#if !defined(MPTC_BIG_ENDIAN) && !defined(MPTC_LITTLE_ENDIAN)
#  if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#    define MPTC_BIG_ENDIAN 1
#    define MPTC_LITTLE_ENDIAN 0
#  else
#    define MPTC_BIG_ENDIAN 0
#    define MPTC_LITTLE_ENDIAN 1
#  endif
#endif

/* ------------------------------------------------------------------ */
/*  Assertions                                                         */
/* ------------------------------------------------------------------ */
/**
 * @def MPTC_ASSERT
 * @brief Override for embedded environments where <assert.h> is undesirable.
 *
 * Default is the standard assert(). On bare-metal targets, you may want:
 *     #define MPTC_ASSERT(x) do { if(!(x)) while(1){} } while(0)
 */
#ifndef MPTC_ASSERT
#  include <assert.h>
#  define MPTC_ASSERT(x) assert(x)
#endif

/* ------------------------------------------------------------------ */
/*  Inline keyword                                                     */
/* ------------------------------------------------------------------ */
/**
 * @def MPTC_INLINE
 * @brief Portable `static inline` with optional forced inlining.
 */
#ifndef MPTC_INLINE
#  if defined(_MSC_VER)
#    define MPTC_INLINE static __forceinline
#  elif defined(__GNUC__) || defined(__clang__)
#    define MPTC_INLINE static inline __attribute__((always_inline))
#  else
#    define MPTC_INLINE static inline
#  endif
#endif

/* ------------------------------------------------------------------ */
/*  Restricted pointer                                                 */
/* ------------------------------------------------------------------ */
/**
 * @def MPTC_RESTRICT
 * @brief Portable `restrict` qualifier for aliasing hints.
 */
#ifndef MPTC_RESTRICT
#  if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#    define MPTC_RESTRICT restrict
#  elif defined(__GNUC__) || defined(__clang__)
#    define MPTC_RESTRICT __restrict
#  else
#    define MPTC_RESTRICT
#  endif
#endif

/* ------------------------------------------------------------------ */
/*  Branchless hint macros (for readability, not enforcement)         */
/* ------------------------------------------------------------------ */
/**
 * @def MPTC_LIKELY
 * @def MPTC_UNLIKELY
 * @brief Branch prediction hints. Used only in non-strict paths.
 */
#if defined(__GNUC__) || defined(__clang__)
#  define MPTC_LIKELY(x)   __builtin_expect(!!(x), 1)
#  define MPTC_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#  define MPTC_LIKELY(x)   (x)
#  define MPTC_UNLIKELY(x) (x)
#endif

#endif /* MPTC_CONFIG_H */
