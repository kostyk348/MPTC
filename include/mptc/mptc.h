/**
 * @file mptc.h
 * @brief MPTC library umbrella header. Includes everything.
 *
 * Include this single header to get the full library API. If you want to
 * minimize code size on a constrained target, include only the specific
 * module headers you need (e.g., `#include <mptc/sintmath.h>`).
 */
#ifndef MPTC_H
#define MPTC_H

#include "mptc_config.h"
#include "mptc_types.h"
#include "mptc/sintmath.h"
#include "mptc/sdft.h"
#include "mptc/triple_pll.h"
#include "mptc/integral_img.h"
#include "mptc/ternary.h"
#include "mptc/bitonic.h"

/** @brief Library version string. */
#define MPTC_VERSION "1.0.0"

/** @brief Library version as a 32-bit integer: 0xMMmmpp. */
#define MPTC_VERSION_HEX 0x010000

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Get the library version string.
 * @return Static C string, e.g. "1.0.0".
 */
const char *mptc_version(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MPTC_H */
