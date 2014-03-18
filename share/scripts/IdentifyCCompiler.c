#ifdef __cplusplus
# error "A C++ compiler has been selected for C."
#endif

/* This file is inspired by the C compiler identification of CMake. */

#if defined(__INTEL_COMPILER) || defined(__ICC)
# define COMPILER_ID Intel
# define CC_ID_TEXT ICC
# define CC_VERSION_TEXT __INTEL_COMPILER __INTEL_COMPILER_BUILD_DATE
#if defined(__GNUC__)
#define CC_GCC_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#endif
#elif defined(__OPEN64__)
# define COMPILER_ID Open64
# define CC_ID_TEXT OPEN64
# define CC_VERSION_TEXT __OPENCC__.__OPENCC_MINOR__.__OPENCC_PATCHLEVEL__
#if defined(__GNUC__)
#define CC_GCC_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#endif
#elif defined(__clang__)
# define COMPILER_ID "Clang"
# define CXX_ID_TEXT CLANG
# define CXX_VERSION_TEXT __clang_major__.__clang_minor__.__clang_patchlevel__
#if defined(__GNUC__)
#define CXX_GCC_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#endif

#elif defined(__BORLANDC__)
# define COMPILER_ID "Borland"

#elif defined(__WATCOMC__)
# define COMPILER_ID "Watcom"

#elif defined(__SUNPRO_C)
# define COMPILER_ID "SunPro"

#elif defined(__HP_cc)
# define COMPILER_ID "HP"

#elif defined(__DECC)
# define COMPILER_ID "Compaq"

#elif defined(__IBMC__)
# define COMPILER_ID "VisualAge"

#elif defined(__PGI)
# define COMPILER_ID "PGI"

#elif defined(__GNUC__)
# define COMPILER_ID GNU
# define CC_ID_TEXT GCC
#define CC_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#define CC_GCC_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__

#elif defined(_MSC_VER)
# define COMPILER_ID "MSVC"
# define CC_ID_TEXT MSVC
# define CC_VERSION_TEXT _MSC_VER
#elif defined(__ADSPBLACKFIN__) || defined(__ADSPTS__) || defined(__ADSP21000__)
/* Analog Devices C++ compiler for Blackfin, TigerSHARC and
   SHARC (21000) DSPs */
# define COMPILER_ID "ADSP"

/* IAR Systems compiler for embedded systems.
   http://www.iar.com
   Not supported yet by CMake
#elif defined(__IAR_SYSTEMS_ICC__)
# define COMPILER_ID "IAR" */

/* sdcc, the small devices C compiler for embedded systems,
   http://sdcc.sourceforge.net  */
#elif defined(SDCC)
# define COMPILER_ID "SDCC"

#elif defined(_COMPILER_VERSION)
# define COMPILER_ID "MIPSpro"

/* This compiler is either not known or is too old to define an
   identification macro.  Try to identify the platform and guess that
   it is the native compiler.  */
#elif defined(__sgi)
# define COMPILER_ID "MIPSpro"

#elif defined(__hpux) || defined(__hpua)
# define COMPILER_ID "HP"

#else /* unknown compiler */
# define COMPILER_ID ""

#endif

CMAKE_COMPILER_ID COMPILER_ID
CC_ID CC_ID_TEXT
CC_VERSION CC_VERSION_TEXT
CC_GCC_VERSION CC_GCC_VERSION_TEXT
