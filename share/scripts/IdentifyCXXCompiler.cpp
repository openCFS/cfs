#ifndef __cplusplus
# error "A C compiler has been selected for C++."
#endif

/* This file is inspired by the C compiler identification of CMake. */

/* This source file must have a .cpp extension so that all C++ compilers
   recognize the extension without flags.  Borland does not know .cxx for
   example.  */

#if defined(__COMO__)
# define COMPILER_ID "Comeau"

#elif defined(__INTEL_COMPILER) || defined(__ICC)
# define COMPILER_ID Intel
# define CXX_ID_TEXT ICC
# define CXX_VERSION_TEXT __INTEL_COMPILER __INTEL_COMPILER_BUILD_DATE
#if defined(__GNUC__)
#define CXX_GCC_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#endif
#elif defined(__OPEN64__)
# define COMPILER_ID Open64
# define CXX_ID_TEXT OPEN64
# define CXX_VERSION_TEXT __OPENCC__.__OPENCC_MINOR__.__OPENCC_PATCHLEVEL__
#if defined(__GNUC__)
#define CXX_GCC_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#endif

#elif defined(__BORLANDC__)
# define COMPILER_ID "Borland"

#elif defined(__WATCOMC__)
# define COMPILER_ID "Watcom"

#elif defined(__SUNPRO_CC)
# define COMPILER_ID "SunPro"

#elif defined(__HP_aCC)
# define COMPILER_ID "HP"

#elif defined(__DECCXX)
# define COMPILER_ID "Compaq"

#elif defined(__IBMCPP__)
# define COMPILER_ID "VisualAge"

#elif defined(__PGI)
# define COMPILER_ID "PGI"

#elif defined(__GNUC__)
# define COMPILER_ID GNU
# define CXX_ID_TEXT GCC
#define CXX_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#define CXX_GCC_VERSION_TEXT __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__

#elif defined(_MSC_VER)
# define COMPILER_ID "MSVC"

#elif defined(__ADSPBLACKFIN__) || defined(__ADSPTS__) || defined(__ADSP21000__)
/* Analog Devices C++ compiler for Blackfin, TigerSHARC and
   SHARC (21000) DSPs */
# define COMPILER_ID "ADSP"

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
CXX_ID CXX_ID_TEXT
CXX_VERSION CXX_VERSION_TEXT
CXX_GCC_VERSION CXX_GCC_VERSION_TEXT
