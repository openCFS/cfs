SET(PARDISO_FOUND 0)

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(LIBEXT "dylib")
ELSE(CFS_DISTRO STREQUAL "MACOSX")
  SET(LIBEXT "so")
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

#-------------------------------------------------------------------------------
# Look for Pardiso license file.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("Pardiso of Olaf Schenk"
  "${CFS_BINARY_DIR}/licenses/pardiso_license.pdf"
  "${CFS_DEPS_ROOT}/pardiso/build_pardiso.pl"
  "build_pardiso.log")

#-------------------------------------------------------------------------------
# Determine paths of PARDISO libraries.
#-------------------------------------------------------------------------------
SET(SCHENK_PARDISO_SERIAL_LIB
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libpardiso.${LIBEXT}")

SET(SCHENK_PARDISO_OPENMP_LIB
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libpardiso_p.${LIBEXT}")

# If a serial lib does not exist and we do not want to use OpenMP make
# sure that we can at least link against the parallel version of the lib.
IF(NOT EXISTS "${SCHENK_PARDISO_SERIAL_LIB}"
    AND NOT USE_OPENMP)
  IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
    SET(SCHENK_PARDISO_SERIAL_LIB "${SCHENK_PARDISO_OPENMP_LIB};-lgomp")
  ELSE(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
    SET(SCHENK_PARDISO_SERIAL_LIB "${SCHENK_PARDISO_OPENMP_LIB};-lguide;-lpthread")
  ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
ENDIF(NOT EXISTS "${SCHENK_PARDISO_SERIAL_LIB}"
  AND NOT USE_OPENMP)

IF(NOT EXISTS "${SCHENK_PARDISO_OPENMP_LIB}")
  SET(SCHENK_PARDISO_OPENMP_LIB "${SCHENK_PARDISO_SERIAL_LIB}")
ENDIF(NOT EXISTS "${SCHENK_PARDISO_OPENMP_LIB}")


SET(SCHENK_PARDISO_SERIAL_LIB
  "${SCHENK_PARDISO_SERIAL_LIB}"
  CACHE FILEPATH "Schenk Pardiso serial library.")

SET(SCHENK_PARDISO_OPENMP_LIB
  "${SCHENK_PARDISO_OPENMP_LIB}"
  CACHE FILEPATH "Schenk Pardiso OpenMP library.")

MARK_AS_ADVANCED(SCHENK_PARDISO_SERIAL_LIB)
MARK_AS_ADVANCED(SCHENK_PARDISO_OPENMP_LIB)

IF(USE_OPENMP)
  IF(SCHENK_PARDISO_OPENMP_LIB)
    SET(PARDISO_LIBRARY "${SCHENK_PARDISO_OPENMP_LIB}")
    SET(PARDISO_FOUND 1)
  ENDIF(SCHENK_PARDISO_OPENMP_LIB)
ELSE(USE_OPENMP)
  IF(SCHENK_PARDISO_SERIAL_LIB)
    SET(PARDISO_LIBRARY "${SCHENK_PARDISO_SERIAL_LIB}")
    #-------------------------------------------------------------------------------
    # Check if we need extra GFortran patch library.
    #-------------------------------------------------------------------------------
    IF(CFS_ARCH STREQUAL I386 AND
       CFS_CXX_COMPILER_NAME STREQUAL GCC AND
       NOT CFS_CXX_COMPILER_VER VERSION_LESS 4.3)
      SET(PARDISO_LIBRARY "${SCHENK_PARDISO_SERIAL_LIB};-lpardiso_gfortran_patch")
    ENDIF(CFS_ARCH STREQUAL I386 AND
          CFS_CXX_COMPILER_NAME STREQUAL GCC AND
          NOT CFS_CXX_COMPILER_VER VERSION_LESS 4.3)
    SET(PARDISO_FOUND 1)
  ENDIF(SCHENK_PARDISO_SERIAL_LIB)
ENDIF(USE_OPENMP)

