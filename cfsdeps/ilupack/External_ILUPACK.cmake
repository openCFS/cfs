# Ilupack is a closed source academic iterative solve by M. Bollhoefer, TU-Braunschweig
# by special agreement we have the source available which must not be redistributed

set(ilupack_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ilupack")
set(ilupack_source  "${ilupack_prefix}/src/ilupack")
set(ilupack_install  "${CMAKE_CURRENT_BINARY_DIR}")

# SET(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -DPRINT_MEM -DPRINT_INLINE")
# SET(CFSDEPS_ILUPACK_C_FLAGS "${CFSDEPS_C_FLAGS} -DPRINT_INLINE")
SET(CFSDEPS_ILUPACK_C_FLAGS "${CFSDEPS_C_FLAGS}")

STRING(REPLACE ";" "^" ILUPACK_AMD_LIBRARY "${AMD_LIBRARY}")
STRING(REPLACE ";" "^" ILUPACK_BLAS_LIBRARY "${BLAS_LIBRARY}")
STRING(REPLACE ";" "^" ILUPACK_LAPACK_LIBRARY "${LAPACK_LIBRARY}")


#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_INSTALL_PREFIX:PATH=${ilupack_install}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_ILUPACK_C_FLAGS}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DCMAKE_Fortran_COMPILER_ID:STRING=${CMAKE_Fortran_COMPILER_ID}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DAMD_LIB:STRING=${ILUPACK_AMD_LIBRARY}
  -DBLAS_LIB:FILEPATH=${ILUPACK_BLAS_LIBRARY}
  -DLAPACK_LIB:FILEPATH=${ILUPACK_LAPACK_LIBRARY}
  -DMETIS_LIB:FILEPATH=${METIS_LIBRARY}
)

IF(CMAKE_CROSSCOMPILING OR MSVC)
  LIST(APPEND CMAKE_ARGS
    -DBUILD_SAMPLES:BOOL=ON
    )
ELSE()
  LIST(APPEND CMAKE_ARGS
    -DBUILD_SAMPLES:BOOL=ON
    )
ENDIF()

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

IF(CFS_DISTRO STREQUAL "MACOSX")
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ilupack/ilupack-patch.cmake.in")
SET(PFN "${ilupack_prefix}/ilupack-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 


#-------------------------------------------------------------------------------
# The ilupack source is closed source! It must not be redistributed!!
# There is a copy at ${CFS_DS_WEBDAV}/ but it has certificate issuses, therefore there
# is also a mirror by the FAU CFS optimization group.
# To obfuscate the download,. the filename is replaced by its md5 sum.
#-------------------------------------------------------------------------------
SET(MIRRORS 
  "${CFS_FAU_MIRROR}/sources/ilupack/${ILUPACK_MD5}"
  "${CFS_DS_WEBDAV}/cfsdeps/sources/ilupack/MD5/${ILUPACK_MD5}")
  
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ilupack/${ILUPACK_GZ}")
SET(MD5_SUM ${ILUPACK_MD5})

SET(DLFN "${ilupack_prefix}/ilupack-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)


PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "ilupack" "${ILUPACK_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${ilupack_prefix}")

SET(ZIPFROMCACHE "${ilupack_prefix}/ilupackzipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${ilupack_prefix}/ilupack-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of ILUPACK libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(ILUPACK_SHARE "${LD}/libclock.a;${LD}/libsparsenew.a;${LD}/libvector.a")  
SET(ILUPACK_OPENMP "${LD}/libspdil.a;${LD}/libtools.a;${LD}/libbasic.a")  

SET(ILUPACK_LIBRARY
  "${ILUPACK_OPENMP};${LD}/libilupack.a;${LD}/libamd.a;${LD}/libcamd.a;${LD}/libblaslike.a;${LD}/libmetis.a;${LD}/libmetisomp.a;${LD}/libmtmetis.a;${LD}/libmumps.a;${LD}/libsparspak.a;${ILUPACK_SHARE};"               
  CACHE FILEPATH "ILUPACK library.")
  
MARK_AS_ADVANCED(ILUPACK_LIBRARY)


SET(ILUPACK_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
MARK_AS_ADVANCED(ILUPACK_INCLUDE_DIR)

