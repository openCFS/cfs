#-------------------------------------------------------------------------------
# Set paths to cgal sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(cgal_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cgal")
set(cgal_source  "${cgal_prefix}/src/cgal")
set(cgal_install  "${CMAKE_CURRENT_BINARY_DIR}")

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
IF(CMAKE_BUILD_TYPE STREQUAL "debug")
  SET(CMAKE_BUILD_TYPE "Debug")
ELSE(CMAKE_BUILD_TYPE STREQUAL "debug")
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF(CMAKE_BUILD_TYPE STREQUAL "debug")

# Make sure that static Boost libs are used under Windows.
SET(CGAL_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -DBOOST_THREAD_USE_LIB")

IF(CFS_CXX_COMPILER_NAME STREQUAL "OPEN64")
  SET(CGAL_CXX_FLAGS "${CGAL_CXX_FLAGS} -mieee-fp -fp-accuracy=strict -DCGAL_DISABLE_ROUNDING_MATH_CHECK")
ENDIF()

SET(CMAKE_ARGS
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_INSTALL_PREFIX:PATH=${cgal_install}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CGAL_CXX_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
)

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

IF(MINGW)
  IF(CFS_ARCH STREQUAL "X86_64")
    LIST(APPEND CMAKE_ARGS
      -C ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/cgal/TryRunResults_GCC45_CENTOS6_WIN64.cmake
      )
  ENDIF()
ENDIF(MINGW)

#-------------------------------------------------------------------------------
# The CGAL external project
#-------------------------------------------------------------------------------
ExternalProject_Add(cgal
  DEPENDS boost zlib gmp mpfr
  PREFIX "${cgal_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/cgal
  URL ${CGAL_URL}/${CGAL_GZ}
  URL_MD5 ${CGAL_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCGAL_INSTALL_BIN_DIR:PATH=bin/${CFS_ARCH_STR}
    -DCGAL_INSTALL_CMAKE_DIR:PATH=${LIB_SUFFIX}/${CFS_ARCH_STR}/CGAL
    -DCGAL_INSTALL_LIB_DIR:PATH=${LIB_SUFFIX}/${CFS_ARCH_STR}
    -DBUILD_SHARED_LIBS:BOOL=OFF
#    -DBoost_INCLUDE_DIR:PATH=${cgal_install}/include
#    -DBoost_LIBRARY_DIRS:PATH=${cgal_install}/${LIB_SUFFIX}/$ARCH_STR
#    -DBoost_THREAD_LIBRARY:FILEPATH=${BOOST_THREAD_LIB}
#    -DBoost_THREAD_LIBRARY_RELEASE:FILEPATH=${BOOST_THREAD_LIB_RELEASE}
#    -DBoost_THREAD_LIBRARY_DEBUG:FILEPATH=${BOOST_THREAD_LIB_DEBUG}
    -DWITH_CGAL_Qt3:BOOL=OFF
    -DWITH_CGAL_Qt4:BOOL=OFF
    -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
    -DZLIB_LIBRARY:PATH=${ZLIB_LIBRARY}
    -DGMPXX_INCLUDE_DIR:PATH=${GMP_INCLUDE_DIR}
    -DGMPXX_LIBRARIES:FILEPATH=${GMPXX_LIBRARY}
    -DGMP_INCLUDE_DIR:PATH=${GMP_INCLUDE_DIR}
    -DGMP_LIBRARIES:FILEPATH=${GMP_LIBRARY}
    -DGMP_LIBRARIES_DIR:PATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    -DMPFR_INCLUDE_DIR:PATH=${MPFR_INCLUDE_DIR}
    -DMPFR_LIBRARIES:FILEPATH=${MPFR_LIBRARY}
)


#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgal/cgal-patch.cmake.in")
SET(PFN "${cgal_prefix}/cgal-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

SET(BOOST_SETUP_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgal/CGAL_SetupBoost.cmake.in")
SET(BOOST_SETUP "${cgal_prefix}/CGAL_SetupBoost.cmake")
CONFIGURE_FILE("${BOOST_SETUP_TEMPL}" "${BOOST_SETUP}" @ONLY) 

#-------------------------------------------------------------------------------
# We do not use the PATCH_COMMAND  of ExternalProject_Add since we do not only
# want to apply the patch script  during configuration time but also if it has
# changed.  Therefore,   we  need  a   dependency  on  the   configured  patch
# script. This can be achieved by  adding an additional build step between the
# download and configure steps.
#
# NOTE: The  patch script should  be designed  in such a  way, that it  can be
# applied to  an already patched  source tree. This  is due to the  fact, that
# ExternalProject_Add only extracts the source if the MD5 sum has has changed.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(cgal custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${cgal_source}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  cgal
)

SET(CGAL_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of CGAL libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(CGAL_LIBRARY
  "${LD}/libCGAL.a"
  CACHE FILEPATH "CGAL library.")

MARK_AS_ADVANCED(CGAL_INCLUDE_DIR)
MARK_AS_ADVANCED(CGAL_LIBRARY)
