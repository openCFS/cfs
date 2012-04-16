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

#-------------------------------------------------------------------------------
# The CGAL external project
#-------------------------------------------------------------------------------
ExternalProject_Add(cgal
  DEPENDS boost zlib-static # gmp mpfr
  PREFIX "${cgal_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/cgal
  URL ${CGAL_URL}/${CGAL_GZ}
  URL_MD5 ${CGAL_MD5}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${cgal_install}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
    -DCGAL_INSTALL_BIN_DIR:PATH=bin/${CFS_ARCH_STR}
    -DCGAL_INSTALL_CMAKE_DIR:PATH=${LIB_SUFFIX}/${CFS_ARCH_STR}/CGAL
    -DCGAL_INSTALL_LIB_DIR:PATH=${LIB_SUFFIX}/${CFS_ARCH_STR}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DBoost_INCLUDE_DIR:PATH=${cgal_install}/include
    -DBoost_LIBRARY_DIRS:PATH=${cgal_install}/${LIB_SUFFIX}/$ARCH_STR
    -DBoost_THREAD_LIBRARY:FILEPATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libboost_thread-mt.a
    -DBoost_THREAD_LIBRARY_RELEASE:FILEPATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libboost_thread-mt.a
    -DBoost_THREAD_LIBRARY_DEBUG:FILEPATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libboost_thread-mt.a
    -DWITH_CGAL_Qt3:BOOL=OFF
    -DWITH_CGAL_Qt4:BOOL=OFF
    -DZLIB_INCLUDE_DIR:PATH=${cgal_install}/include
    -DZLIB_LIBRARY:PATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libz.a
    -DGMPXX_INCLUDE_DIR:PATH=${cgal_install}/include
    -DGMPXX_LIBRARIES:FILEPATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgmpxx.a
    -DGMP_INCLUDE_DIR:PATH=${cgal_install}/include
    -DGMP_LIBRARIES:FILEPATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgmp.a
    -DGMP_LIBRARIES_DIR:PATH=${cgal_install}/${LIB_SUFFIX}
    -DMPFR_INCLUDE_DIR:PATH=${cgal_install}/include
    -DMPFR_LIBRARIES:FILEPATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libmpfr.a
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgal/cgal-patch.cmake.in")
SET(PFN "${cgal_prefix}/cgal-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

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
