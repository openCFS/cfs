#-------------------------------------------------------------------------------
# Set paths to gidpost sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(gidpost_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/gidpost")
set(gidpost_source  "${gidpost_prefix}/src/gidpost")
set(gidpost_install  "${CFS_BINARY_DIR}")

#-------------------------------------------------------------------------------
# The GiDpost external project
#-------------------------------------------------------------------------------
ExternalProject_Add(gidpost
  PREFIX "${gidpost_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/gidpost
  URL ${GIDPOST_URL}/${GIDPOST_ZIP}
  URL_MD5 ${GIDPOST_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/gidpost/CMakeLists.txt ${gidpost_source}/CMakeLists.txt
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${gidpost_install}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DARCH_STR:STRING=${CFS_ARCH_STR}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  )

#-------------------------------------------------------------------------------
# These variables are used to find Gidpost by other projects
#-------------------------------------------------------------------------------
SET(GIDPOST_INCLUDE_DIR "${CFS_BINARY_DIR}/include" CACHE FILEPATH "GiDpost include dir" FORCE)

#-------------------------------------------------------------------------------
# Determine paths of GIDPOST libraries.
#-------------------------------------------------------------------------------
SET(GIDPOST_LIBRARY_DEBUG
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgidpost.a"
  CACHE FILEPATH "GiDpost library" FORCE)
SET(GIDPOST_LIBRARY_RELEASE
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgidpost.a"
  CACHE FILEPATH "GiDpost library" FORCE)

#-------------------------------------------------------------------------------
# Mark paths of GIDPOST libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(GIDPOST_LIBRARY_DEBUG)
MARK_AS_ADVANCED(GIDPOST_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set GIDPOST_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(GIDPOST_LIBRARY "${GIDPOST_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(GIDPOST_LIBRARY "${GIDPOST_LIBRARY_RELEASE}")
ENDIF(DEBUG)

SET(CFS_GIDPOST_VERSION "1.71")