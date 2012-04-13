# The GiDpost external project

set(gidpost_source  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/GiDpost")
set(gidpost_binary  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/GiDpost-build")
set(gidpost_install  "${CMAKE_CURRENT_BINARY_DIR}")

ExternalProject_Add(GiDpost
  DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
  URL ${GIDPOST_URL}/${GIDPOST_ZIP}
  URL_MD5 ${GIDPOST_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/gidpost/CMakeLists.txt ${gidpost_source}/CMakeLists.txt
  SOURCE_DIR ${gidpost_source}
  BINARY_DIR ${gidpost_binary}
  CMAKE_ARGS
#     ${gidpost_lib_args}
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
#    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG:PATH=${gidpost_binary}/lib
#    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=${gidpost_binary}/lib
#    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE:PATH=${gidpost_binary}/lib
#    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=${gidpost_binary}/lib
    -DCMAKE_INSTALL_PREFIX:PATH=${gidpost_install}
    -DARCH_STR:STRING=${CFS_ARCH_STR}
#    -DARCH:STRING=${CFS_ARCH}
#    -DREV:STRING=${CFS_REV}
#    -DDIST:STRING=${CFS_DIST}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
    ${Gidpost_EXTRA_ARGS}
  )

#  ExternalProject_Add_Step(Gidpost InstallUUIDInclude
#    COMMAND ${CMAKE_COMMAND} -E copy_directory ${gidpost_source}/src/uuid/include/gidpost/uuid ${gidpost_install}/include/gidpost/uuid
#    DEPENDEES install
#  )

# These variables are used to find Gidpost by other projects
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

# Safety check for existance of libs and headers
IF(EXISTS "${GIDPOST_LIBRARY}"
    AND EXISTS "${GIDPOST_INCLUDE_DIR}/gidpost.h")
  
  SET(GIDPOST_FOUND 1)
  
ELSE(EXISTS "${GIDPOST_LIBRARY}"
    AND EXISTS "${GIDPOST_INCLUDE_DIR}/gidpost.h")
  
  MESSAGE(FATAL_ERROR "Could not find ${GIDPOST_LIBRARY} or ${GIDPOST_INCLUDE_DIR}/gidpost.h")
  
ENDIF(EXISTS "${GIDPOST_LIBRARY}"
  AND EXISTS "${GIDPOST_INCLUDE_DIR}/gidpost.h")
