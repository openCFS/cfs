# The qtcurve external project for ParaView
set(qtcurve_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/qtcurve")
set(qtcurve_install "${CMAKE_CURRENT_BINARY_DIR}/qtcurve")
set(qtcurve_source "${qtcurve_prefix}/src/qtcurve")

#set(qtcurve_binary "${CMAKE_CURRENT_BINARY_DIR}/qtcurve-build")

SET(QTCURVE_BASE "QtCurve")
SET(QTCURVE_VER "KDE4-1.8.4")
SET(QTCURVE_MD5 "3082f48597c9584b24ba144040f1ac34")
IF(WIN32)
  SET(PRECOMPILED_PCKG_NAME "qtcurve_${QTCURVE_VER}_${CFS_ARCH_STR}_${TOOLSET_ID}.zip")
ELSE(WIN32)
  SET(PRECOMPILED_PCKG_NAME "qtcurve_${QTCURVE_VER}_${CFS_ARCH_STR}_${FC_ID}.zip")
ENDIF(WIN32)
SET(PRECOMPILED_PCKG_FILE "${CFS_DEPS_CACHE_DIR}/precompiled/CFSDEPS/${PRECOMPILED_PCKG_NAME}")
  
SET(PREFIX_DIR "${qtcurve_prefix}")

SET(ZIPFROMCACHE "${qtcurve_prefix}/qtcurve-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${qtcurve_prefix}/qtcurve-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The qtcurve external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(qtcurve
    PREFIX "${qtcurve_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(qtcurve
    DEPENDS ${CFS_PV_DEPENDENCIES}
    PREFIX ${qtcurve_prefix}
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/qt4
    URL http://pkgs.fedoraproject.org/repo/pkgs/qtcurve-kde4/${QTCURVE_BASE}-${QTCURVE_VER}.tar.bz2/${QTCURVE_MD5}/${QTCURVE_BASE}-${QTCURVE_VER}.tar.bz2
    URL_MD5 3082f48597c9584b24ba144040f1ac34
    UPDATE_COMMAND ""
    SOURCE_DIR ${qtcurve_source}
#    BINARY_DIR ${qtcurve_binary}
    INSTALL_DIR ${qtcurve_install}
    CMAKE_ARGS
      -DBUILD_SHARED_LIBS:BOOL=ON
      -DQTC_QT_ONLY:BOOL=ON
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      -DBUILD_TESTING:BOOL=OFF
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
      -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    INSTALL_COMMAND ""
  )
  
  IF("${CFS_DEPS_TOCACHE}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(qtcurve cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
