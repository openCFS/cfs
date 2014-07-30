# The qtcurve external project for ParaView
set(qtcurve_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/qtcurve")
set(qtcurve_install "${CMAKE_CURRENT_BINARY_DIR}/qtcurve")
set(qtcurve_source "${qtcurve_prefix}/src/qtcurve")

#set(qtcurve_binary "${CMAKE_CURRENT_BINARY_DIR}/qtcurve-build")

ExternalProject_Add(qtcurve
  DEPENDS ${CFS_PV_DEPENDENCIES}
  PREFIX ${qtcurve_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/qt4
  URL http://pkgs.fedoraproject.org/repo/pkgs/qtcurve-kde4/QtCurve-KDE4-1.8.4.tar.bz2/3082f48597c9584b24ba144040f1ac34/QtCurve-KDE4-1.8.4.tar.bz2
  URL_MD5 3082f48597c9584b24ba144040f1ac34
  UPDATE_COMMAND ""
  SOURCE_DIR ${qtcurve_source}
#  BINARY_DIR ${qtcurve_binary}
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
