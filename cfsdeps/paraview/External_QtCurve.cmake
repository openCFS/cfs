# The qtcurve external project for ParaView
set(qtcurve_source "${CMAKE_CURRENT_BINARY_DIR}/qtcurve")
set(qtcurve_install "${CMAKE_CURRENT_BINARY_DIR}")

set(qtcurve_binary "${CMAKE_CURRENT_BINARY_DIR}/qtcurve-build")

ExternalProject_Add(qtcurve
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/qt4
  URL http://pkgs.fedoraproject.org/repo/pkgs/qtcurve-kde4/QtCurve-KDE4-1.8.4.tar.bz2/3082f48597c9584b24ba144040f1ac34/QtCurve-KDE4-1.8.4.tar.bz2
  URL_MD5 3082f48597c9584b24ba144040f1ac34
  UPDATE_COMMAND ""
  SOURCE_DIR ${qtcurve_source}
  BINARY_DIR ${qtcurve_binary}
  INSTALL_DIR ${qtcurve_install}
  CMAKE_CACHE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DQTC_QT_ONLY:BOOL=ON
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    -DBUILD_TESTING:BOOL=OFF
    -DCMAKE_CXX_FLAGS:STRING=${pv_tpl_cxx_flags}
    -DCMAKE_C_FLAGS:STRING=${pv_tpl_c_flags}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_CFG_INTDIR}
    ${pv_tpl_compiler_args}
  INSTALL_COMMAND ""
  )
