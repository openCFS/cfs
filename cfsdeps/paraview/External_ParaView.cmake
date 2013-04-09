#-------------------------------------------------------------------------------
# Set prefix path and path to paraview sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(paraview_sb_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/paraview_superbuild")
set(paraview_sb_install  "${CMAKE_CURRENT_BINARY_DIR}/paraview")
set(paraview_sb_source  "${paraview_prefix}/src/paraview_superbuild")

#-------------------------------------------------------------------------------
# Info on SuperBuild
# http://paraview.org/Wiki/ParaView:Superbuild
#
# Info on regular ParaView build
# http://paraview.org/Wiki/ParaView:Build_And_Install
#
# Info on cross compiling
# http://www.paraview.org/Wiki/Cross_compiling_ParaView3_and_VTK#Building_ParaView_for_a_new_target_platform
# http://www.vtk.org/Wiki/CMake_Cross_Compiling#How_to_cross_compile_ParaView.2C_Python
# http://www.vtk.org/Wiki/BuildingPythonWithCMake
#
# Info on extending ParaView
# http://paraview.org/Wiki/ExtendingParaView
# http://paraview.org/Wiki/Plugin_HowTo
# http://paraview.org/Wiki/Extending_ParaView_at_Compile_Time
#
# Info on adding plugins
# http://www.paraview.org/Wiki/ParaView/Plugin_HowTo
# http://www.kitware.com/media/html/WritingAParaViewReaderPlug-in.html
# http://openfoamwiki.net/index.php/Tip_Build_A_Paraview3_Plugin
#-------------------------------------------------------------------------------

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${paraview_sb_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCFS_DEPS_CACHE_DIR:PATH=${CFS_DEPS_CACHE_DIR}
  -DCFS_SOURCE_DIR:PATH=${CFS_SOURCE_DIR}
  -DENABLE_boost:BOOL=ON
  -DENABLE_ffmpeg:BOOL=ON
  -DENABLE_fontconfig:BOOL=ON
  -DENABLE_freetype:BOOL=ON
  -DENABLE_hdf5:BOOL=ON
  -DENABLE_libxml2:BOOL=ON
  -DENABLE_paraview:BOOL=ON
  -DENABLE_png:BOOL=ON
  -DENABLE_szip:BOOL=ON
  -DENABLE_zlib:BOOL=ON
  -DENABLE_qt:BOOL=ON
  -DUSE_SYSTEM_qt=ON
  -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
  -DENABLE_python:BOOL=ON
  -DUSE_SYSTEM_python=ON
  )

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_CXX_FLAGS:FILEPATH=${CFLAGS}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.6
    )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

#-------------------------------------------------------------------------------
# The paraview-download external project
#-------------------------------------------------------------------------------
ExternalProject_Add(paraview-superbuild
  DEPENDS ${CFS_PV_DEPENDENCIES}
  PREFIX ${paraview_sb_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview
  SOURCE_DIR ${paraview_source}
  URL ${PARAVIEW_SB_URL}/${PARAVIEW_SB_GZ}
  URL_MD5 ${PARAVIEW_SB_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
  )

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  paraview-superbuild
)


