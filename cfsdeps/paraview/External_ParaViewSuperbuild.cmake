#-------------------------------------------------------------------------------
# Project Homepages
# http://www.paraview.org
# http://paraview.org/Wiki/ParaView:Superbuild
#
# Building  ParaView with  the HDF5  reader  plugin for  CFS++ is  a two  step
# process.  First  the ParaView Superbuild  gets downloaded and  configured as
# the paraview-superbuild external project from the current CMake script. This
# involves obtaining  a Git clone of  the PV Superbuild repository,  either by
# using Git  itself or by downloading  a zipped archive of  the repository, in
# the case if the Git executable is not available.
#
# In the  patch step of the  paraview-superbuild a number of  files get copied
# from   the   current   source   directory   to  the   Git   repo   via   the
# paraview-superbuild-patch.cmake.in script.  These  files reflect the changes
# which are  necessary, so that  the PV Superbuild  blends in nicely  with the
# CFS++  build.  This  means, that  the Superbuild,  should reuse  the already
# downloaded  archives  in  our  CFS_DEPS_CACHE_DIR,  and  also  put  its  own
# downloads there.  Moreover,  an extra patch step gets added  to the paraview
# external  project  inside  the PV  Superbuild  (Projects/paraview.cmake  and
# Projects/paraview.patch.cmake.in).
# 
# The changes to the ParaView Superbuild reside in the following files:
#
#   paraview-superbuild-srcdir/
#   ├── CMake
#   │   └── ParaViewModules.cmake
#   ├── Projects
#   │   ├── hdf5.cmake
#   │   ├── paraview.cmake
#   │   ├── paraview.patch.cmake.in
#   │   └── zlib.cmake
#   └── versions.cmake
#
# The  additional  patch step  introduced  in  the paraview  external  project
# (Projects/paraview.cmake in PV Superbuild) makes  sure, that the contents of
# the cfsdeps/paraview/paraview-srcdir,  which reflect the changes  needed for
# the CFS++  features, get copied over  to the PV source  directory during the
# patch phase of ParaView.
#
# Our additions to the off-the-shelf ParaView include at the moment:
#   - The CFS++/NACS HDF5 reader
#   - A reader for Geomview .off polyhedral surface geometry files, which 
#     are e.g. produced by CGAL.
#   - Some additional color maps
#   - Fixes to the calculation of derivatives on second order quads and
#     triangles.
#   - A reduced list of VisItBridge readers, to make sure our own reader is
#     the only one, which can read .h5 files. Otherwise the user always has
#     to choose the correct reader, when opening a file
#   - Some small fixes to the VisItBridge CGNS reader for TRIA6 and PYRA13
#     element types.
#
# These additions reside in the following files:
#
#   paraview-srcdir/
#   ├── ParaViewCore
#   │   └── ServerManager
#   │       └── SMApplication
#   │           └── Resources
#   │               └── readers.xml
#   ├── Qt
#   │   └── Components
#   │       └── Resources
#   │           └── XML
#   │               └── ColorMaps.xml
#   ├── Utilities
#   │   └── VisItBridge
#   │       └── databases
#   │           ├── CGNS
#   │           │   └── avtCGNSFileFormat.C
#   │           ├── CMakeLists.txt
#   │           └── visit_readers.xml
#   └── VTK
#       ├── Common
#       │   └── DataModel
#       │       ├── vtkBiQuadraticQuad.cxx
#       │       ├── vtkQuadraticQuad.cxx
#       │       └── vtkQuadraticTriangle.cxx
#       ├── IO
#       │   ├── CFSReader
#       │   │   ├── CMakeLists.txt
#       │   │   ├── hdf5Common.cc
#       │   │   ├── hdf5Common.h
#       │   │   ├── hdf5Reader.cc
#       │   │   ├── hdf5Reader.h
#       │   │   ├── module.cmake
#       │   │   ├── vtkCFSReader.cxx
#       │   │   ├── vtkCFSReader.h
#       │   │   ├── vtkNACSReader_2013-03-18.cxx
#       │   │   ├── vtkNACSReader_2013-03-18.h
#       │   │   ├── vtkNACSReader_2013-08-12.cxx
#       │   │   └── vtkNACSReader_2013-08-12.h
#       │   └── OFFReader
#       │       ├── CMakeLists.txt
#       │       ├── module.cmake
#       │       ├── vtkOFFReader.cxx
#       │       └── vtkOFFReader.h
#       └── ThirdParty
#           └── hdf5
#               └── CMakeLists.txt
#
# The (super-)build of  ParaView has been tested on machines  running CentOS 6
# and Ubuntu 12.04 x86_64, on  which the bootstrap_devel_machine.sh script has
# been run.
#
# Tips for development:
# =====================
#
# As with any external project inside CFS++, the source and binary directories
# for PV Superbuild and ParaView reside under ${CFS_BINARY_DIR}/cfsdeps.
#
# Therefore, the ParaView source and binary directories end up under
#
# ${CFS_BINARY_DIR}/cfsdeps/paraview-superbuild/src/ \
#                  paraview-superbuild-build/paraview/src
#                   ├── paraview        ->  PARAVIEW_SOURCE_DIR
#                   ├── paraview-build  ->  PARAVIEW_BINARY_DIR
#                   └── paraview-stamp
#
#  Once   the   Superbuild   is   finished,    one   can   therefore   go   to
# PARAVIEW_BINARY_DIR,  set   CMAKE_BUILD_TYPE  to   Debug  and   do  ordinary
# development,  since  all  changes  due  to CFS++  are  already  included  in
# PARAVIEW_SOURCE_DIR. Once  development is finished,  the changes need  to be
# just copied back  to CFS_SOURCE_DIR/cfsdeps/paraview/paraview-srcdir and the
# corresponding patch scripts might need to be changed.
#
#
#
# http://www.paraview.org/Wiki/ParaView_Binaries
# https://github.com/Kitware/ParaView/tree/v3.12.0/SuperBuild
# http://www.paraview.org/Wiki/ParaView:Plugin_Deployment_with_Development_Installs
# http://www.kitware.com/products/html/BuildingExternalProjectsWithCMake2.8.html
# http://personal.cscs.ch/~jfavre/Projects/vtkLEA/vtklea.htm   see  also   the
# README_*.txt           files           in          this           directory.
# -------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to paraview sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(paraview_sb_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/paraview-superbuild")
set(paraview_sb_install  "${CMAKE_CURRENT_BINARY_DIR}/paraview")
set(paraview_sb_source  "${paraview_sb_prefix}/src/paraview-superbuild")

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
  -DENABLE_cgns:BOOL=ON
  -DENABLE_visitbridge:BOOL=ON
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

SET(PFN_TEMPL 
  "${CFS_SOURCE_DIR}/cfsdeps/paraview/paraview-superbuild-patch.cmake.in")
SET(PFN "${paraview_sb_prefix}/paraview-superbuild-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Find the Git executable.
#-------------------------------------------------------------------------------
Find_Package(Git)

#-------------------------------------------------------------------------------
# The paraview-superbuild external project
#-------------------------------------------------------------------------------
IF(GIT_FOUND)
  # Clone Git repo for ParaView Superbuild 4.1.
  ExternalProject_Add(paraview-superbuild
    DEPENDS ${CFS_PV_DEPENDENCIES}
    PREFIX ${paraview_sb_prefix}
    GIT_REPOSITORY http://paraview.org/ParaViewSuperbuild.git
    GIT_TAG "5867b1c34b73aee16cc9411007b2d70a4fad73a4"
    SOURCE_DIR ${paraview_source}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
    )
ELSE()
  # If we do not have the Git executable available, fall back
  # to downloading a zipped archive of the Git repo. 
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
ENDIF()


#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  paraview-superbuild
)


