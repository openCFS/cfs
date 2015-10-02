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
#   (1) The CFS++/NACS HDF5 reader
#   (2) A reader for Geomview .off polyhedral surface geometry files, which 
#       are e.g. produced by CGAL.
#   (3) Some additional color maps (GiD, Matlab, etc.).
#   (4) Fixes to the calculation of derivatives on second order quads and
#       triangles (cf. http://www.vtk.org/Bug/view.php?id=13455).
#   (5) A reduced list of VisItBridge readers, to make sure our own reader is
#       the only one, which can read .h5 files. Otherwise the user always has
#       to choose the correct reader, when opening a file
#   (6) Some small fixes to the VisItBridge CGNS reader for TRIA6 and PYRA13
#       element types.
#
# These additions reside in the following files:
#
#   paraview-srcdir/
#   ├── ParaViewCore
#   │   └── ServerManager
#   │       └── SMApplication
#   │           └── Resources
#   │               └── readers.xml (1) & (2)
#   ├── Qt
#   │   └── Components
#   │       └── Resources
#   │           └── XML
#   │               └── ColorMaps.xml (3)
#   ├── Utilities
#   │   └── VisItBridge
#   │       └── databases (5)
#   │           ├── CGNS
#   │           │   └── avtCGNSFileFormat.C (6)
#   │           ├── CMakeLists.txt
#   │           └── visit_readers.xml
#   └── VTK
#       ├── Common
#       │   └── DataModel (4)
#       │       ├── vtkBiQuadraticQuad.cxx
#       │       ├── vtkQuadraticQuad.cxx
#       │       └── vtkQuadraticTriangle.cxx
#       ├── IO 
#       │   ├── CFSReader (1)
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
#       │   └── OFFReader (2)
#       │       ├── CMakeLists.txt
#       │       ├── module.cmake
#       │       ├── vtkOFFReader.cxx
#       │       └── vtkOFFReader.h
#       └── ThirdParty
#           └── hdf5
#               └── CMakeLists.txt (1)
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
# ${CFS_BINARY_DIR}/cfsdeps/pvsb/src/ \
#                  pvsb-build/paraview/src
#                   ├── paraview        ->  PARAVIEW_SOURCE_DIR
#                   ├── paraview-build  ->  PARAVIEW_BINARY_DIR
#                   └── paraview-stamp
#
#  Once   the   Superbuild   is   finished,    one   can   therefore   go   to
# PARAVIEW_BINARY_DIR and  do ordinary development,  since all changes  due to
# CFS++ are  already included  in PARAVIEW_SOURCE_DIR. If  CMAKE_BUILD_TYPE is
# set  to Debug  in the  CFS_BINARY_DIR,  it will  also  be set  to Debug  for
# ParaView.  Once development is finished, the  changes need to be just copied
# back to CFS_SOURCE_DIR/cfsdeps/paraview/paraview-srcdir.
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
set(pvsb_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/pvsb")
set(pvsb_install  "${CMAKE_CURRENT_BINARY_DIR}/paraview")
set(pvsb_source  "${pvsb_prefix}/src/pvsb")

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
  -DCMAKE_INSTALL_PREFIX:PATH=${pvsb_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
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
  # Of course we want to build ParaView
  -DENABLE_paraview:BOOL=ON
  # We need VisItBridge for the CGNS reader
  -DENABLE_visitbridge:BOOL=ON
  -DENABLE_png:BOOL=ON
  -DENABLE_szip:BOOL=OFF
  # ParaView has builtin versions of the following libraries
  -DENABLE_ffmpeg:BOOL=OFF
  -DENABLE_fontconfig:BOOL=OFF
  -DENABLE_freetype:BOOL=OFF
  -DENABLE_libxml2:BOOL=OFF
  # We use the system Qt4 or the one from CFSDEPS, if the system Qt4 is to old.
  -DENABLE_qt:BOOL=ON
  -DUSE_SYSTEM_qt=ON
  -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
  # We use the system Python, to take advantage of its installed packages.
  -DENABLE_python:BOOL=ON
  -DUSE_SYSTEM_python=ON
  # We use Boost, CGNS, HDF5 and zlib from CFSDEPS instead of Superbuild.
  -DENABLE_boost:BOOL=OFF
  -DENABLE_hdf5:BOOL=OFF
  -DENABLE_cgns:BOOL=OFF
  -DENABLE_zlib:BOOL=OFF
  -DBoost_DIR:PATH=Boost_DIR_NOTFOUND
  -DBoost_INCLUDE_DIR:PATH=${Boost_INCLUDE_DIR}
  -DBoost_LIBRARY_DIRS:STRING=${Boost_LIBRARY_DIR}
  -DCGNS_INCLUDE_DIR:PATH=${CGNS_INCLUDE_DIR}
  -DCGNS_LIBRARY:FILEPATH=${CGNS_SHARED_LIBRARY}
  -DHDF5_CXX_LIBRARY:FILEPATH=${HDF5_CPP_SHARED_LIBRARY}
  -DHDF5_C_LIBRARY:FILEPATH=${HDF5_SHARED_LIBRARY}
  -DHDF5_DIR:PATH=${CFS_BINARY_DIR}/share/cmake/hdf5
  -DHDF5_HL_LIBRARY:FILEPATH=${HDF5_LT_SHARED_LIBRARY}
  -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
  -DZLIB_LIBRARY:FILEPATH=${ZLIB_SHARED_LIBRARY}
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
SET(PFN "${pvsb_prefix}/pvsb-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Find the Git executable.
#-------------------------------------------------------------------------------
Find_Package(Git)

PRECOMPILED_ZIP_CXX(PRECOMPILED_PCKG_FILE "paraview" "no_ver")  
  
SET(PREFIX_DIR "${pvsb_prefix}")

SET(ZIPFROMCACHE "${pvsb_prefix}/pvsb-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${pvsb_prefix}/pvsb-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The ParaView Superbuild pvsb external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(pvsb
    PREFIX "${pvsb_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  IF(0) # GIT_FOUND)
    # Clone Git repo for ParaView Superbuild 4.1.
    ExternalProject_Add(pvsb
      DEPENDS ${CFS_PV_DEPENDENCIES}
      PREFIX ${pvsb_prefix}
      GIT_REPOSITORY http://paraview.org/ParaViewSuperbuild.git
      GIT_TAG "5867b1c34b73aee16cc9411007b2d70a4fad73a4"
      SOURCE_DIR ${pvsb_source}
      PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
      CMAKE_ARGS
        ${CMAKE_ARGS}
    )
  ELSE()
    # If we do not have the Git executable available, fall back
    # to downloading a zipped archive of the Git repo.
    MESSAGE("paraview_source ${paraview_source}")
    ExternalProject_Add(pvsb
      DEPENDS ${CFS_PV_DEPENDENCIES}
      PREFIX ${pvsb_prefix}
      DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview
      SOURCE_DIR ${pvsb_source}
      URL ${PARAVIEW_SB_URL}/${PARAVIEW_SB_TGZ}
      URL_MD5 ${PARAVIEW_SB_MD5}
      PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
      CMAKE_ARGS
        ${CMAKE_ARGS}
    )
  ENDIF()
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(pvsb cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  paraview-superbuild
)


