# mv hdf5-1.8.8.tar.gz HDF5-prefix/src/hdf5-1.8.8.tar.gz
# mv libpng-1.5.7.tar.gz png-prefix/src/libpng-1.5.7.tar.gz
# cp szip-2.1.tar.gz szip-prefix/src/szip-2.1.tar.gz
# CGNS-prefix/src/cgnslib_2.5-5.tar.gz
# VRPN-prefix/src/vrpn_07_29.zip

#-------------------------------------------------------------------------------
# Set prefix path and path to paraview sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(paraview_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/paraview")
set(paraview_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(paraview_source  "${paraview_prefix}/src/paraview")

#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/paraview/paraview-patch.cmake.in")
#SET(PFN "${paraview_prefix}/paraview-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${paraview_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCFS_PYTHON_VERSION_STRING:STRING=${PYTHON_VERSION_STRING}
  -DCFS_PYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE}
  -DCFS_PYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}
  -DCFS_PYTHON_INCLUDE_DIR:FILEPATH=${PYTHON_INCLUDE_DIR}
  -DCFS_DEPS_CACHE_DIR:PATH=${CFS_DEPS_CACHE_DIR}
)

#-------------------------------------------------------------------------------
# The paraview-static external project
# We do not need it at the moment since we use the HDF5 inside ParaView for
# our plugin.
#-------------------------------------------------------------------------------
#ExternalProject_Add(paraview-shared
#  DEPENDS zlib-shared
#  PREFIX ${paraview_prefix}
#  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview
#  SOURCE_DIR ${paraview_source}
#  URL ${HDF5_URL}/${HDF5_GZ}
#  URL_MD5 ${HDF5_MD5}
#  CMAKE_ARGS
#     ${CMAKE_ARGS}
#    -DBUILD_SHARED_LIBS:BOOL=ON
#    -DHDF5_BUILD_CPP_LIB:BOOL=ON
#    -DHDF5_BUILD_HL_LIB:BOOL=ON
#    -DHDF5_BUILD_FORTRAN:BOOL=OFF
#    -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
#    -DZLIB_INCLUDE_DIR:PATH=${paraview_install}/include
#    -DZLIB_LIBRARY:FILEPATH=${paraview_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_SHARED_LIBRARY_PREFIX}z${CMAKE_SHARED_LIBRARY_SUFFIX}
#    )

#-------------------------------------------------------------------------------
# The paraview-static external project
#-------------------------------------------------------------------------------
ExternalProject_Add(paraview-download
  PREFIX ${paraview_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview
  SOURCE_DIR ${paraview_source}
  URL ${PARAVIEW_URL}/${PARAVIEW_GZ}
  URL_MD5 ${PARAVIEW_MD5}
  CMAKE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

# MESSAGE(FATAL_ERROR "CFS_PV_CMAKE_COMMAND ${CFS_PV_CMAKE_COMMAND}")
ExternalProject_Add(paraview-superbuild
  DEPENDS ${CFS_PV_DEPENDS} paraview-download
  PREFIX ${paraview_prefix}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR "${paraview_source}/SuperBuild"
  CMAKE_ARGS
    ${CMAKE_ARGS}
    )

#  CMAKE_COMMAND ${CFS_PV_CMAKE_COMMAND}

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/paraview/paraview-patch.cmake.in")
SET(PFN "${paraview_prefix}/paraview-patch.cmake")
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
ExternalProject_Add_Step(paraview-superbuild custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${paraview_source}
)

# ln -s libhdf5_debug.so.1.8.8 libhdf5.so

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  paraview-superbuild
)


