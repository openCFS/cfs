#-------------------------------------------------------------------------------
# Visualization ToolKit Library
#
# Project Homepage
# http://www.vtk.org/
#-------------------------------------------------------------------------------

set(vtk_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/vtk")
set(vtk_source  "${vtk_prefix}/src/vtk")
set(vtk_install  "${CMAKE_CURRENT_BINARY_DIR}/vtk")

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
IF(CMAKE_BUILD_TYPE STREQUAL "debug")
  SET(CMAKE_BUILD_TYPE "Debug")
ELSE(CMAKE_BUILD_TYPE STREQUAL "debug")
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF(CMAKE_BUILD_TYPE STREQUAL "debug")

SET(VTK_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -DBOOST_THREAD_USE_LIB")

SET(CMAKE_ARGS
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_INSTALL_PREFIX:PATH=${vtk_install}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${VTK_CXX_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
)

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/vtk/vtk-patch.cmake.in")
SET(PFN "${vtk_prefix}/vtk-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.vtk.org/files/release/7.1/${VTK_TAR}"
  "${VTK_URL}/${VTK_TAR}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/vtk/${VTK_TAR}")
SET(MD5_SUM ${VTK_MD5})

SET(DLFN "${vtk_prefix}/vtk-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY) 

PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "vtk" "${VTK_VERSION}") 

# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${vtk_prefix}")

SET(ZIPFROMCACHE "${vtk_prefix}/vtk-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${vtk_prefix}/vtk-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Set linking libraries, the exact order of .a files is of highest importance
#-------------------------------------------------------------------------------
SET(LD "${vtk_install}/lib")
SET(VTK_LIBRARY
     ${LD}/libvtkIOParallel-${VTK_VERSION}.a
     ${LD}/libvtkexoIIc-${VTK_VERSION}.a
     ${LD}/libvtkIOXML-${VTK_VERSION}.a
     ${LD}/libvtkIOXMLParser-${VTK_VERSION}.a
     ${LD}/libvtkexpat-${VTK_VERSION}.a
     ${LD}/libvtkIONetCDF-${VTK_VERSION}.a
     ${LD}/libvtkNetCDF_cxx-${VTK_VERSION}.a
     ${LD}/libvtkNetCDF-${VTK_VERSION}.a
     ${LD}/libvtkhdf5_hl-${VTK_VERSION}.a
     ${LD}/libvtkhdf5-${VTK_VERSION}.a
     ${LD}/libvtkIOImage-${VTK_VERSION}.a
     ${LD}/libvtktiff-${VTK_VERSION}.a
     ${LD}/libvtkpng-${VTK_VERSION}.a
     ${LD}/libvtkjpeg-${VTK_VERSION}.a
     ${LD}/libvtkmetaio-${VTK_VERSION}.a
     ${LD}/libvtkIOGeometry-${VTK_VERSION}.a
     ${LD}/libvtkjsoncpp-${VTK_VERSION}.a
     ${LD}/libvtkIOEnSight-${VTK_VERSION}.a
     ${LD}/libvtkFiltersSMP-${VTK_VERSION}.a
     ${LD}/libvtkFiltersParallel-${VTK_VERSION}.a
     ${LD}/libvtkRenderingCore-${VTK_VERSION}.a
     ${LD}/libvtkParallelCore-${VTK_VERSION}.a
     ${LD}/libvtkIOLegacy-${VTK_VERSION}.a
     ${LD}/libvtkIOCore-${VTK_VERSION}.a
     #${LD}/libvtkIOXML-${VTK_VERSION}.a
     ${LD}/libvtkzlib-${VTK_VERSION}.a
     ${LD}/libvtkFiltersModeling-${VTK_VERSION}.a
     ${LD}/libvtkFiltersGeometry-${VTK_VERSION}.a
     ${LD}/libvtkFiltersSources-${VTK_VERSION}.a
     ${LD}/libvtkFiltersExtraction-${VTK_VERSION}.a
     ${LD}/libvtkFiltersStatistics-${VTK_VERSION}.a
     ${LD}/libvtkalglib-${VTK_VERSION}.a
     ${LD}/libvtkImagingFourier-${VTK_VERSION}.a
     ${LD}/libvtkImagingCore-${VTK_VERSION}.a
     ${LD}/libvtkFiltersGeneral-${VTK_VERSION}.a
     ${LD}/libvtkFiltersCore-${VTK_VERSION}.a
     ${LD}/libvtkDICOMParser-${VTK_VERSION}.a
     ${LD}/libvtkCommonExecutionModel-${VTK_VERSION}.a
     ${LD}/libvtkCommonComputationalGeometry-${VTK_VERSION}.a
     ${LD}/libvtkCommonDataModel-${VTK_VERSION}.a
     ${LD}/libvtkCommonTransforms-${VTK_VERSION}.a
     ${LD}/libvtkCommonSystem-${VTK_VERSION}.a
     ${LD}/libvtkCommonMisc-${VTK_VERSION}.a
     ${LD}/libvtkCommonMath-${VTK_VERSION}.a
     ${LD}/libvtkCommonCore-${VTK_VERSION}.a
     ${LD}/libvtksys-${VTK_VERSION}.a
     dl
  CACHE FILEPATH "VTK library.")
MARK_AS_ADVANCED(VTK_LIBRARY)

#-------------------------------------------------------------------------------
# The VTK external project
#-------------------------------------------------------------------------------

IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(vtk
    PREFIX "${vtk_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(vtk
    #BUILD_COMMAND make -j4
    DEPENDS boost zlib 
    PREFIX "${vtk_prefix}"
    URL ${LOCAL_FILE}
    URL_MD5 ${VTK_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DVTK_Group_Rendering:BOOL=OFF
      -DVTK_Group_StandAlone:BOOL=OFF
      -DModule_vtkFiltersParallel:BOOL=ON
      -DModule_vtkFiltersSMP:BOOL=ON
  #if we add more file reader we should make this optional depending on USE_ENSIGHT
      -DModule_vtkIOEnSight:BOOL=ON
      -DModule_vtkIOParallel:BOOL=ON
      -DModule_vtkIOXML:BOOL=ON
      -DVTK_SMP_IMPLEMENTATION_TYPE:STRING="TBB"
     BUILD_BYPRODUCTS ${VTK_LIBRARY}
  #-DVTK_INSTALL_INCLUDE_DIR:PATH=${CFS_BINARY_DIR}/include/vtk
  #    -DVTK_INSTALL_LIBRARY_DIR:PATH=${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/vtk
  #    -DVTK_INSTALL_ARCHIVE_DIR:PATH=${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/vtk
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(vtk cfsdeps_download
     COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
     DEPENDERS download
     DEPENDS "${DLFN}"
     WORKING_DIRECTORY ${vtk_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(vtk cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF()


#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS ${CFSDEPS} vtk)

SET(VTK_INCLUDE_DIR "${vtk_install}/include/vtk-${VTK_VERSION}" CACHE FILEPATH "VTK include directory.")
MARK_AS_ADVANCED(VTK_INCLUDE_DIR)
set(VTK_DIR "${vtk_install}/lib/cmake/vtk-${VTK_VERSION}") # from https://github.com/statismo/statismo/blob/master/superbuild/External-VTK.cmake

#find_package(VTK) # does not work
# this is how it should be done: http://cmake.3232098.n2.nabble.com/How-to-use-VTK-as-an-ExternalProject-td6002193.html
# explantion: http://cmake.3232098.n2.nabble.com/Question-regarding-External-Project-add-and-VTK-td7587557.html
