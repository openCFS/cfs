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
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
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
# number of the FTP server on the openCFS development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.vtk.org/files/release/7.1/${VTK_TAR}"
  "${CFS_DS_SOURCES_DIR}/vtk/${VTK_TAR}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/vtk/${VTK_TAR}")
SET(MD5_SUM ${VTK_MD5})

SET(DLFN "${vtk_prefix}/vtk-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY) 

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/vtk/license/" DESTINATION "${CFS_BINARY_DIR}/license/vtk" )



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
SET(PRF "lib")
SET(EXT "a")
IF(WIN32)
  SET(PRF "")
  SET(EXT "lib")
  SET(VTK_LIBRARY
     ${LD}/${PRF}vtkIOParallel-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkexoIIc-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOXML-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOXMLParser-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkexpat-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIONetCDF-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkNetCDF_cxx-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkNetCDF-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkhdf5_hl-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkhdf5-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOImage-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtktiff-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkpng-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkjpeg-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkmetaio-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOGeometry-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkjsoncpp-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOEnSight-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersSMP-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersParallel-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkRenderingCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkParallelCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOLegacy-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOCore-${VTK_VERSION}.${EXT}
     #${LD}/${PRF}vtkIOXML-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkzlib-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersModeling-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersGeometry-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersSources-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersExtraction-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersStatistics-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkalglib-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkImagingFourier-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkImagingCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersGeneral-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkDICOMParser-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonExecutionModel-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonComputationalGeometry-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonDataModel-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonTransforms-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonSystem-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonMisc-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonMath-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtksys-${VTK_VERSION}.${EXT}
  CACHE FILEPATH "VTK library.")
ELSE(WIN32)
  SET(VTK_LIBRARY
     ${LD}/${PRF}vtkIOParallel-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkexoIIc-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOXML-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOXMLParser-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkexpat-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIONetCDF-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkNetCDF_cxx-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkNetCDF-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkhdf5_hl-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkhdf5-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOImage-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtktiff-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkpng-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkjpeg-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkmetaio-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOGeometry-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkjsoncpp-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOEnSight-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersSMP-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersParallel-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkRenderingCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkParallelCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOLegacy-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkIOCore-${VTK_VERSION}.${EXT}
     #${LD}/${PRF}vtkIOXML-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkzlib-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersModeling-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersGeometry-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersSources-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersExtraction-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersStatistics-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkalglib-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkImagingFourier-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkImagingCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersGeneral-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkFiltersCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkDICOMParser-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonExecutionModel-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonComputationalGeometry-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonDataModel-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonTransforms-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonSystem-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonMisc-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonMath-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtkCommonCore-${VTK_VERSION}.${EXT}
     ${LD}/${PRF}vtksys-${VTK_VERSION}.${EXT}
     dl
  CACHE FILEPATH "VTK library.")
ENDIF(WIN32)
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
    BUILD_BYPRODUCTS ${VTK_LIBRARY}
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
      # enable builds on systems without OpenGL
      -DVTK_RENDERING_BACKEND:STRING=None
  #if we add more file reader we should make this optional depending on USE_ENSIGHT
      -DModule_vtkIOEnSight:BOOL=ON
      -DModule_vtkIOParallel:BOOL=ON
      -DModule_vtkIOXML:BOOL=ON
      -DVTK_SMP_IMPLEMENTATION_TYPE:STRING="TBB"
     BUILD_BYPRODUCTS ${VTK_LIBRARY}
  #-DVTK_INSTALL_INCLUDE_DIR:PATH=${CFS_BINARY_DIR}/include/vtk
  #    -DVTK_INSTALL_LIBRARY_DIR:PATH=${CFS_BINARY_DIR}/${LIB_SUFFIX}/vtk
  #    -DVTK_INSTALL_ARCHIVE_DIR:PATH=${CFS_BINARY_DIR}/${LIB_SUFFIX}/vtk
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
