# Ilupack is a closed source academic iterative solve by M. Bollhoefer, TU-Braunschweig
# by special agreement we have the source available which must not be redistributed

set(ilupack_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ilupack")
set(ilupack_source  "${ilupack_prefix}/src/ilupack")
set(ilupack_install  "${ilupack_prefix}/tmp")



#-------------------------------------------------------------------------------
# The ilupack source is closed source! It must not be redistributed!!
# There is a copy at ${CFS_DS_WEBDAV}/ but it has certificate issuses, therefore there
# is also a mirror by the FAU CFS optimization group.
# To obfuscate the download,. the filename is replaced by its md5 sum.
#-------------------------------------------------------------------------------
#SET(MIRRORS 
#  "${CFS_FAU_MIRROR}/sources/ilupack/${ILUPACK_MD5}"
#  "${CFS_DS_WEBDAV}/cfsdeps/sources/ilupack/MD5/${ILUPACK_MD5}")
  
#SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ilupack/${ILUPACK_GZ}") 

SET(LOCAL_FILE "/home/sri/code/Ilupack_Task_OPENMP/Ilupack_Task_OPENMP.tgz")

#SET(MD5_SUM ${ILUPACK_MD5})
SET(MD5_SUM 17b39bb3fedd5264cacca5a4d7e23d43)

#SET(DLFN "${ilupack_prefix}/ilupack-download.cmake")
#CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)


PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "ilupack" "${ILUPACK_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${ilupack_install}")


IF ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  SET(MAKE_PLATFORM "MAC")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET(MAKE_PLATFORM "GNU64")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
 SET(MAKE_PLATFORM "INTEL64")
endif()



# After the installation we copy to cfs
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ilupack/ilupack-post_install.cmake.in")
SET(PI "${ilupack_prefix}/ilupack-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

SET(ZIPFROMCACHE "${ilupack_prefix}/ilupackzipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${ilupack_prefix}/ilupack-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of ILUPACK libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(ILUPACK_SHARE "${LD}/libclock.a;${LD}/libsparsenew.a;${LD}/libvector.a")  
SET(ILUPACK_OPENMP "${LD}/libspdil.a;${LD}/libtools.a;${LD}/libbasic.a")  

SET(ILUPACK_LIBRARY
  "${ILUPACK_OPENMP};${LD}/libilupack.a;${LD}/libamd.a;${LD}/libcamd.a;${LD}/libblaslike.a;${LD}/libmetis.a;${LD}/libmetisomp.a;${LD}/libmtmetis.a;${LD}/libmumps.a;${LD}/libsparspak.a;${ILUPACK_SHARE};"               
  CACHE FILEPATH "ILUPACK library.")
  
MARK_AS_ADVANCED(ILUPACK_LIBRARY)

SET(ILUPACK_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
MARK_AS_ADVANCED(ILUPACK_INCLUDE_DIR)


#-------------------------------------------------------------------------------
# The ilupack external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ilupack
    PREFIX "${ARPACK_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
 
ELSE()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project (double real)
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ilupack
    #DEPENDS ilupack_external_data lapack metis suitesparsecd 
    PREFIX "${ilupack_prefix}"
    SOURCE_DIR "${ilupack_source}"
    URL ${LOCAL_FILE}
    URL_MD5 ${ILUPACK_MD5}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ./COMPILE_ALL.sh "${MAKE_PLATFORM}"
    BUILD_COMMAND ""
    INSTALL_COMMAND ./COPY_LIB.sh "${ilupack_install}"
   )
    
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(ilupack cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${ilupack_prefix}
  )
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache. complex seems sufficient
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(ilupack cfsdeps_zipToCache
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
  ilupack
)

SET(ILUPACK_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
MARK_AS_ADVANCED(ILUPACK_INCLUDE_DIR)




