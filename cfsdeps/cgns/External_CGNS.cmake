#-------------------------------------------------------------------------------
# CFD General Notation System (CGNS)
# Needed for ADF routines by STARCCM+ reader. Also provides cgnsview, which
# can be used to view .cgns files (HDF5 and ADF) and .ccm files (ADF).
#
# Project Homepage
# http://www.cgns.org
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to CGNS sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(cgns_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cgns")
set(cgns_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(cgns_source  "${cgns_prefix}/src/cgns31")
set(cgns25_source  "${cgns_prefix}/src/cgns25")

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${cgns_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_MAKE_PROGRAM:BOOL=${CMAKE_MAKE_PROGRAM}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DENABLE_HDF5:BOOL=ON
  -DENABLE_LEGACY:BOOL=ON
  -DENABLE_64BIT:BOOL=OFF
  -DENABLE_TESTS:BOOL=ON
  -DHDF5_INCLUDE_PATH:PATH=${cgns_install}/include
  -DHDF5_LIBRARY:FILEPATH=${HDF5_SHARED_LIBRARY}
  -DHDF5_NEED_ZLIB:BOOL=ON
  -DZLIB_LIBRARY:FILEPATH=${ZLIB_SHARED_LIBRARY}
  -DCGNS_BUILD_SHARED:BOOL=ON
  -DCGNS_USE_SHARED:BOOL=OFF
  -DBUILD_CGNSTOOLS:BOOL=ON
  # We do not want to see warning messages from external projects
  -DCMAKE_C_FLAGS:STRING=${CFLAGS}
  -DCMAKE_CXX_FLAGS:STRING=${CFLAGS}
)

#-------------------------------------------------------------------------------
# Let's only build the CGNS tools on platforms, which provide a TCL 
# interpreter. There are certainly TCL interpreters on Windows and Mac, but we
# have not installed them in our nightly test systems.
#-------------------------------------------------------------------------------
SET(BUILD_CGNSTOOLS ON)
IF(MINGW)
  SET(BUILD_CGNSTOOLS OFF)
ELSEIF(CFS_DISTRO STREQUAL "MACOSX")
  IF(CMAKE_CROSSCOMPILING)
    SET(BUILD_CGNSTOOLS OFF)
  ENDIF()
ENDIF()

LIST(APPEND CMAKE_ARGS
  -DBUILD_CGNSTOOLS:BOOL=${BUILD_CGNSTOOLS}
)

IF(CFS_DISTRO STREQUAL "MACOSX")
  # Explicitly set build architectures and  system SDK root dir to match
  # those given in CMake.
  
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    "-DCMAKE_EXE_LINKER_FLAGS:STRING=-L/usr/X11/lib -lX11 -lXmu"
  )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgns/cgns-patch.cmake.in")
SET(PFN "${cgns_prefix}/cgns-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://heanet.dl.sourceforge.net/project/cgns/cgnslib_3.1/cgnslib_3.1.4-2.tar.gz"
  "http://mirror.transact.net.au/pub/sourceforge/c/project/cg/cgns/cgnslib_3.1/cgnslib_3.1.4-2.tar.gz"
  "${CGNS_URL}/${CGNS_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/cgns/${CGNS_GZ}")
SET(MD5_SUM ${CGNS_MD5})

SET(DLFN "${cgns_prefix}/cgns-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  )

#-------------------------------------------------------------------------------
# The CGNS external project
#-------------------------------------------------------------------------------
ExternalProject_Add(cgns31
  DEPENDS hdf5-shared zlib
  PREFIX ${cgns_prefix}
  SOURCE_DIR ${cgns_source}
  URL ${LOCAL_FILE}
  URL_MD5 ${CGNS_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  LIST_SEPARATOR ,
  CMAKE_ARGS
     ${CMAKE_ARGS}
    )

#-------------------------------------------------------------------------------
# Add custom download step to be able to download from a list of mirrors
# instead of just a single URL.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(cgns31 cfsdeps_download
   COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
   DEPENDERS download
   DEPENDS "${DLFN}"
   WORKING_DIRECTORY ${cgns_prefix}
)


SET(PFN25_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgns/cgns25-patch.cmake.in")
SET(PFN25 "${cgns_prefix}/cgns25-patch.cmake")
CONFIGURE_FILE("${PFN25_TEMPL}" "${PFN25}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://ftp.de.cw.net/pub/FreeBSD/ports/distfiles/cgnslib_2.5-5.tar.gz"
  "http://pkgs.fedoraproject.org/repo/pkgs/cgnslib/cgnslib_2.5-5.tar.gz/ae2a2e79b99d41c63e5ed5f661f70fd9/cgnslib_2.5-5.tar.gz"
  "${CGNS_URL}/${CGNS25_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/cgns/${CGNS25_GZ}")
SET(MD5_SUM ${CGNS25_MD5})

SET(DLFN25 "${cgns_prefix}/cgns25-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN25}"
  @ONLY
  )

#-------------------------------------------------------------------------------
# The CGNS 2.5 external project
#-------------------------------------------------------------------------------
ExternalProject_Add(cgns25
  DEPENDS hdf5-shared zlib
  PREFIX ${cgns_prefix}
  SOURCE_DIR ${cgns25_source}
  URL ${LOCAL_FILE}
  URL_MD5 ${CGNS25_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN25}"
  LIST_SEPARATOR ,
  CMAKE_ARGS
     ${CMAKE_ARGS}
    )

#-------------------------------------------------------------------------------
# Add custom download step to be able to download from a list of mirrors
# instead of just a single URL.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(cgns25 cfsdeps_download
   COMMAND ${CMAKE_COMMAND} -P "${DLFN25}"
   DEPENDERS download
   DEPENDS "${DLFN25}"
   WORKING_DIRECTORY ${cgns_prefix}
)


#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  cgns31
  cgns25
)

#-------------------------------------------------------------------------------
# Determine paths of CGNS libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(CGNS_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns.3.1${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "CGNS 3.1 library.")
SET(CGNS25_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns.2.5${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "CGNS 2.5 library.")

IF(CFS_OS STREQUAL LINUX)
  SET(CGNS_SHARED_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns${CMAKE_SHARED_LIBRARY_SUFFIX}.3.1"
    CACHE FILEPATH "CGNS 3.1 shared library.")
  SET(CGNS25_SHARED_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns${CMAKE_SHARED_LIBRARY_SUFFIX}.2.5"
    CACHE FILEPATH "CGNS 2.5 shared library.")
ELSE()
  IF(WIN32)
    SET(CGNS_SHARED_LIBRARY
      "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgnsdll.3.1${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      CACHE FILEPATH "CGNS 3.1 shared library.")
    SET(CGNS25_SHARED_LIBRARY
      "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgnsdll.2.5${CMAKE_IMPORT_LIBRARY_SUFFIX}"
      CACHE FILEPATH "CGNS 2.5 shared library.")
  ELSE()
    SET(CGNS_SHARED_LIBRARY
      "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns.3.1${CMAKE_SHARED_LIBRARY_SUFFIX}"
      CACHE FILEPATH "CGNS 3.1 shared library.")
    SET(CGNS25_SHARED_LIBRARY
      "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns.2.5${CMAKE_SHARED_LIBRARY_SUFFIX}"
      CACHE FILEPATH "CGNS 2.5 shared library.")
  ENDIF()
ENDIF()

IF(EXISTS "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns${CMAKE_SHARED_LIBRARY_SUFFIX}")
  FILE(REMOVE "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns${CMAKE_SHARED_LIBRARY_SUFFIX}")
ENDIF()

SET(CGNS_INCLUDE_DIR ${CFS_BINARY_DIR}/include/cgns-3.1 CACHE PATH "CGNS 3.1 include directory")
SET(CGNS25_INCLUDE_DIR ${CFS_BINARY_DIR}/include/cgns-2.5 CACHE PATH "CGNS 2.5 include directory")

MARK_AS_ADVANCED(CGNS_INCLUDE_DIR)
