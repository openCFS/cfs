#-------------------------------------------------------------------------------
# Computational Geometry Algorithms Library (CGAL)
# As of April 2024, CGAL is used in CFS to enable enhanced element intersection  
# computations for non-conforming interfaces and grid interpolation techniques. 
# It is also used for 3D geometry approximations on surface nodes in the 
# layerGeneration feature and the curvilinear PML formulation.
#
# CGAL is only used as a supplementary build option in CFS. It is not part of the 
# default/distributed build, due to its GPL3 licence.
#
# Project Homepage
# http://www.cgal.org/
#
# Description copied from the CGAL project homepage:
# "CGAL is an open source software project that provides easy access to efficient 
# and reliable geometric algorithms in the form of a C++ library. CGAL is used in
# various areas needing geometric computation, such as geographic information systems, 
# computer aided design, molecular biology, medical imaging, computer graphics, and robotics.
# The library offers data structures and algorithms like triangulations, Voronoi diagrams, 
# Boolean operations on polygons and polyhedra, point set processing, arrangements of curves, 
# surface and volume mesh generation, geometry processing, alpha shapes, convex hull algorithms,
# shape reconstruction, AABB and KD trees... Explore the complete list of features and 
# capabilities by visiting the CGAL Package Overview.
# The CGAL data structures and algorithms are distributed under a dual license, namely under 
# the GPL v3+ and, alternatively, under a commercial license by GeometryFactory."
#-------------------------------------------------------------------------------

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(CFS_VERSION_MONTH "04")
set(PACKAGE_NAME "cgal")
set(PACKAGE_VER "5.6.1")
set(CGAL_VER ${PACKAGE_VER}) # for Dependencies.cc
set(PACKAGE_FILE "CGAL-${PACKAGE_VER}.tar.xz")
set(PACKAGE_MD5 "90ccd8f68894ab2f89933163afd66535")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://github.com/CGAL/cgal/releases/download/v${PACKAGE_VER}/${PACKAGE_FILE}") 

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# we do not need to compile as CGAL is now header only
use_c_and_fortran(OFF OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_header_only_standard_variables()

# this is the standard target for cmake projects. The files to package come from the install_manifest.txt
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON)

# set custom build arguments
set(DEPS_ARGS 
  ${DEPS_ARGS}
  -DWITH_EXAMPLES:BOOL=OFF
  -DWITH_DEMOS:BOOL=OFF
  -DWITH_CGAL_Core:BOOL=ON
  -DWITH_CGAL_Qt5:BOOL=OFF
  -DWITH_CGAL_ImageIO:BOOL=OFF
  -DCGAL_DIR:STRING="${CMAKE_BUILD_TYPE}/include/"
)

# --- it follows generic final block for cmake packages with no patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

# generate package creation script. We get the files from an install_manifest.txt
generate_packing_script_manifest()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

# for debugging purposes, dump all dependency variables
dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  create_external_cmake()

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()
    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()  
  endif()  
endif()

# add dependencies for cgal
add_dependencies(cgal boost gmp mpfr)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})


#set(CGAL_LIBRARY "${CFS_BINARY_DIR}/${LIB_SUFFIX}/libCGAL.a" CACHE FILEPATH "CGAL library.")



#-Old code below this line---------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


# #-------------------------------------------------------------------------------
# # Set paths to cgal sources according to ExternalProject.cmake 
# #-------------------------------------------------------------------------------
# set(cgal_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cgal")
# set(cgal_source  "${cgal_prefix}/src/cgal")
# set(cgal_install  "${CMAKE_CURRENT_BINARY_DIR}")

# string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
# IF(CMAKE_BUILD_TYPE STREQUAL "debug")
#   SET(CMAKE_BUILD_TYPE "Debug")
# ELSE(CMAKE_BUILD_TYPE STREQUAL "debug")
#   SET(CMAKE_BUILD_TYPE "Release")
# ENDIF(CMAKE_BUILD_TYPE STREQUAL "debug")

# # Make sure that static Boost libs are used under Windows.
# SET(CGAL_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -DBOOST_THREAD_USE_LIB")


# SET(CMAKE_ARGS
#   -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
#   -DCMAKE_INSTALL_PREFIX:PATH=${cgal_install}
#   -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
#   -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
#   -DCMAKE_CXX_FLAGS:STRING=${CGAL_CXX_FLAGS}
#   -DCMAKE_BUILD_TYPE:STRING=Release
#   -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
# )

# IF(CMAKE_TOOLCHAIN_FILE)
#   LIST(APPEND CMAKE_ARGS
#     -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
#   )
# ENDIF()

# #-------------------------------------------------------------------------------
# # Set names of patch file and template file.
# #-------------------------------------------------------------------------------
# SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgal/cgal-patch.cmake.in")
# SET(PFN "${cgal_prefix}/cgal-patch.cmake")
# CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

# SET(BOOST_SETUP_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgal/CGAL_SetupBoost.cmake.in")
# SET(BOOST_SETUP "${cgal_prefix}/CGAL_SetupBoost.cmake")
# CONFIGURE_FILE("${BOOST_SETUP_TEMPL}" "${BOOST_SETUP}" @ONLY) 

# #-------------------------------------------------------------------------------
# # Set up a list of publicly available mirrors, since the non-standard port 
# # number of the FTP server on the openCFS development server  may not be
# # accessible from behind firewalls.
# # Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# # used to configure the download CMake file for the library.
# #-------------------------------------------------------------------------------
# SET(MIRRORS
#   "https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-${CGAL_VER}/CGAL-${CGAL_VER}.tar.xz"
#   "http://archive.ubuntu.com/ubuntu/pool/universe/c/cgal/cgal_${CGAL_VER}.orig.tar.xz"
#   "https://gforge.inria.fr/frs/download.php/32361/${CGAL_BZ2}"
#   "${CFS_DS_SOURCES_DIR}/cgal/${CGAL_BZ2}"
# )
# SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/cgal/${CGAL_BZ2}")
# SET(MD5_SUM ${CGAL_MD5})

# SET(DLFN "${cgal_prefix}/cgal-download.cmake")
# CONFIGURE_FILE(
#   "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
#   "${DLFN}"
#   @ONLY
# )

# PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "cgal" "${CGAL_VER}")  
  
# # This should be either PREFIX_DIR (install manifest is used for zipping)
# # or INSTALL_DIR (install directory will be zipped)
# SET(TMP_DIR "${cgal_prefix}")

# SET(ZIPFROMCACHE "${cgal_prefix}/cgal-zipFromCache.cmake")
# CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

# SET(ZIPTOCACHE "${cgal_prefix}/cgal-zipToCache.cmake")
# CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

# #-------------------------------------------------------------------------------
# # Determine paths of CGAL libraries.
# #-------------------------------------------------------------------------------
# SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
# SET(CGAL_LIBRARY
#   "${LD}/libCGAL.a"
#   CACHE FILEPATH "CGAL library.")

# MARK_AS_ADVANCED(CGAL_INCLUDE_DIR)
# MARK_AS_ADVANCED(CGAL_LIBRARY)

# #-------------------------------------------------------------------------------
# # The CGAL external project
# #-------------------------------------------------------------------------------
# IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
#   #-------------------------------------------------------------------------------
#   # If precompiled package exists copy files from cache
#   #-------------------------------------------------------------------------------
#   ExternalProject_Add(cgal
#     PREFIX "${cgal_prefix}"
#     DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
#     PATCH_COMMAND ""
#     UPDATE_COMMAND ""
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND ""
#     INSTALL_COMMAND ""
#   )
# ELSE()
#   #-------------------------------------------------------------------------------
#   # If precompiled package does not exist build external project
#   #-------------------------------------------------------------------------------
#   ExternalProject_Add(cgal
#     DEPENDS boost zlib gmp mpfr
#     PREFIX "${cgal_prefix}"
#     URL ${LOCAL_FILE}
#     URL_MD5 ${CGAL_MD5}
#     PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
#     CMAKE_ARGS
#       ${CMAKE_ARGS}
#       -DCGAL_INSTALL_BIN_DIR:PATH=bin
#       -DCGAL_INSTALL_CMAKE_DIR:PATH=${LIB_SUFFIX}/CGAL
#       -DCGAL_INSTALL_LIB_DIR:PATH=${LIB_SUFFIX}
#       -DBUILD_SHARED_LIBS:BOOL=OFF
# #      -DBoost_INCLUDE_DIR:PATH=${cgal_install}/include
# #      -DBoost_LIBRARY_DIRS:PATH=${cgal_install}/${LIB_SUFFIX}/$ARCH_STR
# #      -DBoost_THREAD_LIBRARY:FILEPATH=${BOOST_THREAD_LIB}
# #      -DBoost_THREAD_LIBRARY_RELEASE:FILEPATH=${BOOST_THREAD_LIB_RELEASE}
# #      -DBoost_THREAD_LIBRARY_DEBUG:FILEPATH=${BOOST_THREAD_LIB_DEBUG}
#       -DWITH_CGAL_Qt3:BOOL=OFF
#       -DWITH_CGAL_Qt3/CMakeLists.txt:BOOL=OFF
#       -DWITH_CGAL_Qt4:BOOL=OFF
#       -DWITH_CGAL_Qt4/CMakeLists.txt:BOOL=OFF
#       -DWITH_CGAL_Core/CMakeLists.txt:BOOL=OFF
#       -DWITH_CGAL_ImageIO/CMakeLists.txt:BOOL=OFF
#       -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
#       -DZLIB_LIBRARY:PATH=${ZLIB_LIBRARY}
#       -DGMP_INCLUDE_DIR:PATH=${GMP_INCLUDE_DIR}
#       -DGMP_LIBRARIES:FILEPATH=${GMP_LIBRARY}
#       -DGMP_LIBRARIES_DIR:PATH=${cgal_install}/${LIB_SUFFIX}
#       -DWITH_GMP:BOOL=ON
#       -DMPFR_INCLUDE_DIR:PATH=${MPFR_INCLUDE_DIR}
#       -DMPFR_LIBRARIES:FILEPATH=${MPFR_LIBRARY}
#       -DWITH_MPFR:BOOL=ON
#     BUILD_BYPRODUCTS ${CGAL_LIBRARY}
#   )
  
#   #-------------------------------------------------------------------------------
#   # Add custom download step to be able to download from a list of mirrors
#   # instead of just a single URL.
#   #-------------------------------------------------------------------------------
#   ExternalProject_Add_Step(cgal cfsdeps_download
#     COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
#     DEPENDERS download
#     DEPENDS "${DLFN}"
#     WORKING_DIRECTORY ${cgal_prefix}
#   )
  
#   IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
#     #-------------------------------------------------------------------------------
#     # Add custom step to zip a precompiled package to the cache.
#     #-------------------------------------------------------------------------------
#     ExternalProject_Add_Step(cgal cfsdeps_zipToCache
#       COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
#       DEPENDEES install
#       DEPENDS "${ZIPTOCACHE}"
#       WORKING_DIRECTORY ${CFS_BINARY_DIR}
#     )
#   ENDIF()
# ENDIF()

# # copy the appropriate special license in the license folder
# file(COPY "${CFS_SOURCE_DIR}/cfsdeps/LICENSE.binary.GPL" DESTINATION "${CFS_BINARY_DIR}/license/")

# #-------------------------------------------------------------------------------
# # Add project to global list of CFSDEPS
# #-------------------------------------------------------------------------------
# SET(CFSDEPS
#   ${CFSDEPS}
#   cgal
# )

# SET(CGAL_INCLUDE_DIR "${CFS_BINARY_DIR}/include" "${CFS_BINARY_DIR}/src/cgal/include" "${CFS_BINARY_DIR}/src/cgal-build/include"  "${CFS_BINARY_DIR}/src/cgal/include/CGAL")
