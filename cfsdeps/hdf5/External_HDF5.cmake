# The native openCFS binary data format (export but also reading) is hdf5
# https://www.hdfgroup.org/solutions/hdf5/
#
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "hdf5")

set(PACKAGE_VER "1.14.6")
set(PACKAGE_FILE "hdf5-${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "63426c8e24086634eaf9179a8c5fe9e5")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

#https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz
set(PACKAGE_MIRRORS "https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a fortran compiler
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# generates PACKAGE_LIBARAY with lib<package>.a/.dll - on Windows also the prefix lib is used, what is uncommon.
if(UNIX OR NOT DEBUG) # Linux, macOS and (Windows in release mode) are the same
  set_package_library_list_lib_prefix("hdf5_hl;hdf5")
else() # WIN32 AND DEBUG: on Windows we have the debug libraries with a D suffix
  set_package_library_list_lib_prefix("hdf5_hl_D;hdf5_D") # the libhdf5... comes from the macro
endif()
# creates HDF5_LIBARAY as CACHE variable, hence it will not be overwritten once in cache!
set_standard_variables()

# we need share/cmake/hdf5-config.cmake for cgns, therefore install_manifest.txt
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags

# dump_depencency_variables()

# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DHDF5_INSTALL_BIN_DIR:PATH=bin
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DBUILD_STATIC_LIBS:BOOL=ON
  -DBUILD_TESTING:BOOL=OFF
  -DHDF5_BUILD_CXX:BOOL=OFF # we got rit of the outdated C++ interface
  -DHDF5_BUILD_HL_LIB:BOOL=ON
  -DHDF5_BUILD_FORTRAN:BOOL=OFF
  -DHDF5_BUILD_EXAMPLES:BOOL=OFF
  -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
  -DH5_ZLIB_HEADER:FILEPATH=${ZLIB_INCLUDE_DIR}/zlib.h # we need to tell hdf5 explicitly which zlib header to use
  -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
  -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
  -DHDF5_BUILD_TOOLS:BOOL=OFF # no binaries wanted, use system hdf5 tools
  -DHDF5_BUILD_UTILS:BOOL=OFF 
  -DHDF5_PACK_MACOS_DMG:BOOL=OFF
  -DTEST_SHELL_SCRIPTS:BOOL=OFF
  -DUSE_SHARED_LIBS:BOOL=OFF)

#if(POLICY CMP0075)
#  list(APPEND DEPS_ARGS -DCMAKE_POLICY_DEFAULT_CMP0075=NEW)# prevent Policy CMP0075 is not set: Include file check macros honor
#endif()

# --- it follows generic final block for cmake packages with no patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/"
     DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

# generate package creation script.
generate_packing_script_manifest()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

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
  endif()
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})

# hdf5 depends on zlib
add_dependencies(hdf5 zlib)
