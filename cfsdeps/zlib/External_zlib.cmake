# A Massively Spiffy Yet Delicately Unobtrusive Compression Library
# boost needs this lib
# https://github.com/madler/zlib
# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "zlib")
set(PACKAGE_VER "1.3.1")
set(ZLIB_VER ${PACKAGE_VER}) # required by boost 
set(PACKAGE_FILE "${PACKAGE_NAME}-${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "9855b6d802d7fe5b7bd5b196a2271655")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.


set(PACKAGE_MIRRORS "https://github.com/madler/zlib/releases/download/v${PACKAGE_VER}/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# pure C
use_c_and_fortran(ON OFF)

# When preCICE is built, boost becomes a shared library (see External_Boost.cmake).
# A shared libboost_iostreams.so cannot embed the non-PIC static libz.a, so it links
# the shared libz.so instead - which this cmake build produces anyway alongside the
# static lib. Pack and ship it, too ("static-dynamic"), so the closed package stays
# self-contained. DEPS_ID gives a distinct precompiled cache name, otherwise an older
# static-only zip (without libz.so) would be reused.
if(CFS_BUILD_PRECICE AND UNIX)
  set(DEPS_ID "shared")
  set(DEPS_LIB_TYPE "static-dynamic")
endif()

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
# on macOS and Linux it is libz.a, on Windows zlibstatic.lib
# originally the windows debug libs have a tailing "d" but we patch it away because of boost and it's annoying bjam :(
if(UNIX)
  set_package_library_list("z")
else()
  set_package_library_list("zlibstatic")
endif()

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# we don't need the share/man stuff and the dynamic lib
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 
# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  # CMAKE_INSTALL_PREFIX is ignored!
  -DINSTALL_BIN_DIR:PATH=${DEPS_INSTALL}/bin
  -DINSTALL_INC_DIR:PATH=${DEPS_INSTALL}/include
  -DINSTALL_LIB_DIR:PATH=${DEPS_INSTALL}/lib
  -DINSTALL_MAN_DIR:PATH=${DEPS_INSTALL}/share/man
  -DINSTALL_PKGCONFIG_DIR:PATH=${DEPS_INSTALL}/share/pkconfig )  
  
# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# Generate ${PACKAGE_NAME}-patch.cmake we use for our external project
generate_patches_script() # sets PATCHES_SCRIPT

# the forked metis insists onb bulding executables to bin. Sort them out automatically
generate_packing_script_install_dir()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()
# if not, build newly and possibly pack the stuff
else()
# standard cmake build with patch
create_external_cmake_patched()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()  
    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()      
  endif() 
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})