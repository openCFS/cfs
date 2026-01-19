# The MPFR library is a C library for multiple-precision floating-point computations with correct rounding. 
# https://www.mpfr.org
# It is based on GMP and used by CGAL but in theory CGAL >= 6 should work without, but with CGAL 6.0.1 this practically does not work yet

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "mpfr")
set(PACKAGE_VER "4.2.2")
set(PACKAGE_FILE "${PACKAGE_NAME}-${PACKAGE_VER}.tar.xz")
set(PACKAGE_MD5 "7c32c39b8b6e3ae85f25156228156061")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://www.mpfr.org/mpfr-current/${PACKAGE_FILE}")

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# MPFR is pure C 
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
set_package_library_default()
# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# this is the standard target for configure projects (builds in source). This directory will be zipped
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

set_configure_default()
list(APPEND DEPS_CONFIGURE "--with-gmp=${CMAKE_BINARY_DIR}")
list(APPEND DEPS_CONFIGURE "--enable-shared=no")

# --- it follows generic final block for cmake packages with a patch and no postinstall ---
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# unpatched
assert_unset(PATCHES_SCRIPT)

# filter from DEPS_INSTALL, zip and copy to binary dir 
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
  # add external project step actually building an cmake package including a patch 
  # also genearate the patch script via generate_patches_script()
  create_external_configure()  
  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()
    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()  
  endif()  
endif()

# we need gmp first
add_dependencies(mpfr gmp)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})