# CGAL: The Computational Geometry Algorithms Library
# https://www.cgal.org
# CGAL consists of different (open source) licenes, most are GPL, therefore we have no CGAL in the default build
# CGAL is header only and requires officially for CGAL >= 6 C++17 and gcc >= 11.4 (but it works with gcc 9)  
# we currently provide gmp and mpfr to prevent cfs link errors, but we don't provide them when configuring CGAL
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "cgal")
set(PACKAGE_VER "6.0.1")
set(PACKAGE_FILE "CGAL-${PACKAGE_VER}-library.tar.xz")
set(PACKAGE_MD5 "ea827f6778063e00554ae41f4c845492")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://github.com/CGAL/cgal/releases/download/v${PACKAGE_VER}/${PACKAGE_FILE}") 
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we do not compile as eigen is header only
use_c_and_fortran(OFF OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# cmake package building directly to cfs build. precompiled package via manifest
set(DEPS_PREFIX  "${CMAKE_BINARY_DIR}/cfsdeps/${PACKAGE_NAME}")
set(DEPS_SOURCE  "${DEPS_PREFIX}/src/${PACKAGE_NAME}")

# the clean-<package> target deletes everything to allow a clean make <package>
add_clean_target()

# Don't use use install_manifest.txt
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project 
set_deps_args_default(ON)
# CGAL 6.0.1 should work without gmp/mpfr, but currently produces link errors without gmp
# see https://github.com/CGAL/cgal/issues/8567 and https://github.com/CGAL/cgal/issues/8606
# see also -DCMAKE_OVERRIDDEN_DEFAULT_ENT_BACKEND=BOOST_BACKEND in compiler.cmake
# currently the following options don't do their job: CGAL_DISABLE_GMP=ON, CGAL_CMAKE_EXACT_NT_BACKEND=BOOST_BACKEND
# therefore we build gmp and mpfr on macOS and Linux.  

# --- it follows generic final block for cmake packages with no patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

# filter from DEPS_INSTALL, zip and copy to binary dir 
generate_packing_script_install_dir()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

# dump_depencency_variables()

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

if(UNIX)
  # we would like to use CGAL with boost only, but currently this doesn't work on Linux/macOS, but on Windows
  add_dependencies(cgal mpfr)
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
