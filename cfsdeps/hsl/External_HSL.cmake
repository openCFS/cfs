# ThirdPary-HSL is a helper for the Ipopt optimizer if IPOPT has no MKL (e.g. ARM)
# It requires the solver sourcode Coin-HSL from https://licences.stfc.ac.uk/product/coin-hsl
# the result is the dynamic hsl libraray (libshl.dylib) loaded by Iptop at runtime
# https://github.com/coin-or-tools/ThirdParty-HSL

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "hsl")
set(PACKAGE_VER "2.2.2")
set(PACKAGE_FILE "${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "9edd0d6e9ab57996d90886e743bbc8fa")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

set(COINHSL_VER "2023.05.26") # this is the coinhsl version

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://codeload.github.com/coin-or-tools/ThirdParty-HSL/tar.gz/refs/tags/releases/${PACKAGE_VER}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# HSL is almost pure Fortran 
use_c_and_fortran(ON ON)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
set_package_library_default()
# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# this is the standard target for configure projects (builds in source). This directory will be zipped
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

set_configure_default()
if(False AND USE_METIS)
  string(CONCAT LAPACK_STR "-L${CMAKE_BINARY_DIR}/lib -lmetis") # handle quote stuff
  set(DEPS_CONFIGURE ${DEPS_CONFIGURE} --with-metis-lflags=${LAPACK_STR})
endif()

# --- it follows generic final block for cmake packages with a patch and no postinstall ---
# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# copy local original coinhsl source into ThirdParty-HSL following it' README.md 
generate_patches_script()

# generate package creation script.
set(DEPS_LIB_TYPE "dynamic")
generate_packing_script_install_dir()

# follow again the ThirdParty-HSL README.md. We need a rename 
generate_postinstall_script()

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

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})