# IPOPT is a general purpose open source optimizer 
# https://coin-or.github.io/Ipopt/
# https://github.com/coin-or/Ipopt

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

if(WIN32)
  message(FATAL_ERROR "Compiling IPOPT from source is currently not implemented for Windows.")
endif()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "ipopt")
set(PACKAGE_VER "3.14.10")
set(PACKAGE_FILE "${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "f22da4b75d9c936e607febe1f7e63815")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://github.com/coin-or/Ipopt/archive/refs/tags/releases/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we'll disable fortran for ipopt as it is not needed
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
set_package_library_default()
# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# this is the standard target for configure projects (builds in source). This directory will be zipped
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# org legacy for 3.14.2 mkl only   CONFIGURE_COMMAND env "CXXFLAGS=${CFS_CXX_FLAGS}" env "LIBRARY_PATH=${MKL_LIB_DIR}" ${IPOPT_SOURCE}/configure --with-lapack-lflags="-L${MKL_LIB_DIR} -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm"  --prefix=${IPOPT_INSTALL} --exec-prefix=${IPOPT_INSTALL} --libdir=${IPOPT_INSTALL}/lib64 --enable-shared=no --enable-static --disable-java --disable-sipopt F77=${CMAKE_Fortran_COMPILER} OPT_FFLAGSS=-O3 CXX=${CMAKE_CXX_COMPILER} OPT_CXXFLAGS=-O3
set_configure_default()
set(DEPS_CONFIGURE ${DEPS_CONFIGURE} --disable-f77 --enable-static=on --enable-shared=off --disable-java --disable-sipopt)    

# it is a little tricky to add quoted blocks to CONFIGURE_COMMAND .. --with-lapack="-L... -l..." as the " confuse cmake a lot :(
# therefore we need to create a string for every quoted block
if(USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  string(CONCAT LAPACK_STR "-L${CMAKE_BINARY_DIR}/lib -lopenblas")
elseif(USE_BLAS_LAPACK STREQUAL "MKL")
  # see https://coin-or.github.io/Ipopt/INSTALL.html howecver --with-lapack seems also to work instead of --with-lapack-lflags
  string(CONCAT LAPACK_STR "-L${MKL_LIB_DIR} -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm")
endif()

# for blas and Accelerate hope for the best, it is just to get over configure
if(NOT ${LAPACK_STR} STREQUAL "")
  set(DEPS_CONFIGURE ${DEPS_CONFIGURE} --with-lapack=${LAPACK_STR})
endif()

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/"
     DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

# generate package ceation script. We get the files from an install_manifest.txt
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
  endif()  
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})