# Xerces-C++ is our default XML reader / schema parser/validator
# https://xerces.apache.org/xerces-c/
#
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "xerces")

set(PACKAGE_VER "3.2.4")
set(PACKAGE_FILE "xerces-c-${PACKAGE_VER}.zip")
set(PACKAGE_MD5 "f84488fc2b8f62c4afca2f9943a42c00")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

set(PACKAGE_MIRRORS "https://archive.apache.org/dist/xerces/c/3/sources/${PACKAGE_FILE}"
                    "https://dlcdn.apache.org/xerces/c/3/sources/${PACKAGE_FILE}") 
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a fortran compiler
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

if(WIN32)
  set_package_library_list("xerces-c_3") # whyever?!
else()
  set_package_library_list("xerces-c") # add -c to package name
endif()

set_standard_variables() 

# we don't want bin and share
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project with standard compiler settings
set_deps_args_default(ON) 

# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -Dnetwork:BOOL=OFF 
  -Dthreads:BOOL=ON # OFF gives error that ICU libraries are not found?!
  -Dtranscoder:STRING=iconv)
  
# --- it follows generic final block for cmake packages with no patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

generate_patches_script() # sets PATCHES_SCRIPT

# generate package creation script. We want to skip bin and share and therefore use no manifest
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
