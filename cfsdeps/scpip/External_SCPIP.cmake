# SCPIP is an efficient optimization for structural optimization problems. 
# it implements the Method of Moving Asymptotes (MMA) from Krister Svanberg
# http://www.mathematik.uni-wuerzburg.de/~zillober
#
# SCPIP is not open source! Users need to obtain a permission to use the
# code for academic purupose! 
# Set the password for the file to CFS_KEY_SCPIP (e.g. in .cfs_platform_defaults.cmake)
# alternatively with CFS_DOWNLOAD_SCPIP the encrypted file can be downloaded
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "scpip")
set(SCPIP_VER "2009.01") # for def_use_scpip.hh.in
set(PACKAGE_VER ${SCPIP_VER})
set(PACKAGE_FILE "scpip-cfsdeps.zip")
set(PACKAGE_MD5 "a2767521596ed925e53b9fea6df4d77e")
set(DEPS_VER "-c") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# with CFS_DOWNLOAD_SCPIP the encrypted file can be found
# when we need to download, the key needs to be set from local .cfs_platform_defaults.cmake, environment or cmake -DCFS_DOWNLOAD_SCPIP=...
set_from_env(CFS_DOWNLOAD_SCPIP)
set(PACKAGE_MIRRORS "${CFS_FAU_MIRROR}/sources/${CFS_DOWNLOAD_SCPIP}/${PACKAGE_FILE}") 
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a fortran compiler
use_c_and_fortran(OFF ON)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
set_package_library_default()
# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# this is the standard target for cmake projects. The files to package come from the install_manifest.txt
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 
# add the specific settings for the packge which comes in cmake style
assert_set(LIB_SUFFIX)
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DLIB_SUFFIX=${LIB_SUFFIX})

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# we need to fix the Fortran code
generate_patches_script()

# for scpip just the static lib
generate_packing_script_manifest()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # variables can be set in .cfs_platform_defaults.cmake, with cmake -D or as environment variables.
  if(NOT EXISTS LOCAL_PACKAGE_FILE AND NOT DEFINED CFS_DOWNLOAD_SCPIP)
    message(WARNING "no local encypted scpip sources exists and CFS_DOWNLOAD_SCPIP not set.")
  endif()

  set_from_env(CFS_KEY_SCPIP)
  if(NOT CFS_KEY_SCPIP)
    message(FATAL_ERROR "to build proprietary scpip you need CFS_KEY_SCPIP to be set or have the original code and remove decryption")
  endif()
  set(PACKAGE_KEY ${CFS_KEY_SCPIP})

  create_external_encrypted_cmake_patched()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()
    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()  
  endif()  
endif()

# copy the appropriate special license in the license folder
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/LICENSE.binary.CLOSED" DESTINATION "${CFS_BINARY_DIR}/license/")

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
