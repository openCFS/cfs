# SNOPT (sparse nonlinear optimizer) is a commercial general purpose optimizer
# which can be used efficiently for structural optimization.
#
# The chair for continuous optimization (FAU) has a license for the code.
#
# You need to request the encrypted file snopt-7.2.8-cfsdeps.zip with the source 
# from the optimization group. 
# Place snopt-7.2.8-cfsdeps.zip to the cfsdepscache source/snopt directory
# and set the password for the file to CFS_KEY_SNOPT (e.g. in .cfs_platform_defaults.cmake)
# alternatively with CFS_DOWNLOAD_SNOPT the encrypted file can be downloaded
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "snopt")
set(SNOPT_VER "7.2.8") # for def_use_snopt.hh.in
set(PACKAGE_VER ${SNOPT_VER})
set(PACKAGE_FILE "snopt-${PACKAGE_VER}-cfsdeps.zip")
set(PACKAGE_MD5 "9e75be8400eb878b9cb3d489084af196")
set(DEPS_VER "-b") # CFSDEPS_Fortran_FLAGS 

# with CFS_DOWNLOAD_SNOPT the encrypted file can be found
# when we need to download, the key needs to be set from local .cfs_platform_defaults.cmake, environment or cmake -DCFS_DOWNLOAD_SNOPT=...
set_from_env(CFS_DOWNLOAD_SNOPT)
set(PACKAGE_MIRRORS "${CFS_FAU_MIRROR}/sources/${CFS_DOWNLOAD_SNOPT}/${PACKAGE_FILE}") 
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
set(DEPS_INSTALL "${DEPS_PREFIX}/install")
# not standard!
set(DEPS_SOURCE  "${DEPS_PREFIX}/src/${PACKAGE_NAME}/snopt7")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 
# add the specific settings for the packge which comes in cmake style
assert_set(LIB_SUFFIX)
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DLIB_SUFFIX=${LIB_SUFFIX})

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/"  DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# generate package creation script. Somehow the install_manifest.txt fails for snopt. It is not unacked to lib. Possibly EOL issue?!
generate_packing_script_install_dir()

# copy CMakeLists.txt
generate_patches_script()

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # variables can be set in .cfs_platform_defaults.cmake, with cmake -D or as environment variables.
  if(NOT EXISTS LOCAL_PACKAGE_FILE AND NOT DEFINED CFS_DOWNLOAD_SNOPT)
    message(WARNING "no local encypted snopt sources exists and CFS_DOWNLOAD_SNOPT not set.")
  endif()

  set_from_env(CFS_KEY_SNOPT)
  if(NOT CFS_KEY_SNOPT)
    message(FATAL_ERROR "to build the commerical snopt7 you need CFS_KEY_SNOPT to be set or have the original code and remove decryption")
  endif()
  set(PACKAGE_KEY ${CFS_KEY_SNOPT})

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
