# SNOPT (sparse nonlinear optimizer) is a commercial general purpose optimizer
# which can be used efficiently for structural optimization, especially with many (linear) constraints.
#
# The chair for continuous optimization AM-KO (FAU) has a license for the code.
# The source code can be licensed for about $600 (academic) via http://www.sbsi-sol-optimize.com/asp/sol_snopt.htm
#
# 1) Place the original vendor zip snopt7.7.7.zip (unencrypted) into the cfsdepscache
#    source/snopt directory. No key is needed at all.
# 2) For AM-KO internal usage, set the cmake variables CFS_DOWNLOAD_SNOPT77 and CFS_KEY_SNOPT (e.g. in .cfs_platform_defaults.cmake).
#
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "snopt")
set(SNOPT_VER "7.7.7") # for def_use_snopt.hh.in
set(PACKAGE_VER ${SNOPT_VER})
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER

# the original vendor zip, manually placed in the cfsdepscache sources. Both zips extract to
# a top level snopt7/ directory, which DEPS_SOURCE relies on.
set(PLAIN_PACKAGE_FILE "snopt${SNOPT_VER}.zip")
if(EXISTS "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}/${PLAIN_PACKAGE_FILE}")
  set(PACKAGE_FILE "${PLAIN_PACKAGE_FILE}")
  set(PACKAGE_MD5 "3e2995e48d83051402204d972029b09b")
  # create_external_encrypted_cmake_patched() requires a non-empty key. Giving a password for an unencrypted archive is harmless.
  set(PACKAGE_KEY "none")
else()
  set(PACKAGE_FILE "snopt-${PACKAGE_VER}-cfsdeps.zip")
  set(PACKAGE_MD5 "41deaaebad8bd1b378b35d2836a3554e") 
  # when we need to download, the key needs to be set from local .cfs_platform_defaults.cmake, environment or cmake -DCFS_DOWNLOAD_SNOPT77=...
  set_from_env(CFS_DOWNLOAD_SNOPT77)
  if(NOT EXISTS "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}/${PACKAGE_FILE}" AND NOT DEFINED CFS_DOWNLOAD_SNOPT77)
    message(FATAL_ERROR "neither ${PLAIN_PACKAGE_FILE} nor encrypted ${PACKAGE_FILE} exist and CFS_DOWNLOAD_SNOPT77 not set.")
  endif()

  set_from_env(CFS_KEY_SNOPT)
  if(NOT CFS_KEY_SNOPT)
    message(FATAL_ERROR "CFS_KEY_SNOPT is missing.")
  endif()
  set(PACKAGE_KEY ${CFS_KEY_SNOPT})
endif()

# in case we have the original zip, the mirror is ignored.
set(PACKAGE_MIRRORS "${CFS_FAU_MIRROR}/sources/${CFS_DOWNLOAD_SNOPT77}/${PACKAGE_FILE}")
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

# copy CMakeLists.txt (we don't patch snopt any more, see snopt-patch.cmake.in)
generate_patches_script()

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # PACKAGE_FILE, PACKAGE_MD5 and PACKAGE_KEY have been determined at the top of this file
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
