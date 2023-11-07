# IPOPT is a general purpose open source optimizer 
# Tis is the precompiled code for Windows - makes no fun to do manually
# https://github.com/coin-or/Ipopt

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "ipopt")
set(PACKAGE_VER "3.14.12") # it makes sense to have this in sync with External_IPOPT.cmake
set(PACKAGE_FILE "Ipopt-${PACKAGE_VER}-win64-msvs2019-md.zip")
set(PACKAGE_MD5 "494932d606bd9eaab4443490cfc1eaf2")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.
# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://github.com/coin-or/Ipopt/releases/download/releases%2F${PACKAGE_VER}/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# we use precompiled
use_c_and_fortran(OFF OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
set_package_library_list("ipopt.dll") # the precompiled lib is ipopt.dll.lib 
# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# the unpacked "source" is the precompiled Ipopt - we pack and copy from there
set(DEPS_INSTALL "${DEPS_PREFIX}/src/${PACKAGE_NAME}")

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}")

assert_unset(PATCHES_SCRIPT)

# generate package creation script. We do not get the files from an install_manifest.txt
set(DEPS_LIB_TYPE "static-dynamic") # dynamic is from bin for WIN32
generate_packing_script_install_dir()

# we have a litte to clean-up from the precompiled stuff
generate_postinstall_script()
#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()
# if not, build newly and possibly pack the stuff
else()
  # we need to build the extract the package
  ExternalProject_Add(${PACKAGE_NAME}
    PREFIX ${DEPS_PREFIX}
    BINARY_DIR ${DEPS_SOURCE}
    URL ${PACKAGE_MIRRORS}
    URL_MD5 ${PACKAGE_MD5}
    # DOWNLOAD_DIR is ignored, if URL contains not the file mirror
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}
    # in case the mirrors have different file names we always store to the same
    DOWNLOAD_NAME ${PACKAGE_FILE}
    DOWNLOAD_NO_PROGRESS ON 
    PATCH_COMMAND ""
    CONFIGURE_COMMAND "" 
    BUILD_COMMAND  ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${PACKAGE_LIBRARY} )

  add_postinstall_step()

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()
    # copy what we need from the precompiled Ipopt  
    add_install_dir_to_binary_step()  
  endif()  
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})