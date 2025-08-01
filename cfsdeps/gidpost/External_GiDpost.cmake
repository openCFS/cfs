#-------------------------------------------------------------------------------
# GiDpost: library to create postprocess files for GiD
#
# Project Homepage
# http://gid.cimne.upc.es/gid-plus/tools/gidpost
#-------------------------------------------------------------------------------

clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "gidpost")
set(PACKAGE_VER "2.13") # note that we hide deprecated warnings in SimOutGiD.cc
set(PACKAGE_FILE "gidpost-${PACKAGE_VER}.zip")
set(PACKAGE_MD5 "3672d46422c1d9987e8985d40bcfa247")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names.

# Define the primary and fallback URLs
set(PACKAGE_MIRRORS "https://downloads.gidsimulation.com/Tools/gidpost/${PACKAGE_FILE}"
                    "https://downloads.gidsimulation.com/Tools/gidpost/old/${PACKAGE_FILE}")
                    
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # C/C++
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# libgidpost.a
set_package_library_default()

set_standard_variables() 

# we will use the manifest.txt
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/gidpost/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/gidpost" )

# set DEPS_ARG with defaults for a cmake project and compiler flags
set_deps_args_default(ON)

set(DEPS_ARGS
  ${DEPS_ARGS}
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
  -DZLIB_INCLUDE_DIR:FILEPATH=${ZLIB_INCLUDE_DIR}
  -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
  -DENABLE_FORTRAN_EXAMPLES=OFF
  -DENABLE_HDF5=OFF
  -DENABLE_PARALLEL_EXAMPLE=OFF
  -DENABLE_SHARED_LIBS=OFF )

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# no patch needed
assert_unset(PATCHES_SCRIPT)

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# we build directory to CMAKE_BINARY_DIR, there is no unnecessary stuff in manifest.txt
generate_packing_script_manifest()

assert_unset(POSTINSTALL_SCRIPT)

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # patched cmake project    
  create_external_cmake()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  endif()
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})

# we don't use the hdf5 dependency - use it, if you want to use binary gidpost
add_dependencies(${PACKAGE_NAME} zlib)
