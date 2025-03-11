#-------------------------------------------------------------------------------
# CFD General Notation System (CGNS)
# Needed for ADF routines by STARCCM+ reader. 
# To build cgnsview which can be used to view .cgns files (HDF5 and ADF) and .ccm files (ADF)
# install it manually or reenable BUILD_CGNSTOOL.
#
# Project Homepage
# http://www.cgns.org
#-------------------------------------------------------------------------------

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(CFS_VERSION_MONTH "04")
set(PACKAGE_NAME "cgns")
set(PACKAGE_VER "4.4.0")
set(CGNS_VER ${PACKAGE_VER}) # for Dependencies.cc
set(PACKAGE_FILE "v${PACKAGE_VER}.zip")
set(PACKAGE_MD5 "9779e1478c4be7c89f096a22bc0c0bef")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://github.com/CGNS/CGNS/archive/refs/tags/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# we only have a c++ compiler
use_c_and_fortran(ON OFF)

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

# add the specific settings for the packge which comes in cmake style...
# CMAKE_INSTALL_LIBDIR is ignored (4.3.0), hence we need to patch
set(DEPS_ARGS ${DEPS_ARGS}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
  -DCGNS_BUILD_CGNSTOOLS:BOOL=OFF # handles only cgnstools, we remove tools by patch
  -DCGNS_BUILD_SHARED:BOOL=OFF
  -DCGNS_BUILD_TESTING:BOOL=OFF
  -DCGNS_ENABLE_64BIT:BOOL=OFF
  -DCGNS_ENABLE_BASE_SCOPE:BOOL=OFF
  -DCGNS_ENABLE_FORTRAN:BOOL=OFF
  -DCGNS_ENABLE_HDF5:BOOL=ON
  -DCGNS_ENABLE_TESTS:BOOL=OFF
  -DCGNS_USE_SHARED:BOOL=OFF
  # we have no hdf5 lib and include dir but need hdf5-config.cmake
  -DHDF5_DIR:FILEPATH=${CMAKE_CURRENT_BINARY_DIR}/share/cmake # for hdf5-config.cmake
)
# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# Generate ${PACKAGE_NAME}-patch.cmake we use for our external project
generate_patches_script() # sets PATCHES_SCRIPT

# generate package creation script. We get the files from an install_manifest.txt
generate_packing_script_manifest()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

# dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # add external project step actually building an cmake package including a patch 
  # also genearate the patch script via generate_patches_script()
  create_external_cmake_patched()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  endif()  
endif()

# cgns depends on hdf5
add_dependencies(cgns hdf5)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})