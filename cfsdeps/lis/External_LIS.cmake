# Library of Iterative Solvers
# https://www.ssisc.org/lis/index.en.html

clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "lis")
set(LIS_VER "2.0.34") # for Dependencies.cc which cannot easily inlcude lis.h (fucks up boost)
set(PACKAGE_VER ${LIS_VER})
set(PACKAGE_FILE "lis-${PACKAGE_VER}.zip")
set(PACKAGE_MD5 "5a666ee5bd8af29d3d171771ead78a36")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://www.ssisc.org/lis/dl/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we'll disable fortran for lis by not using saamg which is fast, but very sensitive to system condition
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
set_package_library_default()
# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# this is the standard target for configure projects (builds in source). This directory will be zipped
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# lis has a configure.bat for Windows. There is also a cmake project for lis which seems to be outdated
set_configure_default()
set(DEPS_CONFIGURE ${DEPS_CONFIGURE} --enable-omp=yes --enable-test=no --enable-fma=yes --enable-complex=yes --enable-saamg=no --enable-static --enable-shared=no )

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/"
     DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

# generate package ceation script. 
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

  assert_set(DEPS_SOURCE)

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()

    # we need include/lis_precon.h but lis does not pack it.
    ExternalProject_Add_Step(${PACKAGE_NAME} pre_packing 
      COMMAND ${CMAKE_COMMAND} -E copy ${DEPS_SOURCE}/include/lis_precon.h ${DEPS_INSTALL}/include
      DEPENDEES install
      DEPENDERS cfsdeps_packaging )
  else()
    ExternalProject_Add_Step(${PACKAGE_NAME} copy_lis_precond
      COMMAND ${CMAKE_COMMAND} -E copy ${DEPS_SOURCE}/include/lis_precon.h ${CMAKE_BINARY_DIR}/include
      DEPENDEES install )
  endif()
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
