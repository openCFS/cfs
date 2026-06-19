# Library of Iterative Solvers
# https://www.ssisc.org/lis/index.en.html

clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "lis")
set(LIS_VER "2.1.11") # for Dependencies.cc which cannot easily inlcude lis.h (fucks up boost)
set(PACKAGE_VER ${LIS_VER})
set(PACKAGE_FILE "lis-${PACKAGE_VER}.zip")
set(PACKAGE_MD5 "ed5cf7cf764ed9574b2497f4dd2f68ab")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

if(USE_OPENMP)
  set(DEPS_ID "OPENMP")
else()
  set(DEPS_ID "NO_OPENMP")
endif()

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
if(UNIX)
  set(DEPS_INSTALL "${DEPS_PREFIX}/install")
else()
  set(DEPS_INSTALL "${DEPS_SOURCE}/win/my_install")
endif()

# these settings are for UNIX only. For configure.bat we have only a reduced set
set_configure_default()
if(USE_OPENMP) # don't combine with setting DEPS_ID - mixes up order of called macros.
  list(APPEND DEPS_CONFIGURE --enable-omp=yes)
  if(APPLE)
    assert_set(OpenMP_LIBDIR) # compiler.cmake, no system path on homebrew >= Oct 2022
    list(APPEND DEPS_CONFIGURE_ENV LDFLAGS=-L${OpenMP_LIBDIR} CPPFLAGS=-I${OpenMP_C_INCLUDE_DIR})
  endif()
else()
  list(APPEND DEPS_CONFIGURE --enable-omp=no)
endif()
# Up form version 2.1.11, algorithm for complex matrix is enabled
list(APPEND DEPS_CONFIGURE --enable-test=no --enable-fma=yes --enable-complex=yes --enable-saamg=no --enable-static --enable-shared=no )

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

if(WIN32)
  generate_patches_script() # for UNIX create_external_configure() asserts this not to be set
endif()

# generate package creation script. We do not get the files from an install_manifest.txt
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
  if(WIN32)
    assert_set(PATCHES_SCRIPT)
    # the standard DEPS_INSTALL has issues with Windows, stay in build dir. enabling complex breaks computations	    
    set(WIN_CONFIGURE --disable-test --prefix my_install ) 
    if(USE_OPENMP)
      list(APPEND WIN_CONFIGURE --enable-omp)
    endif()
    if(CMAKE_C_COMPILER_ID MATCHES "Intel")
      list(APPEND WIN_CONFIGURE --enable-intelc) # for IntelLLVM we need to patch Makefile.in icl -> icx
    endif()
    find_program(NMAKE_PROGRAM nmake)
    ExternalProject_Add("${PACKAGE_NAME}"
      PREFIX "${DEPS_PREFIX}"
      # Windows needs to have condigure.bat and successive nmake exececuted in <source>/win
      BINARY_DIR ${DEPS_SOURCE}/win
      URL "${PACKAGE_MIRRORS}"
      URL_MD5 "${PACKAGE_MD5}"
      DOWNLOAD_DIR "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}"
      DOWNLOAD_NAME "${PACKAGE_FILE}"
      DOWNLOAD_NO_PROGRESS ON 
      PATCH_COMMAND ${CMAKE_COMMAND} -P "${PATCHES_SCRIPT}"
      CONFIGURE_COMMAND ${DEPS_SOURCE}/win/configure.bat ${WIN_CONFIGURE}
      BUILD_COMMAND ${NMAKE_PROGRAM}
      INSTALL_COMMAND ${NMAKE_PROGRAM} install
      BUILD_BYPRODUCTS ${PACKAGE_LIBRARY} )
  else()
    # standard configure works for macOS and Linux
    create_external_configure()
  endif()

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

    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()  
  endif()
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
