# collection of helpers for dependcies

# dump the dependecies variables
macro(dump_depencency_variables)
  # PACKAGE variables to be set/chcked by the user
  cmake_print_variables(PACKAGE_NAME)  # e.g. arpack
  cmake_print_variables(PACKAGE_VER)   # e.g. 3.7.0
  cmake_print_variables(PACKAGE_FILE)  # e.g. ${PACKAGE_VER}.tar.gz as github uses no ${PACKAGE_NAME}
  cmake_print_variables(PACKAGE_MD5)
  cmake_print_variables(PACKAGE_LIBRARY) # e.g. libarpack.a, see set_package_library_default()
  cmake_print_variables(PACKAGE_MIRRORS) # e.g. the original download, see add_standard_mirrors()
  cmake_print_variables(PACKAGE_KEY)   # for encrypted closed source
  # DEPS variables for building the package. Most important is DEPS_ARGS for cmake packages
  cmake_print_variables(DEPS_ARGS)    # for cmake projects, see set_deps_args_default(ON) # set compiler flags, add compiler and specific settings     
  cmake_print_variables(DEPS_CONFIGURE) # for configure projects, see set_configure_default()
  cmake_print_variables(DEPS_VER)     # either "" or "-a", ... to trigger new precompiles when PACKAGE_VER does not change
  cmake_print_variables(DEPS_ID)      # optional package id like "openmp" or "no_openmp"
  cmake_print_variables(DEPS_PREFIX)  # usuallay "${CMAKE_BINARY_DIR}/cfsdeps/${PACKAGE_NAME}"
  cmake_print_variables(DEPS_SOURCE)  # usuallay "${DEPS_PREFIX}/src/${PACKAGE_NAME}"
  cmake_print_variables(DEPS_SUBDIR)  # set when CMakeLists.txt is not directly in DEPS_SOURCE
  cmake_print_variables(DEPS_INSTALL) # with install_manifest.txt directly ${CMAKE_BINARY_DIR} otherwise ${DEPS_PREFIX}/install
  cmake_print_variables(DEPS_LIB_TYPE)# default "static", otherwise "dynamic" or "static-dynamic"
  # kind of internal "class attributes" of DependencyTools. Check with DepsPackaging*.cmake.in
  cmake_print_variables(USE_C_CXX)    # use either C or C++ compiler. See use_c_and_fortran()
  cmake_print_variables(USE_FORTRAN)
  cmake_print_variables(FORCE_C_CXX)  # enforce the C/C++ label to overwrite the cfs compiler. See force_c_cxx_compiler()
  cmake_print_variables(PRECOMPILED_PCKG_FILE) # full patch to precompiled, e.g. /Users/fwein/code/cfsdepscache/precompiled/arpack_3.7.0_MACOSX_12.6_ARM64_F-GNU-11.2.0.zip
  cmake_print_variables(DEPS_BUILD_MESSAGE)    # message to be printed before building
  cmake_print_variables(LOCAL_PACKAGE_FILE)    # full path to local source 
  cmake_print_variables(PATCHES_SCRIPT)
  cmake_print_variables(POSTINSTALL_SCRIPT)
  cmake_print_variables(ZIPTOCACHE_SCRIPT)
endmacro()

# reset dependencies variables such that we less likely use leftovers from last dependency
#
# resetting these variables makes the assert_set() checks more useful.
# Call as early as possible
macro(clear_depencency_variables)
  unset(PACKAGE_NAME)
  unset(PACKAGE_VER)
  unset(PACKAGE_FILE)
  unset(PACKAGE_MD5)
  unset(PACKAGE_LIBRARY)
  unset(PACKAGE_MIRRORS)
  unset(PACKAGE_KEY)
  unset(DEPS_ARGS)
  unset(DEPS_CONFIGURE)
  unset(DEPS_VER)
  unset(DEPS_ID)
  unset(DEPS_PREFIX)
  unset(DEPS_SOURCE)
  unset(DEPS_SUBDIR)
  unset(DEPS_INSTALL)
  set(DEPS_LIB_TYPE "static") # default
  unset(USE_C_CXX)
  unset(USE_FORTRAN)
  unset(FORCE_C_CXX)
  unset(PRECOMPILED_PCKG_FILE)
  unset(DEPS_BUILD_MESSAGE)
  unset(LOCAL_PACKAGE_FILE)
  unset(PATCHES_SCRIPT)
  unset(POSTINSTALL_SCRIPT)
  unset(ZIPTOCACHE_SCRIPT)
endmacro()

# indicate if C/C++ and/or Fortran is used. 
# Triggers package name, DEPS_ARGS compiler settings 
# and if CACHE variable*_INLUDE_DIR is generated
macro(use_c_and_fortran IN_USE_C_CXX IN_USE_FORTRAN)
  set(USE_C_CXX ${IN_USE_C_CXX})
  set(USE_FORTRAN ${IN_USE_FORTRAN})
endmacro()

# variant of use_c_and_fortran where we enforce the given compiler name for C/C++ and don't use
# CFSDEPS_C(XX)_FLAGS. We set USE_C_CXX to ON and USE_FORTRAN to OFF.
# Usage is Win32 to use CL for boost and xerces for an intel cfs and vice versa for ginkgo
macro(force_c_cxx_compiler IN_FORCE_C_CXX)
  assert_unset(USE_C_CXX)
  assert_unset(USE_FORTRAN)
  
  set(FORCE_C_CXX ${IN_FORCE_C_CXX})
  set(USE_C_CXX ON)
  set(USE_FORTRAN OFF) # simply not in the use case of currently boost, xerces, ginkgo
endmacro()

# set os sensitive static lib to given variable in cmake cache
macro(set_static_cache_lib TARGET_LIB IN_NAME IN_LABEL)
  set(${TARGET_LIB} ${${TARGET_LIB}} ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}${IN_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE FILEPATH ${IN_LABEL})
  mark_as_advanced(${TARGET_LIB})
endmacro()

# determine common PACKAGE_LIBRARY content for a standard cmake package.
#
# can only be used for the simple standard case.
# sets PACKAGE_LIBRARY. Note that  ${PACKAGE_NAME}_LIBRARY is set in set_standard_variables()
macro(set_package_library_default)
  assert_set(PACKAGE_NAME)
  assert_set(LIB_SUFFIX)
  assert_unset(PACKAGE_LIBRARY)
  # CMAKE_STATIC_LIBRARY_PREFIX is lib for Linux and macOS and empty for Windows
  set(PACKAGE_LIBRARY ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}${PACKAGE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
endmacro()

# set package library for a given list
macro(set_package_library_list IN_LIST)
  assert_set(PACKAGE_NAME)
  assert_set(LIB_SUFFIX)
  assert_unset(PACKAGE_LIBRARY)
  # CMAKE_STATIC_LIBRARY_PREFIX is lib for Linux and macOS and empty for Windows
  # CMAKE_STATIC_LIBRARY_SUFFIX is .a for Linux and macOS and .lib for Windows
  foreach(ITEM ${IN_LIST})
    list(APPEND PACKAGE_LIBRARY ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}${ITEM}${CMAKE_STATIC_LIBRARY_SUFFIX})
  endforeach()   
endmacro()

# set package library for a given list always with lib as prefix, which is not standard for Windows
macro(set_package_library_list_lib_prefix IN_LIST)
  assert_set(PACKAGE_NAME)
  assert_set(LIB_SUFFIX)
  assert_unset(PACKAGE_LIBRARY)
  
  foreach(ITEM ${IN_LIST})
    list(APPEND PACKAGE_LIBRARY ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/lib${ITEM}${CMAKE_STATIC_LIBRARY_SUFFIX})
  endforeach()   
endmacro()


# set standard variables. Kind of late constructor. Not for header-only deps
#
# also set some standard (hidden) CACHE variables with uppercase package name.
# in rare cases <package>_INCLUDE_DIR is not CMAKE_BINARY_DIR/include. Overwrite the setting
# in this case manually
macro(set_standard_variables)
  assert_set(PACKAGE_NAME)
  assert_set(PACKAGE_LIBRARY)
  assert_set(PRECOMPILED_PCKG_FILE)
  
  string(TOUPPER ${PACKAGE_NAME} _UPPER_PACKAGE_NAME) # arpack -> ARPACK

  # as PACKAGE_LIBRARY can be created even manually before (e.g. a list as in boost),
  # this settings shall always work. We have no header only dependency.
  set(${_UPPER_PACKAGE_NAME}_LIBRARY "${PACKAGE_LIBRARY}" CACHE FILEPATH "${PACKAGE_NAME} library.")
  mark_as_advanced(${_UPPER_PACKAGE_NAME}_LIBRARY)

  # almost all packages go to CMAKE_BINARY_DIR/include. Even fortran packages like arpack have includes
  set(${_UPPER_PACKAGE_NAME}_INCLUDE_DIR "${CMAKE_BINARY_DIR}/include" CACHE FILEPATH "${PACKAGE_NAME} include dir.")
  mark_as_advanced(${_UPPER_PACKAGE_NAME}_INCLUDE_DIR)
 
  # cmake package building directly to cfs build. precompiled package via manifest
  set(DEPS_PREFIX  "${CMAKE_BINARY_DIR}/cfsdeps/${PACKAGE_NAME}")
  set(DEPS_SOURCE  "${DEPS_PREFIX}/src/${PACKAGE_NAME}")

  # the clean-<package> target deletes the precompiled package, local cfsdeps and the lib.
  # to allow a clean make <package>, however, it does not remove 
  # the cached <package>_LIBRARY and <package>_INCLUDE_DIR which cannot be overwritten. 
  # For this use overwrite-<package>  
  add_custom_target(clean-${PACKAGE_NAME} cmake -E remove_directory ${DEPS_PREFIX}
     COMMAND cmake -E remove ${PRECOMPILED_PCKG_FILE}
     COMMAND cmake -E remove ${PACKAGE_LIBRARY}
     COMMENT "delete cfsdeps/${PACKAGE_NAME}, the lib and precompiled")
     
  # deletes <package>_LIBRARY and <package>_INCLUDE_DIR, 
  # note that this will directly trigger setting it with new values as it calls cmake!
  # you might usually also want to call the clean-<package> target   
  add_custom_target(overwrite-${PACKAGE_NAME} cmake . -U ${_UPPER_PACKAGE_NAME}_LIBRARY -U ${_UPPER_PACKAGE_NAME}_INCLUDE_DIR
     COMMENT "delete cached ${PACKAGE_NAME}_LIBRARY and _INCLUDE_DIR")
endmacro()

# adds clean-<package> and target for head-only deps. See set_standard_variables() 
macro(add_clean_target)
  assert_set(PACKAGE_NAME)

  # see set_standard_variables()
  add_custom_target(clean-${PACKAGE_NAME} cmake -E remove_directory ${DEPS_PREFIX}
     COMMAND cmake -E remove ${PRECOMPILED_PCKG_FILE}
     COMMENT "delete cfsdeps/${PACKAGE_NAME} and precompiled")
endmacro()


# set variable DEPS_ARGS with default settings cmake dependencies.
#
# The idea is to not set settings, the package complains about not known.
macro(set_deps_args_default SET_COMPILER_FLAGS)

  assert_unset(DEPS_ARGS) # shall be cleared by clear_depencency_variables()
  assert_set(DEPS_INSTALL)

  set(DEPS_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${DEPS_INSTALL}
    # some packages complain about non-used, but e.g. arpack has lib64 on centos6
    -DCMAKE_INSTALL_LIBDIR:PATH=${DEPS_INSTALL}/${LIB_SUFFIX}
    -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
    # as we usually not debug from cfs into a lib, we build everything as Release
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB} )
   
  if(USE_FORTRAN)
    assert_set(CMAKE_Fortran_COMPILER)
    list(APPEND DEPS_ARGS -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER})
    if(${SET_COMPILER_FLAGS} AND CFSDEPS_Fortran_FLAGS) # not set for ifort
      list(APPEND DEPS_ARGS -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS} )
    endif()
  endif() 

  if(USE_C_CXX AND NOT FORCE_C_CXX)
    list(APPEND DEPS_ARGS -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER})
    if(${SET_COMPILER_FLAGS} AND CFSDEPS_C_FLAGS)
      list(APPEND DEPS_ARGS -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS} )     
    endif()
    list(APPEND DEPS_ARGS -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER} )
    if(${SET_COMPILER_FLAGS} AND CFSDEPS_CXX_FLAGS) 
      list(APPEND DEPS_ARGS -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS})
    endif()
  endif()

  # toolchain is used for crosscompiling. If set, also the (cmake) dependecies are crosscompiled  
  if(CMAKE_TOOLCHAIN_FILE)
    list(APPEND DEPS_ARGS -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE} )
  endif()
endmacro()

# set standard settings for configure builds CONFIGURE_ARGS
macro(set_configure_default)
  assert_set(DEPS_INSTALL)
  assert_unset(DEPS_CONFIGURE)

  # configure of ipopt and lis fail with CC=${CMAKE_C_COMPILER} but without path it works
  # so /Library/Developer/CommandLineTools/usr/bin/cc fails but cc with same path works
  if(APPLE)
    # overwrite the compilers by removing it's path (for ipopt 3.14.10)
    get_filename_component(_CC ${CMAKE_C_COMPILER} NAME)
    get_filename_component(_CXX ${CMAKE_CXX_COMPILER} NAME)
    get_filename_component(_FC ${CMAKE_Fortran_COMPILER} NAME)
  else()
    set(_CC ${CMAKE_C_COMPILER})
    set(_CXX ${CMAKE_CXX_COMPILER})
    set(_FC ${CMAKE_Fortran_COMPILER})
  endif()
  set(DEPS_CONFIGURE_ENV CC=${_CC} CXX=${_CXX} FC=${_FC}) 
  set(DEPS_CONFIGURE --prefix=${DEPS_INSTALL} --exec-prefix=${DEPS_INSTALL} --libdir=${DEPS_INSTALL}/lib)
endmacro()

# add standard mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE
#
# LOCAL_PACKAGE_FILE is where full path of the package file in cfsdepscache sources.
# cmake external projects can download from a list of mirrors, but when there is a file
# (our local chached LOCAL_PACKAGE_FILE) only a single entry is allowd for PACKAGE_MIRRORS
macro(add_standard_mirrors_or_set_local)
  assert_set(PACKAGE_MIRRORS)
  assert_set(CFS_DS_SOURCES_DIR)
  assert_set(CFS_DEPS_CACHE_DIR)
  assert_set(PACKAGE_NAME)
  assert_set(PACKAGE_FILE)

  set(LOCAL_PACKAGE_FILE "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}/${PACKAGE_FILE}")
  
  if(EXISTS ${LOCAL_PACKAGE_FILE})
    set(PACKAGE_MIRRORS "file://${LOCAL_PACKAGE_FILE}")
  else()  
    set(PACKAGE_MIRRORS 
      ${PACKAGE_MIRRORS}
      "${CFS_FAU_MIRROR}/sources/${PACKAGE_NAME}/${PACKAGE_FILE}" )
  endif()   
endmacro()

# generate optional script for external project postinstall. Sets POSTINSTALL_SCRIPT
#
# can be used for cleanup before zipping 
macro(generate_postinstall_script)
  
  assert_set(PACKAGE_NAME)
  assert_set(PACKAGE_FILE)
  assert_set(DEPS_PREFIX)
  
  set(POSTINSTALL_SCRIPT "${DEPS_PREFIX}/${PACKAGE_NAME}-post_install.cmake")
  configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/${PACKAGE_NAME}-post_install.cmake.in"
                 "${POSTINSTALL_SCRIPT}" @ONLY)
endmacro()

# generate ${PACKAGE_NAME}-patch.cmake script. Don't call when you don't want to patch

# sets PATCHES_SCRIPT
macro(generate_patches_script)
 
  assert_set(DEPS_PREFIX)
  assert_set(PACKAGE_NAME)

  set(PATCHES_SCRIPT "${DEPS_PREFIX}/${PACKAGE_NAME}-patch.cmake")
  configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/${PACKAGE_NAME}-patch.cmake.in" 
               "${PATCHES_SCRIPT}" @ONLY)
endmacro()  

# Generate PRECOMPILED_PCKG_FILE with path and name and set DEPS_BUILD_MESSAGE
#
# makes use of USE_C_CXX and USE_FORTRAN from use_c_and_fortran()
macro(set_precompiled_pckg_file)
  assert_set(PACKAGE_NAME)
  assert_set(PACKAGE_VER)
  assert_set(CFS_ARCH_STR)
  # DEPS_VER might be "" which is unset
  assert_unset(PRECOMPILED_PCKG_FILE)
  
  if(FORCE_C_CXX)
    set(_CXX_ID_VER "${FORCE_C_CXX}")
  else()
    set(_CXX_ID_VER "${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")
  endif()
  
  
  # first set variables
  if(USE_FORTRAN)
    if(DEFINED CMAKE_Fortran_COMPILER_VERSION AND NOT "${CMAKE_Fortran_COMPILER_VERSION}" STREQUAL "")
      set(_FORTRAN_COMPILER_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")
    else()
      set(_FORTRAN_COMPILER_VERSION "?")
    endif()
  endif()

  # start with the base
  set(_TMP "${CFS_DEPS_CACHE_DIR}/precompiled/${PACKAGE_NAME}_${PACKAGE_VER}${DEPS_VER}_${CFS_ARCH_STR}")

  if(DEPS_ID)
    set(_TMP "${_TMP}_${DEPS_ID}")
  endif()
    
  if(NOT USE_C_CXX AND USE_FORTRAN)
    set(_TMP "${_TMP}_F-${CMAKE_Fortran_COMPILER_ID}-${_FORTRAN_COMPILER_VERSION}")
  elseif(USE_C_CXX AND NOT USE_FORTRAN)
    set(_TMP "${_TMP}_C-${_CXX_ID_VER}")
  elseif(USE_C_CXX AND USE_FORTRAN)
    assert_set(USE_C_CXX)
    assert_set(USE_FORTRAN)
    # combine C_F if same version
    if("${CMAKE_CXX_COMPILER_VERSION}" STREQUAL "${_FORTRAN_COMPILER_VERSION}")
      assert_unset(FORCE_C_CXX)
      set(_TMP "${_TMP}_C_F-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")
    else()
      set(_TMP "${_TMP}_C-${_CXX_ID_VER}")
      set(_TMP "${_TMP}_F-${CMAKE_Fortran_COMPILER_ID}-${_FORTRAN_COMPILER_VERSION}")
    endif()
  endif()
  
  set(PRECOMPILED_PCKG_FILE "${_TMP}.zip")

  # print a message in case we have precompiled enabled but could not find the package 
  if(${CFS_DEPS_PRECOMPILED} AND NOT EXISTS "${PRECOMPILED_PCKG_FILE}")
    message(STATUS "will build ${PACKAGE_NAME}_${PACKAGE_VER}: cannot find ${PRECOMPILED_PCKG_FILE}")
  endif()
endmacro()


# generate packing script for install_manifest.txt cmake package.
# the Manifinstall_manifest.txt is generated by the package cmake install step
# alternatively use generate_external_unpack_precompiled()
macro(generate_packing_script_manifest)
  assert_set(DEPS_PREFIX)
  
  set(ZIPTOCACHE_SCRIPT "${DEPS_PREFIX}/${PACKAGE_NAME}-DepsPackagingManifest.cmake")
  configure_file("${CMAKE_SOURCE_DIR}/cmake_modules/DepsPackagingManifest.cmake.in"
                 "${ZIPTOCACHE_SCRIPT}" @ONLY)  
endmacro()

# generate packing script for given install_dir
# alternatively use generate_packing_script_manifest()
macro(generate_packing_script_install_dir)
  assert_set(DEPS_INSTALL)
  assert_set(DEPS_PREFIX)
  assert_set(DEPS_LIB_TYPE)
  # with DEPS_INSTALL == CMAKE_BINARY_DIR we would pack the whole lib and include for all current packages
  if(${DEPS_INSTALL} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "either DEPS_INSTALL is wrong or you want to call generate_packing_script_manifest")
  endif()
  
  assert_not(${DEPS_INSTALL} ${CMAKE_BINARY_DIR}) # would pack all and probably install_manifest.txt is meant
  set(ZIPTOCACHE_SCRIPT "${DEPS_PREFIX}/${PACKAGE_NAME}-DepsPackagingInstallDir.cmake")
  configure_file("${CMAKE_SOURCE_DIR}/cmake_modules/DepsPackagingInstallDir.cmake.in"
                 "${ZIPTOCACHE_SCRIPT}" @ONLY)  
endmacro()

# generate ExternalProject_Add() for the unpacking from precompiled case
macro(create_external_unpack_precompiled)

  assert_set(PACKAGE_NAME)
  assert_set(PRECOMPILED_PCKG_FILE)
  assert_set(DEPS_PREFIX)

  ExternalProject_Add(${PACKAGE_NAME}
    PREFIX ${DEPS_PREFIX}
    DOWNLOAD_COMMAND ""
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
   
    BUILD_BYPRODUCTS ${PACKAGE_LIBRARY} )
 
  # we need a step as we need to set WORKING_DIRECTORY
  ExternalProject_Add_Step(${PACKAGE_NAME} pre_download
    COMMAND ${CMAKE_COMMAND} -E echo "unpacking ${PRECOMPILED_PCKG_FILE}"
    COMMAND ${CMAKE_COMMAND} -E tar xzf ${PRECOMPILED_PCKG_FILE}
    DEPENDERS download 
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
endmacro()    

# add cmake external project for building a patched cmake project 
#
# be sure to have a working ${PACKAGE_NAME}-patch.cmake applying the specific patches
# be sure to have generate_patches_script() called
macro(create_external_cmake_patched)

  assert_set(PACKAGE_NAME)
  assert_set(DEPS_ARGS)
  assert_set(DEPS_PREFIX)
  assert_set(PATCHES_SCRIPT)
  assert_set(PACKAGE_LIBRARY)

  # URL can take a list of mirrors but when there is a file, it needs to be the only one.
  # file means, that we have the source already in the cfsdeps cache. If not, we store there
  # see add_standard_mirrors()

  # we need to build the package - here in cmake style
  ExternalProject_Add("${PACKAGE_NAME}"
    PREFIX "${DEPS_PREFIX}"
    URL "${PACKAGE_MIRRORS}"
    URL_MD5 "${PACKAGE_MD5}"
    # DOWNLOAD_DIR is ignored, if URL contains not the file mirror
    DOWNLOAD_DIR "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}"
    # in case the mirrors have different file names we always store to the same
    DOWNLOAD_NAME "${PACKAGE_FILE}"
    DOWNLOAD_NO_PROGRESS ON
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PATCHES_SCRIPT}"
    CMAKE_ARGS ${DEPS_ARGS}
    BUILD_BYPRODUCTS ${PACKAGE_LIBRARY} )

  add_postinstall_step() # only if POSTINSTALL_SCRIPT is set
endmacro()



# set DEPS_SUBDIR to set SOURCE_SUBDIR argument of ExternalProject_Add()
macro(create_external_cmake)

  assert_set(PACKAGE_NAME)
  assert_set(DEPS_ARGS)
  assert_set(DEPS_PREFIX)
  assert_unset(PATCHES_SCRIPT)
  
  if(DEPS_SUBDIR)
    set(_SOURCE_SUBDIR "${DEPS_SUBDIR}")
  else()
    set(_SOURCE_SUBDIR ".")
  endif()    

  # URL can take a list of mirrors but when there is a file, it needs to be the only one.
  # file means, that we have the source already in the cfsdeps cache. If not, we store there
  # see add_standard_mirrors()
  
  # we need to build the package - here in cmake style
  ExternalProject_Add("${PACKAGE_NAME}"
    PREFIX "${DEPS_PREFIX}"
    SOURCE_SUBDIR "${_SOURCE_SUBDIR}"
    URL "${PACKAGE_MIRRORS}"
    URL_MD5 "${PACKAGE_MD5}"
    # DOWNLOAD_DIR is ignored, if URL contains not the file mirror
    DOWNLOAD_DIR "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}"
    # in case the mirrors have different file names we always store to the same
    DOWNLOAD_NAME "${PACKAGE_FILE}"
    DOWNLOAD_NO_PROGRESS ON
    CMAKE_ARGS ${DEPS_ARGS}
    BUILD_BYPRODUCTS ${PACKAGE_LIBRARY} )

  add_postinstall_step() # only if POSTINSTALL_SCRIPT is set
endmacro()


# variant for cmake files with encrypted sources which is snopt, scpip and hsl
# you need licenses to handle these files!
macro(create_external_encrypted_cmake_patched)

  assert_set(PACKAGE_NAME)
  assert_set(DEPS_ARGS)
  assert_set(DEPS_PREFIX)
  assert_set(PACKAGE_LIBRARY)
  assert_set(PACKAGE_KEY)
  assert_set(DEPS_SOURCE) # is non-standard for snopt
  assert_set(PATCHES_SCRIPT)

  # URL can take a list of mirrors but when there is a file, it needs to be the only one.
  # file means, that we have the source already in the cfsdeps cache. If not, we store there
  # see add_standard_mirrors()
  
  # we need to build the package - here in cmake style
  ExternalProject_Add("${PACKAGE_NAME}"
    PREFIX "${DEPS_PREFIX}"
    SOURCE_DIR "${DEPS_SOURCE}"
    
    URL "${PACKAGE_MIRRORS}"
    URL_MD5 "${PACKAGE_MD5}"
    # DOWNLOAD_DIR is ignored, if URL contains not the file mirror
    DOWNLOAD_DIR "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}"
    # in case the mirrors have different file names we always store to the same
    DOWNLOAD_NAME "${PACKAGE_FILE}"
    # disable automatic extract of the zip file by cmake (does not work for encrypted files)
    DOWNLOAD_NO_EXTRACT 1
    # after download comes update, then patch. 
    # usually we unpack into ${DEPS_PREFIX}/src/${PACKAGE_NAME}, however manual unzip adds a /snopt7
    UPDATE_COMMAND unzip -q -u -P ${PACKAGE_KEY} ${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}/${PACKAGE_FILE} -d ${DEPS_PREFIX}/src/${PACKAGE_NAME} 
    # now patch the unpacked source
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PATCHES_SCRIPT}"

    CMAKE_ARGS ${DEPS_ARGS}
    
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1 
    
    BUILD_BYPRODUCTS ${PACKAGE_LIBRARY} )

  add_postinstall_step() # only if POSTINSTALL_SCRIPT is set
endmacro()

# packages witch are note cmake based but which use the configure command
# auto patching when PATCHES_SCRIPT is set
macro(create_external_configure)

  assert_set(PACKAGE_NAME)
  assert_set(DEPS_PREFIX)
  assert_set(PACKAGE_LIBRARY)
  assert_set(DEPS_CONFIGURE)

  if(PATCHES_SCRIPT)
    # we need to build the package - here in configure style
    ExternalProject_Add("${PACKAGE_NAME}"
      PREFIX "${DEPS_PREFIX}"
      URL "${PACKAGE_MIRRORS}"
      URL_MD5 "${PACKAGE_MD5}"
      # DOWNLOAD_DIR is ignored, if URL contains not the file mirror
      DOWNLOAD_DIR "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}"
      # in case the mirrors have different file names we always store to the same
      DOWNLOAD_NAME "${PACKAGE_FILE}"
      DOWNLOAD_NO_PROGRESS ON 
      CONFIGURE_COMMAND ${DEPS_CONFIGURE_ENV} ${DEPS_SOURCE}/configure ${DEPS_CONFIGURE} 
      BUILD_BYPRODUCTS ${PACKAGE_LIBRARY} 
      # now patch the unpacked source
      PATCH_COMMAND ${CMAKE_COMMAND} -P "${PATCHES_SCRIPT}" )
  
  else()
    ExternalProject_Add("${PACKAGE_NAME}"
      PREFIX "${DEPS_PREFIX}"
      URL "${PACKAGE_MIRRORS}"
      URL_MD5 "${PACKAGE_MD5}"
      DOWNLOAD_DIR "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}"
      DOWNLOAD_NAME "${PACKAGE_FILE}"
      DOWNLOAD_NO_PROGRESS ON 
      CONFIGURE_COMMAND ${DEPS_CONFIGURE_ENV} ${DEPS_SOURCE}/configure ${DEPS_CONFIGURE} 
      BUILD_BYPRODUCTS ${PACKAGE_LIBRARY} )
  endif()     
    
  add_postinstall_step() # only if POSTINSTALL_SCRIPT is set
endmacro()

# for packages having a manual cleanup
macro(add_postinstall_step)
  assert_set(PACKAGE_NAME)

  if(POSTINSTALL_SCRIPT)
    ExternalProject_Add_Step(${PACKAGE_NAME} post_install
      COMMAND ${CMAKE_COMMAND} -P "${POSTINSTALL_SCRIPT}"
      # COMMAND ${CMAKE_COMMAND} -E echo  post_install ${POSTINSTALL_SCRIPT}
      DEPENDEES install )
  endif()   
endmacro()

# add custom step to zip a precompiled package to the cache and also to cmake_binary_dir
macro(add_external_storage_step)

  assert_set(PACKAGE_NAME)
  assert_set(ZIPTOCACHE_SCRIPT)
  assert_set(DEPS_INSTALL)

  unset(PI_STEP) # could otherwise have side effects for a package which does not have POSTINSTALL_SCRIPT
  if(POSTINSTALL_SCRIPT)
    set(PI_STEP post_install)
  endif()

  ExternalProject_Add_Step(${PACKAGE_NAME} cfsdeps_packaging
    COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE_SCRIPT}"
    DEPENDEES install ${PI_STEP} 
    DEPENDS "${ZIPTOCACHE_SCRIPT}"
    WORKING_DIRECTORY ${DEPS_INSTALL} )
endmacro()

# add install to to binary_dir copy in case we have no add_external_storage_step()
# only needed for install_dir packages, manifest packages already build into binary_dir
macro(add_install_dir_to_binary_step)
  assert_set(DEPS_INSTALL)

  set(INSTALL_DIR_TO_BINARY_SCRIPT "${DEPS_PREFIX}/${PACKAGE_NAME}-DepsCopyInstallDirToBinary.cmake")
  configure_file("${CMAKE_SOURCE_DIR}/cmake_modules/DepsCopyInstallDirToBinary.cmake.in"
                 "${INSTALL_DIR_TO_BINARY_SCRIPT}" @ONLY)  

  ExternalProject_Add_Step(${PACKAGE_NAME} install_dir_to_binary_dir
    COMMAND ${CMAKE_COMMAND} -P "${INSTALL_DIR_TO_BINARY_SCRIPT}"
    DEPENDEES install
    DEPENDS "${INSTALL_DIR_TO_BINARY_SCRIPT}"
    WORKING_DIRECTORY ${DEPS_INSTALL} )
endmacro()


