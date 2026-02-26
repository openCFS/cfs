# in this script we determine distro stuff like:
# build type (CFS_ARCH), distro name (CFS_DISTRO) and version (CFS_DISTRO_VER)
# this forms the CFS_ARCH_STR for precompiled package identification

#cmake_print_variables(CMAKE_SYSTEM_NAME)
#cmake_print_variables(CMAKE_SYSTEM_VERSION)
#cmake_print_variables(CMAKE_SYSTEM_PROCESSOR)

# first setup stuff to obtain distro info
if(UNIX) # linux and macos
 # not windows means Linux and Apple and we can use a convenient shell script based on uname
 set(DISTRO_SCRIPT "${CFS_SOURCE_DIR}/share/scripts/distro.sh")

 # process stuff to get info into cmake format (linux/mac)
 execute_process(COMMAND "${DISTRO_SCRIPT}" -c
                 OUTPUT_VARIABLE CFS_DISTRO_TEST
                 RESULT_VARIABLE RETVAL)

  execute_process(COMMAND "${CMAKE_COMMAND}" -E make_directory "${CFS_BINARY_DIR}/tmp"
                  WORKING_DIRECTORY "${CFS_BINARY_DIR}"
                  RESULT_VARIABLE RETVAL )

  string(REPLACE ";" "\n" CFS_DISTRO_TEST "${CFS_DISTRO_TEST}")
  file(WRITE "${CFS_BINARY_DIR}/tmp/distro_test.cmake" "${CFS_DISTRO_TEST}")

  include("${CFS_BINARY_DIR}/tmp/distro_test.cmake")
  
  # set it early and allow for change check. This is a cache variable, to set it manually, 
  # call cmake with -DCFS_ARCH=... then it is set and won't be overwritten by ARCH
else() # WIN32
  # we don't do disto.sh magic but use CMAKE
  set(DIST "Win")

  # CMAKE_SYSTEM_VERSION is like 10.0.22621 -> https://stackoverflow.com/questions/77409739/overview-cmake-system-version-values-for-windows/77410234#77410234
  if(CMAKE_SYSTEM_VERSION VERSION_LESS "10.0.20000")
     set(REV "10")
   else()
     set(REV "11") # even if we are on Win11 we (still) build for Win10
   endif()      

   # for Intel i7 cmake reporst "AMD64" but we use X86_64
   if(CMAKE_SYSTEM_PROCESSOR MATCHES "ARM64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
     set(ARCH "ARM64")
   else()
     set(ARCH "X86_64")
   endif()
endif()

set(CFS_ARCH "${ARCH}" CACHE STRING "Architecture CFS is built on.")
set_property(CACHE CFS_ARCH PROPERTY STRINGS X86_64 ARM64)
# hide for standard users
if(CFS_ARCH MATCHES "X86_64")
  mark_as_advanced(CFS_ARCH)
endif()


# Now set some global variables containing information about build/target platform.

# Since some major Linux enterprise distros are binary compatible
# across minor versions, we just set the DIST_FAMILY and the major version
# for them in CFS_DISTRO and CFS_DISTRO_VER and only provide more detailed
# infos in CFS_FULL_DISTRO and CFS_FULL_DISTRO_VER.

# on macOS we have REV "26.3" and MAJOR_REV "26.3" and actually the .3 is a minor
# MAJOR_MAJOR_REV is skipping the last dot stuff if there is a dot
string(FIND "${MAJOR_REV}" "." _LAST_DOT REVERSE)
if(_LAST_DOT GREATER -1)
  string(SUBSTRING "${MAJOR_REV}" 0 "${_LAST_DOT}" MAJOR_MAJOR_REV)
else()
  set(MAJOR_MAJOR_REV "${MAJOR_REV}")
endif()

set(CFS_FULL_DISTRO "${DIST}")
set(CFS_FULL_DISTRO_VER "${REV}")
if(APPLE)
  # DIST and DIST_FAMILY is "MACOSX"
  set(CFS_DISTRO "${DIST}")
  set(CFS_DISTRO_VER "${MAJOR_MAJOR_REV}")
elseif(DIST_FAMILY)
  set(CFS_DISTRO "${DIST_FAMILY}")
  set(CFS_DISTRO_VER "${MAJOR_REV}")
else()
  set(CFS_DISTRO "${DIST}")
  set(CFS_DISTRO_VER "${REV}")
endif()
    
# ARCH will be x86_64 or ARM64. E.g. for Apple M1 we can switch from ARM64 to build x86_64
  
# on macOS we enforce CMAKE_OSX_ARCHITECTURES to the single selected value. 
# We do not support multiple build targets in cfs
  
set(CFS_ARCH_STR "${CFS_DISTRO}_${CFS_DISTRO_VER}_${CFS_ARCH}")
set(CFS_BUILD_DISTRO "${DIST}_${REV}_${CFS_ARCH}")

# on macs there is optionally a CMAKE_OSX_ARCHITECTURES which could even be arm64 and x86_64 together.
# within cfs we use ARM64 instead of arm64
if(CFS_ARCH STREQUAL "ARM64")
  # we need a force otherwise it is not set. Might overwrite external settings?!
  set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE INTERNAL "" FORCE)
else()
  set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE INTERNAL "" FORCE)
endif()

set(LIB_SUFFIX "lib")

set(CFS_DISTRO "${CFS_DISTRO}" CACHE INTERNAL "String specifying the distribution CFS is built on.")
set(CFS_DISTRO_VER "${CFS_DISTRO_VER}" CACHE INTERNAL "Version of the distribution CFS is built on.")
# CFS_ARCH (set in windows and unix) might have been changed
set(CFS_ARCH_STR "${CFS_ARCH_STR}" CACHE INTERNAL "String specifying the architecture CFS is built on." FORCE) 
