# in this script we determine distro stuff like:
# build type (CFS_ARCH), distro name (CFS_DISTRO) and version (CFS_DISTRO_VER)
# this forms the CFS_ARCH_STR for precompiled package identification

# first setup stuff to obtain distro info
if(NOT WIN32)
 # not windows means Linux and Apple and we can use a convenient shell script based on uname
 set(DISTRO_SCRIPT "${CFS_SOURCE_DIR}/share/scripts/distro.sh")
else()
  #---------------------------------------------------------------------------
  # We are building on Windows and want to get the version.
  #---------------------------------------------------------------------------
  SET(DISTRO_SCRIPT "${CFS_SOURCE_DIR}/share/scripts/winver.bat")

  #---------------------------------------------------------------------------
  # If we are on Windows, we assume, that we use the Microsoft or Intel
  # toolchains. We want to determine the target architecture of our compiler
  # by compiling a little test program. 
  #---------------------------------------------------------------------------
  SET(CPP_SRC "int main(int argc, char** argv) {}")
  FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/tmp/test.cpp" "${CPP_SRC}")
  TRY_COMPILE(
    RESULT_VAR "${CMAKE_CURRENT_BINARY_DIR}/tmp"
    "${CMAKE_CURRENT_BINARY_DIR}/tmp/test.cpp"
    COPY_FILE "${CMAKE_CURRENT_BINARY_DIR}/tmp/test.exe"
    OUTPUT_VARIABLE BUILD
    )
  
  EXECUTE_PROCESS(
    COMMAND dumpbin /HEADERS "${CMAKE_CURRENT_BINARY_DIR}/tmp/test.exe"
    WORKING_DIRECTORY "${CFS_BINARY_DIR}"
    OUTPUT_VARIABLE DUMPBIN_ARCH
    RESULT_VARIABLE RETVAL
    )
  
  STRING(REPLACE "\n" ";" DUMPBIN_ARCH "${DUMPBIN_ARCH}")
  foreach(line IN ITEMS ${DUMPBIN_ARCH})
    IF(line MATCHES "machine \\(")
      #	MESSAGE(STATUS "line: ${line}")
      STRING(REPLACE "(" ";" output "${line}")
      LIST(GET output 1 output)
      STRING(REPLACE ")" ";" output "${output}")
      LIST(GET output 0 TARGET_ARCH)
      # MESSAGE(STATUS "ARCH: ${ARCH}")
    ENDIF()
  endforeach()
endif()

# process stuff to get info into cmake format (windows and linux/mac)
execute_process(COMMAND "${DISTRO_SCRIPT}" -c
                OUTPUT_VARIABLE CFS_DISTRO_TEST
                RESULT_VARIABLE RETVAL)

execute_process(COMMAND "${CMAKE_COMMAND}" -E make_directory "${CFS_BINARY_DIR}/tmp"
                WORKING_DIRECTORY "${CFS_BINARY_DIR}"
                RESULT_VARIABLE RETVAL )

string(REPLACE ";" "\n" CFS_DISTRO_TEST "${CFS_DISTRO_TEST}")
file(WRITE "${CFS_BINARY_DIR}/tmp/distro_test.cmake" "${CFS_DISTRO_TEST}")

include("${CFS_BINARY_DIR}/tmp/distro_test.cmake")
#MESSAGE("OS: ${OS}")
#MESSAGE("DIST: ${DIST}")
#MESSAGE("DIST_FAMILY: ${DIST_FAMILY}")
#MESSAGE("REV: ${REV}")
#MESSAGE("ARCH: ${ARCH}")
#MESSAGE("SUBARCH: ${SUBARCH}") # not used

# set it early and allow for change check. This is a cache variable, to set it manually, 
# call cmake with -DCFS_ARCH=... then it is set and won't be overwritten by ARCH
set(CFS_ARCH "${ARCH}" CACHE STRING "Architecture CFS is built on.")
set_property(CACHE CFS_ARCH PROPERTY STRINGS X86_64 ARM64)
# hide for standard users
if(CFS_ARCH MATCHES "X86_64")
  mark_as_advanced(CFS_ARCH)
endif()


#-------------------------------------------------------------------------------
# Now set some global variables containing informations about build/target plat.
#-------------------------------------------------------------------------------
if(NOT WIN32)
  #---------------------------------------------------------------------------
  # We are on Unix/Mac. Since some major enterprise distros are binary compatible
  # across minor versions, we just set the DIST_FAMILY and the major version
  # for them in CFS_DISTRO and CFS_DISTRO_VER and only provide more detailed
  # infos in CFS_FULL_DISTRO and CFS_FULL_DISTRO_VER.
  #---------------------------------------------------------------------------
  set(CFS_FULL_DISTRO "${DIST}")
  set(CFS_FULL_DISTRO_VER "${REV}")
  if(DIST_FAMILY)
    set(CFS_DISTRO "${DIST_FAMILY}")
    set(CFS_DISTRO_VER "${MAJOR_REV}")
  else()
    set(CFS_DISTRO "${DIST}")
    set(CFS_DISTRO_VER "${REV}")
  endif()
    
  # ARCH will be x86_64 or ARM64. E.g. for Apple M1 we can switch from ARM64 to build x86_64
  
  # on macOS we enforce CMAKE_OSX_ARCHITECTURES to the single selected value. 
  # We do not support multiple buid targests in cfs
  
  set(CFS_ARCH_STR "${CFS_DISTRO}_${CFS_DISTRO_VER}_${CFS_ARCH}")
  # todo: CFS_BUILD_DISTRO is almost CFS_ARCH_STR, only windows makes a difference?!
  set(CFS_BUILD_DISTRO "${DIST}_${REV}_${CFS_ARCH}")

  # on macs there is optionally a CMAKE_OSX_ARCHITECTURES which could even be arm64 and x86_64 together.
  # within cfs we use ARM64 instead of arm64
  if(CFS_ARCH STREQUAL "ARM64")
    # we need a force otherwise it is not set. Might overwrite external settings?!
    set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE INTERNAL "" FORCE)
  else()
    set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE INTERNAL "" FORCE)
  endif()

else()
  #---------------------------------------------------------------------------
  # We are on Windows. Since Windows is very compatible across versions we
  # use the C++ compiler toolchain version as CFS_DISTRO.
  #---------------------------------------------------------------------------
  # MSVC_SERVICE_PACK is no longer supported
  #    SET(CFS_DISTRO "${CFS_MSVC_SERVICE_PACK}")
  SET(CFS_FULL_DISTRO "${DIST}")
  SET(CFS_FULL_DISTRO_VER "${REV}")
  #    IF(DIST_FAMILY)
  #  SET(CFS_DISTRO "${DIST_FAMILY}")
  #  SET(CFS_DISTRO_VER "${MAJOR_REV}")
  SET(CFS_DISTRO "${DIST}")
  SET(CFS_DISTRO_VER "${REV}")
  SET(CFS_OS "${OS}")
  IF(TARGET_ARCH STREQUAL "x64")
    SET(CFS_ARCH "X86_64")
  ELSE()
    MESSAGE(FATAL_ERROR "Unsupported machine architecture '${TARGET_ARCH}' on Windows!")
  ENDIF()
      
  set(CFS_ARCH_STR "${CFS_DISTRO}_${CFS_ARCH}")
  
  SET(CFS_FULL_DISTRO "${CFS_DISTRO}")
  SET(CFS_FULL_DISTRO_VER "${CMAKE_C_COMPILER_VERSION}")

  SET(CFS_BUILD_DISTRO "${DIST}_${CMAKE_CXX_COMPILER_ID}_MSVC${CFS_MSVC_VERSION}_${ARCH}")
endif()

#message("CFS_DISTRO: ${CFS_DISTRO}")
#message("CFS_DISTRO_VER: ${CFS_DISTRO_VER}")
#message("CFS_ARCH: ${CFS_ARCH}")
#MESSAGE("CFS_ARCH_STR: ${CFS_ARCH_STR}")
#MESSAGE("CFS_BUILD_DISTRO ${CFS_BUILD_DISTRO}")
#MESSAGE("CMAKE_OSX_ARCHITECTURES ${CMAKE_OSX_ARCHITECTURES}")

set(LIB_SUFFIX "lib")

set(CFS_DISTRO "${CFS_DISTRO}" CACHE INTERNAL "String specifying the distribution CFS is built on.")
set(CFS_DISTRO_VER "${CFS_DISTRO_VER}" CACHE INTERNAL "Version of the distribution CFS is built on.")
# CFS_ARCH (set in windows and unix) might have been changed
set(CFS_ARCH_STR "${CFS_ARCH_STR}" CACHE INTERNAL "String specifying the architecture CFS is built on." FORCE) 
