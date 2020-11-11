IF(APPLE)

  #-----------------------------------------------------------------------------
  # We are building on Unix and want to get informations about the distro.
  #-----------------------------------------------------------------------------
  SET(DISTRO_SCRIPT "${CFS_SOURCE_DIR}/share/scripts/distro.sh")
ELSE()
  IF(WIN32)
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
#	MESSAGE(STATUS "ARCH: ${ARCH}")
      ENDIF()
    endforeach()
    
  ELSE()
    #---------------------------------------------------------------------------
    # We are building on Unix and want to get informations about the distro.
    #---------------------------------------------------------------------------
    SET(DISTRO_SCRIPT "${CFS_SOURCE_DIR}/share/scripts/distro.sh")
  ENDIF()
ENDIF()


#-------------------------------------------------------------------------------
# Get informations about the distro / Windows version in CMake format.
#-------------------------------------------------------------------------------
EXECUTE_PROCESS(
  COMMAND "${DISTRO_SCRIPT}" -c
  OUTPUT_VARIABLE CFS_DISTRO_TEST
  RESULT_VARIABLE RETVAL
  )

EXECUTE_PROCESS(
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${CFS_BINARY_DIR}/tmp"
  WORKING_DIRECTORY "${CFS_BINARY_DIR}"
  RESULT_VARIABLE RETVAL
  )

STRING(REPLACE ";" "\n" CFS_DISTRO_TEST "${CFS_DISTRO_TEST}")
FILE(WRITE "${CFS_BINARY_DIR}/tmp/distro_test.cmake" "${CFS_DISTRO_TEST}")

INCLUDE("${CFS_BINARY_DIR}/tmp/distro_test.cmake")
#MESSAGE("OS: ${OS}")
#MESSAGE("DIST: ${DIST}")
#MESSAGE("DIST_FAMILY: ${DIST_FAMILY}")
#MESSAGE("REV: ${REV}")
#MESSAGE("ARCH: ${ARCH}")
#MESSAGE("SUBARCH: ${SUBARCH}")

#-------------------------------------------------------------------------------
# Now set some global variables containing informations about build/target plat.
#-------------------------------------------------------------------------------
IF(NOT WIN32)
  #---------------------------------------------------------------------------
  # We are on Unix. Since some major enterprise distros are binary compatible
  # across minor versions, we just set the DIST_FAMILY and the major version
  # for them in CFS_DISTRO and CFS_DISTRO_VER and only provide more detailed
  # infos in CFS_FULL_DISTRO and CFS_FULL_DISTRO_VER.
  #---------------------------------------------------------------------------
  SET(CFS_FULL_DISTRO "${DIST}")
  SET(CFS_FULL_DISTRO_VER "${REV}")
  IF(DIST_FAMILY)
    SET(CFS_DISTRO "${DIST_FAMILY}")
    SET(CFS_DISTRO_VER "${MAJOR_REV}")
  ELSE()
    SET(CFS_DISTRO "${DIST}")
    SET(CFS_DISTRO_VER "${REV}")
  ENDIF()

    
  SET(CFS_OS "${OS}")
  SET(CFS_ARCH "${ARCH}")
  SET(CFS_ARCH_STR "${CFS_DISTRO}_${CFS_DISTRO_VER}_${CFS_ARCH}")
  # Determine the subarchitecture of the platform.
  SET(CFS_SUBARCH "${SUBARCH}")
    
  SET(CFS_BUILD_OS "${OS}")
  SET(CFS_BUILD_DISTRO "${DIST}_${REV}_${ARCH}")
  SET(CFS_TARGET_OS "${OS}")

  # MESSAGE("APPLE = ${APPLE}")
  # MESSAGE("CFS_ARCH = ${CFS_ARCH}")
  # MESSAGE("CFS_DISTRO = ${CFS_DISTRO}")
  # MESSAGE("CFS_DISTRO_VER = ${CFS_DISTRO_VER}")
  # MESSAGE("MAJOR_REV = ${MAJOR_REV}")
    
        
ELSE()
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
      
  #SET(CFS_ARCH_STR "${CFS_DISTRO}_${CFS_ARCH}")
  SET(CFS_ARCH_STR "")

  SET(CFS_FULL_DISTRO "${CFS_DISTRO}")
  SET(CFS_FULL_DISTRO_VER "${CMAKE_C_COMPILER_VERSION}")

  SET(CFS_BUILD_OS "${OS}")
  SET(CFS_BUILD_DISTRO "${DIST}_${CFS_CXX_COMPILER_NAME}_MSVC${CFS_MSVC_VERSION}_${ARCH}")
  SET(CFS_TARGET_OS "${OS}")

ENDIF()

# MESSAGE("CFS_DISTRO: ${CFS_DISTRO}")
# MESSAGE("CFS_DISTRO_VER: ${CFS_DISTRO_VER}")
# MESSAGE("CFS_ARCH: ${CFS_ARCH}")
# MESSAGE("CFS_ARCH_STR: ${CFS_ARCH_STR}")
# MESSAGE("CFS_SUBARCH: ${CFS_SUBARCH}")
# MESSAGE("CFS_BUILD_DISTRO ${CFS_BUILD_DISTRO}")


#-----------------------------------------------------------------------------
# Set a few distribution specific variables (FORTRAN libs, ...)
# GFORTRAN_LIBRARY is defined in cmake_modules/FindFortranLibs.cmake
#-----------------------------------------------------------------------------
IF(CFS_DISTRO STREQUAL "SUSE" OR
   CFS_DISTRO STREQUAL "OPENSUSE" OR
   CFS_DISTRO STREQUAL "SLE")

ELSEIF(CFS_DISTRO STREQUAL "DEBIAN")

ELSEIF(CFS_DISTRO STREQUAL "LMDE")

ELSEIF(CFS_DISTRO STREQUAL "UBUNTU")

ELSEIF(CFS_DISTRO STREQUAL "LINUXMINT")

ELSEIF(CFS_DISTRO STREQUAL "FEDORA" OR    
    CFS_DISTRO STREQUAL "REDHAT" OR
    CFS_DISTRO STREQUAL "RED" OR
    CFS_DISTRO STREQUAL "RHEL" OR
    CFS_DISTRO STREQUAL "CENTOS")

ELSEIF(CFS_DISTRO MATCHES "MSVC")
  # Just to support Microsoft toolchain.
ELSEIF(CFS_DISTRO MATCHES "WINDOWS")
ELSEIF(CFS_DISTRO MATCHES "WIN10")

ELSEIF(CFS_DISTRO STREQUAL "MACOSX")

  IF(CMAKE_OSX_ARCHITECTURES)
    LIST(LENGTH CMAKE_OSX_ARCHITECTURES NUM_ARCH)

    IF(NUM_ARCH GREATER 1)
      MESSAGE(FATAL_ERROR "Only binaries for one Mac OS X architecture can be built at a time!")
    ENDIF(NUM_ARCH GREATER 1)

    LIST(GET CMAKE_OSX_ARCHITECTURES 0 OSX_ARCH)
    # MESSAGE("OSX_ARCH ${OSX_ARCH}")

    IF(OSX_ARCH STREQUAL "i386")
      SET(CFS_ARCH "I386")
    ELSE(OSX_ARCH STREQUAL "i386")
      IF(OSX_ARCH STREQUAL "x86_64")
	SET(CFS_ARCH "X86_64")
      ELSE(OSX_ARCH STREQUAL "x86_64")
	MESSAGE(WARNING "Building for architecture ${OSX_ARCH} not supported!")
      ENDIF(OSX_ARCH STREQUAL "x86_64")
    ENDIF(OSX_ARCH STREQUAL "i386")

  ELSE(CMAKE_OSX_ARCHITECTURES)
    IF(CFS_DISTRO_VER VERSION_GREATER 10.5)
      SET(CMAKE_OSX_ARCHITECTURES "x86_64")
    ELSE(CFS_DISTRO_VER VERSION_GREATER 10.5)
      SET(CMAKE_OSX_ARCHITECTURES "i386")
    ENDIF(CFS_DISTRO_VER VERSION_GREATER 10.5)
  ENDIF(CMAKE_OSX_ARCHITECTURES)
ELSE()
  MESSAGE(WARNING "CFS_DISTRO '${CFS_DISTRO}' not supported!")
ENDIF()

#-----------------------------------------------------------------------------
# Determine name of library directory and system bits.
#-----------------------------------------------------------------------------
IF(CFS_ARCH STREQUAL "I386")
  # Set path suffix for system libs
  SET(LIB_SUFFIX "lib")
  SET(SYSTEM_BITS "32")
ELSEIF(CFS_ARCH STREQUAL "X86_64")
  # Set path suffix for system libs
  SET(LIB_SUFFIX "lib64")
  SET(SYSTEM_BITS "64")
ELSEIF(CFS_ARCH STREQUAL "IA64")
  # Set path suffix for system libs
  SET(LIB_SUFFIX "lib")
  SET(SYSTEM_BITS "64")
ENDIF()

SET(SYSTEM_BITS "${SYSTEM_BITS}"
  CACHE INTERNAL "Determins if this is a '32' or '64' bit system.")
SET(CFS_DISTRO "${CFS_DISTRO}"
  CACHE INTERNAL "String specifying the distribution CFS is built on.")
SET(CFS_DISTRO_VER "${CFS_DISTRO_VER}"
  CACHE INTERNAL "Version of the distribution CFS is built on.")
SET(CFS_ARCH "${CFS_ARCH}"
  CACHE INTERNAL "Architecture CFS is built on.")
SET(CFS_ARCH_STR "${CFS_ARCH_STR}"
  CACHE INTERNAL "String specifying the architecture CFS is built on.")
SET(CFS_SUBARCH "${CFS_SUBARCH}"
  CACHE INTERNAL "Sub-Architecture CFS is built on.")
