#-------------------------------------------------------------------------------
# Collect output of compiler version command to determine compiler type and
# version information. 
#-------------------------------------------------------------------------------
EXEC_PROGRAM("${CMAKE_CXX_COMPILER} --version | head -1"
  ARGS
  OUTPUT_VARIABLE CFS_CXX_COMPILER_INFO
  RETURN_VALUE RETVAL)

EXEC_PROGRAM("${CMAKE_Fortran_COMPILER} --version | head -1"
  ARGS
  OUTPUT_VARIABLE CFS_FORTRAN_COMPILER_INFO
  RETURN_VALUE RETVAL)

#-------------------------------------------------------------------------------
# Determine what equivalent GNU version the compiler has, to check if it is
# compatible with the GNU C++ compiler on the system PATH.
#-------------------------------------------------------------------------------
SET(TEST_FILE "${CFS_BINARY_DIR}/include/test_gnu_version.hh")
FILE(WRITE "${TEST_FILE}" "<CFS_GNU_VER>__GNUC__.__GNUC_MINOR__</CFS_GNU_VER>")
EXEC_PROGRAM("${CMAKE_CXX_COMPILER} -E ${TEST_FILE}"
  ARGS
  OUTPUT_VARIABLE CFS_CXX_COMPILER_GNU_VER
  RETURN_VALUE RETVAL)
STRING(REGEX MATCH
  "<CFS_GNU_VER>.*</CFS_GNU_VER>"
  CFS_CXX_COMPILER_GNU_VER "${CFS_CXX_COMPILER_GNU_VER}")
STRING(REPLACE " " "" CFS_CXX_COMPILER_GNU_VER
  "${CFS_CXX_COMPILER_GNU_VER}")
STRING(REGEX MATCH "[0-9]+\\.[0-9]"
  CFS_CXX_COMPILER_GNU_VER "${CFS_CXX_COMPILER_GNU_VER}")

EXEC_PROGRAM("g++ -E ${TEST_FILE}"
  ARGS
  OUTPUT_VARIABLE GNU_CXX_COMPILER_VER
  RETURN_VALUE RETVAL)
STRING(REGEX MATCH
  "<CFS_GNU_VER>.*</CFS_GNU_VER>"
  GNU_CXX_COMPILER_VER "${GNU_CXX_COMPILER_VER}")
STRING(REPLACE " " "" GNU_CXX_COMPILER_VER
  "${GNU_CXX_COMPILER_VER}")
STRING(REGEX MATCH "[0-9]+\\.[0-9]"
  GNU_CXX_COMPILER_VER "${GNU_CXX_COMPILER_VER}")

# MESSAGE("CFS_CXX_COMPILER_GNU_VER: ${CFS_CXX_COMPILER_GNU_VER}")
# MESSAGE("GNU_CXX_COMPILER_VER: ${GNU_CXX_COMPILER_VER}")


# MESSAGE("C++ info: ${CFS_CXX_COMPILER_INFO}")
# MESSAGE("FORTRAN info: ${CFS_FORTRAN_COMPILER_INFO}")

#-------------------------------------------------------------------------------
# Check if we are using the GNU C++ compiler
#-------------------------------------------------------------------------------
IF(CMAKE_COMPILER_IS_GNUCXX)

  # MESSAGE("We are using the GNU C++ compiler. ${CMAKE_CXX_COMPILER}")

  #-----------------------------------------------------------------------------
  # Set the compiler name and compiler version
  #-----------------------------------------------------------------------------
  SET(CFS_CXX_COMPILER_NAME "GCC")
  SET(CFS_CXX_COMPILER_VER "${CFS_CXX_COMPILER_INFO}")

  STRING(REGEX MATCH 
    "[0-9]+\\.[0-9]+\\.[0-9]+"
    CFS_CXX_COMPILER_VER
    ${CFS_CXX_COMPILER_VER})

  #-----------------------------------------------------------------------------
  # Check if compiler has OpenMP support. GCC >= 4.2 has.
  #-----------------------------------------------------------------------------
  IF(USE_OPENMP)
    IF(CFS_CXX_COMPILER_VER GREATER "4.1")
      SET(CFS_C_FLAGS "-fopenmp")
    ELSE(CFS_CXX_COMPILER_VER GREATER "4.1")
      MESSAGE("You chose to use OpenMP but your current GCC does not support it!")
    ENDIF(CFS_CXX_COMPILER_VER GREATER "4.1")
  ENDIF(USE_OPENMP)

  #-----------------------------------------------------------------------------
  # Determine compiler/linker flags according to build type
  #-----------------------------------------------------------------------------
  IF(DEBUG)

    SET(CFS_C_FLAGS "-std=c++98 -Wall -pedantic -Werror -fmessage-length=0 ${CFS_C_FLAGS}")
    # -Wold-style-cast Warnings about old C style casts. Since external libraries
    # make extensive use of it, we switch it off. To filter out the warnings in our own
    # code a command line like the following might be used
    # fgrep 'warning: use of old-style cast' out.txt | grep CFS_SOURCE_DIR | sort -u > old-style-cast.txt
    SET(CFS_CXX_FLAGS "-ftemplate-depth-55")
    SET(CFS_SUPPRESSIONS "-Wno-long-long -Wno-unknown-pragmas -Wno-comment")
    SET(CHECK_MEM_ALLOC 1)
    SET(CHECK_TYPE_CASTS 1)

    IF(CFS_PROFILING)
      SET(CFS_PROF_FLAGS "-pg")
    ENDIF(CFS_PROFILING)

  ELSE(DEBUG)

    SET(CFS_SUPPRESSIONS "-Wno-long-long -Wno-unknown-pragmas -Wno-comment")
    SET(CFS_C_FLAGS "-std=c++98 -Wall -pedantic -Werror -fmessage-length=0 ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-ftemplate-depth-55")

    IF(CFS_ARCH STREQUAL "I386")
      SET(CFS_OPT_FLAGS "-m32 -march=pentium4")
    ENDIF(CFS_ARCH STREQUAL "I386")

    IF(CFS_ARCH STREQUAL "X86_64")
      SET(CFS_OPT_FLAGS "-m64 -march=k8")
    ENDIF(CFS_ARCH STREQUAL "X86_64")

    IF(CFS_ARCH STREQUAL "IA64")
      SET(CFS_OPT_FLAGS "-m64 -mtune-arch=itanium2")
    ENDIF(CFS_ARCH STREQUAL "IA64")

  ENDIF(DEBUG)

ENDIF(CMAKE_COMPILER_IS_GNUCXX)

#-------------------------------------------------------------------------------
# Check if we are using the GNU Fortran (95) compiler
#-------------------------------------------------------------------------------
IF(CFS_FORTRAN_COMPILER_INFO MATCHES "GNU")
  #-----------------------------------------------------------------------------
  # Determine name of Fortran compiler
  #-----------------------------------------------------------------------------
  SET(CFS_FORTRAN_COMPILER_NAME "GNU")

  #-----------------------------------------------------------------------------
  # Determine version of Fortran compiler
  #-----------------------------------------------------------------------------
  SET(CFS_FORTRAN_COMPILER_VER ${CFS_FORTRAN_COMPILER_INFO})
  STRING(REGEX MATCH 
    "[0-9]+\\.[0-9]+\\.[0-9]+"
    CFS_FORTRAN_COMPILER_VER
    ${CFS_FORTRAN_COMPILER_VER})
ENDIF(CFS_FORTRAN_COMPILER_INFO MATCHES "GNU")


#-------------------------------------------------------------------------------
# Check for Intel C++ compiler
#-------------------------------------------------------------------------------
IF(CFS_CXX_COMPILER_INFO MATCHES "ICC")
  #-----------------------------------------------------------------------------
  # Set the compiler name and compiler version
  #-----------------------------------------------------------------------------
  SET(CFS_CXX_COMPILER_NAME "ICC")

  STRING(REGEX MATCH 
    "[0-9]+\\.[0-9]+"
    CFS_CXX_COMPILER_VER
    ${CFS_CXX_COMPILER_INFO})

  #-----------------------------------------------------------------------------
  # Check for the case that Intel compiler is just GCC 4.2 compatible but the
  # system GCC is version 4.3. In this case Intel C++ fails to compile due
  # to incompatibilities.
  #-----------------------------------------------------------------------------
  IF(CFS_CXX_COMPILER_GNU_VER EQUAL "4.2" AND
     GNU_CXX_COMPILER_VER EQUAL "4.3")
    MESSAGE(FATAL_ERROR "Intel C+ + ${CFS_CXX_COMPILER_VER} is known to be broken with g++ ${GNU_CXX_COMPILER_VER}! Intel C++ 11.x should work.")
  ELSE(CFS_CXX_COMPILER_GNU_VER EQUAL "4.2" AND
      GNU_CXX_COMPILER_VER  EQUAL "4.3")
    IF(CFS_CXX_COMPILER_GNU_VER LESS GNU_CXX_COMPILER_VER)
      MESSAGE("Intel C++ ${CFS_CXX_COMPILER_GNU_VER} might not work with g++ ${GNU_CXX_COMPILER_VER}.")
    ENDIF(CFS_CXX_COMPILER_GNU_VER LESS GNU_CXX_COMPILER_VER)
  ENDIF(CFS_CXX_COMPILER_GNU_VER EQUAL "4.2" AND
    GNU_CXX_COMPILER_VER EQUAL "4.3")

  #-----------------------------------------------------------------------------
  # Check for a parallel compiler
  #-----------------------------------------------------------------------------
  IF(USE_OPENMP)
    SET(CFS_C_FLAGS "-openmp")
  ENDIF(USE_OPENMP)

  #-----------------------------------------------------------------------------
  # Determine compiler/linker flags according to build type
  #-----------------------------------------------------------------------------
  IF(DEBUG)

    SET(CFS_C_FLAGS "-g -ansi -w1 -Wcheck -Werror ${CFS_C_FLAGS}")
    #---------------------------------------------------------------------------
    # Disable warnings about hidden overriden functions of base classes,
    # unknown pragmas (openmp, etc.) and multiline comments.
    #---------------------------------------------------------------------------
    SET(CFS_SUPPRESSIONS "-wd1125,654 -Wno-unknown-pragmas -Wno-comment")
    SET(CHECK_MEM_ALLOC 1)
    SET(CHECK_TYPE_CASTS 1)

    IF(CFS_PROFILING)
      SET(CFS_PROF_FLAGS "-pg")
    ENDIF(CFS_PROFILING)

  ELSE(DEBUG)

    SET(CFS_C_FLAGS "-ansi -w0 -Werror ${CFS_C_FLAGS}")
    SET(CFS_SUPPRESSIONS "-Wno-unknown-pragmas -Wno-comment")

  ENDIF(DEBUG)


ENDIF(CFS_CXX_COMPILER_INFO MATCHES "ICC")

#-------------------------------------------------------------------------------
# Check for Intel Fortran compiler
#-------------------------------------------------------------------------------
IF(CFS_FORTRAN_COMPILER_INFO MATCHES "IFORT")
  #-----------------------------------------------------------------------------
  # Set Intel Fortran compiler name and version
  #-----------------------------------------------------------------------------
  SET(CFS_FORTRAN_COMPILER_NAME "IFORT")

  STRING(REGEX MATCH 
    "[0-9]+\\.[0-9]+"
    CFS_FORTRAN_COMPILER_VER
    ${CFS_FORTRAN_COMPILER_INFO})

  #-----------------------------------------------------------------------------
  # Set Intel Fortran library paths in dedicated variables
  #-----------------------------------------------------------------------------
  STRING(REGEX REPLACE "bin/ifort" "lib" IFORT_LIB_PATH
      "${CMAKE_Fortran_COMPILER}")
  LINK_DIRECTORIES(${IFORT_LIB_PATH})
    
  SET(CFS_FORTRAN_DYNRT_LIBS ifcoremt dl)
  IF(NOT CFS_ARCH STREQUAL "IA64")
    SET(CFS_FORTRAN_DYNRT_LIBS 
      svml
      ${CFS_FORTRAN_DYNRT_LIBS})
  ENDIF(NOT CFS_ARCH STREQUAL "IA64")
  SET(CFS_FORTRAN_STATRT_LIBS
    ifcoremt_pic
    irc
    )
    
ENDIF(CFS_FORTRAN_COMPILER_INFO MATCHES "IFORT")


#-------------------------------------------------------------------------------
# Check if we have a valid combination of C++ and Fortran compilers.
#-------------------------------------------------------------------------------
IF(CFS_CXX_COMPILER_NAME STREQUAL "" OR
    CFS_FORTRAN_COMPILER_NAME STREQUAL "")

  MESSAGE(SEND_ERROR "Combination of C++ and Fortran compilers not supported!")

ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "" OR
  CFS_FORTRAN_COMPILER_NAME STREQUAL "")

#-------------------------------------------------------------------------------
# Set compiler/linker flags for all build types
#-------------------------------------------------------------------------------
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFS_C_FLAGS} ${CFS_PROF_FLAGS} ${CFS_OPT_FLAGS} ${CFS_SUPPRESSIONS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CFS_C_FLAGS} ${CFS_CXX_FLAGS} ${CFS_PROF_FLAGS} ${CFS_OPT_FLAGS} ${CFS_SUPPRESSIONS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CFS_PROF_FLAGS}")
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CFS_PROF_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CFS_PROF_FLAGS}")

# MESSAGE("C++ name: ${CFS_CXX_COMPILER_NAME}")
# MESSAGE("C++ version: ${CFS_CXX_COMPILER_VER}")
# MESSAGE("FORTRAN Compiler ${CFS_FORTRAN_COMPILER_NAME}")
# MESSAGE("FORTRAN Compiler version ${CFS_FORTRAN_COMPILER_VER}")
# MESSAGE("CMAKE BUILD TYPE: ${CMAKE_BUILD_TYPE}")
# MESSAGE("CFS_ARCH_STR: ${CFS_ARCH_STR}")
