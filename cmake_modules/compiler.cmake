#-------------------------------------------------------------------------------
# Determine what equivalent GNU version the compiler has, to check if it is
# compatible with the GNU C++ compiler on the system PATH.
#-------------------------------------------------------------------------------
EXEC_PROGRAM("${PERL} ${CFS_SOURCE_DIR}/share/scripts/identify_compiler.pl g++ ${CFS_SOURCE_DIR}/share/scripts/IdentifyCXXCompiler.cpp cmake > ${CFS_BINARY_DIR}/CMakeFiles/out.cmake"
  ARGS
  OUTPUT_VARIABLE CC_COMPILER_INFO
  RETURN_VALUE RETVAL)

INCLUDE(${CFS_BINARY_DIR}/CMakeFiles/out.cmake)

SET(GNU_CXX_COMPILER_VER "${CXX_VERSION}")

#-------------------------------------------------------------------------------
# Collect output of compiler version command to determine compiler type and
# version information. 
#-------------------------------------------------------------------------------
# EXEC_PROGRAM("${PERL} ${CFS_DEPS_ROOT}/utils/perl/identify_compiler.pl ${CMAKE_C_COMPILER} ${CFS_DEPS_ROOT}/utils/compiler/IdentifyCCompiler.c cmake > ${CFS_BINARY_DIR}/CMakeFiles/out.cmake"
#   ARGS
#   OUTPUT_VARIABLE CC_COMPILER_INFO
#   RETURN_VALUE RETVAL)

# INCLUDE(${CFS_BINARY_DIR}/CMakeFiles/out.cmake)


EXEC_PROGRAM("${PERL} ${CFS_SOURCE_DIR}/share/scripts/identify_compiler.pl ${CMAKE_CXX_COMPILER} ${CFS_SOURCE_DIR}/share/scripts/IdentifyCXXCompiler.cpp cmake > ${CFS_BINARY_DIR}/CMakeFiles/out.cmake"
  ARGS
  OUTPUT_VARIABLE CXX_COMPILER_INFO
  RETURN_VALUE RETVAL)

INCLUDE(${CFS_BINARY_DIR}/CMakeFiles/out.cmake)
#-------------------------------------------------------------------------------
# Set the C++ compiler name and compiler version
#-------------------------------------------------------------------------------
SET(CFS_CXX_COMPILER_NAME ${CXX_ID})
SET(CFS_CXX_COMPILER_VER "${CXX_VERSION}")
SET(CFS_CXX_COMPILER_GNU_VER "${CXX_GCC_VERSION}")

EXEC_PROGRAM("${PERL} ${CFS_SOURCE_DIR}/share/scripts/identify_compiler.pl ${CMAKE_Fortran_COMPILER} ${CFS_SOURCE_DIR}/share/scripts/IdentifyFortranCompiler.F90 cmake > ${CFS_BINARY_DIR}/CMakeFiles/out.cmake"
  ARGS
  OUTPUT_VARIABLE FORTRAN_COMPILER_INFO
  RETURN_VALUE RETVAL)

INCLUDE(${CFS_BINARY_DIR}/CMakeFiles/out.cmake)

#-------------------------------------------------------------------------------
# Set the Fortran compiler name and compiler version
#-------------------------------------------------------------------------------
SET(CFS_FORTRAN_COMPILER_NAME ${FC_ID})
SET(CFS_FORTRAN_COMPILER_VER "${FC_VERSION}")

#-------------------------------------------------------------------------------
# Check if we are using the GNU C++ compiler
#-------------------------------------------------------------------------------
IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")

  # MESSAGE("We are using the GNU C++ compiler. ${CMAKE_CXX_COMPILER}")

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
    SET(CFS_C_FLAGS "-std=c++98 -Wall -fmessage-length=0 ${CFS_C_FLAGS}")
    # -Wold-style-cast Warnings about old C style casts. Since external libraries
    # make extensive use of it, we switch it off. To filter out the warnings in our own
    # code a command line like the following might be used
    # fgrep 'warning: use of old-style cast' out.txt | grep CFS_SOURCE_DIR | sort -u > old-style-cast.txt
    # 
    # -frounding-math: is needed for CGAL library
    SET(CFS_CXX_FLAGS "-ftemplate-depth-55 -frounding-math")
    SET(CFS_SUPPRESSIONS "-Wno-long-long -Wno-unknown-pragmas -Wno-comment")
    SET(CHECK_MEM_ALLOC 1)

    IF(CFS_PROFILING)
      SET(CFS_PROF_FLAGS "-pg")
    ENDIF(CFS_PROFILING)

  ELSE(DEBUG)

    SET(CFS_SUPPRESSIONS "-Wno-long-long -Wno-unknown-pragmas -Wno-comment")
    SET(CFS_C_FLAGS "-std=c++98 -Wall -fmessage-length=0 ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-ftemplate-depth-55")

    IF(CFS_ARCH STREQUAL "I386")
      SET(CFS_OPT_FLAGS "-m32 -march=pentium4")
    ENDIF(CFS_ARCH STREQUAL "I386")

    IF(CFS_ARCH STREQUAL "X86_64")
      SET(CFS_OPT_FLAGS "-m64 -march=k8 -msse2")
    ENDIF(CFS_ARCH STREQUAL "X86_64")

    IF(CFS_ARCH STREQUAL "IA64")
      SET(CFS_OPT_FLAGS "-m64 -mtune-arch=itanium2 -msse2")
    ENDIF(CFS_ARCH STREQUAL "IA64")
    
    IF(CFS_CXX_COMPILER_VER EQUAL "4.2" AND
       CFS_GCC42_OPT_SWITCHES)       
     SET(CFS_OPT_FLAGS "${CFS_GCC42_OPT_SWITCHES}")
    ENDIF(CFS_CXX_COMPILER_VER EQUAL "4.2" AND
          CFS_GCC42_OPT_SWITCHES)

    IF(CFS_CXX_COMPILER_VER EQUAL "4.3")
      IF(CFS_GCC43_OPT_SWITCHES)
        SET(CFS_OPT_FLAGS "${CFS_GCC43_OPT_SWITCHES}")
      ENDIF(CFS_GCC43_OPT_SWITCHES)
      
      IF(CFS_BLAS_LAPACK STREQUAL "ACML" AND CFS_ARCH STREQUAL "X86_64")
        SET(CFS_OPT_FLAGS "${CFS_OPT_FLAGS} -mveclibabi=acml")
      ENDIF(CFS_BLAS_LAPACK STREQUAL "ACML" AND CFS_ARCH STREQUAL "X86_64")
    ENDIF(CFS_CXX_COMPILER_VER EQUAL "4.3")

    IF(CFS_CXX_COMPILER_VER EQUAL "4.4")
        SET(CFS_OPT_FLAGS "${CFS_OPT_FLAGS} -mtune=native -march=native")
    ENDIF(CFS_CXX_COMPILER_VER EQUAL "4.4")
      
  ENDIF(DEBUG)
  
  IF(NOT USE_OPENMP)
    SET(CFS_C_FLAGS "-Werror ${CFS_C_FLAGS}")
  ENDIF(NOT USE_OPENMP)

  IF(NOT USE_INTERPOLATION)
    SET(CFS_C_FLAGS "-pedantic ${CFS_C_FLAGS}")
  ENDIF(NOT USE_INTERPOLATION)


ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")

#-------------------------------------------------------------------------------
# Check for Intel C++ compiler
#-------------------------------------------------------------------------------
IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  #-----------------------------------------------------------------------------
  # Check for the case that Intel compiler is just GCC 4.2 compatible but the
  # system GCC is version 4.3. In this case Intel C++ fails to compile due
  # to incompatibilities.
  #-----------------------------------------------------------------------------
  IF(CFS_CXX_COMPILER_GNU_VER EQUAL "4.2" AND
     GNU_CXX_COMPILER_VER EQUAL "4.3")
    MESSAGE(FATAL_ERROR "Intel C++ ${CFS_CXX_COMPILER_VER} is known to be broken with g++ ${GNU_CXX_COMPILER_VER}! Intel C++ 11.x should work.")
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
    SET(CHECK_MEM_ALLOC 1)

    IF(CFS_PROFILING)
      SET(CFS_PROF_FLAGS "-pg")
    ENDIF(CFS_PROFILING)
  ELSE(DEBUG)
    SET(CFS_C_FLAGS "-O3 -ansi -w0 -Werror ${CFS_C_FLAGS}")
    SET(CFS_SUPPRESSIONS "-wd1125,654,980 -Wno-unknown-pragmas -Wno-comment")

    IF(CFS_CXX_COMPILER_VER MATCHES "10\\.")
      SET(CFS_OPT_FLAGS "${CFS_INTEL10_OPT_SWITCHES}")
    ENDIF(CFS_CXX_COMPILER_VER MATCHES "10\\.")

    IF(CFS_CXX_COMPILER_VER MATCHES "11\\.")
      SET(CFS_OPT_FLAGS "${CFS_INTEL11_OPT_SWITCHES}")
    ENDIF(CFS_CXX_COMPILER_VER MATCHES "11\\.")
  ENDIF(DEBUG)

  #---------------------------------------------------------------------------
  # Disable warnings about hidden overriden functions of base classes,
  # unknown pragmas (openmp, etc.) and multiline comments.
  #---------------------------------------------------------------------------
  SET(CFS_SUPPRESSIONS "-wd1125,654")
  SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-unknown-pragmas")
  SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-comment")

  IF(CFS_CXX_COMPILER_VER MATCHES "11\\.")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -fno-builtin-std::basic_istream::get")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -fno-builtin-std::max")
  ENDIF(CFS_CXX_COMPILER_VER MATCHES "11\\.")

ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")

#-------------------------------------------------------------------------------
# Check for Intel Fortran compiler
#-------------------------------------------------------------------------------
IF(CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")
  #-----------------------------------------------------------------------------
  # Set Intel Fortran library paths in dedicated variables
  #-----------------------------------------------------------------------------
  STRING(REGEX REPLACE "bin/ifort" "lib" IFORT_LIB_PATH
      "${CMAKE_Fortran_COMPILER}")

  #-----------------------------------------------------------------------------
  # Paths for libraries and compilers are different for version 11 compilers.
  # There is even a difference between Linux and Mac OS X! Thank you Intel!!!
  #-----------------------------------------------------------------------------
  IF(CFS_FORTRAN_COMPILER_VER MATCHES "11\\." AND
      NOT CFS_DISTRO STREQUAL "MACOSX")
    IF(CFS_ARCH STREQUAL "I386")
      SET(LD "ia32")
    ENDIF(CFS_ARCH STREQUAL "I386")
    IF(CFS_ARCH STREQUAL "X86_64")
      SET(LD "intel64")
    ENDIF(CFS_ARCH STREQUAL "X86_64")
    IF(CFS_ARCH STREQUAL "IA64")
      SET(LD "ia64")
    ENDIF(CFS_ARCH STREQUAL "IA64")
    
    STRING(REGEX REPLACE "bin/${LD}/ifort" "lib/${LD}" IFORT_LIB_PATH
	"${CMAKE_Fortran_COMPILER}")
  ENDIF(CFS_FORTRAN_COMPILER_VER MATCHES "11\\." AND
    NOT CFS_DISTRO STREQUAL "MACOSX")

  LINK_DIRECTORIES(${IFORT_LIB_PATH})
    
  SET(CFS_FORTRAN_DYNRT_LIBS ifcoremt imf dl)
  IF(NOT CFS_ARCH STREQUAL "IA64")
    SET(CFS_FORTRAN_DYNRT_LIBS 
      svml
      ${CFS_FORTRAN_DYNRT_LIBS})
  ENDIF(NOT CFS_ARCH STREQUAL "IA64")
  SET(CFS_FORTRAN_STATRT_LIBS
    ifcoremt_pic
    irc
    )
    
ENDIF(CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")


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
