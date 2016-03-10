SET(IDCOMP_TEMPL "${CFS_SOURCE_DIR}/share/scripts/identify_compiler.cmake.in")
SET(COMPILER_ID_FILE "${CFS_BINARY_DIR}/CMakeFiles/out.cmake")
SET(IDENTIFY_COMPILER_SRC "IdentifyCXXCompiler.cpp")

include(CheckCXXSourceCompiles)

#-------------------------------------------------------------------------------
# Determine what equivalent GNU version the compiler has, to check if it is
# compatible with the GNU C++ compiler on the system PATH.
#-------------------------------------------------------------------------------
IF(UNIX OR MINGW)
  #-----------------------------------------------------------------------------
  # Let's first find out some infos about the system GNU compiler first.
  #-----------------------------------------------------------------------------
  SET(COMPILER "g++")
  
  SET(ID_GXX "${CFS_BINARY_DIR}/share/scripts/identify_gxx.cmake")
  CONFIGURE_FILE("${IDCOMP_TEMPL}" "${ID_GXX}" @ONLY)
  
  EXECUTE_PROCESS(
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${CFS_BINARY_DIR}/tmp"
    WORKING_DIRECTORY "${CFS_BINARY_DIR}"
    RESULT_VARIABLE RETVAL
    )
  
  EXECUTE_PROCESS(
    COMMAND "${CMAKE_COMMAND}" -P "${ID_GXX}"
    WORKING_DIRECTORY "${CFS_BINARY_DIR}/tmp"
    RESULT_VARIABLE RETVAL
    )
  
  INCLUDE("${COMPILER_ID_FILE}")
  
  SET(GNU_CXX_COMPILER_VER "${CXX_VERSION}")
ENDIF(UNIX OR MINGW)



#-------------------------------------------------------------------------------
# Collect output of compiler version command to determine compiler type and
# version information. 
#-------------------------------------------------------------------------------
SET(COMPILER "")

SET(ID_CXX "${CFS_BINARY_DIR}/share/scripts/identify_cxx.cmake")
CONFIGURE_FILE("${IDCOMP_TEMPL}" "${ID_CXX}" @ONLY)

EXECUTE_PROCESS(
  COMMAND "${CMAKE_COMMAND}" -P "${ID_CXX}"
  WORKING_DIRECTORY "${CFS_BINARY_DIR}/tmp"
  RESULT_VARIABLE RETVAL
  )

INCLUDE("${COMPILER_ID_FILE}")

#-------------------------------------------------------------------------------
# Set the C++ compiler name and compiler version
#-------------------------------------------------------------------------------
SET(CFS_CXX_COMPILER_NAME ${CXX_ID})
SET(CFS_CXX_COMPILER_VER "${CXX_VERSION}")
SET(CFS_CXX_COMPILER_GNU_VER "${CXX_GCC_VERSION}")

#-------------------------------------------------------------------------------
# Now let's find out info about Fortran compiler.
#-------------------------------------------------------------------------------
SET(IDENTIFY_COMPILER_SRC "IdentifyFortranCompiler.F90")
SET(ID_FORTRAN "${CFS_BINARY_DIR}/share/scripts/identify_fortran.cmake")
CONFIGURE_FILE("${IDCOMP_TEMPL}" "${ID_FORTRAN}" @ONLY)

EXECUTE_PROCESS(
  COMMAND "${CMAKE_COMMAND}" -P "${ID_FORTRAN}"
  WORKING_DIRECTORY "${CFS_BINARY_DIR}/tmp"
  RESULT_VARIABLE RETVAL
  )

INCLUDE("${COMPILER_ID_FILE}")

#-------------------------------------------------------------------------------
# Set the Fortran compiler name and compiler version
#-------------------------------------------------------------------------------
SET(CFS_FORTRAN_COMPILER_NAME ${FC_ID})
SET(CFS_FORTRAN_COMPILER_VER "${FC_VERSION}")

#-------------------------------------------------------------------------------
# Check if compiler supports OpenMP
#-------------------------------------------------------------------------------
FIND_PACKAGE(OpenMP)

IF(OPENMP_FOUND)

  #-----------------------------------------------------------------------------
  # The USE_OPENMP option triggers the usage of OpenMP versions of external
  # libraries and the compilation of CFS++ with OpenMP compiler switches.
  #-----------------------------------------------------------------------------
  SET(USE_OPENMP "${USE_OPENMP_DEFAULT}" CACHE BOOL "Turn on support for OpenMP. Needs GCC >= 4.2 or Intel C++, Clang 3.7 is not mature yet.")
  
  #-----------------------------------------------------------------------------
  # Check if compiler has OpenMP support. GCC >= 4.2 has.
  #-----------------------------------------------------------------------------
  IF(USE_OPENMP)
    SET(CFS_C_FLAGS "${OpenMP_C_FLAGS}")
    SET(CFS_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
    # MESSAGE("Use OpenMP-Flags ${OpenMP_C_FLAGS}")
    # sets to -qopenmp for icc and -fopenmp for gcc
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
# Check if we are using the GNU C++ compiler
#-------------------------------------------------------------------------------
IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC" OR CFS_CXX_COMPILER_NAME STREQUAL "CLANG")

  # MESSAGE("We are using the GNU C++ compiler. ${CMAKE_CXX_COMPILER}")
  SET(CFS_CXX_FLAGS "-DBOOST_SYSTEM_NO_DEPRECATED=1 ${CFS_CXX_FLAGS}")
  
  # Obtain major version number of GCC or Clang
  STRING(REPLACE "." ";" CFS_CXX_COMPILER_VER_LIST ${CFS_CXX_COMPILER_VER})
  LIST(GET CFS_CXX_COMPILER_VER_LIST 0 CFS_CXX_COMPILER_MAJOR_VER)

  # we assue C++11 for CFS for any compiler
  SET(CFS_CXX_FLAGS "-std=c++11 -Wno-error=unused-variable -DBOOST_NO_AUTO_PTR ${CFS_CXX_FLAGS}")
  SET(CFS_C_FLAGS "-std=c11")

  IF(CFS_CXX_COMPILER_NAME STREQUAL "CLANG")
     # Fix a problem in Boost 1.5.2 bimap.
     SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-redeclared-class-member")
  ENDIF()
  
  #-----------------------------------------------------------------------------
  # Determine compiler/linker flags according to build type
  #-----------------------------------------------------------------------------
  IF(DEBUG)
    SET(CFS_C_FLAGS "-Wall -fmessage-length=0 ${CFS_C_FLAGS}")
    # -Wold-style-cast Warnings about old C style casts. Since external libraries
    # make extensive use of it, we switch it off. To filter out the warnings in our own
    # code a command line like the following might be used
    # fgrep 'warning: use of old-style cast' out.txt | grep CFS_SOURCE_DIR | sort -u > old-style-cast.txt
    # 
    SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wall -ftemplate-depth-100")

    # -frounding-math: is needed for CGAL library
    IF(USE_CGAL)
      SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -frounding-math")
    ENDIF(USE_CGAL)

    SET(CHECK_MEM_ALLOC 1)
  ELSE()
    # release!
    SET(CFS_C_FLAGS "-Wall -fmessage-length=0 ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wall -ftemplate-depth-100")

    IF(CFS_ARCH STREQUAL "X86_64" AND NOT MINGW)
      SET(CFS_OPT_FLAGS "-m64 -march=k8 -msse2")
    ENDIF()
     
  ENDIF() # end debug/release
  
  #-----------------------------------------------------------------------------
  # Disable some annoying warnings.
  #-----------------------------------------------------------------------------
  SET(CFS_SUPPRESSIONS "-Wno-long-long -Wno-unknown-pragmas -Wno-comment -Wno-strict-aliasing -Wno-deprecated")
  IF(CFS_CXX_COMPILER_VER MATCHES "4.8" OR CFS_CXX_COMPILER_VER VERSION_GREATER "4.8")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-unused-local-typedefs") 
  ENDIF()

  SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-attributes")

  IF(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # required for boost:  error: unused typedef 'boost_static_assert_typedef_890
    # also boost: /include/boost/bimap/support/iterator_type_by.hpp:128:1: error: class member cannot be redeclared 
    # ResultHandler.cc: error: expression with side effects will be evaluated despite being used as an operand to 'typeid' "if( typeid(*fct) == typeid(FieldCoefFunctor<Double>"
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-overloaded-virtual -Wno-unused-local-typedefs -Wno-redeclared-class-member -Wno-potentially-evaluated-expression")

    STRING(TOUPPER "${CMAKE_CXX_COMPILER_ID}" CFS_CXX_COMPILER_NAME)
    SET(CFS_CXX_COMPILER_VER ${CMAKE_CXX_COMPILER_VERSION})
  ENDIF()

  IF(APPLE)
    IF(NOT CMAKE_CROSSCOMPILING)
      #-------------------------------------------------------------------------
      # The C++ linker creates compact exception stack unwinding data since
      # MacOS X 10.6. This can cause page long linker warnings, which cannot be
      # deactivated. So we disable it here alltogether (cf. man unwinddump).
      #-------------------------------------------------------------------------
      SET(CFS_LINKER_FLAGS "-Wl,-no_compact_unwind")
    ENDIF()

    # Check wether the compiler has the -sysroot= flag.
    SET(CXX_HAS_SYSROOT_FLAG_SOURCE "
int
main ()
{
  return 0;
}
")
    SET(CMAKE_REQUIRED_FLAGS "-sysroot=${CMAKE_OSX_SYSROOT}")
    CHECK_CXX_SOURCE_COMPILES("${CXX_HAS_SYSROOT_FLAG_SOURCE}" CXX_HAS_SYSROOT_FLAG)
    UNSET(CMAKE_REQUIRED_FLAGS)
  ENDIF()

  IF(COVERAGE)
    SET(CFS_C_FLAGS "-fprofile-arcs -ftest-coverage ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-fprofile-arcs -ftest-coverage ${CFS_CXX_FLAGS}")
    SET(CFS_LINKER_FLAGS "-fprofile-arcs -ftest-coverage ${CFS_LINKER_FLAGS}")
  ENDIF(COVERAGE)
  
  IF(NOT USE_OPENMP)
    SET(CFS_C_FLAGS "-Werror -Wcomment ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-Werror -Wcomment ${CFS_CXX_FLAGS}")
  ENDIF(NOT USE_OPENMP)

  IF(NOT USE_CGAL)
    SET(CFS_C_FLAGS "-pedantic ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-pedantic ${CFS_CXX_FLAGS}")
  ELSE()
    SET(CFS_C_FLAGS "-frounding-math ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-frounding-math ${CFS_CXX_FLAGS}")
  ENDIF()
# end gcc/clang section
ELSEIF(MSVC)
  #-------------------------------------------------------------------------------
  # Check for Visual Studio C++ compiler
  #-------------------------------------------------------------------------------

  include(CMakeDetermineVSServicePack)
  DetermineVSServicePack( CFS_MSVC_SERVICE_PACK )
  string(TOUPPER "${CFS_MSVC_SERVICE_PACK}" CFS_MSVC_SERVICE_PACK)
  SET(CFS_MSVC_SERVICE_PACK "MS${CFS_MSVC_SERVICE_PACK}")
  
  if( CFS_MSVC_SERVICE_PACK )
    message(STATUS "Detected: ${CFS_MSVC_SERVICE_PACK}")
  endif()

  # Disable some warnings. For details google for 'MSDN C/C++ Build Errors'.

  # For details google for 'MSDN Checked Iterators', '_SCL_SECURE_NO_WARNINGS'
  # and 'MSDN Security Enhancements in the CRT', _CRT_SECURE_NO_WARNINGS
  SET(CFS_CXX_FLAGS "/wd4996")

  # 'identifier' : macro redefinition
  # The macro identifier is defined twice. The compiler uses the second
  # macro definition.
  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /wd4005")

  # '%$S': virtual function overrides '%$pS', previous versions of the
  # compiler did not override when parameters only differed by
  # const/volatile qualifiers
  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /wd4373")
 


#  MESSAGE(FATAL_ERROR "CFS_CXX_COMPILER_NAME ${CFS_CXX_COMPILER_NAME} Compiler")

#  MESSAGE(FATAL_ERROR "MSVC Compiler")
#-------------------------------------------------------------------------------
# Check for Intel C++ compiler
#-------------------------------------------------------------------------------
ELSEIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  #-----------------------------------------------------------------------------
  # Determine compiler/linker flags according to build type
  #-----------------------------------------------------------------------------
  IF(DEBUG)
    SET(CFS_C_FLAGS "-g -ansi -w1 -Wcheck -Werror ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-g -ansi -w1 -Wcheck -Werror ${CFS_CXX_FLAGS}")
    SET(CHECK_MEM_ALLOC 1)
  ELSE()
    # release case
    SET(CFS_C_FLAGS "-ansi -w0 -Werror ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-ansi -w0 -Werror ${CFS_CXX_FLAGS}")
    SET(CFS_SUPPRESSIONS "-wd1125,654,980 -Wno-unknown-pragmas -Wno-comment")
  ENDIF()
  
  #---------------------------------------------------------------------------
  # Disable warnings about hidden overriden functions of base classes,
  # unknown pragmas (openmp, etc.) and multiline comments.
  # Additionally the following warnings / remarks are suppressed:
  #    191: type qualifier is meaningless on cast type
  #    279: controlling expression is constant
  #    654: overloaded virtual function 
  #   1125: entity-kind "entity" is hidden by "entity" -- virtual function 
  #         override intended?
  #   1170: invalid redeclaration of nested class
  #   2259: non-pointer conversion from "type" to "type" may lose significant bits
  #---------------------------------------------------------------------------
  SET(CFS_SUPPRESSIONS "-wd191,279,654,1125,1170,2259")
  SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-unknown-pragmas")
  SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-comment")

  #---------------------------------------------------------------------------
  # The  following flags and  defines are  necessary due  to incompatibilities
  # with  some versions  of the  GCC runtime  environment. The  problems often
  # occur in external libraries, most notably Boost. The favored policy should
  # be,  to fix  those problems  locally  by patching  the external  libraries
  # instead of  introducing global  flags and defines  for CFS++,  which might
  # break other stuff.
  #---------------------------------------------------------------------------
  IF(CFS_CXX_COMPILER_VER MATCHES "11\\.")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -fno-builtin-std::basic_istream::get")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -fno-builtin-std::max")
  ENDIF(CFS_CXX_COMPILER_VER MATCHES "11\\.")
  
  #---------------------------------------------------------------------------
  # The  intel  compiler might  not  know  the  function __builtin_isnan  (and
  #  isinf), so redirect that to isnan
  #---------------------------------------------------------------------------
  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -D__builtin_isnan=::isnan -D__builtin_isinf=::isinf")
  # end icc section
ELSEIF(CFS_CXX_COMPILER_NAME STREQUAL "OPEN64")
  IF(NOT USE_CGAL)
    SET(CFS_C_FLAGS "-pedantic ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-pedantic ${CFS_CXX_FLAGS}")
  ELSE(NOT USE_CGAL)
    SET(CFS_C_FLAGS "-mieee-fp -fp-accuracy=strict -DCGAL_DISABLE_ROUNDING_MATH_CHECK ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-mieee-fp -fp-accuracy=strict -DCGAL_DISABLE_ROUNDING_MATH_CHECK ${CFS_CXX_FLAGS}")
  ENDIF()
ENDIF() # close all CXX compiler specific blocks

# common for all compilers
IF(PROFILING)
 SET(CFS_PROF_FLAGS "-g -fno-omit-frame-pointer")
ENDIF()
    


#-------------------------------------------------------------------------------
# Check for Intel Fortran compiler
#-------------------------------------------------------------------------------
IF(CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")
  #-----------------------------------------------------------------------------
  # Set Intel Fortran library paths in dedicated variables
  #-----------------------------------------------------------------------------
  STRING(REGEX REPLACE "bin/ifort" "lib" IFORT_LIB_PATH "${CMAKE_Fortran_COMPILER}")

  #-----------------------------------------------------------------------------
  # Paths for libraries and compilers are different for version 11 compilers.
  # There is even a difference between Linux and Mac OS X! Thank you Intel!!!
  #-----------------------------------------------------------------------------
  IF(CFS_FORTRAN_COMPILER_VER MATCHES "11\\." AND NOT CFS_DISTRO STREQUAL "MACOSX")
    STRING(REGEX REPLACE "bin/(.*)/ifort" "lib/\\1" IFORT_LIB_PATH "${CMAKE_Fortran_COMPILER}")
  ENDIF()

  LINK_DIRECTORIES(${IFORT_LIB_PATH})
    
#  SET(CFS_FORTRAN_DYNRT_LIBS "ifcoremt;imf;dl")
#  IF(NOT CFS_ARCH STREQUAL "IA64")
#    SET(CFS_FORTRAN_DYNRT_LIBS 
#      "svml"
#      ${CFS_FORTRAN_DYNRT_LIBS})
#  ENDIF(NOT CFS_ARCH STREQUAL "IA64")
#  SET(CFS_FORTRAN_STATRT_LIBS
#    "ifcoremt_pic"
#    "irc"
#    )
ENDIF(CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")


#-------------------------------------------------------------------------------
# Check if we have a valid combination of C++ and Fortran compilers.
#-------------------------------------------------------------------------------
IF(CFS_CXX_COMPILER_NAME STREQUAL "" OR  CFS_FORTRAN_COMPILER_NAME STREQUAL "")
  MESSAGE(SEND_ERROR "Combination of C++ and Fortran compilers not supported!")
ENDIF()

#-------------------------------------------------------------------------------
# Set compiler/linker flags for all build types
#-------------------------------------------------------------------------------
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFS_C_FLAGS} ${CFS_PROF_FLAGS} ${CFS_OPT_FLAGS} ${CFS_SUPPRESSIONS}")
STRING(STRIP "${CMAKE_C_FLAGS}" CMAKE_C_FLAGS)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CFS_CXX_FLAGS} ${CFS_PROF_FLAGS} ${CFS_OPT_FLAGS} ${CFS_SUPPRESSIONS}")
STRING(STRIP "${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CFS_PROF_FLAGS} ${CFS_LINKER_FLAGS}")
STRING(STRIP "${CMAKE_EXE_LINKER_FLAGS}" CMAKE_EXE_LINKER_FLAGS)
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CFS_PROF_FLAGS} ${CFS_LINKER_FLAGS}")
STRING(STRIP "${CMAKE_MODULE_LINKER_FLAGS}" CMAKE_MODULE_LINKER_FLAGS)
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CFS_PROF_FLAGS} ${CFS_LINKER_FLAGS}")
STRING(STRIP "${CMAKE_SHARED_LINKER_FLAGS}" CMAKE_SHARED_LINKER_FLAGS)

# MESSAGE("C++ name: ${CFS_CXX_COMPILER_NAME}")
# MESSAGE("C++ version: ${CFS_CXX_COMPILER_VER}")
# MESSAGE("FORTRAN Compiler ${CFS_FORTRAN_COMPILER_NAME}")
# MESSAGE("FORTRAN Compiler version ${CFS_FORTRAN_COMPILER_VER}")
# MESSAGE("CMAKE BUILD TYPE: ${CMAKE_BUILD_TYPE}")
# MESSAGE("CFS_ARCH_STR: ${CFS_ARCH_STR}")

