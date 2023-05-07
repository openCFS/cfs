include(CheckCXXSourceCompiles) # exotic cmake service function for macOS

# for compiler id and version we use CMAKE_X_COMPILER_ID and CMAKE_X_COMPILER_VERSION which is (for copy & paste)
# CMAKE_C_COMPILER_ID, CMAKE_C_COMPILER_VERSION, 
# CMAKE_CXX_COMPILER_ID, CMAKE_CXX_COMPILER_VERSION
# CMAKE_Fortran_COMPILER_ID, CMAKE_Fortran_COMPILER_VERSION

# see https://cmake.org/cmake/help/v3.25/variable/CMAKE_LANG_COMPILER_ID.html
# note that there is "Clang" and "AppleClang", so use CMAKE_CCX_COMPILER_ID MATCHES "Clang" which identifies both


IF(WIN32)
  # we currently support 
  #    Visual Studio 2017 or 2019 plus Intel Fortran (2019/2020)
  # or
  #    Intel 2020 C++ and Fortran
  # the later requires either VC2017 or VC2019 environment
  # in order to build boost correctly, we need to know which env is used
  SET(CFS_MSVC_VERSION $ENV{VisualStudioVersion}) # TODO: shall be replaced by CMAKE_X_COMPILER_VERSION
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Prepare for macOS
#-------------------------------------------------------------------------------
# pracically one needs brew to compile cfs on macOS. The base is system dependent
if(CFS_ARCH MATCHES "ARM64") # don't care for non-Apple - no harm
  set(CFS_BREW_BASE "/opt/homebrew")
else()
  set(CFS_BREW_BASE "/usr/local")
endif()    

#-------------------------------------------------------------------------------
# Check if compiler supports OpenMP
#-------------------------------------------------------------------------------
find_package(OpenMP)

if(USE_OPENMP)
  set(CFS_C_FLAGS "${OpenMP_C_FLAGS}")
  if(WIN32)
    cmake_print_variables(OpenMP_CXX_FLAGS)
  endif()  
  set(CFS_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CFS_CXX_FLAGS}")
  if(APPLE)
    # best is to use homebrew llvm: https://stackoverflow.com/questions/43555410/enable-openmp-support-in-clang-in-mac-os-x-sierra-mojave
    include_directories(AFTER SYSTEM "${CFS_BREW_BASE}/include")
    set(CFS_LINKER_FLAGS "${CFS_LINKER_FLAGS} -lomp -L${CFS_BREW_BASE}/lib ")
  endif()   
endif() # USE_OPENMP

# Clang can be UNIX (macOS as AppleClang or Linux) or Windows (MSVC bundled).
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # we assue C++14 for CFS for any compiler (including icc below)
  set(CFS_CXX_FLAGS "-std=c++14 -Wuninitialized -Wno-error=unused-variable -Wno-error=maybe-uninitialized -DBOOST_NO_AUTO_PTR ${CFS_CXX_FLAGS}")
  set(CFS_C_FLAGS "-std=c11")

  if(WIN32)
   set(CFS_CXX_FLAGS " -D_WIN32_WINNT=0x0A00")
  endif()

  IF(CFS_FSANITIZE)
    SET(CFS_CXX_FLAGS " -fsanitize=address ${CFS_CXX_FLAGS}")
    SET(CFS_C_FLAGS " -fsanitize=address ${CFS_C_FLAGS}")
  ENDIF()
  
  #-----------------------------------------------------------------------------
  # Determine compiler/linker flags according to build type
  #-----------------------------------------------------------------------------
  if(DEBUG)
    set(CFS_C_FLAGS "-Wall -fmessage-length=0 ${CFS_C_FLAGS}")
    # -Wold-style-cast Warnings about old C style casts. Since external libraries
    # make extensive use of it, we switch it off. To filter out the warnings in our own
    # code a command line like the following might be used
    # fgrep 'warning: use of old-style cast' out.txt | grep CFS_SOURCE_DIR | sort -u > old-style-cast.txt
    #
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wall -ftemplate-depth-100")

    # -frounding-math: is needed for CGAL library
    if(USE_CGAL)
      set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -frounding-math")
    endif()  

    set(CHECK_MEM_ALLOC 1)
  else()
    # release!
    set(CFS_C_FLAGS "-Wall -fmessage-length=0 ${CFS_C_FLAGS}")
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wall -ftemplate-depth-100")

    # the CFS_OPT_FLAGS are also used in FindCFSDEPS.cmake for CFSDEPS_*_FLAGS
    if(CFS_NATIVE)
      # -m64 -> 32 bit int and 64 bit pointers and long
      # further candidates: https://developer.amd.com/wordpress/media/2020/04/Compiler%20Options%20Quick%20Ref%20Guide%20for%20AMD%20EPYC%207xx2%20Series%20Processors.pdf
      # -march=native - has up to 20% boost against nonsense k8 and slight different results on some tests
      # -Ofast Maximize performance - almost no additional effect
      # -funroll-all-loops Enable unrolling - almost no additional effect
      # -flto Link time optimization - extremely slow linking on gcc with almost no effect
      # -param prefetch-latency=300 - Generate memory preload in structions - negative effect
      # some more for AMD Clang in the link above
      set(CFS_OPT_FLAGS "-march=native -Ofast ") 
    else()
      set(CFS_OPT_FLAGS "-m64 -O3 ")
    endif()
  endif() # end debug/release

  # Disable some annoying warnings.
  # note we have at least gcc 4.8
  # -Wno-overflow because of /boost/iostreams/filter/gzip.hpp:674:13: error: overflow in implicit constant conversion [-Werror=overflow]
  # -Wno-parentheses because of boost 1.67: warning: unnecessary parentheses in declaration of 'assert_not_arg' [-Wparentheses]
  SET(CFS_SUPPRESSIONS "-Wno-long-long -Wno-unknown-pragmas -Wno-comment -Wno-parentheses -Wno-unused-function ")
  # stuff that was removed with boost 1.67: -Wno-strict-aliasing -Wno-deprecated -Wno-attributes -Wno-unused-local-typedefs -Wno-overflow

  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0) # there is no >= and also there is no 5.0.0.0
    # /home/fwein/code/trunk/cfs/debug/include/boost/archive/detail/iserializer.hpp:65:1: error: this use of "defined" may not be portable [-Werror=expansion-to-defined]
     #if ! DONT_USE_HAS_NEW_OPERATOR
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-address -Wno-error=address -Wno-expansion-to-defined ")
  ENDIF()

  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0)
    # -Wno-misleading-indentation -Wno-error=placement-new are for gcc 6.1.1 and boost 1.61 maybe remove when a newer boost ist available!
    # however cfs has linkin problems with boost 1.61 hence we have also -Wno-address for boost 1.58
    # we must not set this to CFS_SUPRESSIONS because these also become CMAKE_C_FLAGS and then the following happens:
    # CheckFortranRuntime.cmake (CFS) ->  FortranCInterface.cmake (system) -> Detect.cmake (system) -> try_compile(FortranCInterface_COMPILED
    # this calls a C test (cfs/BUILD/CMakeFiles/FortranCInterface -> make) and reports "command line option *** is valid for C++ but not for C"
    # for debug with -Werror this fails and as a result Fortran name mangling does not work (BUILD/include/def_cfs_fortran_interface.hh is empty)
    SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-misleading-indentation -Wno-placement-new")
  ENDIF()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7.0)
    # on macOS with gcc-7.1 
    # /include/boost/archive/detail/iserializer.hpp:208:9: error: this use of "defined" may not be portable  #if DONT_USE_HAS_NEW_OPERATOR
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-expansion-to-defined ")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0)
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-stringop-truncation ")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 10.0)
    # to prevent boost error with gcc 11:  ‘this’ pointer is null [-Werror=nonnull]
    # try to remove when boost ist updated beyond 1.73
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-nonnull ")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 11.0)
    # to prevent boost error with gcc 11:  thread_data.hpp: error: comparison of integer expressions of different signedness:
    # try to remove when boost ist updated beyond 1.73
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-sign-compare ")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 12.0)
    # to prevent error with gcc 12: struct std::unary_function’ is deprecated
    # unary_function is used in cfs and in boost. Shall be removed in cfs when boost is updated
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-deprecated-declarations ")
    # gcc13 complains this about CGAL
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-catch-value -Wno-dangling-reference ")
    if(DEBUG)
      # gcc13 spams the console with this output but it might be worse to check first if some
      # warnings make sense - but disable for debug
      set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-overloaded-virtual ")
    endif() 
  endif()

  # most specific -Wno-error= are for plain old boost and gcc >= 6. Check to skip them for newer boost than 1.58
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # required for boost:  error: unused typedef 'boost_static_assert_typedef_890
    # also boost: /include/boost/bimap/support/iterator_type_by.hpp:128:1: error: class member cannot be redeclared
    # ResultHandler.cc: error: expression with side effects will be evaluated despite being used as an operand to 'typeid' "if( typeid(*fct) == typeid(FieldCoefFunctor<Double>"
    # include/boost/smart_ptr/detail/sp_counted_base_clang.hpp: warning: '_Atomic' is a C11 extension
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-overloaded-virtual -Wno-redeclared-class-member -Wno-potentially-evaluated-expression -Wno-c11-extensions")
    # -Wno-constant-conversion: boost/iostreams/filter/gzip.hpp:674:16: error: implicit conversion from 'const int' to 'char' changes value from 139 to -117
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-constant-conversion")
    # include/muParserBytecode.h:51:7: error: anonymous types declared in an anonymous union are an extension
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-nested-anon-types")
    # not all gcc options are compatible with clang (mac)
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-unknown-warning-option")
    # boost 1.73 'unary_function<bool, unsigned long>' is deprecated
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Wno-deprecated-declarations")
  ENDIF()

  if(UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # no need for AppleClang
    # we are in the strict Linux case
    # clang 15 on openSUSE would otherwise not compile sniot
    set(CFS_LINKER_FLAGS "${CFS_LINKER_FLAGS$} -no-pie")
  endif()

  if(APPLE)
    # Check wether the compiler has the -sysroot= flag.
    # note that in the first run this is executed before distro.cmake where the architecture might be defined!
    set(CXX_HAS_SYSROOT_FLAG_SOURCE "int main() { return 0; }")
    set(CMAKE_REQUIRED_FLAGS "-sysroot=${CMAKE_OSX_SYSROOT}")
    CHECK_CXX_SOURCE_COMPILES("${CXX_HAS_SYSROOT_FLAG_SOURCE}" CXX_HAS_SYSROOT_FLAG)
    unset(CMAKE_REQUIRED_FLAGS)

    # practically we need homebrew for libquadmath (on x86_64) we might already have this from openmp
    if(NOT USE_OPENMP)
      set(CFS_LINKER_FLAGS "${CFS_LINKER_FLAGS} -L${CFS_BREW_BASE}/lib")  
    endif()

    # warning! don't do -Wl,-nowarn_compact_unwind to prevent unwind warnings! This kills exception catching!
  endif() # apple

  # see https://en.wikipedia.org/wiki/Gcov
  IF(CFS_COVERAGE)
    SET(CFS_C_FLAGS "-fprofile-arcs -ftest-coverage ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-fprofile-arcs -ftest-coverage ${CFS_CXX_FLAGS}")
    SET(CFS_LINKER_FLAGS "-fprofile-arcs -ftest-coverage ${CFS_LINKER_FLAGS}")
  ENDIF()

  IF(NOT USE_OPENMP)
    IF(NOT USE_PHIST_EV OR USE_PHIST_CG)
      SET(CFS_C_FLAGS "-Werror ${CFS_C_FLAGS}")
      SET(CFS_CXX_FLAGS "-Werror ${CFS_CXX_FLAGS}")
    ENDIF()
  ENDIF()


  IF(USE_CGAL) # does CGAL really use C?
    SET(CFS_C_FLAGS "-frounding-math ${CFS_C_FLAGS}")
    SET(CFS_CXX_FLAGS "-frounding-math ${CFS_CXX_FLAGS}")
  ENDIF()
# end gcc/clang section
ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  #-------------------------------------------------------------------------------
  # Check for Visual Studio C++ compiler
  #-------------------------------------------------------------------------------

  # Disable some warnings. For details google for 'MSDN C/C++ Build Errors'.

  # For details google for 'MSDN Checked Iterators', '_SCL_SECURE_NO_WARNINGS'
  # and 'MSDN Security Enhancements in the CRT', _CRT_SECURE_NO_WARNINGS
  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /wd4996")

  # 'identifier' : macro redefinition
  # The macro identifier is defined twice. The compiler uses the second
  # macro definition.
  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /wd4005")

  # '%$S': virtual function overrides '%$pS', previous versions of the
  # compiler did not override when parameters only differed by
  # const/volatile qualifiers
  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /wd4373")

  # prevent warning about data loss when converting size_t
  set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /wd4267")

  # this is for WINDOWS 10. It is essential to have the same version for the boost libs see External_Boost.cmake
  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /D_WIN32_WINNT=0x0A00")
  
  # support for alternative tokens requires the following
  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /permissive- /Zc:twoPhase-")

  # This flag is mandatory on any Winwdows variant, otherwise you get cannot find "libboost_serialization-clangw16-mt-x64-1_81.lib"
  # see https://www.boost.org/doc/libs/1_81_0/libs/config/doc/html/index.html
  # issue: https://discourse.cmake.org/t/linking-on-windows-requires-versioned-boost-libraries/7119/3
  set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /DBOOST_ALL_NO_LIB")

# check for Intel oneAPI llvm based compiler
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM") # Windows (icx) or UNIX (icpx). Interface seems compatible
  # mandatoy for Windows, correct but not necessary on UNIX
  set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -DBOOST_ALL_NO_LIB")

  if(WIN32)
   set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -D_WIN32_WINNT=0x0A00")
  endif()

  # also icx on Windows with MSVC command line interface seems to understand gcc style
  set(CFS_SUPPRESSIONS "-Wno-overloaded-virtual -Wno-deprecated-declarations -Wno-comment ")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 2023)
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-enum-constexpr-conversion -Wno-deprecated-builtins ")
  endif()

# Check for Intel C++ compiler (classic compiler) - to be depreciated mid 2023
ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  #-----------------------------------------------------------------------------
  # Determine compiler/linker flags according to build type
  #-----------------------------------------------------------------------------
  IF(UNIX)
    IF(DEBUG)
      SET(CFS_C_FLAGS "-g -c99 -w1 -Wcheck -Werror ${CFS_C_FLAGS}")
      SET(CFS_CXX_FLAGS "-std=c++14 -g -w1 -Wcheck -Werror ${CFS_CXX_FLAGS}")
      SET(CFSDEPS_CXX_FLAGS "-std=c++14 -g ${CFSDEPS_CXX_FLAGS}") # TODO move this to FindCFSDEPS.cmake
      SET(CHECK_MEM_ALLOC 1)
    ELSE()
      # release case
      SET(CFS_C_FLAGS "-c99 -w0 -Werror ${CFS_C_FLAGS}")
      SET(CFS_CXX_FLAGS "-std=c++14 -w0 -Werror ${CFS_CXX_FLAGS}")
      SET(CFSDEPS_CXX_FLAGS "-std=c++14 -w0 ${CFSDEPS_CXX_FLAGS}")
      SET(CFS_SUPPRESSIONS "-wd1125,654,980 -Wno-unknown-pragmas -Wno-comment")
    ENDIF()
  ELSE()
    # this is for WINDOWS 10
    SET(CFS_C_FLAGS "${CFS_C_FLAGS} /D_WIN32_WINNT=0x0A00 /Qstd=c99")
    SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /D_WIN32_WINNT=0x0A00 /DBOOST_ALL_NO_LIB /Qstd=c++14")
    SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} /D_WIN32_WINNT=0x0A00 /Qstd=c++14")
    IF(DEBUG)
      SET(CFS_C_FLAGS "/Z7 /w1 /Wcheck /Werror ${CFS_C_FLAGS}")
      SET(CFS_CXX_FLAGS "/Z7 /W1 /Wcheck /Werror ${CFS_CXX_FLAGS}")
      SET(CFSDEPS_CXX_FLAGS "/Z7 ${CFSDEPS_CXX_FLAGS}")
      SET(CHECK_MEM_ALLOC 1)
    ELSE()
      # release case
      SET(CFS_CXX_FLAGS "/W0 ${CFS_CXX_FLAGS}")
      SET(CFSDEPS_CXX_FLAGS "/W0 ${CFSDEPS_CXX_FLAGS}")
      # remark #10441: The Intel(R) C++ Compiler Classic (ICC) is deprecated and will be removed from product release in the second half of 2023.
      SET(CFS_SUPPRESSIONS "/Qdiag-disable1125,654,980,10441")
    ENDIF()
  ENDIF()

  # on windows deactivate all MS special features
  # should be an easy way to activate use of alternative tokens for intel
  # unfortunaley, this leads to many, many, many, issues using boost
  # --> deactivated......
  #  IF(WIN32)
	  #    SET(CFS_C_FLAGS "${CFS_C_FLAGS} /Za -wd2358")
	  #  SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /Za -wd2358")
	  #    SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} /Za -wd2358")
  #  ENDIF()

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
  IF(UNIX)
    SET(CFS_SUPPRESSIONS "-wd191,279,654,1125,1170,2259")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-unknown-pragmas -Wno-comment")
  ENDIF()

  #---------------------------------------------------------------------------
  # The  following flags and  defines are  necessary due  to incompatibilities
  # with  some versions  of the  GCC runtime  environment. The  problems often
  # occur in external libraries, most notably Boost. The favored policy should
  # be,  to fix  those problems  locally  by patching  the external  libraries
  # instead of  introducing global  flags and defines  for openCFS,  which might
  # break other stuff.
  #---------------------------------------------------------------------------
  IF(CMAKE_CXX_COMPILER_VERSION MATCHES "11\\.")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -fno-builtin-std::basic_istream::get")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -fno-builtin-std::max")
  ENDIF(CMAKE_CXX_COMPILER_VERSION MATCHES "11\\.")
  # Open64 compiler support removed. See svn version 15997
ENDIF() # close all CXX compiler specific blocks

# common for all compilers
# adds debug information to the code such that vtune, valgrind, ... can show the lines of the hotspots
# this is different from adding gprof support by -pg wich adds changes the code to generate an output file
IF(CFS_PROFILING)
 SET(CFS_PROF_FLAGS "-g -fno-omit-frame-pointer")
ENDIF()



#-------------------------------------------------------------------------------
# Check for Intel Fortran compiler
#-------------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  #-----------------------------------------------------------------------------
  # Set Intel Fortran library paths in dedicated variables
  #-----------------------------------------------------------------------------
  STRING(REGEX REPLACE "bin/ifort" "lib" IFORT_LIB_PATH "${CMAKE_Fortran_COMPILER}")

  #-----------------------------------------------------------------------------
  # Paths for libraries and compilers are different for version 11 compilers.
  # There is even a difference between Linux and Mac OS X! Thank you Intel!!!
  #-----------------------------------------------------------------------------
  IF(CMAKE_Fortran_COMPILER_VERSION MATCHES "11\\." AND NOT CFS_DISTRO STREQUAL "MACOSX")
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

  IF(WIN32)
   
    # TODO: this block could be doubled with FindIntelMKL.cmake ! 
    
    #-----------------------------------------------------------------------------
    # Copy compiler runtime .dlls to bin directory
    #-----------------------------------------------------------------------------
    SET(INTEL_DLLS
        libifcoremd.dll
        libifportmd.dll
        libiomp5md.dll
        libmmd.dll
        svml_dispmd.dll
    )
    IF(DEBUG)
      SET(INTEL_DLLS 
        ${INTEL_DLLS}
        libifcoremdd.dll
        libmmdd.dll
      )
    ENDIF()
    
    SET(LIB_DEST_DIR "${CFS_BINARY_DIR}/bin/")
    GET_FILENAME_COMPONENT(INTEL_COMPILER_DIR ${CMAKE_Fortran_COMPILER} PATH)    

	  set(ICC_REDIST_DIR "${INTEL_COMPILER_DIR}/../../redist/intel64_win/compiler/")
    # for  parallel studio pre oneApi was set(ICC_REDIST_DIR "${INTEL_COMPILER_DIR}/../../redist/intel64/compiler/")
    
    MESSAGE(STATUS "Copying INTEL redistributable files from ${ICC_REDIST_DIR} to ${LIB_DEST_DIR}")
    FOREACH(lib IN LISTS INTEL_DLLS)
      FILE(COPY ${ICC_REDIST_DIR}/${lib} DESTINATION ${LIB_DEST_DIR})
    ENDFOREACH()
  ENDIF(WIN32)
endif() # Intel

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

if(WIN32)
  cmake_print_variables(CMAKE_CXX_FLAGS)
  cmake_print_variables(CFS_SUPPRESSIONS)
endif()
