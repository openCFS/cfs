include(CheckCXXSourceCompiles) # exotic cmake service function for macOS

# for compiler id and version we use CMAKE_X_COMPILER_ID and CMAKE_X_COMPILER_VERSION which is (for copy & paste)
# CMAKE_C_COMPILER_ID, CMAKE_C_COMPILER_VERSION, 
# CMAKE_CXX_COMPILER_ID, CMAKE_CXX_COMPILER_VERSION
# CMAKE_Fortran_COMPILER_ID, CMAKE_Fortran_COMPILER_VERSION

# see https://cmake.org/cmake/help/v3.25/variable/CMAKE_LANG_COMPILER_ID.html
# note that there is "Clang" and "AppleClang", so use CMAKE_CCX_COMPILER_ID MATCHES "Clang" which identifies both
# there is classic "Intel" (deprecated) and "IntelLLVM" which is also a clang variant
# "GNU" for gcc, g++ and gfortan, "MSVC" for Microsoft Visual Studio  

# we use the following variables: 
# CFS_CXX_FLAGS for openCFS itself without warnings and CFS_OPT_FLAGS. There is no need for CFS_C_FLAGS and CFS_Fortran_FLAGS as we are pure C++
# CFS_SUPPRESSIONS enable all warnings (-Wall) and the suppress what cannot be modified by code or what is from a lib (e.g. boost)
# CFS_OPT_FLAGS this is what we set for release build to CFS_CXX_FLAGS and CFSDEPS_C/CXX_FLAGS
# CFSDEPS_C/CXX/Fortran_FLAGS what we give most cfsdeps. Always release and the openCFS compilers
# CFS_LINKER_FLAGS for openCFS only
# note that CMAKE_CXX_FLAGS in ccmake (often empty) can be used to add own options for a build!

if(DEBUG)
  set(CHECK_MEM_ALLOC 1)
endif()

# make sure openmp_blas.cmake run already
if(USE_OPENMP)
  assert_set(OpenMP_FOUND)
endif()

# Clang can be UNIX (macOS as AppleClang or Linux) or Windows (MSVC bundled).
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CFS_CXX_FLAGS "-std=c++14 -ftemplate-depth-100 -DBOOST_NO_AUTO_PTR ${CFS_CXX_FLAGS}")

  if(USE_CGAL)
    # CGAL seems to require -frounding-math 
    # warning gcc 12+13 crash creating a hdf5 file using hdf5 1.8.20 
    # probably -frounding-math needs to be in sync with all cfsdeps
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -frounding-math") 
  endif() 

  # this is set via ccmake (advanded) and produces very slow code which can check during runtime for memory issues   
  if(CFS_FSANITIZE)
    set(CFS_CXX_FLAGS " -fsanitize=address ${CFS_CXX_FLAGS}")
  endif()
  
  # see https://en.wikipedia.org/wiki/Gcov
  if(CFS_COVERAGE)
    set(CFS_CXX_FLAGS "-fprofile-arcs -ftest-coverage ${CFS_CXX_FLAGS}")
    set(CFS_LINKER_FLAGS "-fprofile-arcs -ftest-coverage ${CFS_LINKER_FLAGS}")
  endif()
  
  # adds debug information to the code such that vtune, valgrind, ... can show the lines of the hotspots
  # this is different from adding gprof support by -pg wich adds changes the code to generate an output file
  if(CFS_PROFILING)
   set(CFS_PROF_FLAGS "-g -fno-omit-frame-pointer")
  endif()  
  
  # we set CFS_OPT_FLAGS for release and CFSDEPS 
  if(CFS_NATIVE)
    # -m64 -> 32 bit int and 64 bit pointers and long
    # further candidates: https://developer.amd.com/wordpress/media/2020/04/Compiler%20Options%20Quick%20Ref%20Guide%20for%20AMD%20EPYC%207xx2%20Series%20Processors.pdf
    # -march=native - has up to 20% boost against nonsense k8 and slight different results on some tests
    # -Ofast Maximize performance - almost no additional effect
    # -funroll-all-loops Enable unrolling - almost no additional effect
    # -flto Link time optimization - extremely slow linking on gcc with almost no effect
    # -param prefetch-latency=300 - Generate memory preload in structions - negative effect
    # some more for AMD Clang in the link above
    set(CFS_OPT_FLAGS "-m64 -march=native -Ofast ") 
  else()
    set(CFS_OPT_FLAGS "-m64 -O3 ")
  endif()
   
  if(UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # no need for AppleClang
    # we are in the strict Linux case
    # clang 15 on openSUSE would otherwise not compile sniot
    set(CFS_LINKER_FLAGS "${CFS_LINKER_FLAGS$} -no-pie")
  endif()

  if(APPLE)
    # linker issue with xcode 15 on mac
    # https://developer.apple.com/forums/thread/735426
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0")
      set(CFS_LINKER_FLAGS "${CFS_LINKER_FLAGS} -ld64") # was -ld_classic and for >= 16 -ld64 is preferred and shall also work with clang 15 
    endif()

    # warning! don't do -Wl,-nowarn_compact_unwind to prevent unwind warnings! This kills exception catching!
    set(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -arch ${CMAKE_OSX_ARCHITECTURES} -isysroot ${CMAKE_OSX_SYSROOT}")
    set(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -arch ${CMAKE_OSX_ARCHITECTURES} -isysroot ${CMAKE_OSX_SYSROOT}")
  endif() # apple

  # some compilers ignore -w, others ignore -Wno-everything - both have the same purpose to disable all warning
  set(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} ${CFS_OPT_FLAGS} -w -Wno-everything")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "14")
    set(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -Wno-int-conversion -Wno-implicit-function-declaration")
  endif()
  set(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} ${CFS_OPT_FLAGS} -w -Wno-everything")
  if(USE_CGAL) # see comment about -frounding-math for gcc 12+13 above
    set(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -frounding-math ")
    set(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -frounding-math ")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "14")
    set(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -Wno-implicit-int") # needed for SuperLU (and maybe others)
  endif()

  # enable all warnings, then disable the ones we cannot prevent (e.g. from lib includes).
  # better is always changes to code, updateing/patching libs, guarding includes by pragmas
  set(CFS_SUPPRESSIONS "-Wall -Wuninitialized -Wno-error=unused-variable -Wno-error=maybe-uninitialized")
  if(debug)
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Werror") # does not allow most warnings -> clang debug pipeline
  endif()
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-long-long -Wno-unknown-pragmas -Wno-comment -Wno-address -Wno-error=address -Wno-unused-function ")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0)
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-stringop-truncation ")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12.0)
    # -Wno-deprecated-declarations: boost 1.78 hash.hpp: struct std::unary_function’ is deprecated 
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-deprecated-declarations -Wno-stringop-overflow -Wno-array-bounds")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13.0)
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-overloaded-virtual ")
    if(USE_CGAL)
      set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-catch-value -Wno-dangling-reference ")
    endif()  
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-overloaded-virtual -Wno-redeclared-class-member -Wno-potentially-evaluated-expression -Wno-c11-extensions")
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-overloaded-virtual -Wno-potentially-evaluated-expression")
    # -Wno-constant-conversion: boost/iostreams/filter/gzip.hpp:674:16: error: implicit conversion from 'const int' to 'char' changes value from 139 to -117
    #set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-constant-conversion")
    # include/muParserBytecode.h:51:7: error: anonymous types declared in an anonymous union are an extension
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-nested-anon-types")
    # from boost 1.78 hash.hpp 'unary_function<const std::error_category *, unsigned long>' is deprecated
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-deprecated-declarations")
    # not all gcc options are compatible with clang (mac)
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-unknown-warning-option")
    
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
      set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-deprecated-builtins") # at least AppleClang
      set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-enum-constexpr-conversion")
    endif()

    # clang >= 18 is missing some standard c++ functions without that flag
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 18)
      set(CFS_LINKER_FLAGS "${CFS_LINKER_FLAGS} -static-libstdc++")
    endif()
  endif() # end Clang   
endif() # CXX GNU+Clang  

# check for Intel oneAPI llvm based compiler
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM") # Windows (icx) or UNIX (icpx). Interface seems compatible
  # BOOST_ALL_NO_LIB is mandatoy for Windows, correct but not necessary on UNIX
  if(WIN32)
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} -Qstd=c++14 -D_WIN32_WINNT=0x0A00  -DBOOST_ALL_NO_LIB /fp:precise")
    set(CFSDEPS_C_FLAGS " -D_WIN32_WINNT=0x0A00 /fp:precise")
    set(CFSDEPS_CXX_FLAGS " -D_WIN32_WINNT=0x0A00 /fp:precise")
  else()
    set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS}  -std=c++14 -fp-model=precise")
    set(CFSDEPS_C_FLAGS " -fp-model=precise")
    set(CFSDEPS_CXX_FLAGS " -fp-model=precise")
  endif()
  
  set(CFS_OPT_FLAGS "-O3")
  
  set(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} ${CFS_OPT_FLAGS} -w -Wno-everything")
  set(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} ${CFS_OPT_FLAGS} -w -Wno-everything")
  if(WIN32)
    # error: cannot use 'throw' with exceptions disabled
    set(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} /EHsc")
  endif() 
  if(USE_CGAL) # remove when we use header only CGAL
    set(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -frounding-math ")
    set(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -frounding-math ")
  endif()

  # also icx on Windows with MSVC command line interface seems to understand gcc style
  set(CFS_SUPPRESSIONS "-Wno-overloaded-virtual -Wno-deprecated-declarations -Wno-comment ")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 2023)
    # to allow typeid(*fct) we need -Wno-potentially-evaluated-expression
    set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-enum-constexpr-conversion -Wno-deprecated-builtins -Wno-potentially-evaluated-expression")
    if(WIN32)
      set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-unused-variable -Wno-unused-private-field -Wno-microsoft-unqualified-friend -Wno-macro-redefined")
    endif()
  endif()
endif() # IntelLLVM

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # this is for WINDOWS 10. It is essential to have the same version for the boost libs see External_Boost.cmake
  set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /D_WIN32_WINNT=0x0A00")
  
  # support for alternative tokens requires the following
  set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /permissive- /Zc:twoPhase-")

  # This flag is mandatory on any Winwdows variant, otherwise you get cannot find "libboost_serialization-clangw16-mt-x64-1_81.lib"
  # see https://www.boost.org/doc/libs/1_81_0/libs/config/doc/html/index.html
  # issue: https://discourse.cmake.org/t/linking-on-windows-requires-versioned-boost-libraries/7119/3
  set(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /DBOOST_ALL_NO_LIB")

  # Disable some warnings. For details google for 'MSDN C/C++ Build Errors'.

  # For details google for 'MSDN Checked Iterators', '_SCL_SECURE_NO_WARNINGS'
  # and 'MSDN Security Enhancements in the CRT', _CRT_SECURE_NO_WARNINGS
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4996")
  # 'identifier' : macro redefinition
  # The macro identifier is defined twice. The compiler uses the second
  # macro definition.
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4005")
  # '%$S': virtual function overrides '%$pS', previous versions of the
  # compiler did not override when parameters only differed by
  # const/volatile qualifiers
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4373")
  # prevent warning about data loss when converting size_t
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4267")
  # prevent warning about data loss when converting size_t
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4267")
  # muparser stuff
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4251")
  # type conversion
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4244")
  # template stuff
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4661")
  # Unreferenzierte lokale Variable (e.g. catch)
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4101")
  # Unbekanntes Pragma "GCC".
  set(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} /wd4068")

  # silencing all warnings with /w for cfsdeps does not work as cmake sets /W4 which overrides /w
  set(CFSDEPS_CXX_FLAGS "${CFS_CXX_FLAGS} /EHsc ${CFS_SUPPRESSIONS} /wd4310")

endif() # MSVC

# Check for Intel C++ compiler (classic compiler) - to be depreciated mid 2023
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  #-----------------------------------------------------------------------------
  # Determine compiler/linker flags according to build type
  #-----------------------------------------------------------------------------
  IF(UNIX)
    IF(DEBUG)
      SET(CFS_CXX_FLAGS "-std=c++14 -g -w0 ${CFS_CXX_FLAGS}")
      SET(CFSDEPS_CXX_FLAGS "-std=c++14 -g ${CFSDEPS_CXX_FLAGS}") 
    ELSE()
      # release case
      SET(CFS_CXX_FLAGS "-std=c++14 -w0 ${CFS_CXX_FLAGS}")
      SET(CFSDEPS_CXX_FLAGS "-std=c++14 -w0 ${CFSDEPS_CXX_FLAGS}")
      SET(CFS_SUPPRESSIONS "-wd1125,654,980 -Wno-unknown-pragmas -Wno-comment")
    ENDIF()
  ELSE()
    # this is for WINDOWS 10
    SET(CFS_CXX_FLAGS "${CFS_CXX_FLAGS} /D_WIN32_WINNT=0x0A00 /DBOOST_ALL_NO_LIB /Qstd=c++14")
    SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} /D_WIN32_WINNT=0x0A00 /Qstd=c++14")
    IF(DEBUG)
      SET(CFS_CXX_FLAGS "/Z7 /W0 ${CFS_CXX_FLAGS}")
      SET(CFSDEPS_CXX_FLAGS "/Z7 ${CFSDEPS_CXX_FLAGS}")
    ELSE()
      # release case
      SET(CFS_CXX_FLAGS "/W0 ${CFS_CXX_FLAGS}")
      SET(CFSDEPS_CXX_FLAGS "/W0 ${CFSDEPS_CXX_FLAGS}")
      # remark #10441: The Intel(R) C++ Compiler Classic (ICC) is deprecated and will be removed from product release in the second half of 2023.
      SET(CFS_SUPPRESSIONS "/Qdiag-disable1125,654,980,10441")
    ENDIF()
  ENDIF()

  IF(UNIX)
    SET(CFS_SUPPRESSIONS "-wd191,279,654,1125,1170,2259")
    SET(CFS_SUPPRESSIONS "${CFS_SUPPRESSIONS} -Wno-unknown-pragmas -Wno-comment")
  ENDIF()
endif() # ends classic Intel C++

# Fortran compilers
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  set(CFSDEPS_Fortran_FLAGS "${CFS_OPT_FLAGS} -w") # -w is like -Wno-everything but recognized by gfortan
endif()  

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU" AND (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)) # gfortran version >= 10 (gcc10+)
  set(CFSDEPS_Fortran_FLAGS "${CFSDEPS_Fortran_FLAGS} -fallow-argument-mismatch") # was once --std=legacy # see https://github.com/Reference-LAPACK/lapack/issues/353
endif()  

if(WIN32 AND ${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel") # ifx and ifort 
   # prevent the following error on arpack - but shall not harm for the other Fortran cfsdeps, too 
   # ucrt.lib(api-ms-win-crt-math-l1-1-0.dll) : error LNK2005: ldexp already defined in libmmt.lib(ldexp_iface_c99.obj)
   set(CFSDEPS_Fortran_FLAGS "/fpp /libs:dll")
   if(USE_OPENMP)                                                                                                                                                                         
      set(CFSDEPS_Fortran_FLAGS "${CFSDEPS_Fortran_FLAGS} /threads ")                                                                                                                  
   endif() 
endif()

# in CheckFortanRuntime.cmake and redistributables.cmake we copy redistributable intel libs for deployment 

#cmake_print_variables(CFSDEPS_C_FLAGS)
#cmake_print_variables(CFSDEPS_CXX_FLAGS)
#cmake_print_variables(CFSDEPS_Fortran_FLAGS)
#cmake_print_variables(CFS_OPT_FLAGS)
#cmake_print_variables(CFS_SUPPRESSIONS)
#cmake_print_variables(CFS_CXX_FLAGS)
#cmake_print_variables(CFS_LINKER_FLAGS)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CFS_CXX_FLAGS} ${CFS_SUPPRESSIONS}") # does not overwrite the cache variable
if(NOT DEBUG) # note that CFS_OPT_FLAGS are set for cfsdeps even for debug as we store the zip without build information
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CFS_OPT_FLAGS}")
endif()
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CFS_LINKER_FLAGS}")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CFS_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CFS_LINKER_FLAGS}")

#cmake_print_variables(CMAKE_CXX_FLAGS)
#cmake_print_variables(CMAKE_EXE_LINKER_FLAGS)



