#-------------------------------------------------------------------------------
# This script  is responsible for the  determination of the correct  include and
# linker  parameters for  the Intel  Math Kernel  Library (MKL).  We distinguish
# between Windows, MacOS and Linux.
#
# On Linux  we determine the  required linker flags by compiling one of the MKL
# examples and reading the flags from the Makefile output (cf. below).
# another option: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#-------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Function to read MKL version from header file mkl.h or mkl_version.h
#-----------------------------------------------------------------------------
function(MKL_VERSION_FROM_HEADER)
  if(EXISTS "${MKL_ROOT_DIR}/include/mkl_version.h")
    file(STRINGS "${MKL_ROOT_DIR}/include/mkl_version.h" MKL_HEADER)
  elseif(EXISTS "${MKL_ROOT_DIR}/include/mkl.h")
    file(STRINGS "${MKL_ROOT_DIR}/include/mkl.h" MKL_HEADER)
  else()
    message(FATAL_ERROR "MKL header could not be found for MKL_ROOT_DIR=${MKL_ROOT_DIR} at include/mkl_version.h or include/mkl.h")
  endif()

  foreach(line IN LISTS MKL_HEADER)
    IF(line MATCHES "#define __INTEL_MKL")
      STRING(REPLACE "_" ";" MKL_LIST "${line}")
      LIST(GET MKL_LIST 4 ITEM_FOUR)

      IF(ITEM_FOUR STREQUAL "")
        LIST(GET MKL_LIST 5 MKL_MAJOR_VERSION)
        STRING(STRIP "${MKL_MAJOR_VERSION}" MKL_MAJOR_VERSION)
      ELSEIF(ITEM_FOUR STREQUAL "MINOR")
        LIST(GET MKL_LIST 6 MKL_MINOR_VERSION)
        STRING(STRIP "${MKL_MINOR_VERSION}" MKL_MINOR_VERSION)
      ELSEIF(ITEM_FOUR STREQUAL "UPDATE")
        LIST(GET MKL_LIST 6 MKL_UPDATE)
        STRING(STRIP "${MKL_UPDATE}" MKL_UPDATE)
      ENDIF()
    ENDIF()
  endforeach()

  SET(MKL_MAJOR_VERSION ${MKL_MAJOR_VERSION} PARENT_SCOPE)
  SET(MKL_MINOR_VERSION ${MKL_MINOR_VERSION} PARENT_SCOPE)
  SET(MKL_UPDATE ${MKL_UPDATE} PARENT_SCOPE)
  IF(${MKL_MAJOR_VERSION} LESS 2019)
    MESSAGE(FATAL_ERROR "MKL version has to be at least 2019.0.0 MKL_VERSION: ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}")
  ENDIF()

endfunction(MKL_VERSION_FROM_HEADER)

if(MSVC)
  #-----------------------------------------------------------------------------
  # If not specified by the user, try to determine proper MKL root directory.
  #-----------------------------------------------------------------------------
  if(NOT MKL_ROOT_DIR)
    set(MKL_POSSIBLE_ROOT_DIRS
      "e:/dev/intel/MKL/composer_xe_2013"
      "e:/dev/intel/MKL/10.0.5.025"
      "$ENV{MKLROOT}"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.4.228/windows/mkl" )

    find_file(MKL_H
      "include/mkl.h"
      PATHS ${MKL_POSSIBLE_ROOT_DIRS}
      DOC "MKL root dir"
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH )

    STRING(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_H}")
    #  MESSAGE(FATAL_ERROR "MKL_ROOT_DIR ${MKL_ROOT_DIR}")
  endif()
 
  # e.g. MKL_ROOT_DIR=C:/Program Files (x86)/Intel/oneAPI/mkl/latest
  assert_dir_exists(${MKL_ROOT_DIR})

  #---------------------------------------------------------------------------
  # Read mkl version from header files
  #---------------------------------------------------------------------------
  MKL_VERSION_FROM_HEADER()
  set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include" CACHE PATH "mkl include dir")
  
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 19.1)
    set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/intel64_win" CACHE PATH "here we assume the mkl libs")
  else()
    set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/intel64" CACHE PATH "here we assume the mkl libs")
  endif()
  assert_dir_exists(${MKL_LIB_DIR})
  set(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL.")

  # mkl_solver* libs are deprecated as of v10
  # see: https://web.archive.org/web/20091027212420/https://software.intel.com/en-us/articles/mkl_solver_libraries_are_deprecated_libraries_since_version_10_2_Update_2/

  set(MKL_BLAS_LIB
    ${MKL_LIB_DIR}/mkl_intel_lp64_dll.lib
    ${MKL_LIB_DIR}/mkl_intel_thread_dll.lib
    ${MKL_LIB_DIR}/mkl_core_dll.lib )
	
  if(USE_OPENMP)	
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 19.1)
      set(LIBOMP ${MKL_ROOT_DIR}/../compiler/lib/intel64_win/libiomp5md.lib) 
      assert_file_exists(${LIBOMP})
    else() 
      # there seems to be bo static openmp-lib ?!
      # set(LIBOMP "${MKL_ROOT_DIR}/../../compiler/latest/windows/redist/intel64_win/compiler/libiomp5md.dll")
    endif()
    list(APPEND MKL_BLAS_LIB ${LIBOMP})
  endif() 	

  set(MKL_LAPACK_LIB ${MKL_BLAS_LIB})

  set(DEPS_SEQUENTIAL
    ${MKL_LIB_DIR}/mkl_intel_lp64.lib
    ${MKL_LIB_DIR}/mkl_intel_thread.lib
    ${MKL_LIB_DIR}/mkl_core.lib
    ${MKL_LIB_DIR}/mkl_sequential.lib
    ${LIBOMP} ) # LIBOMP is empty without openMP and with too new MSVC

  # Copy MKL redistributable dlls to bin/ directory
  
  # redistributable directory
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 19.1)
    set(MKL_REDIST_DIR "${MKL_ROOT_DIR}/../redist/intel64/mkl/")
  else()
    # C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler>
    set(MKL_REDIST_DIR "${MKL_ROOT_DIR}/redist/intel64/")
  endif()
  assert_dir_exists(${MKL_REDIST_DIR})

  set(LIB_DEST_DIR "${CFS_BINARY_DIR}/bin/")
  # copy all dlls to binary dir
  message(STATUS "Copying MKL redistributable files from ${MKL_REDIST_DIR} to ${LIB_DEST_DIR}")
  file(COPY ${MKL_REDIST_DIR} DESTINATION ${LIB_DEST_DIR} FILES_MATCHING PATTERN "*.dll")

# END MSVC
elseif(APPLE) # note tha APPLE is ALSO UNIX! 
  # for APPLE we do it as with MSVC, we forget about the complex LINUX compile_mkl_test.sh.in stuff
  # this is possibly also an option for Linux to get rid of the ugly stuff. 
  # base is the interactive configurator: https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html
  
  if(NOT MKL_ROOT_DIR)
    set(MKL_POSSIBLE_ROOT_DIRS 
       "/opt/intel/mkl"
       $ENV{MKLROOT}) # set by compilervars.sh intel64 

    find_file(MKL_H
      "include/mkl.h"
      PATHS ${MKL_POSSIBLE_ROOT_DIRS}
      DOC "MKL root dir"
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH)

    string(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_H}")
    #message(STATUS "MKL_ROOT_DIR set from include/mkl.h: ${MKL_ROOT_DIR}")
  endif()


  #---------------------------------------------------------------------------
  # Read mkl version from header files
  #---------------------------------------------------------------------------
  MKL_VERSION_FROM_HEADER()
  set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include" CACHE PATH "mkl include directory")
  
  set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib" CACHE PATH "here we assume the mkl libs")

  # https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html
  if(USE_OPENMP)
    set(MKL_THREAD_LIB ${MKL_LIB_DIR}/libmkl_intel_thread.a)
  else()
    set(MKL_THREAD_LIB ${MKL_LIB_DIR}/libmkl_sequential.a)
  endif()
  
  # lp64 is 32bit integer, ilp64 is 64bit integer 

  # -liomp5 is only required in the USE_OPENMP case but it should not harm
  set(MKL_BLAS_LIB 
    ${MKL_LIB_DIR}/libmkl_intel_lp64.a
    ${MKL_THREAD_LIB}
    ${MKL_LIB_DIR}/libmkl_core.a
    -liomp5 
    -lpthread 
    -lm 
    -ldl)

  # for macOS and clang and MKL we need to link -lgcc_s.1.dylib from the libgfortran.a location
  # to prevent undefined references ___divtf3, ___eqtf2, ___gttf2, ..., ___subtf3, ___unordtf2
  # according to nm there is no static lib to provide this (used by libgfortran.a)
  # Edit: at least on macOS there is a static libgcc.a, see GFORTAN_LIBGCC_STATIC in CheckFortranRuntime.cmake
  if(NOT GFORTRAN_LIBRARY)
    message(FATAL_ERROR "neeed path to libgcc_s.1.dylib from GFORTRAN_LIBRARY")
  endif()
  
  get_filename_component(GFORTRAN_LIB_DIR ${GFORTRAN_LIBRARY} DIRECTORY)
  set(MKL_BLAS_LIB
    ${MKL_BLAS_LIB}
    -L${GFORTRAN_LIB_DIR}
    -lgcc_s.1)

  # the path for libimp5 is not set by defeault: LD_LIBRARY_PATH=$MKLROOT/../compiler/lib/ works, 
  # but it is easier to copy the file to the lib-dir.
  file(COPY ${MKL_ROOT_DIR}/../compiler/lib/libiomp5${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
     
  set(MKL_LAPACK_LIB ${MKL_BLAS_LIB})

  set(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL." FORCE)


  mark_as_advanced(MKL_H)
  mark_as_advanced(MKL_ROOT_DIR)

elseif(UNIX AND NOT APPLE) # neither MSVC and neither APPLE. Hence UNIX and Linux. Note that APPLE is also UNIX
  if(NOT MKL_ROOT_DIR)
    set(MKL_POSSIBLE_ROOT_DIRS
       "/opt/intel/mkl"
       "/opt/intel/oneapi/mkl/latest"
       "/share/programs/intel-mkl/latest/mkl"
       $ENV{MKLROOT}) # set by compilervars.sh intel64

    find_file(MKL_H
      "include/mkl.h"
      PATHS ${MKL_POSSIBLE_ROOT_DIRS}
      DOC "MKL root dir"
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH)
    mark_as_advanced(MKL_H)

    string(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_H}")
    message(STATUS "MKL_ROOT_DIR set from include/mkl.h: ${MKL_ROOT_DIR}")
  endif()
  set(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL.")
  mark_as_advanced(MKL_ROOT_DIR)
  
  if(NOT EXISTS ${MKL_ROOT_DIR}/lib/intel64/libmkl_core.a)
    message(FATAL_ERROR "MKL_ROOT_DIR=${MKL_ROOT_DIR} but ${MKL_ROOT_DIR}/lib/intel64/libmkl_core.a is not found") 
  endif()
  set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/intel64" CACHE PATH "here we assume the mkl libs")
  message(STATUS "setting MKL_LIB_DIR=${MKL_LIB_DIR} (validated)")

  #---------------------------------------------------------------------------
  # Read mkl version from header files
  #---------------------------------------------------------------------------
  MKL_VERSION_FROM_HEADER()

  # set include dir
  set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include")

  if(USE_OPENMP)
    set(MKL_OMP_LIB "${MKL_LIB_DIR}/libiomp5.so")
    # strange that we search here for alternative locations of libiomp5.so but not for the other libs later ...
    if(NOT EXISTS "${MKL_OMP_LIB}")
      set(MKL_OMP_LIB "${MKL_ROOT_DIR}/../compiler/lib/intel64_lin/libiomp5.so")
    endif()
    # path for oneAPI 2020.2
    if(NOT EXISTS "${MKL_OMP_LIB}")
      set(MKL_OMP_LIB "${MKL_ROOT_DIR}/../../compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so")
    endif()
    # copy over
    file(COPY ${MKL_OMP_LIB} DESTINATION ${LIBRARY_OUTPUT_PATH})
    set(MKL_OMP_LIB_LINE "-L${LIBRARY_OUTPUT_PATH} -liomp5")
    set(MKL_THREADING_LIB "${MKL_LIB_DIR}/libmkl_intel_thread.a")
  else()
    set(MKL_THREADING_LIB "${MKL_LIB_DIR}/libmkl_sequential.a")
  endif()

  # set the link line, essentailly from
  # see https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html
  # regarting the correct Fortran interface library (libmkl_gf_lp64):
  # see: https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/560573
  # Select interface layer: Fortran API with 32 bit integer (lp64, not ilp64)
  set(MKL_FORTRAN_INTERFACE_LIB "${MKL_LIB_DIR}/libmkl_gf_lp64.a")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message(STATUS "Using Intel Fortran compiler (ifort or ifx): use libmkl_intel_lp64.a")
    set(MKL_FORTRAN_INTERFACE_LIB "${MKL_LIB_DIR}/libmkl_intel_lp64.a")
  endif()
  set(MKL_BLAS_LIB
     ${MKL_LIB_DIR}/libmkl_blas95_lp64.a
     ${MKL_LIB_DIR}/libmkl_lapack95_lp64.a 
     -Wl,--start-group
     ${MKL_FORTRAN_INTERFACE_LIB}
     ${MKL_THREADING_LIB}
     ${MKL_LIB_DIR}/libmkl_core.a
     -Wl,--end-group
     ${MKL_OMP_LIB_LINE}
     -lpthread
     -lm
     -ldl)
  set(MKL_LAPACK_LIB ${MKL_BLAS_LIB})
else() # end UNIX
  message(FATAL_ERROR "unhandled system type")
endif()  


#-------------------------------------------------------------------------------
# Set BLAS, LAPACK and PARDISO libraries depending on the MKL version.
#-------------------------------------------------------------------------------

# see also External_OpenBLAS and External_LAPACK, where the setting is quite different:
# e.g. BLAS_LIBRARY=-Wl,--start-group;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_intel_lp64.a;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_gnu_thread.a;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_core.a;-Wl,--end-group;-L/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/../compiler/lib/intel64;-liomp5;-lpthread;-lm;-ldl
# e.g. LAPACK_LIBRARY=
if(USE_BLAS_LAPACK STREQUAL "MKL")
  set(BLAS_LIBRARY "${MKL_BLAS_LIB}")
  if(MKL_MAJOR_VERSION LESS 10)
    set(LAPACK_LIBRARY "${MKL_LAPACK_LIB}")
  endif(MKL_MAJOR_VERSION LESS 10)  
endif()
set(PARDISO_LIBRARY "${MKL_PARDISO_LIB}")

#-------------------------------------------------------------------------------
# Status message of found MKL
#-------------------------------------------------------------------------------
if(DEFINED MKL_UPDATE)
  message(STATUS "Using Intel MKL version ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}.")
else()
  message(STATUS "Using Intel MKL version ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.")
endif()

message(STATUS "defining MKL link-line via MKL_BLAS_LIB=${MKL_BLAS_LIB}")# VERBOSE-TRACE are only supported from cmake 3.15

mark_as_advanced(MKL_LIB_DIR)
mark_as_advanced(MKL_INCLUDE_DIR)
