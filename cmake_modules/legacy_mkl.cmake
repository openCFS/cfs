# this is for systems open MKL/oneAPI before 2021 and shall be removed once we get rid of the centos6 pipeline

message(WARNING "MKL was not found via MKLConfig.cmake! 
Either your MKL installation is too old (< 2021.3) or the possible paths (MKL_POSSIBLE_PATHS) above need to be refined.
please install a recent MKL >= 2021.3 and/or set the correct hints in MKL_POSSIBLE_PATHS to find MKLConfig.cmake
Trying the old way instead - not guaraneed to work and not maintained any more ...")
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
  set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include" CACHE PATH "here we assume the mkl includes")

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
    # path for oneAPI 2024.0: There is a new "Unified Directory Layout", https://www.intel.com/content/www/us/en/developer/articles/release-notes/onemkl-release-notes.html
    if(NOT EXISTS "${MKL_OMP_LIB}")
      set(MKL_OMP_LIB "${MKL_ROOT_DIR}/../../compiler/${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}/lib/libiomp5.so")
    endif()
    # this needs to be a copy, otherwise centos6-gcc7 does not link as it has no setvars
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
  set(MKL_LIBS "${MKL_FORTRAN_INTERFACE_LIB},${MKL_THREADING_LIB},${MKL_LIB_DIR}/libmkl_core.a")


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
  set(MKL_VERSION "${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}")
else()
  set(MKL_VERSION "${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}")
endif()

# some debug output
# cmake_print_variables(CMAKE_LINK_GROUP_USING_RESCAN)
# cmake_print_variables(CMAKE_CXX_LINK_GROUP_USING_RESCAN_SUPPORTED)

message(STATUS "defining MKL link-line via MKL_BLAS_LIB=${MKL_BLAS_LIB}")# VERBOSE-TRACE are only supported from cmake 3.15
mark_as_advanced(MKL_LIB_DIR)
mark_as_advanced(MKL_INCLUDE_DIR)
message(STATUS "Using Intel MKL version ${MKL_VERSION}")