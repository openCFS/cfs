cmake_minimum_required (VERSION 3.6)

include(ExternalProject)
include(CMakePrintHelpers) # cmake_print_variables()

project (SuiteSparse LANGUAGES C )

option(BUILD_STATIC "build static libraries" ON)
option(USE_OPENMP "enable OpenMP" ON)
# inspired by https://github.com/sergiud/SuiteSparse
option(ALLOW_GPL_EXTENSIONS "let the solvers use optional GPL features" ON)
option(ALLOW_64BIT_BLAS "look for 32-bit BLAS only or concurrent 32- and 64-bit BLAS" OFF)
set(SUGGEST_BLAS_LIBRARIES "" CACHE STRING "override autodetection by suggeting lib(s) for BLAS_LIBRARIES")

option(ENABLE_AMD "routines for pre-ordering sparse matrices prior to factorization" ON)
option(ENABLE_CAMD "routines for permuting sparse matrices prior to factorization" ON)
option(ENABLE_CCOLAMD "constrained column approximate minimum degree ordering" ON)
option(ENABLE_COLAMD "column approximate minimum degree ordering" ON)
option(ENABLE_CHOLMOD "sparse CHOLesky MODification package" ON)
option(ENABLE_UMFPACK "routines solving sparse linear systems via LU factorization" ON)

# this are the common arguments for all individual cmake projects
set(CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX})
if(BUILD_STATIC)
  list(APPEND CMAKE_ARGS -DNSTATIC:BOOL=OFF)
else()
  list(APPEND CMAKE_ARGS -DNSTATIC:BOOL=ON)
endif()

if(USE_OPENMP)
  set(OPENMP -DNOPENMP:BOOL=OFF)
else()
  set(OPENMP -DNOPENMP:BOOL=ON)
endif()

# overwriting these values is actually wrong, as the first find( BLAS x) is true
if(NOT SUGGEST_BLAS_LIBRARIES STREQUAL "")
  set(BLAS_LIB -DBLAS_LIBRARIES=${SUGGEST_BLAS_LIBRARIES} -DLAPACK_LIBRARIES=${SUGGEST_BLAS_LIBRARIES})
endif()


# the CMAKE_SOURCE_DIR points to the root directory, independent from where it's CMakeLists.txt is called
# the add the desired subpackages. Via DEPENDS the proper order is guaranteed
set(ARGS ${CMAKE_ARGS} ${OPENMP})
if(NOT SUGGEST_BLAS_LIBRARIES STREQUAL "")
  list(APPEND ARGS ${BLAS_LIB})
endif()
ExternalProject_Add(SuiteSparse_config
  PREFIX "SuiteSparse_config"
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/SuiteSparse_config"
  BINARY_DIR "${CMAKE_SOURCE_DIR}/SuiteSparse_config/build"
  INSTALL_COMMAND "" # we disable instant install and perform install later via INSTALL_LIST 
  CMAKE_ARGS ${ARGS} -DALLOW_64BIT_BLAS=${ALLOW_64BIT_BLAS})

# here we collect install targets. Via INSTALL_COMMAND "" in the ExternalProject_Add() we prevent
# instant installation. Via this list a make install install in the build directory performes this
# by calling make install in all modules build directories. The install target is set via CMAKE_INSTALL_PREFIX
# in ExternalProject_Add()
set(INSTALL_LIST ${CMAKE_SOURCE_DIR}/SuiteSparse_config/build)

if (ENABLE_AMD OR ENABLE_CHOLMOD OR ENABLE_UMFPACK)
  if (NOT ENABLE_AMD)
    message(FATAL_ERROR "enable ENABLE_AMD, required by other package")
    # unfortunately set(ENABLE_AMD ON FORCE) is ignored :(
  endif()
  ExternalProject_Add(AMD
    PREFIX "AMD"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/AMD"
    BINARY_DIR "${CMAKE_SOURCE_DIR}/AMD/build"
    DEPENDS SuiteSparse_config
    INSTALL_COMMAND ""
    CMAKE_ARGS ${CMAKE_ARGS} )
  list(APPEND INSTALL_LIST ${CMAKE_SOURCE_DIR}/AMD/build)   
endif()

if (ENABLE_CAMD OR ENABLE_CHOLMOD)
  if (NOT ENABLE_CAMD)
    message(FATAL_ERROR "enable ENABLE_CAMD, required by CHOLMOD")
  endif()
  ExternalProject_Add(CAMD
    PREFIX "CAMD"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/CAMD"
    BINARY_DIR "${CMAKE_SOURCE_DIR}/CAMD/build"
    DEPENDS SuiteSparse_config
    INSTALL_COMMAND ""    
    CMAKE_ARGS ${CMAKE_ARGS} )
  list(APPEND INSTALL_LIST ${CMAKE_SOURCE_DIR}/CAMD/build)    
endif()

if (ENABLE_CCOLAMD OR ENABLE_CHOLMOD)
  if (NOT ENABLE_CCOLAMD)
    message(FATAL_ERROR "enable ENABLE_CCOLAMD, required by CHOLMOD")
  endif()
  ExternalProject_Add(CCOLAMD
    PREFIX "CCOLAMD"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/CCOLAMD"
    BINARY_DIR "${CMAKE_SOURCE_DIR}/CCOLAMD/build"
    DEPENDS SuiteSparse_config
    INSTALL_COMMAND ""    
    CMAKE_ARGS ${CMAKE_ARGS} )
  list(APPEND INSTALL_LIST ${CMAKE_SOURCE_DIR}/CCOLAMD/build)    
endif()

if (ENABLE_COLAMD OR ENABLE_CHOLMOD)
  if (NOT ENABLE_COLAMD)
    message(FATAL_ERROR "enable ENABLE_COLAMD, required by CHOLMOD")
  endif()
  ExternalProject_Add(COLAMD
    PREFIX "COLAMD"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/COLAMD"
    BINARY_DIR "${CMAKE_SOURCE_DIR}/COLAMD/build"
    DEPENDS SuiteSparse_config
    INSTALL_COMMAND ""    
    CMAKE_ARGS ${CMAKE_ARGS} )
  list(APPEND INSTALL_LIST ${CMAKE_SOURCE_DIR}/COLAMD/build)    
endif()

if (ENABLE_CHOLMOD)
  set(ARGS ${CMAKE_ARGS} ${OPENMP})
  if(NOT SUGGEST_BLAS_LIBRARIES STREQUAL "")
    list(APPEND ARGS ${BLAS_LIB})
  endif()
  if(NOT ALLOW_GPL_EXTENSIONS) # switch off gpl features
    list(APPEND ARGS -DNMATRIXOPS:BOOL=ON -DNMODIFY:BOOL=ON -DNSUPERNODAL:BOOL=ON)
  endif()
  ExternalProject_Add(CHOLMOD
    PREFIX "CHOLMOD"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/CHOLMOD"
    BINARY_DIR "${CMAKE_SOURCE_DIR}/CHOLMOD/build"
    DEPENDS SuiteSparse_config COLAMD CCOLAMD AMD
    INSTALL_COMMAND ""    
    CMAKE_ARGS ${ARGS} -DALLOW_64BIT_BLAS=${ALLOW_64BIT_BLAS} )
  list(APPEND INSTALL_LIST ${CMAKE_SOURCE_DIR}/CHOLMOD/build)    
endif()

if (ENABLE_UMFPACK)
  set(ARGS ${CMAKE_ARGS} ${OPENMP})
  if(NOT SUGGEST_BLAS_LIBRARIES STREQUAL "")
    list(APPEND ARGS ${BLAS_LIB})
  endif()

  ExternalProject_Add(UMFPACK
    PREFIX "UMFPACK"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/UMFPACK"
    BINARY_DIR "${CMAKE_SOURCE_DIR}/UMFPACK/build"
    DEPENDS SuiteSparse_config 
    INSTALL_COMMAND ""    
    CMAKE_ARGS ${ARGS} -DALLOW_64BIT_BLAS=${ALLOW_64BIT_BLAS} )
  list(APPEND INSTALL_LIST ${CMAKE_SOURCE_DIR}/UMFPACK/build)      
endif()

cmake_print_variables(CMAKE_MAKE_PROGRAM)

# generate suitsparse_install.cmake for root make install with INSTALL_LIST
configure_file("${CMAKE_SOURCE_DIR}/suitesparse_install.cmake.in" "${CMAKE_BINARY_DIR}/suitesparse_install.cmake" @ONLY)
install(SCRIPT "${CMAKE_BINARY_DIR}/suitesparse_install.cmake")
