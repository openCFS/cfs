# Lis: a Library of Iterative Solvers for Linear Systems (http://www.ssisc.org/lis/)
# https://github.com/norihiro-w/lis-cmake

cmake_minimum_required(VERSION 2.8.7)
project(Lis C Fortran)

INCLUDE(CMakeFortranInformation)
INCLUDE(FortranCInterface)
INCLUDE(cmake/CheckTypeSizes.cmake)
INCLUDE(cmake/FindIncludeHeader.cmake)

######################################################################
# Options

OPTION(LIS_ENABLE_COMPLEX "Enable complex" ON)
OPTION(LIS_ENABLE_OPENMP "Use OpenMP" ON)
OPTION(LIS_ENABLE_MPI "Use MPI" OFF)
OPTION(LIS_ENABLE_FORTRAN "build Fortran 77 interfaces" OFF)
OPTION(LIS_ENABLE_SAAMG "build SA-AMG preconditioner" OFF)
OPTION(LIS_ENABLE_QUAD "enable quadruple precision operations" OFF)
OPTION(LIS_ENABLE_SSE2 "use Intel Streaming SIMD Extensions" ON)
OPTION(LIS_BUILD_TEST "Build tests" OFF)

######################################################################
# Install directories

set(INSTALL_BIN_DIR
  "${CMAKE_INSTALL_PREFIX}/bin"
  CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR
  "${CMAKE_INSTALL_PREFIX}/${LIB_SUFFIX}"
  CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR
  "${CMAKE_INSTALL_PREFIX}/include"
  CACHE PATH "Installation directory for headers")

file(READ ${CMAKE_CURRENT_SOURCE_DIR}/include/lis.h _lis_h_contents)
string(REGEX REPLACE ".*#define[ \t]+LIS_VERSION[ \t]+\"([0-9A-Za-z.]+).*"
  "\\1" LIS_VERSION ${_lis_h_contents})

######################################################################
# General build setting

# Set build directories
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )
# Enable Visual Studio project folder grouping
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

######################################################################
# Lis build settting

ADD_DEFINITIONS(-DHAVE_CONFIG_H)
IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF()
# Create Configure.h
CONFIGURE_FILE (include/lis_config.h.in ${PROJECT_BINARY_DIR}/lis_config.h)
INCLUDE_DIRECTORIES (${PROJECT_BINARY_DIR}/)

IF(LIS_ENABLE_OPENMP)
	FIND_PACKAGE( OpenMP REQUIRED )
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	IF(NOT ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC"))
		SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
	ENDIF()
ENDIF()

IF(LIS_ENABLE_MPI)
	FIND_PACKAGE(MPI)
	IF(MPI_FOUND)		
		SET(CMAKE_C_COMPILER ${MPI_COMPILER})
		SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
		ADD_DEFINITIONS(-DUSE_MPI)
	ELSE(MPI_FOUND)
		MESSAGE (FATAL_ERROR "Aborting: MPI implementation is not found!")
	ENDIF(MPI_FOUND)			
ENDIF()

IF(LIS_ENABLE_COMPLEX)
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DHAVE_COMPLEX /DUSE_COMPLEX /DCOMPLEX" )
	SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /DHAVE_COMPLEX /DUSE_COMPLEX /DCOMPLEX" )
ENDIF()

######################################################################
# Add subdirectories

ADD_SUBDIRECTORY( src )
IF (BUILD_TEST)
	ADD_SUBDIRECTORY( test )
ENDIF()
