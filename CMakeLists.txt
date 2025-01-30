# cmake 3.16 is required for essential options of the new FindPython.cmake 
cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)
message(STATUS "Using CMake ${CMAKE_VERSION} with generator ${CMAKE_GENERATOR}")

# Include global definitions for openCFS development server (host name, etc.)
include("cmake_modules/DevelopmentServer.cmake")
include(CMakePrintHelpers)

if(POLICY CMP0135)
  #>= cmake 3.24 wants a policy for DOWNLOAD_EXTRACT_TIMESTAMP with external projects. is demanded
  # howerver, setting this to NEW makes make always build the complete source
  cmake_policy(SET CMP0135 OLD)
endif()

# Note: "Policy CMP0077 is not set: option() honors normal variables." 
# is likely when the USE_* option is set, e.g. in .cfs_platform_defaults.cmake
# then the variable is preused before option() - better set USE_*_DEFAULT!

# output some information if we are on a CI server
if(DEFINED ENV{CI_SERVER})
  message("GitLab CI pipeline: $ENV{CI_PIPELINE_URL}")
  message("GitLab CI job: $ENV{CI_JOB_URL}")
endif()

# The name of our project is "CFS".  CMakeLists files in this
# project can refer to the root source directory of the project as
# ${CFS_SOURCE_DIR} and to the root binary directory of the project
# as ${CFS_BINARY_DIR}. We use C, C++ and Fortran as languages.
if("$ENV{CI_JOB_STAGE}" MATCHES "test(_extra)?")
  message(STATUS "detected running in CI pipeline 'test' or 'test_extra' stage via environment variable CI_JOB_STAGE=$ENV{CI_JOB_STAGE}")
  # just some debug/info output
  if(DEFINED PIPELINE_TEST_STAGE_CONFIG) # check if the variable is already set
    message(WARNING "PIPELINE_TEST_STAGE_CONFIG does not need to be set any more, we use environment variables to do it automatically")
    if(PIPELINE_TEST_STAGE_CONFIG)
      message(STATUS "variable PIPELINE_TEST_STAGE_CONFIG set to ON in the test-stage")
    else()
      message(WARNING "why is PIPELINE_TEST_STAGE_CONFIG=OFF in the test stage? maybe the closed build")
    endif()
  else() # set the variable to get easy if statements
    set(PIPELINE_TEST_STAGE_CONFIG ON)  
  endif()  
endif()
if(PIPELINE_TEST_STAGE_CONFIG)
  message(STATUS "No compiler check due to PIPELINE_TEST_STAGE_CONFIG variable (=${PIPELINE_TEST_STAGE_CONFIG})")
  if(NOT DEFINED ENV{CI_SERVER}) # only issue a warning if we're not on a CI server
    message(WARNING "No compiler check - this should only occur in the testing stage of the gitlab pipeline")
  endif()
  set(CMAKE_MAKE_PROGRAM "none")
  project(CFS NONE)
else()
# this else part is has no indent and goes for almost 400 lines. Search for PIPELINE_TEST_STAGE_CONFIG
project(CFS C CXX Fortran) 

# most early identification of detected compilers with their variables to support coding cmake files
# cmake prints this information also for the initial configure but not for following runs, therefore explicitly 
message(STATUS "CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID} CMAKE_C_COMPILER_VERSION=${CMAKE_C_COMPILER_VERSION}")
message(STATUS "CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")

if(WIN32)
  cmake_print_variables(CMAKE_LINKER)	
endif()

# Make sure that output of started programs is in English, so that we can
# grep for defined messages.
set(ENV{LC_MESSAGES} "en_EN")
set(ENV{LANG} C)

# adds the cmake feature cmake_dependent_option
include(CMakeDependentOption)

# Include some utility macros (for building ext. libs, etc...)
include("cmake_modules/CFS_macros.cmake")

# CFS version number.  Simply the last modification date. 
# It is good practice to keep this up to date, independent on the significance of the change (CFS_NAME can stay)
set(CFS_VERSION_YEAR "24")
set(CFS_VERSION_MONTH "11")
set(CFS_VERSION "${CFS_VERSION_YEAR}.${CFS_VERSION_MONTH}")

# Version Name. Every major release gets a label/name
# shall be updated less frequently than VERSION
set(CFS_NAME "Ideal Implementation")

colormsg(HIBLUE "Configuration for openCFS ${CFS_VERSION}: ${CFS_NAME} starting")

# Set hostname of build host
site_name(CFS_BUILD_HOST)
set(CFS_BUILD_HOST "${CFS_BUILD_HOST}" CACHE INTERNAL "CFS_BUILD_HOST")

set(CMAKE_MODULE_PATH ${CFS_SOURCE_DIR}/cmake_modules)

#-----------------------------------------------------------------------------
# Now come the options. Sanity checks only below all options
#-----------------------------------------------------------------------------

# Set default values according to environment or platform. The user can
# specify his own defaults in $HOME/.cfs_platform_defaults.cmake
include("cmake_modules/platform_defaults.cmake")

# this are general options for configuring cfs and also provide support options for code developers
set(CFS_PARDISO ${CFS_PARDISO_DEFAULT} CACHE STRING "Pardiso implementation for CFS. May be one of MKL or SCHENK. Irrelevant if disabled via USE_PARDISO")
set_property(CACHE CFS_PARDISO PROPERTY STRINGS MKL SCHENK)
# Switch on/off usage of cfsdepscache
option(CFS_DEPS_PRECOMPILED "Use (create) precompiled cfsdeps packages" ${CFS_DEPS_PRECOMPILED_DEFAULT})
option(DEBUG "Build debug version of the code which allows debugging, logging and has asserts - up to 10 times slower." ${DEBUG_DEFAULT})
option(CFS_NATIVE "Compile for native CPU for speed, might not be transferable" ${CFS_NATIVE_DEFAULT})
option(CFS_PROFILING "Set profiling compiler options for Release for vtune and valgrind (-g)." OFF)
option(CFS_COVERAGE "Turn on coverage compiler options (gcc, clang) to be analyzed by gcov." OFF)
option(CFS_FSANITIZE "Use -fsanitize to check for memory leaks and general errors. Very costly during runtime!" OFF)
# background for reordering is that on arm64 metis fails but sloand and no reordering works
set(CFS_REORDERING ${CFS_REORDERING_DEFAULT} CACHE STRING "Overwrite implementation of default for reordering")
set_property(CACHE CFS_REORDERING PROPERTY STRINGS default Metis Sloan noReordering)

# Do we want to use expression templates
option(USE_EXPRESSION_TEMPLATES "Use expression templates for matrix and vector operators." ON)
# Set a sane default value for the openCFS schema path
set(XMLSCHEMA "${CFS_SOURCE_DIR}/share/xml" CACHE STRING "Path where the XML schema for the input files is located.")
# this is by default the submodule
set(TESTSUITE_DIR "${TESTSUITE_DIR_DEFAULT}" CACHE PATH "Path where the openCFS testsuite is located.")

# BUILD_* generates optional artefacts. When you generate a file, consider to trigger it via BUILD_*
option(BUILD_CFSTOOL "Build openCFS converter and testing tool (cfstool)."  ${BUILD_CFSTOOL_DEFAULT})
option(BUILD_CFSDAT "Build openCFS data analysis tool (cfsdat)."  ${BUILD_CFSDAT_DEFAULT})
option(BUILD_UNIT_TESTS "Build openCFS unit tests (cfstest)."  ${BUILD_UNIT_TESTS_DEFAULT})

option(BUILD_TESTING "Build the testsuite tree, required BUILD_CFSTOOL" ${BUILD_TESTING_DEFAULT}) 

# USE_* triggers more or less optional external dependencies.
set(USE_BLAS_LAPACK ${USE_BLAS_LAPACK_DEFAULT} CACHE STRING "BLAS/LAPACK implementation. May be one of OPENBLAS, MKL, NETLIB, (Apple) ACCELERATE")
# we set convenience USE_MKL, USE_OPENBLAS, USE_NETLIB, USE_ACCELERATE in openmp_blas.cmake
set_property(CACHE USE_BLAS_LAPACK PROPERTY STRINGS OPENBLAS MKL NETLIB ACCELERATE)
option(USE_OPENMP "Enable support for OpenMP" ${USE_OPENMP_DEFAULT})
option(USE_GIDPOST "Use the GiD library for simulation output."  ${USE_GIDPOST_DEFAULT})
# USE_ENSIGHT will we mirrored to USE_VTK which is only used for Ensight
option(USE_ENSIGHT "Support for Ensight Gold file format for input.." ${USE_ENSIGHT_DEFAULT})
option(USE_EMBEDDED_PYTHON "External Python Optimizer" ${USE_EMBEDDED_PYTHON_DEFAULT})
option(USE_SGPP    "Enable support for sparse grid interpolation."   ${USE_SGPP_DEFAULT})
option(USE_PARDISO "Use Pardiso, see also advanced option CFS_PARDISO (MKL or SCHENK)" ${USE_PARDISO_DEFAULT})
option(USE_ARPACK "Use the ARPACK eigensolver library. Required for Pade Approximation Linearization Method (PALM) solver." ${USE_ARPACK_DEFAULT})
option(USE_FEAST "Use the FEAST eigenvalue solver (commumity edition)." ${USE_FEAST_DEFAULT})
option(USE_SUITESPARSE "CHOLMOD/UMPFACK direct solvers with GPL configuration" ${USE_SUITESPARSE_DEFAULT} )
option(USE_LIS "Use (L)ibrary of (I)terative (S)olvers (OpenMP parallel)" ${USE_LIS_DEFAULT})
option(USE_SUPERLU "Use SuperLU direct solver and ILUs. Required for Pade Approximation Linearization Method (PALM) solver." ${USE_SUPERLU_DEFAULT})
option(USE_PETSC "Suite of MPI parallel linear solver, requires local mpi installation" ${USE_PETSC_DEFAULT})
option(USE_PHIST_EV "phist eigensolver library from sppexa." ${USE_PHIST_EV_DEFAULT})
option(USE_PHIST_CG "phist linear library from sppexa." ${USE_PHIST_CG_DEFAULT})
option(USE_METIS "Use Metis Graph partitioning." ${USE_METIS_DEFAULT})
# CGNS - CFD General Notation System
option(USE_CGNS "Enable support for the CFD General Notation System." ${USE_CGNS_DEFAULT})
# SCPIP optimizer library.
option(USE_SCPIP "Use the external academic SCPIP optimizer." ${USE_SCPIP_DEFAULT})
# MMA/GCMMA optimizer library.
option(USE_DUMAS "Use the MMA/GCMMA C++ implementations by Jérémie Dumas." ${USE_DUMAS_DEFAULT})
# commercial SNOPT optimizer library.
option(USE_SNOPT "Use the external commercial SNOPT optimizer." ${USE_SNOPT_DEFAULT})
# open source IPOPT optimizer library.
option(USE_IPOPT "Use the external open source IPOPT optimizer." ${USE_IPOPT_DEFAULT})
# open source SGP optimizer library.
option(USE_SGP "Use the external open source SGP optimizer." ${USE_SGP_DEFAULT})
option(USE_CGAL "Activate fast element mapping using CGAL library." ${USE_CGAL_DEFAULT})
option(USE_EIGEN "Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms." ${USE_EIGEN_DEFAULT})
# Use libfbi toggles the inclusion of the Fast Box Intersection library.
# Please note, that either CGAL or libfbi can be used at a time, but not both!
option(USE_LIBFBI "Use fast box intersection library." ${USE_LIBFBI_DEFAULT})
option(USE_FLANN "Enable support for FLANN - Fast Library for Approximate Nearest Neighbors." ${USE_FLANN_DEFAULT})
# we choose here and hande USE_XERCES and USE_LIBXML2 in the conditional options
set(USE_XML_READER ${USE_XML_READER_DEFAULT} CACHE STRING "Choose either 'xerces' or 'libxml2' as xml reader.")
set_property(CACHE USE_XML_READER PROPERTY STRINGS xerces libxml2)

option(USE_PRECICE "Enable preCICE support" ${USE_PRECICE_DEFAULT} )


# icc has issues with boost asio, so simply disable this rare feature
if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(USE_STREAMING OFF)
else()
  set(USE_STREAMING ON)
endif()    

#-----------------------------------------------------------------------------
# Do sanitiy checks after all the settings are done
# Note, that dependend options come next
#-----------------------------------------------------------------------------

# currently, we use eigen only in CGAL for geometry approximation.
if (NOT USE_CGAL AND USE_EIGEN)
  message(WARNING "\nYour configuration USE_EIGEN=ON and USE_CGAL=OFF has no effect. \nEigen is currently only used by CGAL, for surface-geometry approximation in the automatic layer generation and curvilinear PML.")
endif()

if(USE_CGAL AND NOT USE_EIGEN)
  message(WARNING "\nUsing USE_CGAL=ON and USE_EIGEN=OFF will likely fail since CGAL needs Eigen.")
endif()
  
if(USE_CGAL AND USE_LIBFBI)
  message(FATAL_ERROR "Either USE_CGAL or USE_LIBFBI can be actived, but not both!")
endif()

if(USE_PARDISO)
  # PARDISO form MKL also requires BLAS/LAPACK from MKL
  # with CFS_PARDISO=SCHENK can be used with OPENBLAS or NETLIB
  if(CFS_PARDISO STREQUAL "MKL" AND NOT USE_BLAS_LAPACK STREQUAL "MKL")
    message(FATAL_ERROR "CFS_PARDISO=MKL requires USE_BLAS_LAPACK=MKL. Disable USE_PARDISO or set CFS_PARDSIO to SCHENK if not using USE_BLAS_LAPACK=MKL.")
  endif()
endif(USE_PARDISO)

# Since MKL Pardiso only works with MKL BLAS/LAPACK, we force that here.
if(USE_PARDISO AND CFS_PARDISO STREQUAL "MKL")
  if(NOT USE_BLAS_LAPACK STREQUAL "MKL")
    message(FATAL_ERROR "with USE_PARDISO and CFS_PARDISO=MKL one also needs USE_BLAS_LAPACK=MKL") 
  endif()
endif()

if(BUILD_TESTING AND NOT EXISTS ${TESTSUITE_DIR}/TESTSUIT AND NOT DEFINED ENV{CI_PIPELINE_NO_TESTSUITE})
  message(WARNING "BUILD_TESTING is on but TESTSUITE_DIR seems invalid, deactivate BUILD_TESTING")
endif()

if(BUILD_TESTING AND NOT BUILD_CFSTOOL)
  message(FATAL_ERROR "BUILD_CFSTOOL required for BUILD_TESTING")
endif()

if(USE_SGP AND NOT USE_BLAS_LAPACK STREQUAL "MKL")
  message(FATAL_ERROR "USE_SGP requires USE_BLAS_LAPACK=MKL")    
endif()

#-----------------------------------------------------------------------------
# Do conditional settings after sanity checks have passed
#-----------------------------------------------------------------------------

# Do we have mpi? There are multiple mpi installations, one should stick with one
# but we don't build mpi by cfsdeps to allow usage of tuned hpc mpi instances.
# phist requires ghost, ghost requires hwloc. One can set hwloc also for petsc
if(USE_PETSC)
  set(USE_MPI "ON")
  include("cmake_modules/mpi.cmake")
else()
  set(USE_MPI "OFF") 
endif()

# phist requires ghost
if(USE_PHIST_EV OR USE_PHIST_CG)
  set(BUILD_GHOST ON CACHE BOOL "ghost kernel library (optionally cuda based) from sppexa/FAU. Required for phist" FORCE)
endif()

if(BUILD_GHOST)
  set(BUILD_HWLOC ON CACHE BOOL "toggled by USE_PHIST and BUILD_GHOST" FORCE)  
else() 
  set(BUILD_HWLOC OFF CACHE BOOL "toggled by USE_PHIST and BUILD_GHOST" FORCE)  
endif()

if(USE_XML_READER STREQUAL "xerces")
  set(USE_XERCES ON CACHE BOOL "toggled by USE_XML_READER" FORCE)
  set(USE_LIBXML2 OFF CACHE BOOL "toggled by USE_XML_READER" FORCE)
else()
  set(USE_XERCES OFF CACHE BOOL "toggled by USE_XML_READER" FORCE)
  set(USE_LIBXML2 ON CACHE BOOL "toggled by USE_XML_READER" FORCE)
endif()  

if(USE_ENSIGHT)
  set(USE_VTK ON  CACHE INTERNAL "Surplus comment for cached variable") 
else()
  set(USE_VTK OFF CACHE INTERNAL "Surplus comment for cached variable")      
endif()

if(DEBUG)
  set(CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "Only DEBUG and RELEASE types are valid." FORCE)
  
  # Check for certain error conditions (in DEBUG mode)
  set(CHECK_INDEX 1)
  set(CHECK_INITIALIZED 1)
  set(CHECK_MEM_ALLOC 1)
  set(CHECK_MEMORY 1)
else()
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Only DEBUG and RELEASE types are valid." FORCE)
endif()

#-----------------------------------------------------------------------------
# Mark some of the options as advanced to clear up the standard selection
#-----------------------------------------------------------------------------
mark_as_advanced(CFS_PARDISO)
mark_as_advanced(CFS_REORDERING)
mark_as_advanced(CMAKE_BUILD_TYPE)
mark_as_advanced(CMAKE_BACKWARDS_COMPATIBILITY)
mark_as_advanced(CMAKE_INSTALL_PREFIX)
mark_as_advanced(BUILD_GHOST)
mark_as_advanced(BUILD_HWLOC)
mark_as_advanced(XMLSCHEMA)
mark_as_advanced(USE_XERCES)
mark_as_advanced(USE_LIBXML2)
mark_as_advanced(USE_PHIST_CG)

#-----------------------------------------------------------------------------
# Prepare further build
#-----------------------------------------------------------------------------

# Find executables of a few required programs
include("cmake_modules/FindPrograms.cmake")

# Set path where shared/static libraries should go. We only have lib
set(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for building all libraries.")

# set RPATH to a relative path (Linux and theoretically but not pracically macOS)
# this ensures dynamically linked libraries are found
# check the set value with `chrpath -l bin/cfs` (Linux) or 'otool -l bin/cfs' (macOS) 
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
# this is no typo with $ORIGIN: https://stackoverflow.com/tags/rpath/info
set(CMAKE_BUILD_RPATH "$ORIGIN/../lib") # this is needed such that origin is first
set(CMAKE_INSTALL_REMOVE_ENVIRONMENT_RPATH TRUE) # should remove oneapi stuff, but does not seem to work
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib") # for some reason this seems to be ignored because rpath remains the same after install

# Set path where executables should go
set(EXECUTABLE_OUTPUT_PATH ${CFS_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")

# find conditonally OpenMP and load the proper BLAS (especially MKL)
include("cmake_modules/openmp_blas.cmake")

# Set compiler flags according to output target and debug setting
# TODO: for in distro.cmake the architecture is set, hence distro.cmake
#       shall run before compiler.cmake. However vor win32 compiler settings are used in distr.sh 
include("cmake_modules/compiler.cmake")

# Determine on what UNIX/Linux/Windows distribution we are building on
include("cmake_modules/distro.cmake")

# Determine if some Fortran system libs are installed (blas, lapack, g2c...)
include("cmake_modules/CheckFortranRuntime.cmake")

# Perform some configure checks
include("cmake_modules/ConfigureChecks.cmake")

# Create and set a path where temporaries should go
set(CFS_TEMP_DIR "${CFS_BINARY_DIR}/tmp")
file(MAKE_DIRECTORY "${CFS_TEMP_DIR}")

# Find CFSDEPS, i.e. gidpost, HDF5, pardiso, arpack...
include("cmake_modules/FindCFSDEPS.cmake")

# Set path of XML schema. If env variable XMLSCHEMA is set, use it!
set(MYXMLSCHEMA $ENV{XMLSCHEMA})
if(NOT ${MYXMLSCHEMA} STREQUAL "")
  set(XMLSCHEMA ${MYXMLSCHEMA} CACHE STRING "Path where the XML schema for the input files is located." FORCE)
endif()

# Set date of configure run.
TODAY(CFS_CONF_DATE)

# Determine URL of subversion repository, global revision of working copy
# and other informations from working copy
include("cmake_modules/DetermineWorkingCopyInfos.cmake")

endif(PIPELINE_TEST_STAGE_CONFIG) # this is the end of the non indent else() for the if(PIPELINE_TEST_STAGE_CONFIG) case

# Configure Dart testing support.  This should be done before any
# message(FATAL_ERROR ...) commands are invoked (???) and after the languages
# for the project have been defined since CMake needs to know about the 
# language-specific build tools.
# Switch support for testing ON if the default testsuite dir is used, i.e.
# the git submodule 'Testing'
if(BUILD_TESTING OR DEFINED ENV{CI_PIPELINE_USE_BUILDTESTS})
  # check if TESTSUITE_DIR is correct (contains TESTSUIT)
  if(BUILD_TESTING AND NOT DEFINED ENV{CI_PIPELINE_NO_TESTSUITE})
    if(EXISTS ${TESTSUITE_DIR}/TESTSUIT)
      add_subdirectory(${TESTSUITE_DIR} "${CMAKE_CURRENT_BINARY_DIR}/testsuite")
      message(STATUS "including testsuite via ${CMAKE_CURRENT_BINARY_DIR}/testsuite")
    else()
      set(MSG "TESTSUITE_DIR '${TESTSUITE_DIR}' does not contain a directoy 'TESTSUIT'.\n")
      set(MSG "${MSG}Please specify a proper location!\n")
      set(MSG "${MSG}Did you clone with 'git clone --recurse-submodules'?\n")
      set(MSG "${MSG}If not go to your source directory and do\n")
      set(MSG "${MSG}   'git submodule update --init --recursive'\n")
      set(MSG "${MSG}   'git pull --recurse-submodules'\n")
      message(FATAL_ERROR "${MSG}")
    endif()
  endif()
  # set a build name
  set(BUILDNAME "${CFS_WC_REVISION} on ${DIST} ${REV} with ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}" CACHE STRING "Name of build on the dashboard")
  mark_as_advanced(BUILDNAME)
  enable_testing()
endif(BUILD_TESTING OR DEFINED ENV{CI_PIPELINE_USE_BUILDTESTS})
if(NOT PIPELINE_TEST_STAGE_CONFIG)

# add packing and build test
#  - only used in CI pipeline
if(DEFINED ENV{CI_PIPELINE_USE_BUILDTESTS} OR "$ENV{CI_JOB_STAGE}" MATCHES "build")
  message(STATUS "running in CI pipeline 'build' stage detected -> creating build and pack tests")
  include(CTest)
  # Create distributable binary via cpack running ctest, alternatively use create_binary_archive.py
  if(UNIX) # pack to a tar.gz archive
    add_test(Pipeline_Build_cpack-TGZ cpack -G TGZ --verbose)
  else() # probably windows, we do zip
    add_test(Pipeline_Build_cpack-ZIP cpack -G ZIP --verbose)
  endif()
  # Use an additional test during the build stage to grep the compiler errors
  add_test(Pipeline_Build_check-compiler-errors grep "0 Compiler errors" build.log)
  enable_testing()
  message(STATUS "submitting tests to: ${CTEST_SUBMIT_URL}")
endif()

if(NOT DEFINED ENV{CI_PIPELINE_USE_BUILDTESTS} AND NOT BUILD_TESTING)
  # remove Ctest config if testing is off
  # will probably not be necessay any mode with https://gitlab.kitware.com/cmake/cmake/-/merge_requests/258
  file(REMOVE "${CFS_BINARY_DIR}/CTestTestfile.cmake")
endif()

# ------------------------------------------------------------------------
# Add stuff according to the settings
# ------------------------------------------------------------------------

# Make sure the compiler can find include files from CFS.
include_directories(
  ${BOOST_INCLUDE_DIR} 
  ${CFS_BINARY_DIR}/include
  ${CFS_SOURCE_DIR}/source)

if(USE_CGAL)
  include_directories(${CGAL_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})
endif(USE_CGAL)

if(USE_VTK)
  include_directories(${VTK_INCLUDE_DIR})
endif(USE_VTK)


# Make sure the linker can find the CFS base libraries once it is built.
link_directories(${CFSDEPS_LIBRARY_DIR} ${LIBRARY_OUTPUT_PATH})

# Recurse into the "source" subdirectory.  This does not actually
# cause another cmake executable to run.  The same process will
# walk through the project's entire directory structure.
add_subdirectory(source/OLAS/algsys)
add_subdirectory(source/OLAS/graph)
add_subdirectory(source/OLAS/precond)
add_subdirectory(source/OLAS/solver)
add_subdirectory(source/OLAS/utils)
if(USE_ARPACK)
  add_subdirectory(source/OLAS/external/arpack)
endif()
if(USE_PHIST_CG OR USE_PHIST_EV)
  add_subdirectory(source/OLAS/external/phist)
endif()
add_subdirectory(source/OLAS/external/lapack)
if(USE_PARDISO)
  add_subdirectory(source/OLAS/external/pardiso)
endif(USE_PARDISO)
if(USE_FEAST)
  add_subdirectory(source/OLAS/external/feast)
endif()
if(USE_SUITESPARSE)
  add_subdirectory(source/OLAS/external/umfpack)
  add_subdirectory(source/OLAS/external/cholmod)
endif(USE_SUITESPARSE)
if(USE_SUPERLU)
  add_subdirectory(source/OLAS/external/superlu)
endif(USE_SUPERLU)
if(USE_LIS)
  add_subdirectory(source/OLAS/external/lis)
endif(USE_LIS)
if(USE_PETSC)
  add_subdirectory(source/OLAS/external/petsc)
endif(USE_PETSC)
add_subdirectory(source/CoupledPDE)
add_subdirectory(source/DataInOut)
add_subdirectory(source/DataInOut/Logging)
add_subdirectory(source/DataInOut/ParamHandling)
add_subdirectory(source/DataInOut/SimInOut/AnsysCDB)
add_subdirectory(source/DataInOut/SimInOut/AnsysFile)
add_subdirectory(source/DataInOut/SimInOut/internalMesh)
add_subdirectory(source/DataInOut/SimInOut/gmsh)
add_subdirectory(source/DataInOut/SimInOut/RefElems)
add_subdirectory(source/DataInOut/SimInOut/hdf5)
if(USE_GIDPOST)
  add_subdirectory(source/DataInOut/SimInOut/GiD)
endif(USE_GIDPOST)
add_subdirectory(source/DataInOut/SimInOut/Unverg)
if(USE_ENSIGHT)
  add_subdirectory(source/DataInOut/SimInOut/VTKBased/Ensight)
endif()
if(USE_CGNS)
  add_subdirectory(source/DataInOut/SimInOut/CGNS)
endif(USE_CGNS)
if(USE_EMBEDDED_PYTHON)
  add_subdirectory(source/DataInOut/SimInOut/python)
endif(USE_EMBEDDED_PYTHON)

add_subdirectory(source/DataInOut/SimInOut/TextOutput)
add_subdirectory(source/DataInOut/SimInOut/InfoResultOutput)
if(USE_STREAMING)
  add_subdirectory(source/DataInOut/SimInOut/Streaming)
endif()
add_subdirectory(source/DataInOut/ScatteredDataInOut)
add_subdirectory(source/Domain)
add_subdirectory(source/Domain/CoefFunction)
add_subdirectory(source/Domain/CoordinateSystems)
add_subdirectory(source/Domain/ElemMapping)
add_subdirectory(source/Domain/Mesh)
add_subdirectory(source/Domain/Results)
add_subdirectory(source/Driver)
add_subdirectory(source/Driver/TimeSchemes)
add_subdirectory(source/FeBasis)
add_subdirectory(source/Forms)
add_subdirectory(source/Forms/LinForms)
add_subdirectory(source/Forms/BiLinForms)
add_subdirectory(source/General)
add_subdirectory(source/Materials)
add_subdirectory(source/MatVec)
add_subdirectory(source/ODESolve)
add_subdirectory(source/PDE)
add_subdirectory(source/PDE/LatticeBoltzmannSolver)
add_subdirectory(source/Utils)
add_subdirectory(source/Utils/mathParser)
add_subdirectory(source/Utils/EvalIntegrals)
add_subdirectory(source/Optimization)
add_subdirectory(source/Optimization/Design)
add_subdirectory(source/Optimization/Optimizer)
add_subdirectory(source/main)
if(BUILD_CFSTOOL)
  add_subdirectory(source/cfstool)
endif()
if(BUILD_CFSDAT)
  add_subdirectory(source/cfsdat)
endif()
if(BUILD_UNIT_TESTS)
  add_subdirectory(source/unittests)
endif()

# Build the actual headers with defines from templates.
include("cmake_modules/configure_def_headers.cmake")
endif(NOT PIPELINE_TEST_STAGE_CONFIG)
# all copy/configure only targets
add_subdirectory(share/doc) #building documentation (see 'make help')
add_subdirectory(share/scripts)
add_subdirectory(share/xml)
add_subdirectory(share/xsl)
add_subdirectory(share/python)

# generate compile_commands.json in the build dir which is used e.g. by vscode
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "") # since project() sets a cache entry for it

# set cpack config options
# see https://cmake.org/Wiki/CMake:CPackConfiguration
string(STRIP "${CFS_GIT_COMMIT}" CPACK_PACKAGE_VERSION)
if(CFS_WC_MODIFIED)
  set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}-modified")
endif()
set(CPACK_CMAKE_GENERATOR "Ninja") # this disables a rebuild i.e. "CPack: - Run preinstall target for..." which seems to be only done for "Unix Makefiles"

# make sure that on cpack -G ZIP or make install everything is copied. make install is uncommon
include("cmake_modules/redistributables.cmake")
