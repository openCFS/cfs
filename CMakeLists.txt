cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)

# Include global definitions for openCFS development server (host name, etc.)
include("cmake_modules/DevelopmentServer.cmake")

# Set CMake policies for newer versions of CMake.
# Enter 'cmake --help-policy CMP000X' to get information about the meaning
# of a certain policy.
if(COMMAND CMAKE_POLICY)
  if(POLICY CMP0038)
    cmake_policy(SET CMP0038 NEW)
  endif(POLICY CMP0038)
endif(COMMAND CMAKE_POLICY)
# Note: "Policy CMP0077 is not set: option() honors normal variables." 
# is likely when the USE_* option is set, e.g. in .cfs_platform_defaults.cmake
# then the variable is preused before option() - better set USE_*_DEFAULT!


# The name of our project is "CFS".  CMakeLists files in this
# project can refer to the root source directory of the project as
# ${CFS_SOURCE_DIR} and to the root binary directory of the project
# as ${CFS_BINARY_DIR}. We use C, C++ and Fortran as languages.
project(CFS C CXX Fortran)

# Make sure that output of started programs is in English, so that we can
# grep for defined messages.
set(ENV{LC_MESSAGES} "en_EN")
set(ENV{LANG} C)

# Include some standard macros
include(CMakeDependentOption)

# Include some utility macros (for building ext. libs, etc...)
include("cmake_modules/CFS_macros.cmake")

# CFS version number.  Simply the last modification date. 
# It is good practice to keep this up to date (CFS_NAME can stay)
set(CFS_VERSION_YEAR "20")
set(CFS_VERSION_MONTH "09")
set(CFS_VERSION "${CFS_VERSION_YEAR}.${CFS_VERSION_MONTH}")

# Version Name. Every major release gets a label/name
# shall be updated less frequently than VERSION
set(CFS_NAME "Frohes Fensterln")

colormsg(HIBLUE "Configuration for openCFS ${CFS_VERSION}: ${CFS_NAME} starting")

# Set hostname of build host
site_name(CFS_BUILD_HOST)
set(CFS_BUILD_HOST "${CFS_BUILD_HOST}" CACHE INTERNAL "CFS_BUILD_HOST")

# Set user who builds openCFS, in order to be able to track back who built a
# given binary.
if(WIN32)
  set(REG_KEY "HKEY_CURRENT_USER\\Software\\Microsoft\\Windows")
  set(REG_KEY "${REG_KEY}\\CurrentVersion\\Explorer")
  execute_process(
    COMMAND reg query "${REG_KEY}" /v "Logon User Name"
    WORKING_DIRECTORY "."
    OUTPUT_VARIABLE OUTPUT)

  string(REPLACE "\n" ";" OUTPUT "${OUTPUT}")

  foreach(line IN LISTS OUTPUT)
    if(line MATCHES "Logon User Name")
      string(REPLACE "Logon User Name" "" line "${line}")
      string(REPLACE "REG_SZ" "" line "${line}")
      string(STRIP "${line}" CFS_BUILD_USER)
    endif()
  endforeach()

elseif(UNIX) # includes APPLE
  execute_process(
    COMMAND id -unr # get user name (possibly DOMAIN\user)
    COMMAND sed s/[A-Z]\\+\\\\// # remove 'DOMAIN\' from 'DOMAIN\user' ('\' must be escaped in cmake)
    WORKING_DIRECTORY "."
    RESULT_VARIABLE res OUTPUT_VARIABLE out)
  if(res)
    message("could not determine user: ${res}, ${out}")
  endif(res)
  string(STRIP "${out}" CFS_BUILD_USER)
endif()

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
option(CFS_PROFILING "Set profiling compiler options for Release for vtune and valgrind (-g)." OFF)
option(CFS_COVERAGE "Turn on coverage compiler options (gcc, clang) to be analyzed by gcov." OFF)
option(CFS_FSANITIZE "Use -fsanitize to check for memory leaks and general errors. Very costly during runtime!" OFF)
# Do we want to use explicit or implicit template instantiation
option(EXPLICIT_TEMPLATE_INSTANTIATION "Use explicit template instantiation." ${EXPLICIT_TEMPLATE_INSTANTIATION_DEFAULT})
# Do we want to use expression templates
option(USE_EXPRESSION_TEMPLATES "Use expression templates for matrix and vector operators." ON)
# ICC needs a gcc to get header from. Newest ICC might not be fully compatible with gcc 8.* yet
# will be set as advanced in non-intel case
set(CFS_ICC_GCC_VERSION ${CFS_ICC_GCC_VERSION_DEFAULT} CACHE STRING "GCC version where ICC takes the headers from")
# Set a sane default value for the openCFS schema path
set(XMLSCHEMA "${CFS_SOURCE_DIR}/share/xml" CACHE STRING "Path where the XML schema for the input files is located.")
# this is by default the submodule
set(TESTSUITE_DIR "${TESTSUITE_DIR_DEFAULT}" CACHE PATH "Path where the openCFS testsuite is located.")

# BUILD_* generates optional artefacts. When you generate a file, consider to trigger it via BUILD_*
option(BUILD_CFSTOOL "Build openCFS converter and testing tool (cfstool)."  ${BUILD_CFSTOOL_DEFAULT})
option(BUILD_CFSDAT "Build openCFS data analysis tool (cfsdat)."  ${BUILD_CFSDAT_DEFAULT})
option(BUILD_UNIT_TESTS "Build openCFS unit tests (cfstest)."  ${BUILD_UNIT_TESTS_DEFAULT})
option(BUILD_HDFVIEW "Build the HDF5 file viewer from the hdf5lib" ${BUILD_HDFVIEW_DEFAULT}) # advanced
option(BUILD_H5TOOL "Build openCFS own's h5tool" ON) #andvanced
# shall not work for WIN32
if(UNIX)
  option(BUILD_CPACK_BINARY "Create distributable binary via cpack running ctest, alternatively use create_binary_archive.py" ${BUILD_CPACK_BINARY_DEFAULT})
endif()
option(BUILD_TESTING "Build the testsuite tree, required BUILD_CFSTOOL" ${BUILD_TESTING_DEFAULT}) 

# USE_* triggers more or less optional external dependencies.
set(USE_BLAS_LAPACK ${USE_BLAS_LAPACK_DEFAULT} CACHE STRING "BLAS/LAPACK implementation. May be one of OPENBLAS, MKL, NETLIB")
set_property(CACHE USE_BLAS_LAPACK PROPERTY STRINGS OPENBLAS MKL NETLIB)
option(USE_OPENMP "Enable support for OpenMP" ${USE_OPENMP_DEFAULT})
option(USE_GIDPOST "Use the GiD library for simulation output."  ${USE_GIDPOST_DEFAULT})
option(USE_GMSH "Use Gmsh meshes for simulation input." ${USE_GMSH_DEFAULT})
option(USE_GMV "Use General Mesh Viewer (GMV) for simulation in/output."  ${USE_GMV_DEFAULT})
# USE_ENSIGHT will we mirrored to USE_VTK which is only used for Ensight
option(USE_ENSIGHT "Support for Ensight Gold file format for input.." ${USE_ENSIGHT_DEFAULT})
option(USE_EMBEDDED_PYTHON "External Python Optimizer" ${USE_EMBEDDED_PYTHON_DEFAULT})
option(USE_SGPP    "Enable support for sparse grid interpolation."   ${USE_SGPP_DEFAULT})
option(USE_PARDISO "Use Pardiso, see also advanced option CFS_PARDISO (MKL or SCHENK)" ${USE_PARDISO_DEFAULT})
option(USE_ARPACK "Use the ARPACK eigensolver library. Required for Pade Approximation Linearization Method (PALM) solver." ${USE_ARPACK_DEFAULT})
option(USE_FEAST "Use the FEAST eigenvalue solver." ${USE_FEAST_DEFAULT})
option(USE_ILUPACK "Iterative solver package with own algebraic multilevel prconditioners." ${USE_ILUPACK_DEFAULT})
option(USE_ILUPACK_PARALLEL "Ilupack with openmp. Experimental version." OFF) # advanced
option(USE_SUITESPARSE "Use CholMod and UMFPACK direct solvers from SuiteSparse." ${USE_SUITESPARSE_DEFAULT})
option(USE_LIS "Use (L)ibrary of (I)terative (S)olvers (OpenMP parallel)" ${USE_LIS_DEFAULT})
option(USE_SUPERLU "Use SuperLU direct solver and ILUs. Required for Pade Approximation Linearization Method (PALM) solver." ${USE_SUPERLU_DEFAULT})
option(USE_PETSC "Suite of MPI parallel linear solver, requires local mpi installation" ${USE_PETSC_DEFAULT})
option(USE_PHIST_EV "phist eigensolver library from sppexa." ${USE_PHIST_EV_DEFAULT})
option(USE_PHIST_CG "phist linear library from sppexa." ${USE_PHIST_CG_DEFAULT})
option(USE_METIS "Use Metis Graph partitioning." ${USE_METIS_DEFAULT})
# FFTW - Fastest Fourier Transform in the West (C subroutine library) 
option(USE_FFTW "Activate DFT and filtering capabilities of cfstool." ${USE_FFTW_DEFAULT})
# CGNS - CFD General Notation System
option(USE_CGNS "Enable support for the CFD General Notation System." ${USE_CGNS_DEFAULT})
# SCPIP optimizer library.
option(USE_SCPIP "Use the external academic SCPIP optimizer." ${USE_SCPIP_DEFAULT})
# commercial SNOPT optimizer library.
option(USE_SNOPT "Use the external commercial SNOPT optimizer." ${USE_SNOPT_DEFAULT})
# open source IPOPT optimizer library.
option(USE_IPOPT "Use the external open source IPOPT optimizer." ${USE_IPOPT_DEFAULT})
option(USE_CGAL "Activate fast element mapping using CGAL library." ${USE_CGAL_DEFAULT})
# Use libfbi toggles the inclusion of the Fast Box Intersection library.
# Please note, that either CGAL or libfbi can be used at a time, but not both!
option(USE_LIBFBI "Use fast box intersection library." ${USE_LIBFBI_DEFAULT})
option(USE_FLANN "Enable support for FLANN - Fast Library for Approximate Nearest Neighbors." ${USE_FLANN_DEFAULT})
# we choose here and hande USE_XERCES and USE_LIBXML2 in the conditional options
set(USE_XML_READER ${USE_XML_READER_DEFAULT} CACHE STRING "Choose either 'xerces' or 'libxml2' as xml reader.")
set_property(CACHE USE_XML_READER PROPERTY STRINGS xerces libxml2)

#-----------------------------------------------------------------------------
# Do sanitiy checks after all the settings are done
# Note, that dependend options come next
#-----------------------------------------------------------------------------

if(USE_FEAST AND NOT USE_BLAS_LAPACK STREQUAL "MKL" )
  message(FATAL_ERROR "USE_FEAST requires Intel MKL! Set USE_BLAS_LAPACK=MKL")
endif()

if(USE_ILUPACK AND USE_ILUPACK_PARALLEL)
  message(FATAL_ERROR "Either USE_ILUPACK or USE_ILUPACK_PARALLEL. Both cannot be swithed on together")
endif()

if(BUILD_CFSDAT AND NOT USE_CGAL)
  message(FATAL_ERROR "BUILD_CFSDAT requires USE_CGAL, correct your configuration.")
endif()

if(USE_CGAL AND USE_LIBFBI)
  message(FATAL_ERROR "Either USE_CGAL or USE_LIBFBI can be actived, but not both!")
endif()

if(USE_PETSC AND USE_ILUPACK_PARALLEL )
  message(FATAL_ERROR "Using PETSc with ilupack will lead to build error (to be fixed in future), switch either one to  OFF")
endif()

if(USE_PARDISO)
  # PARDISO form MKL also requires BLAS/LAPACK from MKL
  # with CFS_PARDISO=SCHENK can be used with OPENBLAS or NETLIB
  if(CFS_PARDISO STREQUAL "MKL" AND NOT USE_BLAS_LAPACK STREQUAL "MKL")
    message(FATAL_ERROR "CFS_PARDISO=MKL requires USE_BLAS_LAPACK=MKL. Disable USE_PARDISO or set CFS_PARDSIO to SCHENK if not using USE_BLAS_LAPACK=MKL.")
  endif()
endif(USE_PARDISO)

if(USE_FEAST)
  # the current FEAST implementation can be from MKL (only for symmetric versions)
  # or from the origninal online source (COMMUNITY).
  # Both require MKL. To get rid of this dependency one would have to implement
  # the reverse communication interface to include another lineare solver (default=pardiso).
  # note that feast has an option to be built w/o mkl
  if(NOT USE_BLAS_LAPACK STREQUAL "MKL")
    message(FATAL_ERROR "USE_FEAST requires Intel MKL! Set USE_BLAS_LAPACK=MKL")
  endif()
endif()

# Since MKL Pardiso only works with MKL BLAS/LAPACK, we force that here.
if(USE_PARDISO AND CFS_PARDISO STREQUAL "MKL")
  if(NOT USE_BLAS_LAPACK STREQUAL "MKL")
    message(FATAL_ERROR "with USE_PARDISO and CFS_PARDISO=MKL one also needs USE_BLAS_LAPACK=MKL") 
  endif()
endif()

if(BUILD_TESTING AND NOT EXISTS ${TESTSUITE_DIR}/TESTSUIT)
  message(WARNING "BUILD_TESTING is on but TESTSUITE_DIR seems invalid, deactivate BUILD_TESTING")
  set(BUILD_TESTING OFF FORCE) # automatism is legaciy behavior
endif()

if(BUILD_TESTING AND NOT BUILD_CFSTOOL)
  message(FATAL_ERROR "BUILD_CFSTOOL required for BUILD_TESTING")
endif()

if(NOT USE_METIS)
  if(USE_IPOPT OR USE_ILUPACK OR USE_SUITESPARSE)
    message(FATAL_ERROR "USE_METIS required for IPOPT, ILUPACK and SUITESPARSE")
  endif() 
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

if(USE_GMV)
  set(USE_GMV_INPUT ON CACHE INTERNAL "Use GMV files as input for simulation. Debugging may not work!")
  set(USE_GMV_OUTPUT ON CACHE INTERNAL "Use GMV files as output for simulation.")
endif(USE_GMV)

cmake_dependent_option(USE_FEAST_COMMUNITY "Use the community version of FEAST" ON USE_FEAST OFF)
cmake_dependent_option(USE_FEAST_COMMUNITY_PRECOMPILED "Use the precompiled FEAST libraries" OFF USE_FEAST_COMMUNITY OFF)

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
mark_as_advanced(CMAKE_BUILD_TYPE)
mark_as_advanced(CMAKE_BACKWARDS_COMPATIBILITY)
mark_as_advanced(CMAKE_INSTALL_PREFIX)
mark_as_advanced(EXPLICIT_TEMPLATE_INSTANTIATION)
mark_as_advanced(BUILD_HDFVIEW)
mark_as_advanced(BUILD_H5TOOL)
mark_as_advanced(BUILD_GHOST)
mark_as_advanced(BUILD_HWLOC)
mark_as_advanced(XMLSCHEMA)
mark_as_advanced(USE_XERCES)
mark_as_advanced(USE_LIBXML2)
mark_as_advanced(USE_ILUPACK_PARALLEL)
mark_as_advanced(USE_PHIST_CG)
mark_as_advanced(USE_FEAST_COMMUNITY)
mark_as_advanced(USE_FEAST_COMMUNITY_PRECOMPILED)
if(NOT CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  mark_as_advanced(CFS_ICC_GCC_VERSION) # shall be seen in intel case
endif()

#-----------------------------------------------------------------------------
# Prepare further build
#-----------------------------------------------------------------------------

# Find executables of a few required programs
include("cmake_modules/FindPrograms.cmake")

# Set compiler flags according to output target and debug setting
include("cmake_modules/compiler.cmake")

# Determine on what UNIX/Linux/Windows distribution we are building on
include("cmake_modules/distro.cmake")

# Determine if some Fortran system libs are installed (blas, lapack, g2c...)
include("cmake_modules/CheckFortranRuntime.cmake")

# Perform some configure checks
include("cmake_modules/ConfigureChecks.cmake")

# Set path where shared/static libraries should go
if(CFS_ARCH STREQUAL "X86_64")
  set(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib64/${CFS_ARCH_STR} CACHE INTERNAL "Single output directory for building all libraries.")
else()
  set(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR} CACHE INTERNAL "Single output directory for building all libraries.")
endif()

# Set path where executables should go
set(EXECUTABLE_OUTPUT_PATH ${CFS_BINARY_DIR}/bin/${CFS_ARCH_STR}
  CACHE INTERNAL "Single output directory for building all executables.")

# Create symlinks to lib & bin dirs allowing the build to run on (newer) linux
if(UNIX)
  # lib dir
  get_filename_component(wd ${LIBRARY_OUTPUT_PATH} DIRECTORY)
  add_custom_target(symlink_lib_dir ALL
     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CFS_ARCH_STR} "LINUX_${CFS_ARCH}"
     WORKING_DIRECTORY ${wd})
  # bin dir
  add_custom_target(symlink_bin_dir ALL
     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CFS_ARCH_STR} "LINUX_${CFS_ARCH}"
     WORKING_DIRECTORY "${CFS_BINARY_DIR}/bin")
endif(UNIX)

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

# Set date of configure run. TODO: make platform independent!!!
TODAY(CFS_CONF_DATE)

# Determine URL of subversion repository, global revision of working copy
# and other informations from working copy
include("cmake_modules/DetermineWorkingCopyInfos.cmake")

# Configure Dart testing support.  This should be done before any
# message(FATAL_ERROR ...) commands are invoked (???) and after the languages
# for the project have been defined since CMake needs to know about the 
# language-specific build tools.
# Switch support for testing ON if the default testsuite dir is used, i.e.
# the git submodule 'Testing'
if(BUILD_TESTING)
  # check if TESTSUITE_DIR is correct (contains TESTSUIT)
    if(EXISTS ${TESTSUITE_DIR}/TESTSUIT)
      add_subdirectory(${TESTSUITE_DIR} "${CMAKE_CURRENT_BINARY_DIR}/testsuite")
    else()
      set(MSG "TESTSUITE_DIR '${TESTSUITE_DIR}' does not contain a directoy 'TESTSUIT'.\n")
      set(MSG "${MSG}Please specify a proper location!\n")
      set(MSG "${MSG}Did you clone with 'git clone --recurse-submodules'?\n")
      set(MSG "${MSG}If not go to your source directory and do\n")
      set(MSG "${MSG}   'git submodule update --init --recursive'\n")
      set(MSG "${MSG}   'git pull --recurse-submodules'\n")
      message(FATAL_ERROR "${MSG}")
    endif()
  # set a build name
  set(BUILDNAME "${CFS_WC_REVISION} on ${DIST} ${REV} with ${CXX_ID} ${CXX_VERSION}" CACHE STRING "Name of build on the dashboard")
  mark_as_advanced(BUILDNAME)
    if(GITLAB_TEST_TAG) # configure with -DGITLAB_TEST_TAG:STRING=`date +%Y%m%d-%H%M` to activate
    # get the build host
    execute_process(COMMAND ${CFS_BINARY_DIR}/bin/cfs --version --noColor
                    COMMAND grep CFS_BUILD_HOST
                    COMMAND awk "{print $2}"
                    RESULT_VARIABLE var
                    OUTPUT_VARIABLE ORIGINAL_BUILD_HOST)
    string(STRIP "${ORIGINAL_BUILD_HOST}" ORIGINAL_BUILD_HOST) # important for comparison !
    #message("CFS_BUILD_HOST='${CFS_BUILD_HOST}' ORIGINAL_BUILD_HOST='${ORIGINAL_BUILD_HOST}'")
    if(NOT (CFS_BUILD_HOST STREQUAL "${ORIGINAL_BUILD_HOST}") )
      # extract the short name (first 2 domain name components)
      string(REGEX MATCH "[^.]+([.][^.]+)?" ORIGINAL_BUILD_HOST "${ORIGINAL_BUILD_HOST}")
      message("current build host CFS_BUILD_HOST='${CFS_BUILD_HOST}'")
      message("build host where the binary is from ORIGINAL_BUILD_HOST='${ORIGINAL_BUILD_HOST}'")
    endif()
  endif(GITLAB_TEST_TAG)
  enable_testing()
endif(BUILD_TESTING)
# add a packing test
if(BUILD_CPACK_BINARY)
  include(CTest)
  add_test(Install_Pack_Linux_TGZ cpack -G TGZ --verbose)
  enable_testing()
endif()

if(NOT BUILD_CPACK_BINARY AND NOT BUILD_TESTING)
  # remove Ctest config if testing is off
  # will probably not be necessay any mode with https://gitlab.kitware.com/cmake/cmake/-/merge_requests/258
  file(REMOVE "${CFS_BINARY_DIR}/CTestTestfile.cmake")
endif()

# ------------------------------------------------------------------------
# Add stuff according to the settings
# ------------------------------------------------------------------------

# Make sure the compiler can find include files from CFS.
include_directories(
  ${Boost_INCLUDE_DIR} 
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
if(USE_ILUPACK OR USE_ILUPACK_PARALLEL)
  add_subdirectory(source/OLAS/external/ilupack)
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
if(USE_GMV)
  add_subdirectory(source/DataInOut/SimInOut/gmv)
endif(USE_GMV)
if(USE_GMSH)
  add_subdirectory(source/DataInOut/SimInOut/gmsh)
endif(USE_GMSH)
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
add_subdirectory(source/DataInOut/SimInOut/TextOutput)
add_subdirectory(source/DataInOut/SimInOut/InfoResultOutput)
add_subdirectory(source/DataInOut/SimInOut/Streaming)
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
add_subdirectory(share/doc) #building documentation (see 'make help')
add_subdirectory(share/scripts)
add_subdirectory(share/mkmesh)
add_subdirectory(share/xml)
add_subdirectory(share/xsl)
add_subdirectory(share/python)

# set cpack config options
# see https://cmake.org/Wiki/CMake:CPackConfiguration
string(STRIP "${CFS_GIT_COMMIT}" CPACK_PACKAGE_VERSION)
if(CFS_WC_MODIFIED)
  set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}-modified")
endif(CFS_WC_MODIFIED)

# choose what to install
# see https://cmake.org/cmake/help/v3.9/command/install.html#installing-directories

if(UNIX)
  # tailing / in source leads to rename to destionation
  #install(PROGRAMS "${CFS_BINARY_DIR}/bin" DESTINATION "." PATTERN "[A-Z]*" EXCLUDE)
  install(DIRECTORY "${EXECUTABLE_OUTPUT_PATH}/" DESTINATION "bin/LINUX_${CFS_ARCH}" USE_SOURCE_PERMISSIONS)
  install(DIRECTORY "${LIBRARY_OUTPUT_PATH}/" DESTINATION "lib64/LINUX_${CFS_ARCH}" USE_SOURCE_PERMISSIONS PATTERN "*.a" EXCLUDE) # lib64 only works for 64bit ...
else(UNIX)
  message("Install currently restricted to UNIX")
endif(UNIX)
install(DIRECTORY "${CFS_BINARY_DIR}/share/" DESTINATION "share" USE_SOURCE_PERMISSIONS PATTERN "CMakeFiles" EXCLUDE PATTERN "Makefile" EXCLUDE PATTERN "*.dir" EXCLUDE)

include(CPack)
