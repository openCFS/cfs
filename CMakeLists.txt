CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#-------------------------------------------------------------------------------
# The name of our project is "CFS".  CMakeLists files in this
# project can refer to the root source directory of the project as
# ${CFS_SOURCE_DIR} and to the root binary directory of the project
# as ${CFS_BINARY_DIR}. We use C, C++ and Fortran as languages.
#-------------------------------------------------------------------------------
PROJECT(CFS C CXX Fortran)

#-------------------------------------------------------------------------------
# Don't build anything unless the version of CMake is high enough.
# The matching ELSE/ENDIF should be the last lines in the file.
#-------------------------------------------------------------------------------
IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.3)


  #-----------------------------------------------------------------------------
  # CFS version number.  An even minor number corresponds to releases.
  #-----------------------------------------------------------------------------
  SET(CFS_MAJOR_VERSION 0)
  SET(CFS_MINOR_VERSION 6)
  SET(CFS_BUILD_VERSION 0)
  SET(CFS_VERSION
    "${CFS_MAJOR_VERSION}.${CFS_MINOR_VERSION}.${CFS_BUILD_VERSION}")



  #-----------------------------------------------------------------------------
  # Switch on/off external libraries and other options for compile
  #-----------------------------------------------------------------------------
  OPTION(USE_LAPACK "Use the LAPACK linear algebra library." ON)
  OPTION(USE_ARPACK "Use the ARPACK eigensolver library." ON)
  OPTION(USE_PARDISO "Use Pardiso" ON)
  OPTION(USE_METIS "Use Metis Graph partitioning." ON)
  OPTION(USE_ILUPACK "Use Preconditioners from ILUPACK." ON)
  OPTION(USE_BLAS "Use the basic linear algebra systems library." ON)
  OPTION(USE_XERCES "Use the Xerces-C XML parser library." ON)
  OPTION(USE_TCL "Use Tool Command Language (TCL) for scripting." ON)
  OPTION(USE_PYTHON "Use Python for scripting." OFF)
  OPTION(USE_MESH "Use the mesh file reader for input from ANSYS/GiD." ON)
  OPTION(USE_GIDPOST "Use the GiD library for simulation output." ON)
  OPTION(USE_HDF5 "Use HDF5 for simulation in/output." OFF)
  OPTION(USE_GMV "Use General Mesh Viewer (GMV) for simulation in/output." ON)
  OPTION(USE_UNV "Use Universal File for simulation output." ON)
  OPTION(USE_ANSYSRST "Use ANSYS .RST format for simulation output." OFF)
  OPTION(USE_MPCCI "Enable support for the mesh parallel coupling interface."
    OFF)

  OPTION(DEBUG "Build debug version of the code." ON)
  OPTION(PROFILING "Turn on profiling (DEBUG/RELEASE)." OFF)
  OPTION(CONVERTER "Build CFS++ converter tool (cfstool)." OFF)

  # Set a sane default value for the CFS++ schema path
  SET(XMLSCHEMA "${CFS_SOURCE_DIR}/XML" CACHE STRING
    "Path where the XML schema for the input files is located.")


  #-----------------------------------------------------------------------------
  # Mark some of the options as advanced to clear up the standard selection
  #-----------------------------------------------------------------------------
  MARK_AS_ADVANCED(USE_UNV)
  MARK_AS_ADVANCED(USE_GMV)
  MARK_AS_ADVANCED(USE_MESH)
  MARK_AS_ADVANCED(USE_XERCES)
  MARK_AS_ADVANCED(CMAKE_BUILD_TYPE)
  MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
  MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)



  #-----------------------------------------------------------------------------
  # Set some variables according to build type DEBUG / RELEASE
  #-----------------------------------------------------------------------------
  IF(DEBUG)
    SET(CMAKE_BUILD_TYPE "DEBUG" CACHE STRING
      "Only DEBUG and RELEASE types are valid." FORCE)
    SET(_DEBUG 1)

    #---------------------------------------------------------------------------
    # Defines related to tracing of functions and memory usage
    #---------------------------------------------------------------------------
    SET(TRACE 1)
    SET(TRACE_ALL 1)
    SET(MEMTRACE 1)

    #---------------------------------------------------------------------------
    # Check for certain error conditions (in DEBUG mode)
    #---------------------------------------------------------------------------
    SET(CHECK_INDEX 1)
    SET(CHECK_INITIALIZED 1)
    SET(CHECK_MEM_ALLOC 1)
    SET(CHECK_MEMORY 1)
    SET(CHECK_TYPE_CASTS 1)

  ELSE(DEBUG)
    SET(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
      "Only DEBUG and RELEASE types are valid." FORCE)
    SET(ENABLE_PROFILING OFF CACHE BOOL
      "Turn on profiling (DEBUG only)."FORCE)

    #---------------------------------------------------------------------------
    # Activate expression templates (in RELEASE mode)
    #---------------------------------------------------------------------------
    SET(EXPR_TEMPLATES 1)
  ENDIF(DEBUG)

  #-----------------------------------------------------------------------------
  # Set compiler flags according to output target and debug setting
  # ( if necessesary )
  #-----------------------------------------------------------------------------
  INCLUDE("cmake_modules/compiler.cmake")

  #-----------------------------------------------------------------------------
  # Determine if some Fortran system libs are installed (blas, lapack, g2c...)
  #-----------------------------------------------------------------------------
  INCLUDE("cmake_modules/FindFortranLibs.cmake")

  #-----------------------------------------------------------------------------
  # Determine on what UNIX/Linux/Windows distribution we are building on
  #-----------------------------------------------------------------------------
  INCLUDE("cmake_modules/distro.cmake")

  #-----------------------------------------------------------------------------
  # Find executables of a few required programs
  #-----------------------------------------------------------------------------
  INCLUDE("cmake_modules/FindPrograms.cmake")


  #-----------------------------------------------------------------------------
  # Set path where shared/static libraries should go
  #-----------------------------------------------------------------------------
  IF(CFS_ARCH STREQUAL "X86_64")
    SET(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib64 CACHE INTERNAL
      "Single output directory for building all libraries.")
  ELSE(CFS_ARCH STREQUAL "X86_64")
    SET(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib CACHE INTERNAL
      "Single output directory for building all libraries.")
  ENDIF(CFS_ARCH STREQUAL "X86_64")

  #-----------------------------------------------------------------------------
  # Set path where executables should go
  #-----------------------------------------------------------------------------
  SET(EXECUTABLE_OUTPUT_PATH ${CFS_BINARY_DIR}/bin CACHE INTERNAL
    "Single output directory for building all executables.")


  #-----------------------------------------------------------------------------
  # Find CFSDEPS, i.e. gidpost, HDF5, pardiso, arpack...
  #-----------------------------------------------------------------------------
  INCLUDE("cmake_modules/FindCFSDEPS.cmake")
  IF(NOT CFSDEPS_FOUND)
    MESSAGE(FATAL_ERROR "CFSDEPS could not be found!")
  ENDIF(NOT CFSDEPS_FOUND)

  #-----------------------------------------------------------------------------
  # Find Boost
  #-----------------------------------------------------------------------------
  FIND_PACKAGE(Boost)
  MARK_AS_ADVANCED(Boost_INCLUDE_DIR)
  IF(NOT Boost_FOUND)
    MESSAGE(SEND_ERROR "Boost could not be found!")
    SUGGEST_INSTALL_PCKG()
  ENDIF(NOT Boost_FOUND)

  #-----------------------------------------------------------------------------
  # If USE_PYTHON option is defined find Python library
  #-----------------------------------------------------------------------------
  IF(USE_PYTHON)
    FIND_PACKAGE(PythonLibs)
    IF(NOT PYTHON_LIBRARY)
      MESSAGE(SEND_ERROR "Python library could not be found!")
      SUGGEST_INSTALL_PCKG()
    ENDIF(NOT PYTHON_LIBRARY)

    SET(USE_SCRIPTING 1)
    SET(USE_SCRIPTING_PYTHON 1)

  ENDIF(USE_PYTHON)
  
  #-----------------------------------------------------------------------------
  # If USE_TCL option is defined find TCL library
  #-----------------------------------------------------------------------------
  IF(USE_TCL)
    FIND_PACKAGE(TCL)
    IF(NOT TCL_LIBRARY)
      MESSAGE(SEND_ERROR "TCL library could not be found!")
      SUGGEST_INSTALL_PCKG()
    ENDIF(NOT TCL_LIBRARY)

    SET(USE_SCRIPTING 1)
    SET(USE_SCRIPTING_TCL 1)

    MARK_AS_ADVANCED(
      TCL_INCLUDE_PATH
      TCL_LIBRARY
      TCL_LIBRARY_DEBUG
      TCL_TCLSH
      TK_INCLUDE_PATH
      TK_LIBRARY
      TK_LIBRARY_DEBUG
      TK_WISH
      )
  ENDIF(USE_TCL)

  #-----------------------------------------------------------------------------
  # If USE_BLAS option is defined find BLAS library
  #-----------------------------------------------------------------------------
  IF(USE_BLAS)
    IF(NOT BLAS_FOUND)
      MESSAGE(SEND_ERROR "BLAS library could not be found!")
      SUGGEST_INSTALL_PCKG()
    ENDIF(NOT BLAS_FOUND)
  ENDIF(USE_BLAS)

  #-----------------------------------------------------------------------------
  # If USE_LAPACK option is defined find LAPACK library
  #-----------------------------------------------------------------------------
  IF(USE_LAPACK)
    IF(NOT BLAS_FOUND)
      MESSAGE(SEND_ERROR "BLAS library could not be found!")
      SUGGEST_INSTALL_PCKG()
    ENDIF(NOT BLAS_FOUND)
    IF(NOT LAPACK_FOUND)
      MESSAGE(SEND_ERROR "LAPACK library could not be found!")
      SUGGEST_INSTALL_PCKG()
    ENDIF(NOT LAPACK_FOUND)
  ENDIF(USE_LAPACK)

  #-----------------------------------------------------------------------------
  # If USE_XERCES option is defined find Xerces library
  #-----------------------------------------------------------------------------
  IF(USE_XERCES)
    INCLUDE("cmake_modules/FindXerces.cmake")
    IF(NOT XERCES_FOUND)
      MESSAGE(SEND_ERROR "Xerces library could not be found!")
      SUGGEST_INSTALL_PCKG()
    ENDIF(NOT XERCES_FOUND)
  ENDIF(USE_XERCES)

  #-----------------------------------------------------------------------------
  # Set path of XML schema. If env variable XMLSCHEMA is set, use it!
  #-----------------------------------------------------------------------------
  SET(MYXMLSCHEMA $ENV{XMLSCHEMA})
  IF(NOT ${MYXMLSCHEMA} STREQUAL "")
    SET(XMLSCHEMA ${MYXMLSCHEMA} CACHE STRING
      "Path where the XML schema for the input files is located."
      FORCE)
  ENDIF(NOT ${MYXMLSCHEMA} STREQUAL "")

  #-----------------------------------------------------------------------------
  # Find MpCCI / MPICH
  #-----------------------------------------------------------------------------
  IF(USE_MPCCI)
    INCLUDE("cmake_modules/FindMpCCI.cmake")    
    IF(NOT MPCCI_FOUND)
      MESSAGE(FATAL_ERROR "MpCCI could not be found!")
    ENDIF(NOT MPCCI_FOUND)
  ENDIF(USE_MPCCI)



  #-----------------------------------------------------------------------------
  # Set date of configure run. TODO: make platform independent!!!
  #-----------------------------------------------------------------------------
  EXEC_PROGRAM("date"
    ARGS +%D
    OUTPUT_VARIABLE CFS_CONF_DATE
    RETURN_VALUE RETVAL)

  #-----------------------------------------------------------------------------
  # Set hostname of build host
  #-----------------------------------------------------------------------------
  EXEC_PROGRAM("hostname"
    ARGS
    OUTPUT_VARIABLE CFS_BUILD_HOST
    RETURN_VALUE RETVAL)

  #-----------------------------------------------------------------------------
  # Determine URL of subversion repository
  #-----------------------------------------------------------------------------
  EXEC_PROGRAM("${SVN} info ${CFS_SOURCE_DIR} | head -2 | tail -1 | cut -d' ' -f2"
    ARGS
    OUTPUT_VARIABLE CFS_SUBVERSION_REPOS
    RETURN_VALUE RETVAL)

  #-----------------------------------------------------------------------------
  # Determine global revision of working copy
  #-----------------------------------------------------------------------------
  EXEC_PROGRAM(${SVNVERSION}
    ARGS ${CFS_SOURCE_DIR}
    OUTPUT_VARIABLE CFS_SUBVERSION_REV
    RETURN_VALUE RETVAL)



  #-----------------------------------------------------------------------------
  # Make sure the compiler can find include files from CFS.
  #-----------------------------------------------------------------------------
  INCLUDE_DIRECTORIES(
    ${Boost_INCLUDE_DIR} 
    ${CFSDEPS_INCLUDE_DIR} 
    ${CFS_BINARY_DIR}/include
    ${CFS_SOURCE_DIR}/source
    ${CFS_SOURCE_DIR}/source/OLAS)

  #-----------------------------------------------------------------------------
  # Make sure the linker can find the CFS base libraries once it is built.
  #-----------------------------------------------------------------------------
  LINK_DIRECTORIES(${CFSDEPS_LIBRARY_DIR} ${CFS_BINARY_DIR}/lib)



  #-----------------------------------------------------------------------------
  # Recurse into the "source" subdirectory.  This does not actually
  # cause another cmake executable to run.  The same process will
  # walk through the project's entire directory structure.
  #-----------------------------------------------------------------------------
  SUBDIRS(source/OLAS/algsys)
  SUBDIRS(source/OLAS/graph)
  SUBDIRS(source/OLAS/matvec)
  #  SUBDIRS(source/multigrid)
  SUBDIRS(source/OLAS/precond)
  SUBDIRS(source/OLAS/solver)
  SUBDIRS(source/OLAS/utils)
  #  SUBDIRS(source/main)

  IF(USE_ARPACK)
    SUBDIRS(source/OLAS/external/arpack)
  ENDIF(USE_ARPACK)

  IF(USE_ILUPACK)
    SUBDIRS(source/OLAS/external/ilupack)
  ENDIF(USE_ILUPACK)

  IF(USE_LAPACK)
    SUBDIRS(source/OLAS/external/lapack)
  ENDIF(USE_LAPACK)

  IF(USE_PARDISO)
    SUBDIRS(source/OLAS/external/pardiso)
  ENDIF(USE_PARDISO)

  SUBDIRS(source/CoupledPDE)
  SUBDIRS(source/DataInOut)
  SUBDIRS(source/DataInOut/CommandLine)
  SUBDIRS(source/DataInOut/DefineFiles)
  SUBDIRS(source/DataInOut/Logging)
  SUBDIRS(source/DataInOut/ParamHandling)

  IF(USE_SCRIPTING)
    SUBDIRS(source/DataInOut/Scripting)

    IF(USE_SCRIPTING_PYTHON)
      SUBDIRS(source/DataInOut/Scripting/python)
    ENDIF(USE_SCRIPTING_PYTHON)

    IF(USE_SCRIPTING_TCL)
      SUBDIRS(source/DataInOut/Scripting/tcl)
    ENDIF(USE_SCRIPTING_TCL)
  ENDIF(USE_SCRIPTING)


  IF(USE_MESH)
    SUBDIRS(source/DataInOut/SimInOut/AnsysFile)
  ENDIF(USE_MESH)

  IF(USE_GMV)
    OPTION(USE_GMV_INPUT
      "Use GMV files as input for simulation. Attention! Debugging may not work!"
      OFF)

    OPTION(USE_GMV_OUTPUT
      "Use GMV files as output for simulation."
      ON)

    MARK_AS_ADVANCED(USE_GMV_INPUT USE_GMV_OUTPUT)
    SUBDIRS(source/DataInOut/SimInOut/gmv)
  ENDIF(USE_GMV)

  IF(USE_HDF5)
    SUBDIRS(source/DataInOut/SimInOut/xmdf)
  ENDIF(USE_HDF5)

  IF(USE_GIDPOST)
    SUBDIRS(source/DataInOut/SimInOut/GiD)
  ENDIF(USE_GIDPOST)

  IF(USE_UNV)
    SUBDIRS(source/DataInOut/SimInOut/Unverg)
  ENDIF(USE_UNV)

  IF(USE_ANSYSRST)
    SUBDIRS(source/DataInOut/SimInOut/AnsysRST)
  ENDIF(USE_ANSYSRST)

  SUBDIRS(source/DataInOut/SimInOut/TextOutput)
  SUBDIRS(source/Domain)
  SUBDIRS(source/Domain/GridCFS)
  SUBDIRS(source/Driver)
  SUBDIRS(source/Elements)
  #  SUBDIRS(source/Estimator)
  SUBDIRS(source/Forms)
  SUBDIRS(source/General)
  SUBDIRS(source/Materials)
  SUBDIRS(source/Matrix)

  IF(USE_MPCCI)
    SUBDIRS(source/MpCCIcpl)
  ENDIF(USE_MPCCI)

  SUBDIRS(source/ODEDescr)
  SUBDIRS(source/ODESolve)
  SUBDIRS(source/ParamIdent)
  SUBDIRS(source/PDE)
  SUBDIRS(source/Utils)
  SUBDIRS(source/Utils/mathParser)
  SUBDIRS(source/main)

  IF(CONVERTER)
    SUBDIRS(source/cfstool)
  ENDIF(CONVERTER)

  #-----------------------------------------------------------------------------
  # Build the actual headers with defines from templates.
  #-----------------------------------------------------------------------------
  INCLUDE("cmake_modules/configure_def_headers.cmake")

  #-----------------------------------------------------------------------------
  # Add targets for building documentation (see 'make help')
  #-----------------------------------------------------------------------------
  SUBDIRS(Doc/developer)
  SUBDIRS(Doc/user)

  #-----------------------------------------------------------------------------
  # If the version of CMake was too old, complain and build nothing.
  # These should be the last lines in this file.
  #-----------------------------------------------------------------------------
ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.3)
  MESSAGE(SEND_ERROR
    "This version of CMake is too old to build CFS. "
    "Please upgrade to CMake 2.4.x")
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.3)
