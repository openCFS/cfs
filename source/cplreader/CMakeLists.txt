# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: LSE.
# Date:   05.09.2003.
# --------------------------------------------------------------------------
#

#-------------------------------------------------------------------------------
# Don't build anything unless the version of CMake is high enough.
#-------------------------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

#-------------------------------------------------------------------------------
# This is project cplreader
#-------------------------------------------------------------------------------
PROJECT(cplreader C CXX Fortran)

#-------------------------------------------------------------------------------
# Determine if we part of CFS build or standalone
#-------------------------------------------------------------------------------
IF(CFS_VERSION)
  SET(CPLREADER_STANDALONE 0)
ELSE(CFS_VERSION)
  SET(CPLREADER_STANDALONE 1)
ENDIF(CFS_VERSION)

IF(CPLREADER_STANDALONE)
  OPTION(CPLREADER_CFX "" OFF)
ELSE(CPLREADER_STANDALONE)
  OPTION(CPLREADER_CFX "" ON)  
ENDIF(CPLREADER_STANDALONE)

OPTION(CPLREADER_FASTEST "" ON)  
OPTION(CPLREADER_OPENFOAM "" ON)  

IF(NOT CPLREADER_STANDALONE)
  MARK_AS_ADVANCED(CPLREADER_CFX)
  MARK_AS_ADVANCED(CPLREADER_FASTEST)
  MARK_AS_ADVANCED(CPLREADER_OPENFOAM)
ELSE(NOT CPLREADER_STANDALONE)
  SET(BOOST_PROGRAM_OPTIONS_LIB "/usr/lib64/libboost_program_options.a"
    CACHE FILE_PATH "BOOST_PROGRAM_OPTIONS_LIB")
  SET(BOOST_SERIALIZATION_LIB "/usr/lib64/libboost_serialization.a"
    CACHE FILE_PATH "BOOST_SERIALIZATION_LIB")
  SET(BOOST_FILESYSTEM_LIB "/usr/lib64/libboost_filesystem.a"
    CACHE FILE_PATH "BOOST_FILESYSTEM_LIB")
  
  SET(EXECUTABLE_OUTPUT_PATH ${cplreader_BINARY_DIR}/bin CACHE INTERNAL
    "Single output directory for building all executables.")
 
ENDIF(NOT CPLREADER_STANDALONE)


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${MPCCI_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR})

CONFIGURE_FILE(cplreaderdefs.hh.in ${CMAKE_CURRENT_BINARY_DIR}/cplreaderdefs.hh)

#-------------------------------------------------------------------------------
# CFX specific settings
#-------------------------------------------------------------------------------
IF(CPLREADER_CFX)
  # Set a sane default value for the CFS++ schema path
  SET(CFX_ROOT "/home/data/programs/ansys/v110/CFX" CACHE STRING
    "Path where CFX is located.")
  #SET(CFX_ROOT "/home/data/programs/ansys/CFX-10.0")
  #SET(CFX_ROOT "/home/simon/CFS_BUILD/CFX")
  MARK_AS_ADVANCED(CFX_ROOT)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CFX_ROOT}/include)

  IF(CFS_ARCH STREQUAL "I386")
    LINK_DIRECTORIES(${CFX_ROOT}/lib/linux)
  ELSE(CFS_ARCH STREQUAL "I386")
    LINK_DIRECTORIES(${CFX_ROOT}/lib/linux-amd64)
  ENDIF(CFS_ARCH STREQUAL "I386")
ENDIF(CPLREADER_CFX)

SET(CPPSRC 
  cplreader.cc
  filereader.cc
  MpCCIexch.cc
  params.cc
)

IF(NOT CPLREADER_STANDALONE)
  SET(CPPSRC
    ${CPPSRC}
    integlib/elemIntegr.cc
    integlib/integr.cc
    )
ENDIF(NOT CPLREADER_STANDALONE)

IF(CPLREADER_CFX)
  SET(CPPSRC
    ${CPPSRC}
    CFX/glibchack.cc
    CFX/filereader_CFX.cc
    CFXexport/FileReader_CFXexport.cc
    )

  SET(FSRC
    CFX/OpenFile.F
    CFX/CloseFile.F
    CFX/iocntf_wrap.F
    CFX/RedHdr.F
    CFX/ReadLong.F
    CFX/RedSht.F
    CFX/tools.F
    )
ENDIF(CPLREADER_CFX)

IF(CPLREADER_FASTEST)
  SET(CPPSRC
    ${CPPSRC}
    FASTEST/filereader_FASTEST.cc
    )
    
    CONFIGURE_FILE(FASTEST/FASTEST-FileRenamer.pl.in ${CFS_BINARY_DIR}/bin/FASTEST-FileRenamer.pl @ONLY)
ENDIF(CPLREADER_FASTEST)

IF(CPLREADER_OPENFOAM)
  SET(CPPSRC
    ${CPPSRC}
    OPENFOAM/filereader_OPENFOAM.cc
    )
ENDIF(CPLREADER_OPENFOAM)

ADD_EXECUTABLE(cplreader ${CPPSRC} ${FSRC})

IF(NOT CPLREADER_STANDALONE)
  SET(TARGET_LL ${TARGET_LL}
    siminouthdf5
    forms
    elements
    matrix
    utils
    general
    )
ENDIF(NOT CPLREADER_STANDALONE)

IF(CPLREADER_CFX)
  SET(TARGET_LL
    ${TARGET_LL}
    io # CFX I/O library
    )
ENDIF(CPLREADER_CFX)

IF(NOT CPLREADER_STANDALONE)
  SET(TARGET_LL ${TARGET_LL}
    ${BOOST_PROGRAM_OPTIONS_LIB}
    # DIRTY, DIRTYYYY HACK!
    datainout
    simouttext
    datainout
    )
ELSE(NOT CPLREADER_STANDALONE)
  SET(TARGET_LL ${TARGET_LL}
    ${BOOST_PROGRAM_OPTIONS_LIB}
    ${BOOST_SERIALIZATION_LIB}
    ${BOOST_FILESYSTEM_LIB}
  )
ENDIF(NOT CPLREADER_STANDALONE)

TARGET_LINK_LIBRARIES(cplreader
  ${TARGET_LL}
  )

IF(CPLREADER_CFX)
  SET_TARGET_PROPERTIES(cplreader PROPERTIES
    LINK_FLAGS -Wl,--wrap,__ctype_b
    LINKER_LANGUAGE CXX)
ENDIF(CPLREADER_CFX)
