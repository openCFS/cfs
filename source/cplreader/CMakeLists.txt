# ------------- MAKEFILE FOR FILEREADER FOR COUPLING -----------------------
# Author: LSE.
# Date:   05.09.2003.
# --------------------------------------------------------------------------
#

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


SET(CPPSRC 
  cplreader.cc
  filereader.cc
  MpCCIexch.cc
  params.cc
  CFX/glibchack.cc
  CFX/filereader_CFX.cc
  FASTEST/filereader_FASTEST.cc
  integlib/elemIntegr.cc
  integlib/integr.cc
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

ADD_EXECUTABLE(cplreader ${CPPSRC} ${FSRC})

SET(TARGET_LL ${TARGET_LL}
  forms
  elements
  matrix
  utils
  general
  io # CFX I/O library
  boost_program_options
# DIRTY, DIRTYYYY HACK!
  datainout
  simouttext
  datainout
)

TARGET_LINK_LIBRARIES(cplreader
  ${TARGET_LL}
#  boost_serialization
#  boost_filesystem
  )

SET_TARGET_PROPERTIES(cplreader PROPERTIES
  LINK_FLAGS -Wl,--wrap,__ctype_b
  LINKER_LANGUAGE CXX)

CONFIGURE_FILE(FASTEST/FASTEST-FileRenamer.pl.in ${CFS_BINARY_DIR}/bin/FASTEST-FileRenamer.pl @ONLY)