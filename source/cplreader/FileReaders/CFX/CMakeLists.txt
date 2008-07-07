#-------------------------------------------------------------------------------
# CFX specific settings
#-------------------------------------------------------------------------------

# Set a sane default value for the CFS++ schema path
SET(CFX_ROOT "/home/data/programs/ansys/v110/CFX" CACHE STRING
  "Path where CFX is located.")
#SET(CFX_ROOT "/home/data/programs/ansys/CFX-10.0")
#SET(CFX_ROOT "/home/simon/CFS_BUILD/CFX")
MARK_AS_ADVANCED(CFX_ROOT)

INCLUDE_DIRECTORIES(${CFX_ROOT}/include)

IF(CFS_ARCH STREQUAL "I386")
  SET(CFX_IO_LIB ${CFX_ROOT}/lib/linux/libio.a)
ELSE(CFS_ARCH STREQUAL "I386")
  SET(CFX_IO_LIB ${CFX_ROOT}/lib/linux-amd64/libio.a)
ENDIF(CFS_ARCH STREQUAL "I386")

SET(CPLREADER_CFX_SRCS
  FileReader_CFX.cc
  glibchack.cc
  OpenFile.F
  CloseFile.F
  iocntf_wrap.F
  RedHdr.F
  ReadLong.F
  RedSht.F
  tools.F
)

ADD_LIBRARY(cplreader_cfx STATIC ${CPLREADER_CFX_SRCS})

TARGET_LINK_LIBRARIES(cplreader_cfx
  ${CFX_IO_LIB} # CFX I/O library
)