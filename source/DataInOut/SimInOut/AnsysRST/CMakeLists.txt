
IF(CFS_ARCH STREQUAL "I386")
  SET(CMAKE_Fortran_FLAGS "-safe-cray-ptr -fcray-pointer -DLINUXIA32_SYS")
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  SET(CMAKE_Fortran_FLAGS "-safe-cray-ptr -fcray-pointer -DLINUXOP64_SYS")
ENDIF(CFS_ARCH STREQUAL "X86_64")

IF(DEBUG)
  SET(CMAKE_Fortran_FLAGS "-g ${CMAKE_Fortran_FLAGS}")
ENDIF(DEBUG)

SET(SIMOUTANSYSRST_SRCS simOutputRST.cc)

ADD_LIBRARY(simoutansysrst STATIC ${SIMOUTANSYSRST_SRCS})

TARGET_LINK_LIBRARIES(simoutansysrst
  CFSAnsysBinlibIface
  ${BOOST_FILESYSTEM_LIB}
  )

ADD_LIBRARY(CFSAnsysBinlibIface SHARED AnsysBinlibIface.cc)
TARGET_LINK_LIBRARIES(CFSAnsysBinlibIface
  CFSDynLd
  )

# Dynamic library interface for ANSYS 10.0 libbin
SUBDIRS(AnsysBinlibIfaceRev100)

# Dynamic library interface for ANSYS 11.0 libbin
SUBDIRS(AnsysBinlibIfaceRev110)