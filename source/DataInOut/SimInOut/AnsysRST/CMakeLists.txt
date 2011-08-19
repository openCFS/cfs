IF(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")
  IF(CFS_FORTRAN_COMPILER_VER LESS 4.1.0)
    MESSAGE(SEND_ERROR
      "The ANSYS RST writer requires a Fortran 95 compiler!")
  ENDIF(CFS_FORTRAN_COMPILER_VER LESS 4.1.0)
ENDIF(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")

IF(CFS_ARCH STREQUAL "I386")
  SET(CMAKE_Fortran_FLAGS "-safe-cray-ptr -fcray-pointer -DLINUXIA32_SYS")
  
  # Set possible paths for Intel Fortran libraries at LSE and Klagenfurt
  IF(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")  
    LINK_DIRECTORIES(
    /opt/pckg/intel/fc/10.1.015/lib
    /opt/pckg/intel/Compiler/11.0/081/lib/ia32
    /opt/pckg/intel/Compiler/11.1/069/lib/ia32
    /home/data/programs/intel/Compiler/fc/10.0.025/lib
    /home/data/programs/intel/Compiler/fc/10.1.018/lib
    /home/data/programs/intel/Compiler/11.0/081/lib/ia32
    )
  ENDIF(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")  
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  SET(CMAKE_Fortran_FLAGS "-safe-cray-ptr -fcray-pointer -DLINUXOP64_SYS")

  # Set possible paths for Intel Fortran libraries at LSE and Klagenfurt
  IF(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")  
    LINK_DIRECTORIES(
    /opt/pckg/intel/fce/10.0.025/lib
    /opt/pckg/intel/fce/10.1.015/lib
    /opt/pckg/intel/Compiler/11.0/081/lib/intel64
    /opt/pckg/intel/Compiler/11.1/069/lib/intel64
    /home/data/programs/intel/Compiler/fce/10.0.025/lib
    /home/data/programs/intel/Compiler/fce/10.1.018/lib
    /home/data/programs/intel/Compiler/11.0/081/lib/intel64
    )
  ENDIF(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")  
    
ENDIF(CFS_ARCH STREQUAL "X86_64")

IF(DEBUG)
  SET(CMAKE_Fortran_FLAGS "-g ${CMAKE_Fortran_FLAGS}")
ENDIF(DEBUG)

SET(SIMOUTANSYSRST_SRCS simOutputRST.cc)

ADD_LIBRARY(simoutansysrst STATIC ${SIMOUTANSYSRST_SRCS})

TARGET_LINK_LIBRARIES(simoutansysrst
  CFSAnsysBinlibIface
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_SYSTEM_LIB}
  )

ADD_LIBRARY(CFSAnsysBinlibIface SHARED AnsysBinlibIface.cc)
TARGET_LINK_LIBRARIES(CFSAnsysBinlibIface
  CFSDynLd
  )

# Dynamic library interface for ANSYS 10.0 libbin
SUBDIRS(AnsysBinlibIfaceRev100)

# Dynamic library interface for ANSYS 11.0 libbin
SUBDIRS(AnsysBinlibIfaceRev110)

# At the moment we only have libs for x86_64
IF(CFS_ARCH STREQUAL "X86_64")
  # Dynamic library interface for ANSYS 12.0 libbin
  SUBDIRS(AnsysBinlibIfaceRev120)

  # Dynamic library interface for ANSYS 13.0 libbin
  SUBDIRS(AnsysBinlibIfaceRev130)
ENDIF(CFS_ARCH STREQUAL "X86_64")

