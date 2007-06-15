SET(OPTIMIZATION_SRCS DesignElement.cc DesignSpace.cc Optimization.cc PiezoSIMP.cc SIMP.cc 
                      TransferFunction.cc Condition.cc)

IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
  LINK_DIRECTORIES(${IPOPT_LINK_DIR})

   SET(TARGET_LL
    ipopt 
    gfortranbegin 
    gcc_s)
   
   SET( OPTIMIZATION_SRCS
    ${OPTIMIZATION_SRCS}
    IPOPT.cc)
ENDIF(USE_IPOPT)

IF(USE_SCPIP)
  INCLUDE_DIRECTORIES(${SCPIP_INCLUDE_DIR})
  LINK_DIRECTORIES(${SCIP_LINK_DIR})

   SET(TARGET_LL 
    ${TARGET_LL}
    scpip)
    
   SET( OPTIMIZATION_SRCS
    ${OPTIMIZATION_SRCS}
    SCPIP.cc SCPIPBase.cc)
ENDIF(USE_SCPIP)


ADD_LIBRARY(optimization STATIC ${OPTIMIZATION_SRCS})

IF(CFS_DISTRO STREQUAL "SUSE" OR CFS_DISTRO STREQUAL "OPENSUSE")
  IF(CFS_DISTRO_VER GREATER 9.3)
    SET(TARGET_LL ${TARGET_LL} gfortran)
  ENDIF(CFS_DISTRO_VER GREATER 9.3)
ENDIF(CFS_DISTRO STREQUAL "SUSE" OR CFS_DISTRO STREQUAL "OPENSUSE")

SET(TARGET_LL ${TARGET_LL} g2c)

TARGET_LINK_LIBRARIES(optimization ${TARGET_LL})

#-----------------------
# SET(SCPIP_SRCS scpip_test.cc)
SET(SCPIP_SRCS SCPIP_HS071.cc SCPIPBase.cc)
ADD_EXECUTABLE(scpip ${SCPIP_SRCS})
TARGET_LINK_LIBRARIES(scpip general utils logging scripting utils matrix datainout domain materials utils
                     lapack blas 
                     /space/fwein/packages/boost_dev/lib/libboost_serialization-gcc.so
                     /space/fwein/packages/scpip/libscpip.a g2c gfortran m)
                           