INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR})

IF(MKL_INCLUDE_DIR)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
ENDIF(MKL_INCLUDE_DIR)

# List of source codes for current target.
SET(UNIT_TESTS_SRC
  unittests.cc
#  bsplinetests.cc
#  bicubicinttests.cc
)

# Main CFS++ Dat executable target.
ADD_EXECUTABLE(cfstestbin ${UNIT_TESTS_SRC})

SET(TARGET_LL 
  matvec
  design
  cfsgeneral
  utils
  datainout
  paramh
)

IF(MKL_BLAS_LIB)
  LIST(APPEND TARGET_LL ${MKL_BLAS_LIB})
ENDIF(MKL_BLAS_LIB)

TARGET_LINK_LIBRARIES(cfstestbin ${TARGET_LL})

SET(CFS_LINK_FLAGS "")

if(CFS_CXX_COMPILER_NAME STREQUAL "GCC" AND CFS_CXX_COMPILER_VER VERSION_GREATER "6.0") # there is no >= and also there is no 5.0.0.0
  # the current boost test case is crap :(
  set(CMAKE_CXX_FLAGS "${CFS_SUPPRESSIONS} -Wno-terminate")
endif() 

SET_TARGET_PROPERTIES(cfstestbin PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")
