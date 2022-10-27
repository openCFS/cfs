include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CFS_BINARY_DIR}/include)


if(MKL_INCLUDE_DIR)
  include_directories(${MKL_INCLUDE_DIR})
endif()

# List of source codes for current target.
set(UNIT_TESTS_SRC
  unittests.cc
  muparsertests.cc
  testbed.cc
#  bsplinetests.cc
  #cubicinttests.cc
)

set(TARGET_LL 
  matvec
  design
  cfsgeneral
  utils
  datainout)
#  paramh
#  utils)

if(USE_EMBEDDED_PYTHON)
  # see FindPrograms.cmake on how the PYTHON_ stuff is set
  set(TARGET_LL ${TARGET_LL} ${PYTHON_LIBRARY})
  include_directories(${PYTHON_INCLUDE_DIR})
  include_directories(${PYTHON_SITE_PACKAGES_DIR})
  set(UNIT_TESTS_SRC ${UNIT_TESTS_SRC} embeddedpython.cc)
endif()

# Main CFS++ Dat executable target.
add_executable(cfstest ${UNIT_TESTS_SRC})


IF(MKL_BLAS_LIB)
  LIST(APPEND TARGET_LL ${MKL_BLAS_LIB})
ENDIF(MKL_BLAS_LIB)

TARGET_LINK_LIBRARIES(cfstest ${TARGET_LL})

# this is only necessary for legacy system with glibc <= 2.17 (check with ldd --version)
# new apple system have no more librt dll, probably same for Windows
if(NOT APPLE AND NOT WIN32)
  set(CFS_LINK_FLAGS "-lrt")
endif()

if(CFS_CXX_COMPILER_NAME STREQUAL "GCC" AND CFS_CXX_COMPILER_VER VERSION_GREATER "6.0") # there is no >= and also there is no 5.0.0.0
  # the current boost test case is crap :(
  set(CMAKE_CXX_FLAGS "${CFS_SUPPRESSIONS} -Wno-terminate")
endif() 

SET_TARGET_PROPERTIES(cfstest PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")
