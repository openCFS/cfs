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
#  cubicinttests.cc
)

set(TARGET_LL
  matvec
  design
  cfsgeneral
  utils
  datainout)

if(USE_EMBEDDED_PYTHON)
  list(APPEND TARGET_LL ${Python_LIBRARIES})
  set(UNIT_TESTS_SRC ${UNIT_TESTS_SRC} embeddedpython.cc)
endif()

# Main openCFS Dat executable target.
add_executable(cfstest ${UNIT_TESTS_SRC})


IF(MKL_BLAS_LIB)
  LIST(APPEND TARGET_LL ${MKL_BLAS_LIB})
ENDIF(MKL_BLAS_LIB)


if(USE_BLAS_LAPACK STREQUAL "NETLIB")
  assert_set(LAPACK_LIBRARY)
  list(APPEND TARGET_LL ${LAPACK_LIBRARY};${BLAS_LIBRARY})
endif()

TARGET_LINK_LIBRARIES(cfstest ${TARGET_LL})

# Linux-specific link flags:
# - librt is needed for glibc <= 2.17
# - --start-group fixes circular dependencies between static libraries
#   (e.g., materials->utils for BiLinInterpolate/TriLinInterpolate).
#   The linker auto-adds --end-group at the end if missing.
if(NOT APPLE AND NOT WIN32)
  set(CFS_LINK_FLAGS "-Wl,--start-group -lrt")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0) # there is no >= and also there is no 5.0.0.0
  # the current boost test case is crap :(
  set(CMAKE_CXX_FLAGS "${CFS_SUPPRESSIONS} -Wno-terminate")
endif()

# Even if we did not compile openCFS itself with OpenMP enabled, we may have
# compiled external libraries like LIS with OpenMP support.
if(OPENMP_FOUND)
  set(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

SET_TARGET_PROPERTIES(cfstest PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")
