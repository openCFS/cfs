
IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
  LINK_DIRECTORIES(${IPOPT_LINK_DIR})
ENDIF(USE_IPOPT)

# List of source codes for current target.
SET(CFS_SRCS CFS.cc)

# Main openCFS executable target.
ADD_EXECUTABLE(cfs ${CFS_SRCS})

add_dependencies(cfs cfsdeps)

# Define a list of libraries openCFS depends upon. At the moment there
# are some nasty circular dependencies which should be resolved in
# the future (see e.g. pde <-> domain, domain <-> gridcfs !!!)
  
SET(TARGET_LL
  cfsgeneral
  utils
  datainout
  paramh
  utils
  domain
  elemmapping
  ${CFS_FORTRAN_DYNRT_LIBS} )

IF(USE_XERCES)
  LIST(APPEND TARGET_LL ${XERCES_LIBRARY})
ENDIF(USE_XERCES)

IF(USE_LIBXML2)
  LIST(APPEND TARGET_LL ${LIBXML2_LIBRARY})
ENDIF(USE_LIBXML2)  

IF(USE_LIBFBI)
  LIST(APPEND TARGET_LL ${BOOST_THREAD_LIB} )
ENDIF()

IF(USE_PHIST_EV OR USE_PHIST_CG)
  LIST(APPEND TARGET_LL ${CUDA_LIBS} ${HWLOC_LIBRARY})
ENDIF()

if(USE_BLAS_LAPACK STREQUAL "NETLIB")
  assert_set(LAPACK_LIBRARY)
  list(APPEND TARGET_LL ${LAPACK_LIBRARY})
  list(APPEND TARGET_LL ${BLAS_LIBRARY})
elseif(USE_BLAS_LAPACK STREQUAL "MKL")
  list(APPEND TARGET_LL ${MKL_BLAS_LIB})
endif()

IF(DEPS_SEQUENTIAL)
  LIST(APPEND TARGET_LL ${DEPS_SEQUENTIAL})
ENDIF(DEPS_SEQUENTIAL)

# lapack is in NETLIB, OPENBLAS and MKL
# LIST(APPEND TARGET_LL -lpthread -lquadmath )

SET(CFS_LINK_FLAGS "")

# Even if we did not compile openCFS itself with OpenMP enabled, we may have
# compiled external libraries like LIS with OpenMP support.
IF(OPENMP_FOUND)
  SET(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} ${OpenMP_CXX_FLAGS}  ")
ENDIF()

IF(CMAKE_GENERATOR STREQUAL "Visual Studio 16 2019")
  SET(INTEL_LIB_DIR "${INTEL_COMPILER_DIR}/../../compiler/lib/intel64_win")
  TARGET_LINK_DIRECTORIES(cfs PUBLIC ${INTEL_LIB_DIR})
ENDIF()

IF(USE_PRECICE)
  find_package(precice CONFIG REQUIRED)

  # Retrieve the interface include directories from the imported target
  get_target_property(precice_includes precice::precice INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Using preCICE include directories: ${precice_includes}")

  # Append the 'precice' subdirectory to the existing include paths
  # I really don't know why precice config misses the /precice path...
  foreach(dir ${precice_includes})
    list(APPEND corrected_includes "${dir}/precice")
  endforeach()

  message(STATUS "Corrected preCICE include directories: ${corrected_includes}")

  # Add these include directories to the target
  target_include_directories(cfs PRIVATE ${corrected_includes})
  LIST(APPEND TARGET_LL precice::precice)
ENDIF()

# Target link libraries for openCFS executable target.
#message(STATUS "cfs TARGET_LL=${TARGET_LL}")
TARGET_LINK_LIBRARIES(cfs ${TARGET_LL})
if(MKL_LIBS AND CMAKE_VERSION VERSION_GREATER "3.24") # TODO: possibly remove this once legacy_mkl.cmake is removed
  message(STATUS "Setting link group")
  target_link_libraries(cfs "$<LINK_GROUP:RESCAN,${MKL_LIBS}>") # keept the group together
endif()
SET_TARGET_PROPERTIES(cfs PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")

# debug
#cmake_print_properties(TARGETS cfs PROPERTIES BUILD_RPATH INSTALL_RPATH INSTALL_REMOVE_ENVIRONMENT_RPATH)