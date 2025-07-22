# Locate MPI for PETSc wrapper
include(FindMPI)
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_DIRS})

SET(PETSC_SRC PETSCSolver.cc PETScCommon.cc)

SET_SOURCE_FILES_PROPERTIES(${PETSC_SRC}
  PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing")

ADD_LIBRARY(petsc-olas STATIC ${PETSC_SRC})

# Ensure PETSc wrapper target uses MPI includes
target_include_directories(petsc-olas PUBLIC ${MPI_CXX_INCLUDE_DIRS})

SET(TARGET_LL
  ${PETSC_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS})

TARGET_LINK_LIBRARIES(petsc-olas ${TARGET_LL} ${MPI_C_LIBRARIES})
