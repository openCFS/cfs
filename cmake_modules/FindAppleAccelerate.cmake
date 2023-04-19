# Apple has the Accelerate Framework as BLAS/LAPACK implementation
#
# as it is system provided we don't bundle with openCFS

# control with environment variable VECLIB_MAXIMUM_THREADS
# has on M1 with single core double speed compared to openblas
# complex valued eigenvalue problems fail/are not supported

if(NOT APPLE)
  message(FATAL_ERROR "The Apple Accelerate Framework BLAS/LAPACK can only be search on macOS")
endif()

set(BLA_VENDOR Apple)
find_package(BLAS REQUIRED)

if(NOT BLAS_FOUND)
  message(FATAL_ERROR "CMake cannot find Apple Accelerate Framework")
endif()

# cmake uses BLAS_LIBRARIES and openCFS uses BLAS_LIBRARY, how stupid :(
set(BLAS_LIBRARY ${BLAS_LIBRARIES})
# BLAS_LINKER_FLAGS are empty 

