# FEAST Eigenvalue Solver
# http://www.feast-solver.org/
# Originally has only Makefiles, we provide a CMakeLists.txt for Windows and use it always
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "feast")
set(PACKAGE_VER "4.0") # note that the Version is hardcoded in CMakeLists.txt
set(FEAST_VER ${PACKAGE_VER}) # for --version
set(PACKAGE_FILE "feast_${PACKAGE_VER}.tgz")
set(PACKAGE_MD5 "e4e6b47de276c203de2c0e9e7d9e5a65")
set(PACKAGE_MIRRORS "https://gitlab.com/api/v4/projects/12930334/packages/generic/cfsdeps/sources/${PACKAGE_FILE}")  
set(DEPS_VER "-d") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# we cannot link a parallel compiled suitesparse with debug without openmp
# we need to set DEPS_ID before calling set_precompiled_pckg_file()
if(USE_OPENMP)
  set(DEPS_ID "OPENMP")
else()
  set(DEPS_ID "NO_OPENMP")
endif()
# MKL is a special code feature 
if(USE_BLAS_LAPACK STREQUAL "MKL")
  set(DEPS_ID "${DEPS_ID}-MKL")
else()  
  set(DEPS_ID "${DEPS_ID}-NO_MKL")
endif()

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a fortran compiler
use_c_and_fortran(OFF ON)

set_precompiled_pckg_file()

# libfeast.a
set_package_library_default()

set_standard_variables() 

# this is the standard target for cmake projects. The files to package come from the install_manifest.txt
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(OFF) # don't set compiler flags, we need to change

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU") 
  # this is the original call via Makefile, where -fallow-argument-mismatch is already set in FindCFSDEPS.cmake to CFSDEPS_Fortran_FLAGS
  # gfortran -O3 -fopenmp -ffree-line-length-none -ffixed-line-length-none -cpp  -DMKL -m64 -O3  -w -fallow-argument-mismatch -c .../dzfeast.f90 -o .../dzfeast.o
  set(ADDITIONAL_FORTRAN_ARGS "-ffree-line-length-none -ffixed-line-length-none -cpp ")
  if(USE_OPENMP)
    set(ADDITIONAL_FORTRAN_ARGS "-fopenmp ${ADDITIONAL_FORTRAN_ARGS}")
  endif() 
else()
  # this is the original call via Makefile
  # ifort -O3 -qopenmp -cpp -fPIC -DMKL -c .../dzfeast.f90 -o .../dzfeast.o
  set(ADDITIONAL_FORTRAN_ARGS "-cpp -fPIC ")
  if(USE_OPENMP)
    set(ADDITIONAL_FORTRAN_ARGS "-qopenmp ${ADDITIONAL_FORTRAN_ARGS}")
  endif()  
endif()
if(USE_BLAS_LAPACK STREQUAL "MKL")
  set(ADDITIONAL_FORTRAN_ARGS "${ADDITIONAL_FORTRAN_ARGS} -DMKL")
endif() 
list(APPEND DEPS_ARGS "-DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS} ${ADDITIONAL_FORTRAN_ARGS}") 

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# copies your CMakeLists.txt a
generate_patches_script()

assert_unset(POSTINSTALL_SCRIPT)

# generate package creation script. We get the files from an install_manifest.txt
generate_packing_script_manifest()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()
# if not, build newly and possibly pack the stuff
else()
  create_external_cmake_patched()  
  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  endif()
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
