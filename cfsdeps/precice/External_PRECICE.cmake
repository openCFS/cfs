# External_PRECICE.cmake - Integrate a system provided preCICE package
#
# preCICE has a rather involved dependency chain (Boost, Eigen, MPI, PETSc, …)
# and we currently rely on a user / package manager provided installation.  This
# external therefore just validates that such an installation is available and
# exposes the usual *_INCLUDE_DIR / *_LIBRARY helper variables so that the
# regular openCFS targets can link against the imported target.

clear_depencency_variables()

set(PACKAGE_NAME "precice")
set(PACKAGE_VER "3.1.2")

find_package(precice ${PACKAGE_VER} CONFIG REQUIRED)

if(NOT TARGET precice::precice)
  message(FATAL_ERROR
    "preCICE CMake package was found but the imported target precice::precice "
    "is missing. Please check your installation.")
endif()

# Provide compatibility variables that match other CFS dependency modules
set(PRECICE_INCLUDE_DIR "${precice_INCLUDE_DIRS}")
set(PRECICE_LIBRARY "${precice_LIBRARIES}")
set(PRECICE_LIBRARIES "${precice_LIBRARIES}")

# Provide a no-op build target so the cfsdeps aggregate target can depend on it
if(NOT TARGET precice)
  add_custom_target(precice
    COMMENT "System-installed preCICE detected - nothing to build in cfsdeps")
endif()

set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
