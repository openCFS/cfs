# LibTorch Dependency (PyTorch C++ API)
# https://pytorch.org

clear_depencency_variables()

if(NOT CFS_ARCH STREQUAL "ARM64")
  message(FATAL_ERROR "LibTorch is only supported on macOS for this project.")
endif()


set(PACKAGE_NAME "libtorch")
set(LIBTORCH_VER "2.6.0.dev20241018")
set(PACKAGE_VER ${LIBTORCH_VER})
set(PACKAGE_FILE "libtorch-macos-arm64-${LIBTORCH_VER}.zip")
set(PACKAGE_MD5 "0032fc91729ea256f724c113684582cf")  
# set(PACKAGE_FILE "libtorch-macos-arm64-${LIBTORCH_VER}.zip")
# set(PACKAGE_MD5 "c2d460b9493f41fef6b1251ef6d567b0")  

# Set download mirrors or use local package file
set(PACKAGE_MIRRORS "https://download.pytorch.org/libtorch/nightly/cpu/${PACKAGE_FILE}")
add_standard_mirrors_or_set_local()

# Sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

set(DEPS_INSTALL "${DEPS_PREFIX}/install")
# Determine paths of libraries and make them visible via ccmake
set_package_library_default()
# Set hidden cache variables *_LIBRARY, *_INCLUDE, and some defaults
set_standard_variables()

set(LIBTORCH_INCLUDE_DIR "${DEPS_PREFIX}/src/libtorch/include")
# set(LIBTORCH_LIBRARY "${DEPS_PREFIX}/libtorch/lib/libtorch.dylib;${DEPS_PREFIX}/libtorch/lib/libc10.dylib")

# If we have precompiled binaries, unpack them
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #create_external_unpack_precompiled()

# Otherwise, download and extract the precompiled package directly
else()
  ExternalProject_Add(
    ${PACKAGE_NAME}
    PREFIX "${DEPS_PREFIX}"
    URL "${PACKAGE_MIRRORS}"
    URL_MD5 "${PACKAGE_MD5}"
    DOWNLOAD_DIR "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}"
    DOWNLOAD_NAME "${PACKAGE_FILE}"
    CONFIGURE_COMMAND ""  # Skip configure
    BUILD_COMMAND ""  # Skip build
    INSTALL_COMMAND ""  # Precompiled, no install needed
  )
endif()

# Set C++17 for the libtorch module only
# set_target_properties(libtorch PROPERTIES CXX_STANDARD 17)
# set_target_properties(libtorch PROPERTIES CXX_STANDARD_REQUIRED ON)


# Add this project to the global list of dependencies
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})