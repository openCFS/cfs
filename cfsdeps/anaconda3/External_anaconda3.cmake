# anaconda by continuum analytics

#-------------------------------------------------------------------------------
# Set prefix path and path to sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(anaconda3_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/anaconda3")
set(anaconda3_install "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/anaconda3/install")

SET(MIRRORS
  "file:///home/ftoth/Downloads/Anaconda3-4.2.0-Linux-x86_64.sh" # local for testing
  "https://repo.continuum.io/archive/${ANACONDA3_SH}"
  "${ANACONDA3_URL}/${ANACONDA3_SH}"
  )
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/anaconda3/${ANACONDA3_SH}")
SET(MD5_SUM "${ANACONDA3_MD5}")

SET(DLFN "${anaconda3_prefix}/anaconda3-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

#-------------------------------------------------------------------------------
# The ANACONDA3 external project
#-------------------------------------------------------------------------------
find_program(BASH bash)
MARK_AS_ADVANCED(BASH)
ExternalProject_Add(anaconda3
    PREFIX ${anaconda3_prefix}
    SOURCE_DIR ${anaconda3_source}
    URL ${LOCAL_FILE}
    URL_MD5 "${ANACONDA3_MD5}"
    DOWNLOAD_NO_EXTRACT 1
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "${BASH}" "${LOCAL_FILE}" "-b" "-p" "${anaconda3_install}"
)

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
ExternalProject_Add_Step(anaconda3 cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
     DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${anaconda3_prefix}
)

ExternalProject_Add(anaconda3-vtk
    DEPENDS anaconda3
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    DOWNLOAD_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "${anaconda3_install}/bin/conda" "install" "--yes" "-c" "menpo" "vtk"
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  anaconda3
  anaconda3-vtk
)

set(PYTHON_LIBRARY "${anaconda3_install}/lib" CACHE PATH
	"Path to Python library")


