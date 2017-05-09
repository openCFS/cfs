# This script sets varaibles appropriate to this site.
# It is should be included at the start of every test
# on this site.

# include site specific vars
include("${CMAKE_CURRENT_LIST_DIR}/site_specific_vars.cmake")

# load data on development server
include("${NIGHTLY_SOURCE_DIR}/cmake_modules/DevelopmentServer.cmake")

# This initializes the test and must be invoked in the actual test file
macro(INIT_TEST)
# set source and binary dirs
SET(CTEST_SOURCE_DIRECTORY "${NIGHTLY_SOURCE_DIR}")
get_filename_component(FNAME "${CMAKE_CURRENT_LIST_FILE}" NAME_WE)
set(CTEST_BINARY_DIRECTORY "${NIGHTLY_BUILD_DIR}/${FNAME}")

# include common macros
include("${CTEST_SOURCE_DIRECTORY}/ctest_scripts/shared/test_macros.cmake")

# Identify distro.
IDENTIFY_DISTRO()
# empty binary directory
#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY}) # this fails if there is no "CMakeLists.txt"
file(REMOVE_RECURSE "${CTEST_BINARY_DIRECTORY}")
file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")

# write ctest config
WRITE_CTEST_CONFIG()

# set ctest & cmake commands to the executables corresponding to the one used to run this script.
SET(CTEST_COMMAND  "\"${CTEST_EXECUTABLE_NAME}\"")
SET(CTEST_CMAKE_COMMAND  "\"${CMAKE_EXECUTABLE_NAME}\"")
endmacro(INIT_TEST)
