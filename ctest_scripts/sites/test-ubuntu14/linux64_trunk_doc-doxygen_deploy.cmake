# This script is can be included by the deploy scipts
# It packs the binaries, copies them to draco and extracts them there
# to test the file run it via cmake -DTEST_NAME=linux64_trunk_doc-doxygen -P filename
message("=============================================================================")
message(" START DEPLOY")
message("=============================================================================\n")
# include site variables & set directories
include("${CMAKE_CURRENT_LIST_DIR}/site_specific_vars.cmake")
set(CTEST_SOURCE_DIRECTORY "${NIGHTLY_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY "${NIGHTLY_BUILD_DIR}/${TEST_NAME}")
set(CTEST_COMMAND  "\"${CTEST_EXECUTABLE_NAME}\"")

message("${CMAKE_CURRENT_LIST_DIR},  ${CTEST_SOURCE_DIRECTORY}, ${CTEST_BINARY_DIRECTORY}")

# create an archive with all necessary files
set(DOC_PATH "${CTEST_BINARY_DIRECTORY}/share/doc/developer/doxygen/html")
set(ARCHIVE_PATH "${CTEST_BINARY_DIRECTORY}/share/doc/developer/doxygen/doxygen.tgz")
message("  Packing binaries to '${ARCHIVE_PATH}' from ${DOC_PATH}... ")
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND} -E tar czf "doxygen.tgz" "html"
    WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}/share/doc/developer/doxygen"
  )
# the created archive is picked up by the CFS server and published online
