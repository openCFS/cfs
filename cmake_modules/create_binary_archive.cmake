# this script can be used to create an archive of the binary directory suitable for distribution of CFS.
#
# Run this script by: cmake -P create_binary_archive.cmake

if(CTEST_BINARY_DIRECTORY) # we are probably in a nighlty test
  set(ARCHIVE_NAME "${HOSTNAME}_${TEST_NAME}")
  set(ARCHIVE_PATH "${NIGHTLY_ARCHIVES_DIR}/${HOSTNAME}_${TEST_NAME}.tgz")
  message("CTEST_BINARY_DIRECTORY=${CTEST_BINARY_DIRECTORY}")
else(CTEST_BINARY_DIRECTORY)
  message("CTEST_BINARY_DIRECTORY not defined: set to current dir '${CMAKE_CURRENT_SOURCE_DIR}'")
  set(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  # determine distro
  exec_program("${CTEST_BINARY_DIRECTORY}/share/scripts/distro.sh -c" OUTPUT_VARIABLE DISTRO_OUT)
  file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeFiles/distro_out.cmake" "${DISTRO_OUT}")
  include("${CTEST_BINARY_DIRECTORY}/CMakeFiles/distro_out.cmake")
  set(ARCHIVE_NAME "CFS_${DIST_FAMILY}_${MAJOR_REV}_${ARCH}")
  set(ARCHIVE_PATH "${CTEST_BINARY_DIRECTORY}/${ARCHIVE_NAME}.tgz")
endif(CTEST_BINARY_DIRECTORY)

set(CFS_BIN_DIR "${CTEST_BINARY_DIRECTORY}/bin/${DIST_FAMILY}_${MAJOR_REV}_${ARCH}")
message("  Searching for cfs executables in ${CFS_BIN_DIR}")
set(CFS_BINARY "${CFS_BIN_DIR}/cfsbin")
if(EXISTS "${CFS_BINARY}")
  
  # define which files to keep
  set(FILES_TO_KEEP 
    "${CTEST_BINARY_DIRECTORY}/bin"
    "${CTEST_BINARY_DIRECTORY}/lib64"
    "${CTEST_BINARY_DIRECTORY}/share"
  )
  # copy to subfolder
  message("  Copy to subfolder '${ARCHIVE_NAME}' ...")
  file(COPY ${FILES_TO_KEEP} DESTINATION ${ARCHIVE_NAME})
  # find static archives and remove them
  file(GLOB_RECURSE STATIC_ARCHIVES "${CTEST_BINARY_DIRECTORY}/${ARCHIVE_NAME}/lib64/*.a")
  file(REMOVE ${STATIC_ARCHIVES})  

  # pack the binaries
  message("  Packing binaries to '${ARCHIVE_PATH}' ... ")
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND} -E tar czf ${ARCHIVE_PATH} ${ARCHIVE_NAME}
    WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}"
  )

  # remove created subfolder
  message("  Removing '${ARCHIVE_NAME}' ...")
  file(REMOVE_RECURSE "${CTEST_BINARY_DIRECTORY}/${ARCHIVE_NAME}")

else()

  message("  CFS binary not found. Was the build successful?")

endif()
