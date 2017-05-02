# this script can be used to create an archive of the binary directory suitable for distribution of CFS.
#
# Run this script by: cmake -P create_binary_archive.cmake

if(NOT CTEST_BINARY_DIRECTORY)
  message("CTEST_BINARY_DIRECTORY not defined: set to current dir '${CMAKE_CURRENT_SOURCE_DIR}'")
  set(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif(NOT CTEST_BINARY_DIRECTORY)

if(NOT DIST_FAMILY)
  # determine distro
  exec_program("${CTEST_BINARY_DIRECTORY}/share/scripts/distro.sh -c" OUTPUT_VARIABLE DISTRO_OUT)
  file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeFiles/distro_out.cmake" "${DISTRO_OUT}")
  include("${CTEST_BINARY_DIRECTORY}/CMakeFiles/distro_out.cmake")
  if(DIST_FAMILY STREQUAL "")# its probably ubuntu ...
    set(DIST_FAMILY ${DIST})
  endif(DIST_FAMILY STREQUAL "")
  if(MAJOR_REV STREQUAL "")# its probably ubuntu ...
    set(MAJOR_REV ${REV})
  endif(MAJOR_REV STREQUAL "")
endif(NOT DIST_FAMILY)

if(HOSTNAME AND TEST_NAME)
  # the script is probably included by a deploy from a nightly test
  set(ARCHIVE_NAME "${HOSTNAME}_${TEST_NAME}")
else()
  set(ARCHIVE_NAME "CFS_${DIST_FAMILY}_${MAJOR_REV}_${ARCH}")
endif(HOSTNAME AND TEST_NAME)

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
  message("  Copy to subfolder '${CTEST_BINARY_DIRECTORY}/${ARCHIVE_NAME}' ...")
  file(COPY ${FILES_TO_KEEP} DESTINATION "${CTEST_BINARY_DIRECTORY}/${ARCHIVE_NAME}")
  # find static archives and remove them
  file(GLOB_RECURSE STATIC_ARCHIVES "${CTEST_BINARY_DIRECTORY}/${ARCHIVE_NAME}/lib64/*.a")
  if(STATIC_ARCHIVES)
    file(REMOVE ${STATIC_ARCHIVES})
  endif(STATIC_ARCHIVES)

  # pack the binaries
  message("  Packing binaries to '${ARCHIVE_NAME}.tgz' ... ")
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND} -E tar czf "${ARCHIVE_NAME}.tgz" ${ARCHIVE_NAME}
    WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}"
  )

  # remove created subfolder
  message("  Removing '${ARCHIVE_NAME}' ...")
  file(REMOVE_RECURSE "${CTEST_BINARY_DIRECTORY}/${ARCHIVE_NAME}")

else()

  message("  CFS binary not found. Was the build successful?")

endif()
