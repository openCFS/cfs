# This script is can be included by the deploy scipts
# It packs the binaries, copies them to draco and extracts them there
message("=============================================================================")
message(" START DEPLOY")
message("=============================================================================\n")
# include site variables & set directories
include("${CMAKE_CURRENT_LIST_DIR}/site_specific_vars.cmake")
set(CTEST_SOURCE_DIRECTORY "${NIGHTLY_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY "${NIGHTLY_BUILD_DIR}/${TEST_NAME}")
set(CTEST_COMMAND  "\"${CTEST_EXECUTABLE_NAME}\"")

# create an archive with all necessary binary files
message("packing binaries")
include("${CTEST_SOURCE_DIRECTORY}/cmake_modules/create_binary_archive.cmake")

# copy the archive to draco
message("copy to draco using scp")
set(NIGHTLY_PATH_DRACO "/opt/programs/cfs_nightly")
set(NIGHTLY_ARCHIVES_PATH_DRACO "${NIGHTLY_PATH_DRACO}/archives")
find_program(SCP_EXECUTABLE scp)
execute_process(COMMAND ${SCP_EXECUTABLE} "${CTEST_BINARY_DIRECTORY}/${ARCHIVE_NAME}.tgz" "draco:${NIGHTLY_ARCHIVES_PATH_DRACO}/." RESULT_VARIABLE RETVAL)
if(RETVAL)
  message("  FAIL")
endif(RETVAL)

# use ssh to extract them on draco
message("use ssh to extract on draco")
find_program(SSH_EXECUTABLE ssh)
string(REPLACE "_" ";" TEST_NAME_LIST "${TEST_NAME}")
message("${TEST_NAME}, ${TEST_NAME_LIST}")
list(GET TEST_NAME_LIST 1 NT_BRANCH)
list(GET TEST_NAME_LIST 2 NT_COMPILER)
set(NIGHTLY_BIN_PATH_DRACO "${NIGHTLY_PATH_DRACO}/${NT_BRANCH}_${NT_COMPILER}")
# set the ssh command
set(SSHCMD "${SSH_EXECUTABLE} draco 'tar -xzf ${NIGHTLY_ARCHIVES_PATH_DRACO}/${ARCHIVE_NAME}.tgz --strip 1 -C ${NIGHTLY_BIN_PATH_DRACO}'")
message("  ${SSHCMD}")
# write it to a temporary file
write_file("/tmp/sshcmd.sh" ${SSHCMD})
# execute it using sh
execute_process(COMMAND sh /tmp/sshcmd.sh RESULT_VARIABLE RETVAL)
if(RETVAL)
  message("  FAIL: ${RETVAL}, ${OUT}, ${ERR}")
endif(RETVAL)
