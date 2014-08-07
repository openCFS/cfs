CONFIGURE_FILE("/vagrant/readme_nightly_binaries.md"
               "/vagrant/README.md"
               @ONLY)

# Create readme HTML for nightly binaries
EXECUTE_PROCESS(
  COMMAND pandoc --highlight-style tango -s -S -H "/vagrant/readme_macosx.css" "/vagrant/README.md" -o README.html
  WORKING_DIRECTORY "${NIGHTLY_ARCHIVES_DIR}"
  RESULT_VARIABLE RETVAL
)
