# Determine informations in case of a Git working copy. Used e.g. by cpack -G ZIP for the package name and structure
set(CFS_GIT_COMMIT "0000000") # might end with -modified in case
set(CFS_GIT_BRANCH "not_available")

if(EXISTS "${CFS_SOURCE_DIR}/.git")
  find_package("Git" 1.6.3 REQUIRED)
  
  # determine long git SHA1, e.g. "c06b95f58bd62caecbcb21f6ebc45cee6a5a0847" 
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD OUTPUT_VARIABLE CFS_GIT_COMMIT)
  string(STRIP "${CFS_GIT_COMMIT}" CFS_GIT_COMMIT)
  
  # Check if the working copy has been modified - empty when there is no uncommited change
  execute_process(COMMAND ${GIT_EXECUTABLE} diff OUTPUT_VARIABLE _MODIFIED)
  # change the commit name if modified, see cfs --version
  if(_MODIFIED)
    set(CFS_GIT_COMMIT "${CFS_GIT_COMMIT}-modified")
  endif()

  # Determine name of local git branch
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD OUTPUT_VARIABLE CFS_GIT_BRANCH)
  string(STRIP "${CFS_GIT_BRANCH}" CFS_GIT_BRANCH)
endif()
