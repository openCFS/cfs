# Determine informations in case of a Git working copy. Used e.g. by cpack -G ZIP for the package name and structure
set(CFS_GIT_COMMIT "0000000") # might end with -modified in case
set(CFS_GIT_BRANCH "not_available")

if(EXISTS "${CFS_SOURCE_DIR}/.git")
  find_package("Git" 1.4 REQUIRED)

  # e.g. "commit c06b95f58bd62caecbcb21f6ebc45cee6a5a0847 (HEAD -> fwein_small_fixes)"
  execute_process(COMMAND ${GIT_EXECUTABLE} log -1 OUTPUT_VARIABLE _GIT_LOG)

  # extract the git SHA1 -> c06b95f58bd62caecbcb21f6ebc45cee6a5a0847
  string(REGEX REPLACE "^(.*\n)?commit ([^\n]+).*" "\\2" CFS_GIT_COMMIT "${_GIT_LOG}")

  # shorten the commit SHA to something short and still unique -> c06b95f58
  # execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short ${_CFS_GIT_COMMIT} OUTPUT_VARIABLE CFS_GIT_COMMIT)
  # string(STRIP "${CFS_GIT_COMMIT}" CFS_GIT_COMMIT)

  # Check if the working copy has been modified - empty when there is no uncommited change
  execute_process(COMMAND ${GIT_EXECUTABLE} diff OUTPUT_VARIABLE _MODIFIED)

  # change the commit name if modified, see cfs --version
  if(_MODIFIED)
    set(CFS_GIT_COMMIT "${CFS_GIT_COMMIT}-modified")
  endif()

  # Determine name of local git branch - note that every commit has a different commit hash
  execute_process(COMMAND ${GIT_EXECUTABLE} branch OUTPUT_VARIABLE _BRANCH_LIST)

  # identify CFS_GIT_BRANCH by the trailing asterix
  string(REGEX REPLACE "^(.*\n)?\\* ([^\n]+).*" "\\2" CFS_GIT_BRANCH "${_BRANCH_LIST}")
endif()
