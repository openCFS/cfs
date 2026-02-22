# This file sets common CTest settings for the Testsuite (used standalone or from within CFS)
# for some reason sometimes the CTEST_* version is used and sometimes ONLY the version WITHOUT CTEST works

set(CTEST_PROJECT_NAME "CFS")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "am-ko.mi.uni-erlangen.de")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=CFS")
set(SUBMIT_URL "${CTEST_DROP_METHOD}://${CTEST_DROP_SITE}${CTEST_DROP_LOCATION}") # the Testsuite needs this version
set(CTEST_DROP_SITE_CDASH TRUE)

# cmake ctest_update only report version information, not change commit
set(CTEST_UPDATE_VERSION_ONLY ON)

# set update type explicitly
set(UPDATE_TYPE "git")

# probably not necessary
find_package(Git)
set(CTEST_GIT_COMMAND ${GIT_EXECUTABLE})

# set a build name
if(DEFINED ENV{CI_JOB_NAME}) # we are in a pipeline
  string(REPLACE ":" ";" BUILDNAME "$ENV{CI_JOB_NAME}") # transform CI_JOB_NAME into a cmake list
  # expand abbreviations used in the ci job names
  string(REPLACE "_ser" " serial" BUILDNAME "${BUILDNAME}")
  string(REPLACE "_par" " parallel" BUILDNAME "${BUILDNAME}")
  if("$ENV{CI_JOB_STAGE}" MATCHES "test(_extra)?") # in a test-stage the name starts with test-type for better sorting in the gitlab UI
    list(REVERSE BUILDNAME) # reverse to be consistant in CDash
  endif()
  list(JOIN BUILDNAME " / " BUILDNAME) # join again with colon
  set(BUILDNAME "$ENV{CI_COMMIT_REF_SLUG} // ${BUILDNAME}") # put branch name in front
else() # default case -> set chache varaible that the user could change
  set(BUILDNAME "${CFS_GIT_BRANCH} ${CFS_GIT_COMMIT} on ${DIST} ${REV} with ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}" CACHE STRING "Name of build on the dashboard")
  mark_as_advanced(BUILDNAME)
endif()
message(STATUS "set BUILDNAME ${BUILDNAME}")