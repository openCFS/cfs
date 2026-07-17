## Description

(Summarise the changes introduced in this merge request)

## Check before assigning to a maintainer for review

* [ ] The history is clean (squashed into meaningful commits and commit messages)
* [ ] New Code is well-documented, understandable and according to our standard (see our [contributing guide](https://gitlab.com/openCFS/cfs/-/blob/master/CONTRIBUTING.md))
* [ ] In case of new functionality, there is a fast test case (also mentioned in the description, <1 second, <1 MB, see [here](https://gitlab.com/openCFS/cfs/-/wikis/Testsuite-Ctest-CDash-and-GitLab) for more information)
* [ ] Every commit passes the pipeline
* [ ] There are no new build warnings on cdash/pipeline on any build (see pipeline jobs)
* [ ] In case of new functionality link userdocu MR

## Maintainer checks before merge

* [ ] Above mentioned points are completed
* [ ] Check for performance regressions in the performance testsuite and determine if they are acceptable
* [ ] Review is approved, and all comments are resolved
* [ ] New tests are actually running (e.g. check if they appear on CDash)


## Additional information

For additional information (labels, when does a pipeline run, etc.) look [here](https://gitlab.com/openCFS/cfs/-/wikis/Pipeline)