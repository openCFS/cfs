Contributing to openCFS
=======================

Looking to __contribute code__ to openCFS? Awesome! Go for it (after reading this guide) ...
You can also contribute by _suggesting features_ or _giving feedback_ which is best done by opening an [issue](https://gitlab.com/openCFS/cfs/-/issues).

Getting Started
---------------

The codebase of openCFS is large, so it can be hard to find your way around at first.
To help you to get started, we provide the following resources
* short build instruction in the projects top level [README](README.md),
* information about [build dependencies for various platforms](share/doc/developer/build-dependencies),
* [getting started guide](https://github.com/openCFS/cfs/wiki/getting-started) in the wiki,
* introduction to the [git workflow](https://github.com/openCFS/cfs/wiki/git-workflow) in the wiki.

Merge Request Guidelines
------------------------

* __Create topic branches__.
  Choose a descriptive branch name for your change.
  Once you're done open a _merge request_ from _your branch/fork_ to the _master_ branch.
  Do not mix unrelated changes in a single branch.
  Branch names of the _CFS_ and _Testsuite_ projects must match (this will be enforced strictly in the CI pipelines).
* __Code against the master branch__.
  Always base your work on top of the master branch.
* __Single topic per change__. 
  Often one finds something that needs fixing but is unrelated to the current work (can happen in the best codebase).
  Please fix it in a separate commit, or (for larger issues) create another topic branch and associated _merge request_.
* __Use a coherent commit history__.
  Make sure your merge request is organized in logical chunks of meaningful individual commits.
  Tidy up your history before submitting a merge request.
* __Keep commits small__. 
  Ideally, the changes bundled into a commit can be described by a one-liner, e.g.
  *added postprocessing result for magnetic energy density in MagEdgePDE*.
  The commit message adds context to the code changes (why this change?) and commits become instructive for new developers.
* __Discuss large changes__ in an issue to make sure nobody is already working on a similar topic. 
  This will also make sure you can get relevant input of senior developers and avoids possible conflicts with related changes happening at the same time.
* Run the __Testsuite__ regularly on your local development machine,
  and make sure every fast-forward commit in your merge request passes the __Pipeline__ with all relevant build configurations.
  We only include contributions that pass all tests in the pipeline. 
* Make sure your contribution __merges fast forward__ into the master branch.
  This is the only way to ensure the code is tested before ending up in the master branch, keeping it stable.
  If, for _good reason_ not every commit in your merge request passes the pipeline, you can create a (fast-forward) merge commit on your feature branch (which must pass all tests).
* Respect the __Assignee__ of the merge request on gitlab.
  To prevent conflicts when cleaning up the history, only a single person should be __assigned__ to a merge request:
  This person is allowed to __work__ (commit, force push, ...) on the branch of the merge request.
  Change the assignment from yourself to someone else (which will get notified) to hand over!
  
Bringing your contribution in shape for a merge adhering to the guidelines above can be hard.
If you're unsure of how to do something please study our guide on [git workflow](https://github.com/openCFS/cfs/wiki/git-workflow),
and do not hesitate to ask for further guidance.

Coding Guidelines
-----------------
* __Follow the Style Guide__ by adapting to the existing code.
  The rules are actually quite relaxed.
  In any case: no tabs, use 2 spaces as indentation.
  For details see https://github.com/openCFS/cfs/wiki/StyleGuide.
* __Separate code from cosmetics__. 
  Keep code and cosmetics in separate commits, e.g. do not combine _white space change_ like adjusting line breaks or indentation with _functional change_.
  This makes code review harder!
  Of course corrections in _areas related to the change_ are welcome.
  Please be mindful of the reviewer going through the diffs and keep white space changes to areas where you also do functional changes.
  If you want to correct indentation, style, etc. we are happy to receive your equally important contribution in a separate topic.
* __Document the code__.
  There is hardly too much documentation in source code! 
  Usually one is very happy about it years later, and it tremendously helps others to understand your intention.
  You're also welcome to contribute to the existing `README.md` files in the source tree, like [source/README.md](source/README.md).

Workflow for Contributing
-------------------------

Please use the following workflow to contribute

1. Discuss your planned contribution: open an **issue**.
   This step is optional for trivial changes, but very useful for larger changes.
   Describe the planned changes and ask for feedback to make sure
   * you are not re-implementing already existing things, and
   * nobody else is working on a topic with possible conflicts.
2. Create a topic-**branch** and start coding.
  As soon as possible, **commit** and **push** your work and
3. open a **merge request**.
   * Mark the merge request as *draft* or *work in progress*. This disables automatically running CI pipelines.
   * Add a description to the merge request (e.g. to-do list, discussion outcome of the issue, ...).
   * Comments, changes, and pipeline results will be collected within the merge request page on gitlab.
   * If you changed the _Testsuite_, open a merge request (with matching branch names) in the _Testsuite_ project too.
4. Once you're done (contribution adheres to guidelines, pipeline passes), start the review: **assign** the merge request to a _maintainer_.
   * The maintainer might assign reviewers to comment on the contribution.
   * If there is still an issue the reviewer/maintainer will tell you about it and assign the issue back to you.
   * Once you're done resolving everything, change the assignment of the MR again.
   * To allow maintainers to directly fix things, you can [allow commits from upstream members](https://docs.gitlab.com/ee/user/project/merge_requests/allow_collaboration.html#allow-commits-from-upstream-members) in your fork.
     By selecting the respective checkbox when creating/editing the MR you can give openCFS maintainers write permission on the MR-source-branch in your fork.
   * Please do not open a new merge request!

> If you contribute from a _fork_ please do not rename the project, i.e. keep `your-namespace/cfs` and `your-namespace/Testsuite`, because renaming will break the CI pipelines.

Notes for Maintainers
---------------------

* Adapt the `CFS_NAME` and `PrintHistory`, bump up version number.
* Generate a __Release__ by pushing an annotated tag `git tag -a`. 
  The annotation message should summarize the changes since the last release.
  Usually, we do two releases per year, tagged as `yyyyS` (summer) and `yyyyW` (winter).
  Do not forget to replace `:release` in the URLs below with the tag, e.g. `2022W`.
* A binary __Package__ can be distributed by using [GitLab's Generic Packages Repository](https://docs.gitlab.com/ee/user/packages/generic_packages/).
  This happens automatically via the Piepline but can be done manually too.
  Push the file using the API: `curl --header "PRIVATE-TOKEN: <personal_access_token>" --upload-file CFS-2022W-Linux.tar.gz "https://gitlab.com/api/v4/projects/12930334/packages/generic/openCFS/:release/CFS-2022W-Linux.tar.gz"`.
  Use a _personal access token_ with _api_ privileges created via the gitlab web-ui.
* Releases and packages are connected by __links__ in the release.
  Links are set up automatically in the pipeline, but can also be added manually.
  A _permalink_ to the latest release is `https://gitlab.com/openCFS/cfs/-/releases/permalink/latest/downloads/<FILEPATH>`.
  In order to make it work, one has to define a `filepath` for the link.
  This is only possible via the API: use `curl --request POST --header "PRIVATE-TOKEN: <personal_access_token>" --data name="<NAME>" --data url="<URL>" --data filepath="/<FILEPATH>" --data link_type="package" "https://gitlab.com/api/v4/projects/12930334/releases/:release/assets/links"`
   and adapt the following:
   - `<NAME>` is the name visible in the web-ui, e.g. `Linux binary (tag.gz)` or `Windows binary (zip)`
   - `<URL>` is the url of the package defined above
   - `<FILEPATH>` is the filename used for _latest_ links, e.g. `CFS-Linux.tar.gz` or `CFS-Windows.zip`
   - `:release` is the release you want to place the link in
