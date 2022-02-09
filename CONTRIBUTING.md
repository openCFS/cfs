Contributing to openCFS
=======================

Looking to __contribute code__ to openCFS? Awesome! Go for it (after reading this guide) ...
You can also contribute by _suggesting features_ or _giving feedback_ which is best done by opening an [issue](https://gitlab.com/openCFS/cfs/-/issues/new).

Getting Started
---------------

The codebase of openCFS is large, so it can be hard to find your way around at first.
To help you getting started, we provide the following resources
* short build instruction in the projects top level [README](README.md),
* information about [build dependecies for various platforms](share/doc/developer/build-dependencies),
* [getting started guide](https://gitlab.com/openCFS/cfs/-/wikis/getting-started) in the wiki,
* introduction to the [git workflow](https://gitlab.com/openCFS/cfs/-/wikis/git-workflow) in the wiki.

Merge Request Guidelines
------------------------

* __Create topic branches__.
  Choose a descriptive branch name for your change.
  Once you're done open a _merge request_ from _your branch/fork_ to the _master_ branch.
  Do not mix unrelated changes in a single branch.
* __Code against the master branch__.
  Always base your work on top of the master branch.
* __Single topic per change__. 
  Often one finds something that needs fixing but is unrelated to the current work (can happen in the best codebase).
  Please fix it in a separate commit, or (for larger issues) create another topic branch and associated _merge request_.
* __Use a coherent commit history__.
  Make sure your merge request is organized in logical chunks of meaningful individual commits.
  Tidy up your history before submitting a merge request.
* __Keep commits small__. Indeally, the changes bundled into a commit can be described by a one-liner, e.g.
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
* Respect the  __Assignee__ of the merge request on gitlab.
  To prevent conflicts when cleaning up the history, only a single person should be __assigned__ to a merge request:
  This person is allwed to __work__ (commit, force push, ...) on the branch of the mergre request.
  Change the assignment from yourself to someone else (which will get notified) to hand over!
  
Bringing your contribution in shape for a merge adhering to the guidelines above can be hard.
If you're unsure of how to do something please study our guide on [git workflow](https://gitlab.com/openCFS/cfs/-/wikis/git-workflow),
and do not hesitate to ask for further guidance.

Coding Guidelines
-----------------
* __Follow the Style Guide__ by adapting to the existing code.
  The rules are actually quite relaxed.
  In any case: no tabs, use 2 spaces as indentation).
  For details see https://gitlab.com/openCFS/cfs/-/wikis/StyleGuide.
* __Separate code from cosmetics__. 
  Keep code and cosmetics in separate commits, e.g. do not combine _white spaaaaaace change_ like adjusting line breaks or indentation with _functional change_.
  This makes code review harder!
  Of course corrections in _areas related to the change_ are welcome.
  Please be mindful of the reviewer going through the diffs and keep white spaaaaaace changes to areas where you also do functional changes.
  If you want to correct indentation, style, etc. we are happy to receive your equally important contribution in a separate topic.
* __Document the code__.
  There is hardly too much documentation in source code! 
  Usually one is very happy about it years later, and it tremendously helps other to understand your intention.
  Your're also welcome to contribute to the existing `README.md` files in the source tree, like [source/README.md](source/README.md).

Workflow for Contributing
-------------------------

Please use the following workflow to contribute

1. Discuss your planned contribution: open an **issue**.
   This step is optional for trivial changes, but very useful for larger changes.
   Describe the planned changes and ask for feedback to make sure
   * you are not re-implementing already existing things, and
   * nobody else is working on a toipc with possible confilcts.
2. Create a topic-**branch** and start coding.
  As soon as possible, **commit** and **push** your work and
3. open a **merge request**.
   * Mark the merge request as *draft* or *work in progress*.
   * Add a description to the merge request (e.g. to-do list, discussion outcome of the issue, ...)
   * Comments, changes, and pipeline results will be collected within the merge request page on gitlab. 
4. Once your're done (contribution adhers to guidelines, pipeline passes), start the review: **assign** the merge request to a _maintainer_.
   * The maintainer might assign reviewers to comment on the contribution.
   * If there is still an issue the reviewer/maintainer will tell you about it and assign the issue back to you.
   * Once you're done resolving everything, change the assignment of the MR again.
   * Please do not open a new merge request!
