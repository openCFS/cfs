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
  Please fix it in a separate topic branch and open another _merge request_.
* __Use a coherent commit history__.
  Make sure each individual commit in your pull request is meaningful and organized in logical chunks.
  Tidy up and squash commits before submitting a merge request.
  Most contributions will fit in a single commit.
* __Discuss large changes__ in an issue to make sure nobody is already working on a similar topic. 
  This will also make sure you can get relevant input of senior developers and avoids possible conflicts with related changes happening at the same time.
* Run the __Testsuite__ regularly on your local development machine,
  and make sure your merge request passes the __Pipeline__ with all relevant build configurations.
  We only include contributions that pass all tests in the pipeline. 
* Make sure your contribution __merges fast forward__ into the master branch.
  This is the only way to ensure the code is tested before ending up in the master branch, keeping it stable.

Coding Guidelines
-----------------
* __Follow the Style Guide__ by adapting to the existing code.
  The rules are actually quite relaxed.
  In any case: no tabs, use 2 spaces as indentation).
  For details see https://gitlab.com/openCFS/cfs/-/wikis/StyleGuide.
* __Separate code from cosmetics__. 
  Keep code and cosmetics in separate commits, e.g. do not combine _white space change_ like adjusting line breaks or indentation with _functional change_.
  This makes code review harder!
  Of course corrections in _areas related to the change_ are welcome.
  Please be mindful of the reviewer going through the diffs and keep white space changes to areas where you also do functional changes.
  If you want to correct indentation, style, etc. we are happy to receive your equally important contribution in a separate topic.
* __Document the code__.
  There is hardly too much documentation in source code! 
  Usually one is very happy about it years later, and it tremendously helps other to understand your intention.
  Your're also welcome to contribute to the existing `README.md` files in the source tree, like [source/README.md](source/README.md).

Handling and Updating a Merge Request
-------------------------------------

Bringing your contribution in shape for a merge adhering to the guidelines above can be hard.
When your work is ready please _assign_ your merge request (MR) to a _maintainer_.
If there is still an issue the maintainer will tell you about it and assign the issue back to you.
Once you're done resolving everything, change the assignment of the MR again.
Please do not open a new merge request.
If you're unsure of how to do that please study our guide on [git workflow](https://gitlab.com/openCFS/cfs/-/wikis/git-workflow),
and do not hesitate to ask for further guidance.
