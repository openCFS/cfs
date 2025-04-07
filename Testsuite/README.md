Testsuite for openCFS
=====================

This repo contains the testsuite for CFS.
It is used as
* a *submodule* in the main CFS repository, and
* a standalone project to test certain CFS builds.

Additionally it provides examples illustrating the input of CFS.

Submodule
---------

The *Testsuite* is referenced in the main CFS repository, thereby tracking with revision of CFS is associated with revision of the testsuite.
Please name associated branches in testsuite and CFS __exactly the same__.

When used as a submodule in CFS, the testsuite is included into the configuration of CFS.


Standalone Mode
---------------

Just configure with setting `CFS_BINARY_DIR` to the path of an CFS installation, e.g.
```shell
cmake -DCFS_BINARY_DIR=build .
```
to use the CFS build in `build`.
This will execute `build/bin/cfs -c config.cmake` to determine the used config options, which will activate the correct tests.

Now you can run
```shell
ctest
```
to run all activated tests.
