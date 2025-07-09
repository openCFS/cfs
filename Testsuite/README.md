Testsuite for openCFS
=====================

This directory contains the testsuite for CFS.
It can be used as a
* _part of CFS_ to verify changes in the main CFS code base during development, and
* _standalone project_ to test certain CFS builds.

Additionally it provides examples illustrating the input of CFS.

Regular Use
-----------

The use of the *Testsuite* in CFS is governed by the configure option `CFS_TESTSUITE`.
When switched on, some additional config options `TESTSUITE_*` will become available.

Standalone Mode
---------------

Just configure with setting `CFS_BINARY_DIR` to the path of an CFS installation and using the `Testsuite` directory as the source, e.g.
```shell
cmake -DCFS_BINARY_DIR=/path/to/cfs/build /path/to/CFS/Testsuite
```
to use the CFS build in `/path/to/cfs/build` using the Testsuite source in `/path/to/CFS/Testsuite`.
This will execute `/path/to/cfs/build/bin/cfs -c config.cmake` to determine the used config options, which will activate the correct tests.

Now you can run
```shell
ctest
```
to run all activated tests.
