CFS++
==============================

[![pipeline status](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/pipeline.svg)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master) [![Testsuite](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/coverage.svg?job=test_stable:gcc4:linux)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master)

* `CodeScore` [![build gcc4](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/coverage.svg?job=build:gcc4:linux)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master) GCC4 build +++ [![build gcc6](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/coverage.svg?job=build:gcc6:linux)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master) GCC6 build +++ [![build gcc7](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/coverage.svg?job=build:gcc7:linux)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master) GCC7 build
* `Testsuite` [![test_all gcc4](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/coverage.svg?job=test_all:gcc4:linux)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master) GCC4 build +++ [![test_all gcc6](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/coverage.svg?job=test_all:gcc6:linux)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master) GCC6 build +++ [![test_all gcc7](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/coverage.svg?job=test_all:gcc7:linux)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master) GCC7 build




This is the git repository of CFS++.

Building
--------

For building you need
* cmake
* C++ and Fortran compilers

We use CMake so it might be as simple as
```
mkdir build && cd build
cmake ..
ccmake . # to review the options
make
```

To get started be sure to visit the [wiki](/../wikis).

Running
------------

After a sucessful build you can execute
```
./bin/cfs
```
For details visit the [user documentation](https://cfs-doc.mdmt.tuwien.ac.at) with plenty of examples.


[Developer Documentation](/share/doc/developer/README.md)
---------------------------------------------------------
Is done by markdown files at appropriate locations direclty in the source tree.
They should be all linked in [`share/doc/developer/README.md`](/share/doc/developer/README.md).

Bugs
----------------
Did you _really_ find a bug? Please create an [issue](/../issues), and then fix it ...