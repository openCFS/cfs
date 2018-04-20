CFS++
==============================

[![pipeline status](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/badges/master/pipeline.svg)](https://cfs-dev.mdmt.tuwien.ac.at/cfs/CFS/commits/master)

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