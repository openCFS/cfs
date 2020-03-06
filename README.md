CFS++
==============================

This is the git repository of CFS++.

Building
--------

For building you need
* [CMake](https://cmake.org)
* C++ and Fortran compilers
* some additional dependecies depending on the build configuration

There are installation instructions for required dependencies for
[Ubuntu](share/doc/developer/build-dependencies/ubuntu.md), 
[Debian](share/doc/developer/build-dependencies/debian.md),  
[Centos](share/doc/developer/build-dependencies/centos.md) and
[Fedora](share/doc/developer/build-dependencies/fedora.md).
They are regularly tested on the newest releases using [docker containers](share/docker/README.md), via our CI/CD [pipeline](.gitlab-ci.yml). 

We use CMake so it might be as simple as
```shell
mkdir build && cd build
cmake ..
cmake .
make
```

To get started be sure to visit the [wiki](/../wikis).

Running
------------

After a sucessful build you can execute
```shell
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
