openCFS
==============================

Multiphysics and more based on the finite element method:
openCFS (Coupled Field Simulation) is a finite element-based multi-physics modelling and simulation tool.
With about 20 years of research-driven development, the core of openCFS is used in scientific research and industrial application.
The modelling strategy focuses on physical fields and their respective couplings.

openCFS is under **[MIT License](LICENSE)**. You may use, distribute and copy it under the license terms.
Any compiled openCFS binaries will include a number of third party libraries with different licenses.
See the output of e.g. `cfs --version` or the `license` directory for details.

Building
--------

Building openCFS is supported under macOS, Windows, and Linux.

For building you need
* [CMake](https://cmake.org),
* C++ and Fortran compilers, and
* some additional dependencies depending on the build configuration.

There are installation instructions for required dependencies for
[Ubuntu](share/doc/developer/build-dependencies/ubuntu.md), 
[Debian](share/doc/developer/build-dependencies/debian.md),  
[Fedora](share/doc/developer/build-dependencies/fedora.md) and
[Red Hat Enterprise Linux (RHEL) and derivatives](share/doc/developer/build-dependencies/rhel.md).
They are regularly tested on the newest releases using [docker containers](share/docker/README.md), via our CI/CD [pipeline](.gitlab-ci.yml). 

We use CMake so it might be as simple as
```shell
mkdir build && cd build
cmake ..
make
```

To get started be sure to visit the [wiki](https://gitlab.com/openCFS/cfs/-/wikis/home).

Running
------------

After a successful build, you can execute
```shell
./bin/cfs
```
For details visit the [user documentation](https://opencfs.gitlab.io/userdocu/) with plenty of examples.

Contributing
------------

Interested in contributing to the project?
We're happy about
* constructive feedback
* code contributions
* scientific collaboration

Please consult our [contributing guide](CONTRIBUTING.md) for details.

Want to suggest an awesome feature? 
Did you _really_ find a bug? 
Please create an [issue](https://gitlab.com/openCFS/cfs/-/issues)!
You can also reach the issue tracker by [email](mailto:contact-project+opencfs-cfs-support@incoming.gitlab.com).  


Resources
---------

* [project homepage](https://www.opencfs.org)
* [user documentation](https://opencfs.gitlab.io/userdocu)
* [related projects](https://gitlab.com/openCFS)
* [developer wiki](https://gitlab.com/openCFS/cfs/-/wikis/home)
* [issue tracker](https://gitlab.com/openCFS/cfs/-/issues)
* [Doxygen](https://opencfs.gitlab.io/cfs/doxygen/) documentation
* [Developer documentation](share/doc/developer/README.md) in the source tree
