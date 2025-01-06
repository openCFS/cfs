CFS++ build dependencies for Centos 7
=====================================

> since Centos 7 is not supported any more the repos are offline.
> Thus, these instructions do not work any more.
> Perhaps they could be reactivated by using the arviced repos at https://vault.centos.org

This should also work (with slight adaptions) on derivates like 
* Fedora
* Centos
* Scientific Linux

For the default build config, we need some tools
```shell
yum install -y make patch m4 file which findutils diffutils
```
Furthermore, we need a c++14-ready compiler,
we enable the [software-collections repo](https://www.softwarecollections.org/en/docs/)
and install GCC7 compilers from devtoolset-7
```shell
yum install -y centos-release-scl
yum install -y devtoolset-7-gcc devtoolset-7-gcc-c++ devtoolset-7-gcc-gfortran
```
One needs to activate the use of devtoolset-7 before building by running 
`scl enable devtoolset-7 bash`
or to make it persistent when using bash
```shell
echo "source /opt/rh/devtoolset-7/enable" >> ~/.bashrc
```

To configure we need CMake 3 available via the EPEL repository
```shell
yum install -y epel-release
yum install -y cmake3 
```
Now we make cmake3 the default, i.e. running `cmake` will execute `cmake3` including the
```shell
alternatives --install /usr/local/bin/cmake cmake /usr/bin/cmake3 20 \
  --slave /usr/local/bin/ctest ctest /usr/bin/ctest3 \
  --slave /usr/local/bin/cpack cpack /usr/bin/cpack3 \
  --slave /usr/local/bin/ccmake ccmake /usr/bin/ccmake3 \
  --family cmake
```
To switch to the original `cmake` one can install an equivalent alternative and then switch by `alternatives --config cmake` [1](https://stackoverflow.com/questions/48831131/cmake-on-linux-centos-7-how-to-force-the-system-to-use-cmake3).

For the Intel MKL installer (`USE_MKL=ON`) we need `cpio`
```shell
yum install -y cpio
```

Intel MKL can be installed using [Intel's YUM repositories](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-yum-repo).

First, add the repo and install the GPG key
```shell
yum-config-manager --add-repo https://yum.repos.intel.com/mkl/setup/intel-mkl.repo
rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
```

Then install the 64-bit version of the most recent MKL release
```shell
LATEST_MKL=$(yum search -y intel-mkl-64bit | grep intel-mkl-64bit | tail -n 1 | awk '{print $1}')
yum install -y $LATEST_MKL
```

Additionally, every developer should have git installed
```shell
yum install -y git
```
