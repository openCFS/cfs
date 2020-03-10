CFS++ build dependencies for Centos
===================================

This sould also work (with slight adaptions) on derivates like 
* [Fedora](./fedora.md)
* Scientific Linux

For older CentOS versions see
* [CentOS 7](./centos:7.md)

To configure we need CMake
```shell
yum install -y cmake
```

For the typical build config we need
```shell
yum install -y make gcc gcc-c++ gcc-gfortran patch m4
```

Intel MKL can be installed using [Intel's YUM repositories](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-yum-repo).

First add the repo (with required dnf-plugin) and install the GPG key
```shell
dnf install -y dnf-plugins-core
dnf config-manager --add-repo https://yum.repos.intel.com/mkl/setup/intel-mkl.repo
rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
```

Run `yum search` once to trigger GPG import.
Since this returns with a non-zero exit code we add and OR for the testing pipeline.
```shell
yum search -y intel || echo "GPG imported ..."
```

Then install the 64-bit version of the most recent MKL release (first `yum search` necessary to import GPG key)
```shell
yum search -y intel-mkl-64bit
LATEST_MKL=$(yum search -y intel-mkl-64bit | grep intel-mkl-64bit | tail -n 1 | awk '{print $1}')
yum install -y $LATEST_MKL
```

Additionally every developer should have git installed
```shell
yum install -y git
```
