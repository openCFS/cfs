CFS++ build dependencies for Fedora Core
========================================

For the minimal build config we need
```shell
dnf install -y make gcc gcc-c++ gcc-gfortran patch m4 cmake
```

Intel MKL can be installed using [Intel's YUM repositories](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-yum-repo).

First add the repo and install the GPG key
```shell
dnf install -y dnf-plugins-core
dnf config-manager --add-repo https://yum.repos.intel.com/mkl/setup/intel-mkl.repo
rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
```

Then install the 64-bit version of the most recent MKL release
```shell
LATEST_MKL=$(dnf search -y intel-mkl-64bit | sed -r -e :a -e '$!N;s/\n[[:space:]]+://;ta' -e 'P;D' | grep ^intel | sort | tail -n 1 | awk '{print $1}')
dnf install -y $LATEST_MKL
```

Additionally every developer should have git installed
```shell
dnf install -y git
```

Further more having python including some handy libraries is recommended
(and required for some parts of the testsuite).
```shell
dnf install -y python3 python3-lxml python3-numpy python3-h5py python3-scipy
```
