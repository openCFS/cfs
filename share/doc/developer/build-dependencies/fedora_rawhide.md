CFS build-dependencies for Fedora
=================================

For the default build config, we need
```shell
dnf install -y make gcc gcc-c++ gcc-gfortran cmake patch m4 findutils diffutils gawk
```

Intel MKL can be installed using the package manager from [Intel's oneAPI repositories](https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/yum-dnf-zypper.html#yum-dnf-zypper).

First, create the repo file
```shell
tee > /tmp/oneAPI.repo << EOF
[oneAPI]
name=Intel® oneAPI repository
baseurl=https://yum.repos.intel.com/oneapi
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
EOF
```
then copy it to the correct place and import the GPG key
```shell
mv /tmp/oneAPI.repo /etc/yum.repos.d
rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
```
Finally, install the most recent MKL release
```shell
yum install -y intel-oneapi-mkl-devel
```

Additionally, every developer should have Git installed
```shell
dnf install -y git
```

Furthermore, having Python including some handy libraries is recommended
(and required for some parts of the Testsuite).
```shell
dnf install -y python3 python3-lxml python3-numpy python3-h5py python3-scipy python3-vtk
```
