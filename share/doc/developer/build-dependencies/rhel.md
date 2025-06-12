CFS build dependencies for RHEL and derivatives
===============================================

These instructions are tested on [rockylinux docker containers](https://hub.docker.com/_/rockylinux), but might also work on
[AlmaLinux](https://almalinux.org/),
[Oracle Linux](https://www.oracle.com/linux/),
[CentOS Stream](https://www.centos.org/centos-stream/) and other RHEL derivatives.

For instructions on older systems see:
* [CentOS 7](./centos_7.md)
* [RockLinux 8](./rockylinux_8.md)

For the default build configuration, we need
```shell
dnf install -y make gcc gcc-c++ gcc-gfortran cmake patch m4 findutils diffutils
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

Additionally, every developer should have git installed
```shell
dnf install -y git
```