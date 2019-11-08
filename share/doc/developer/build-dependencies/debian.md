CFS++ build dependencies for Debian
===================================

For the minimal build config we need
```shell
apt-get update
apt-get install -y gcc g++ gfortran cmake patch m4
```

For the Intel MKL installer (`USE_MKL=ON`) we need `cpio`
```shell
apt-get install -y cpio
```

Intel MKL can also be installed from [Intel's APT repos](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo).
First make sure we have the dependecies for the instructions
```shell
apt-get install -y wget gnupg apt-transport-https ca-certificates
```

Then get the GPG-key and add it
```shell
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
```
Add the repo and update
```shell
sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
apt-get update
```

Determine the latest version of the virtual package `intel-mkl-64bit` and install it
```shell
LATEST_MKL=$(apt-cache search intel-mkl-64bit | tail -n 1 | awk '{print $1}')
apt-get install -y $LATEST_MKL
```

Additionally, every developer should have git installed
```shell
apt-get install -y git
```
