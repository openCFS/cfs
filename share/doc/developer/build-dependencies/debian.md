CFS++ build dependencies for Debian
===================================

For the typical build configuration, we need
```shell
apt-get update
apt-get install -y gcc g++ gfortran cmake patch m4
```

Intel MKL can be installed via the official [Intel oneAPI repository](https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/apt.html).

First, make sure we have the dependencies for the instructions
```shell
apt-get install -y wget gnupg apt-transport-https ca-certificates
```
Then get the PGP-key and add it
```shell
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
```
Add the repo and update
```shell
echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
apt-get update
```
Install the latest MKL version
```shell
apt-get install -y intel-oneapi-mkl-devel
```

Additionally, every developer should have Git installed
```shell
apt-get install -y git
```
