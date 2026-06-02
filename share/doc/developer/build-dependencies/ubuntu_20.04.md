CFS build-dependencies for Ubuntu 20.04 (broken)
================================================

> you will need to install a recent CMake version, the one from the packages is too old

To make the package installation non-interactive we set
```shell
export DEBIAN_FRONTEND=noninteractive
```
which is probably only needed in automated installs (and for testing this instruction via docker).

For the typical build config, we need
```shell
apt-get update
apt-get install -y gcc g++ gfortran cmake patch m4
```

Intel MKL can be installed from [Intel's oenAPI repos](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo).
First, make sure we have the dependencies for the instructions
```shell
apt-get install -y wget gnupg apt-transport-https sudo
```

Then get the GPG key and save it in the keyrings directory
```shell
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
| gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
```
Add the repo and update
```shell
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
apt-get update
```

Install the latest MKL, the output will tell you which MKL version the meta-package refers to
```shell
apt-get install -y intel-oneapi-mkl-devel
```

Additionally, every developer should have Git installed
```shell
apt-get install -y git
```