% CFS++ Nightly Binaries README
% Simon Triebenbacher
% 05-10-2013

## Description

This   directory   contains  the   nightly   binaries   of  CFS++   built   on
`wiki.mdmt.tuwien.ac.at`. Each archive  is named according to the  site it has
been built on and the name of its corresponding `.ctest` script. For the exact
build settings used  for each of the  tests cf. `CFS_SOURCE_DIR/ctest_scripts`
and the directories  below. You can also  check the features of  a given CFS++
binary by issuing  `CFS_BINARY_DIR/bin/cfs -v`. In the  following sections the
individual binary archives get described in more detail.

## Linux

### General Linux Binaries: `oracle5_linux[32|64]_fespace_icc14_release.zip`

The [32-bit](oracle5_linux32_fespace_icc14_release.zip) and [64-bit FeSpace
binaries](oracle5_linux64_fespace_icc14_release.zip) should  run on  every
Linux distribution dating back as far as 2006. This is possible  by
compiling the binaries on Oracle Enterprise Linux 5,  which is a quite old
distribution, but using the latest Intel compiler.   For these binaries
`USE_OPENMP=ON` has been set,  enabling all  parallel features  of CFS++ in
addition to  those of Pardiso  from Intel  MKL.  On  `sedici` and  `dodici`
the  same binaries  are available under
`/opt/pckg/cfs_nightly/[fespace|trunk]_icc`.


### Ubuntu 12.04 LTS: `precise_linux64_[fespace|trunk]_gcc_release.zip`

We   provide  64-bit   [FeSpace](precise_linux64_fespace_gcc_release.zip)  and
[Trunk binaries](precise_linux64_trunk_gcc_release.zip)  for Ubuntu  12.04 LTS
which is used  on the Linux Laptops at the  department in Vienna. `USE_OPENMP`
is `OFF` for these binaries but Pardiso can be used in parallel nonetheless by
setting `OMP_NUM_THREADS=N`.

### CentOS 6: `oracle6_linux64_[fespace|trunk]_gcc_release.zip`

We   provide  64-bit   [FeSpace](oracle6_linux64_fespace_gcc_release.zip)  and
[Trunk binaries](oracle6_linux64_trunk_gcc_release.zip) binaries  for CentOS 6
which  is used  on the  Linux Cluster  machines `sedici`  and `dodici`  at the
department in Vienna. On these machines  the same binaries are available under
`/opt/pckg/cfs_nightly/[fespace|trunk]_gcc_release.zip`.  The   binaries  will
also run  on the  other Redhat  Enterprise Linux (RHEL)  6 clones  like Oracle
Enterprise Linux  6 or  Scientific Linux  6. `USE_OPENMP`  is `OFF`  for these
binaries  but  Pardiso  can  be   used  in  parallel  nonetheless  by  setting
`OMP_NUM_THREADS=N`.

## Windows

The
[32-bit binaries](winxp32_fespace_mingw32_release.zip)  in `winxp32_fespace_mingw32_release.zip`
should run  on all Windows releases  starting with XP.  Due  to limitations in
the              MinGW              runtime             libs,              the
[64-bit binaries](oracle6_linux64_fespace_mingw64_release.zip) in `oracle6_linux64_fespace_mingw64_release.zip`
only run on releases from Windows Vista onwards.

## Mac OS X: `linux32_fespace_macosx[32|64]_release.dmg`

The            [32-bit](oracle6_linux32_fespace_macosx32_release.dmg)            and
[64-bit   OS  X   binaries](oracle6_linux32_fespace_macosx64_release.dmg)  have   been
compiled using  a [custom  built Linux-to-OSX cross-compiler  toolchain](http://www.cmake.org/pipermail/cmake/2013-September/055881.html) and
run on all OS  X versions starting from Snow Leopard (10.6)  but may even work
on  older releases.  Please follow  the  instructions in  the README  document
inside the disk images.

## Mirrors

The     nightly     binaries     are    available     via     two     mirrors.
[`wiki.mdmt`](http://wiki.mdmt.tuwien.ac.at/cfs_nightly)  is   only  available
from                                 within                                 TU
Vienna. The [CFS++ development server](@CFS_DS_WEBDAV@/nightly-builds)
is available to any registered CFS++ developer.

## Documentation

Along  with the  binaries  also  the documentation  gets  rebuilt.  It can  be
accessed    via   [`wiki.mdmt`](http://wiki.mdmt.tuwien.ac.at/cfs_docu)    and
the [development server](@CFS_DS_WEBDAV@/doc/cfs_docu_nightly.zip).

