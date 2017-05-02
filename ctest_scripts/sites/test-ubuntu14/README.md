Nightly builds & tests on test-ubuntu14
=======================================

Gerneral Layout
---------------
Varaibles common to all tests are set in `init_test.cmake` which is incuded by every `*.ctest` file.
Stuff common for all tests should be defined in `ctest_scripts/shared/test_macros.cmake`.
The deploy is done by `deploy_draco.cmake`, which is included from `*_deploy.cmake`.

Workfow for Testing
-------------------
Start working on single tests first running the separately e.g.
```
ctest -VV -S /home/testuser/cfs/trunk/ctest_scripts/sites/test-ubuntu14/linux64_trunk_icc_release.ctest
```
Then continue running a single test using the nightly-script. 
To run only one test adapt `generate_test_name_list.cmake` and run
```
ctest -VV -S /home/testuser/cfs/trunk/ctest_scripts/nightly_test.cmake
```
This way the deploy scripts can be tested too.

You can work on the git repo directly. Git refuses to update if you have local modified content.
In order to commit your changes you need to put them into your own repo using a patch, because 'testuser-wien' is not allowed to commit.


Setup
-----
The machine is a standard ubuntu vm.
A user 'testuser' with password 'CFStester2016' exists.
Sudo is disabled for this user.

### Installed packages
```
subversion git git-svn cmake gfortran vim vim-gtk cmake-curses-gui
```

### Directory Structure
```
mkdir /home/testuser/cfs
mkdir /home/testuser/cfs/build
```

### Checkout Source

Clone the repos

```shell
cd /home/testuser/cfs
git svn clone --stdlayout --prefix origin/ https://cfs.mdmt.tuwien.ac.at/svn/CFS++ trunk
git svn clone --stdlayout --prefix origin/ https://cfs.mdmt.tuwien.ac.at/svn/CFS++_TEST testsuite
```

You have to use the user `testuser-wien` with password `Ti.i8aa8` for svn.

### Schedule Nightly builds
Run `crontab -e` and enter
```
 0 2    *   *   * /share/programs/bin/ctest -S /home/testuser/cfs/trunk/ctest_scripts/nightly_test.cmake -VV > /tmp/nightly_tests.log 2>&1
 0 0    *   *   7 /bin/rm -rf /home/testuser/cfs/build/* > /tmp/nightly_clear.log 2>&1
```

This runs nightly tests at 02:00 every night and deletes the build directory (including cache) once a week (on Sunday).

### ssh keys for deploy
To deploy the tests to draco we use ssh and key-based login.
To setup:
1. create users on both machines
2. allow key based login from test-ubuntu14 to draco
