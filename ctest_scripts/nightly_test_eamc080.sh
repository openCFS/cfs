#!/usr/bin/zsh

# this small script can be used to run nightly builds on eamc080.eam
# it defines some helper functions to call 'svn update' in the appropriate
# directories and then the ctest scripts which build cfs and execute
# the tests in the testsuite

homedir=/home/eamnightly
cfsdir=$homedir/cfs-trunk
masterdir=$homedir/cfs-master
cfsdepsdir=$homedir/cfsdeps
mastertestdir=$homedir/cfstest-master
cfstestdir=$homedir/cfstest-trunk

export PATH=/usr/local/bin:$PATH

function preparecfsdeps() {
  # update cfsdeps
  print "preparing cfsdeps..."
  cd $cfsdepsdir
  svn update
  cd $homedir
}

function preparecfstest() {
  # update cfstest
  print "preparing cfstest..."
  cd $cfstestdir
  svn update
  cd $homedir
  print "preparing mastertest..."
  cd $mastertestdir
  svn update
  cd $homedir
}

function preparecfs() {
  # update cfs
  print "preparing cfs..."
  cd $cfsdir
  rm CMakeLists.txt
  rm CTestConfig.cmake

  svn update

  # copy CTestConfig.cmake in root-directory of cfs
  cp $cfsdir/ctest_scripts/CTestConfig.cmake $cfsdir

  # 'patch' CMakeLists.txt for nightly testing
  {
    print ""
    print "ENABLE_TESTING()"
    print "INCLUDE(CTest)" 
  } >> CMakeLists.txt
  

  # here, we need to delete the build directories before we can proceed
  # they must be empty, since configuring and building is done with ctest
  cd $homedir
  if [ -d $homedir/release ]; then
    rm -rf $homedir/release
  fi
  if [ -d $homedir/debug ]; then
    rm -rf $homedir/debug
  fi
  mkdir $homedir/release
  mkdir $homedir/debug
  cd $homedir
}

function preparemaster() {
  # update cfs
  print "preparing master..."
  cd $masterdir
  rm CMakeLists.txt
  rm CTestConfig.cmake

  svn update

  # copy CTestConfig.cmake in root-directory of cfs
  cp $masterdir/ctest_scripts/CTestConfig.cmake $masterdir

  # 'patch' CMakeLists.txt for nightly testing
  {
    print ""
    print "ENABLE_TESTING()"
    print "INCLUDE(CTest)" 
  } >> CMakeLists.txt
  

  # here, we need to delete the build directories before we can proceed
  # they must be empty, since configuring and building is done with ctest
  cd $homedir
  if [ -d $homedir/release-master ]; then
    rm -rf $homedir/release-master
  fi
  if [ -d $homedir/debug-master ]; then
    rm -rf $homedir/debug-master
  fi
  mkdir $homedir/release-master
  mkdir $homedir/debug-master
  cd $homedir
}

# update all directories
preparecfsdeps
preparecfstest
preparecfs
preparemaster

cd $homedir

# date | mailx -s "starting nightly builds on eam80" schury@am.uni-erlangen.de

# perform the tests for release and debug build of the master
print "running release tests..."
ctest -S $masterdir/ctest_scripts/ctest_eam_gcc_nightly_master_release.cmake -V > $homedir/tests.master.release.log 2>&1
print "running debug tests..."
ctest -S $masterdir/ctest_scripts/ctest_eam_gcc_nightly_master_debug.cmake -V   > $homedir/tests.master.debug.log   2>&1

# perform the tests for release and debug build of the trunk
print "running release tests..."
ctest -S $cfsdir/ctest_scripts/ctest_eam_gcc_nightly_release.cmake -V > $homedir/tests.release.log 2>&1
print "running debug tests..."
ctest -S $cfsdir/ctest_scripts/ctest_eam_gcc_nightly_debug.cmake -V   > $homedir/tests.debug.log   2>&1

# date | mailx -s "finished nightly builds on eam80" schury@am.uni-erlangen.de
