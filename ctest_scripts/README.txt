The scripts in this directory are used at the Chair of Applied Mechatronics in
Klagenfurt   to   produce   nightly   test   reports  for   a   CDash   server
(http://www.cdash.org)  running on rom.ist.uni-klu.ac.at  (short: rom).  It is
accesible through the address http://rom/cdash.


We now use  a password-less ssh key to log into  the Subversion server (lse10)
at  Erlangen. Therefore  we can  now run  automatical updates  of  the working
copies (WCs)  on rom and  sedici. The Ubuntu  8.04 Virtual Box running  on rom
does not run 'svn up' by itself, but just gets the updated WCs from rom.


Here is a summary of the steps  required on the build host rom. For sedici the
steps are quite similar, except that it does not start a VBox.

The exact order of operations on rom is as follows:

- start nightly_test_rom.sh from /etc/crontab

- the   script   checks   if    proper   working   copies   are   present   in
  /home/strieben/Documents/dev/NIGHTLY   and  performs   checkouts/updates  if
  necessary.

- the  scripts sets up  the build  directories and  deletes the  CFSDEPS cache
  directory

- the  contents  of  /home/strieben/Documents/dev/NIGHTLY/CFSDEPS_NIGHTLY  get
  copied to /opt/CFSDEPS so that the global CFSDEPS WC is allways up-to-date

- the script  starts the Ubuntu 8.04  VirtualBox in headless  mode. The Ubuntu
  machine  has in  turn crontab  entries (->  nightly_test_ubuntu804.sh) which
  perform  similar  steps.  At  a   certain  time  the  Ubuntu  machine  halts
  automatically.

- after having started  the Ubuntu VBox in background  the actual test scripts
  for rom get executed.

- if  the builds  returned with  no errors  the CFS++  binaries get  copied to
  /opt/pckg/cfs_nightly/trunk_[gi]cc


For  more  information  concerning  the  setup  of  a  CDash  server  and  the
interaction please browse the references within the ctest_*.cmake files.

