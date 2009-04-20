The scripts in this directory are used at the Chair of Applied Mechatronics in
Klagenfurt   to   produce   nightly   test   reports  for   a   CDash   server
(http://www.cdash.org) running on rom.ist.uni-klu.ac.at (short: rom).


The difficulty  we face  in Klagenfurt  is that we  cannot use  the Subversion
update facility provided  by CTest to update our working  copy for the nightly
build   since   that   would   require   an   automatic   svn+ssh   login   to
lse10.e-technik.uni-erlangen.de. Such a  login without password authentication
is not possible at the moment!


Therefore we use the MacBook Pro of Simon Triebenbacher to run a cronjob which
checks  out the  CFS++  trunk at  a certain  time.  Since MacOS  X caches  the
Subversion password this works automatically.


The exact order of operations on the MacBook is as follows:

- unpack previous working copy from a Samba share

- update  Subversion working  copies  of CFS++  (trunk),  CFSDEPS (trunk)  and
  Testsuite (trunk)

- pack updated working copies to Samba share

The  script nightly_svn_update.sh  is responsible  for this  and is  run  as a
cronjob.


As  a next  step the  build hosts  (rom and  ubuntu804 which  is  a VirtualBox
running on rom) unpack the updated archives and start the test procedure.

The exact order of operations on rom is as follows:

- start nightly_test_rom.sh from /etc/crontab

- the script extracts  the contents of the updated  working copy archives from
  the Samba share.

- the  scripts sets up  the build  directories and  deletes the  CFSDEPS cache
  directory

- the script  starts the Ubuntu 8.04  VirtualBox in headless  mode. The Ubuntu
  machine  has in  turn crontab  entries (->  nightly_test_ubuntu804.sh) which
  perform  similar  steps.  At  a   certain  time  the  Ubuntu  machine  halts
  automatically.

- after having started  the Ubuntu VBox in background  the actual test scripts
  for rom get executed.


For  more  information  concerning  the  setup  of  a  CDash  server  and  the
interaction browse the references within the ctest_*.cmake files.
