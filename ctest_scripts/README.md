Nightly Build/Test scripts for CFS++
====================================

Description
-----------

The scripts in  this directory are used to produce nightly  test reports for a
[CDash](http://www.cdash.org)            server           running           at
[lse17](https://lse17.e-technik.uni-erlangen.de:2001/cdash).

Furthermore, the scripts are used  to build nightly binaries for certain Linux
 and  Windows platforms  at [Measurement  and Actuator  Division,  Institute of
 Mechanics and Mechatronics/E325-A4, TU Wien](http://www.mec.tuwien.ac.at).

Overview
--------

We  use   a  password-less  SSH  key   to  log  into   the  Subversion  server
[lse17](https://lse17.e-technik.uni-erlangen.de:2001/svn)                    at
Erlangen. Therefore we can run automatical updates of the working copies (WCs)
on          `wiki.mdmt.tuwien.ac.at`          from         the          script
`sites/wiki/site_specific_init.cmake`.

A  number of  [Virtual Boxes](http://virtualbox.org)  get started  through the
same     script     by     [Vagrant](http://www.vagrantup.com)     and     the
`nightly_test.cmake` script  gets executed inside  them by the  [Vagrant shell
provisioning mechanism](http://docs.vagrantup.com/v2/provisioning/shell.html).

The  Vagrant base  `.box` image  files are  stored locally  on  `wiki.mdmt` at
`/var/www/html/boxes`.   They are  therefore, accessible  through  the [Apache
server](http://wiki.mdmt.tuwien.ac.at/boxes) from within  TU Wien.  The images
are        also        mirrored         on        the        [lse17        FTP
server](ftp://lse17.e-technik.uni-erlangen.de:40065/boxes).

The  VBoxes get  configured by  the `Vagrantfile`s  in the  corresponding site
directories     to      mount     the     directories      `/opt/pckg`     and
`/home/testuser/Documents/dev/NIGHTLY`  as  shared folders.   This  way it  is
possible for them to access the same working copies.


Summary of Steps on _wiki.mdmt_
-------------------------------


Here is a summary of the steps required on the build host _wiki.mdmt_:

- Set  proper   permissions  for  folders  involved  in   nightly  build  from
  `/etc/crontab` This may  be necessary in case some other  user (esp. root) has
  tempered with them.

        45 23  *  *  *  root /opt/pckg/cmake-2.8.9/bin/cmake -P \
        /home/testuser/Documents/dev/NIGHTLY/CFS_FESPACE_NIGHTLY/ctest_scripts/sites/wiki/set_permissions.cmake

- Create directory for log files for the site _wiki_ from `/etc/crontab`

        20  1  *  *  *  testuser mkdir -p \
        /home/testuser/Documents/dev/NIGHTLY/CFS_FESPACE_NIGHTLY/ctest_scripts/sites/wiki/logs

- Start `nightly_test.cmake` at 1:30  from `/etc/crontab` and write the output
  to a  log file.  The log can  be found  in `nightly_test.log` in  the `logs`
  subdirectory.

        30  1  *  *  *  testuser /opt/pckg/cmake-2.8.9/bin/ctest -S \
        /home/testuser/Documents/dev/NIGHTLY/CFS_FESPACE_NIGHTLY/ctest_scripts/nightly_test.cmake > \
        /home/testuser/Documents/dev/NIGHTLY/CFS_FESPACE_NIGHTLY/ctest_scripts/sites/wiki/logs/nightly_test.log > 2>&1

- The script sets  the proper global variables. Especially  the site directory
  gets determined  and the scripts  in the `sites/wiki` directory  get sourced
  and executed.

- The sourced  scripts are responsible  for checking if proper  working copies
  are present  in `/home/testuser/Documents/dev/NIGHTLY` and checkouts/updates
  get performed if necessary.

- The  nightly VBoxes  (_oracle6_ and  _precise_)  get started  by Vagrant  in
  non-blocking  mode and  nightly tests  are started  inside them  through the
  shell  provisioning mechanism. Log  files are  stored for  all sites  in the
  corresponding `logs` subdirectories.   These include `nightly_test.log`, the
  CTest logs  from `${CFS_BINARY_DIRECTORY}/Testing` and  the `CMakeCache.txt`
  files.

- The tests for _wiki.mdmt_ (i.e. Intel compiler 12.1 and 13.1) get executed.

- Wait until 7:30.

- The  contents of  `/home/testuser/Documents/dev/NIGHTLY/CFSDEPS_NIGHTLY` get
  copied to  `/opt/pckg/CFSDEPS` so  that the global  CFSDEPS WC for  Trunk is
  allways up-to-date.

- The  contents   of  `/home/testuser/Documents/dev/NIGHTLY/CFSDEPSCACHE`  get
  copied  to `/opt/pckg/CFSDEPSCACHE`  so that  the global  CFSDEPS  cache for
  Trunk and FeSpace is allways up-to-date.

- If  the builds  inside  the _oracle6_  (binaries  for `RHEL_6_X86_64`  which
  includes  CentOS, Scientific  and  Oracle Linux  and `MINGW_X86_64`  Windows
  64-bit) and  _precise_ (binaries for  `UBUNTU_12.04_X86_64`) VBoxes returned
  with     no    errors,     the    CFS++     binaries    get     packed    to
  `/opt/pckg/cfs_nightly/archives`.

- At   7:30  the   packed   binaries  and   documentation   get  copied   from
  `/opt/pckg/cfs_nightly/archives`  to  `/var/www/html/cfs_[nightly|docu]` and
  are    therefore   also    accessible    from   within    TU   Wien    under
  [binaries](http://wiki.mdmt.tuwien.ac.at/cfs_nightly)                    and
  [documentation](http://wiki.mdmt.tuwien.ac.at/cfs_|docu).

- The      binaries     get      also     unpacked      and      copied     to
  `/opt/pckg/cfs_nightly/trunk_[gi]cc`  and  may  therefore be  accessed  from
  _wiki_, _sedici_ or _dodici_.

- Since hard  disk drive  space is  very limited on  the `/home`  partition of
  _wiki_,  we start the  remaining VBoxes  one after  another. The  Windows XP
  32-bit VBox gets  started every day. Afterwards, depending  on the week day,
  either the VBoxes  from the SuSE, RedHat or  Debian/Ubuntu family of distros
  get started in blocking mode.

- Usually, the  test schedule  should be finished  by 14:00. Therefore,  it is
  best to wait  until afternoon, if changes to the  nightly test schedule have
  to be made or VBoxes have to started manually.

    To start  a VBox manually, perform the following steps:

    +  Check  if  VBoxes  are  still running  by  issueing:  `vboxmanage  list
    runningvms`

    + Shutdown or kill VBoxes. E.g `vboxmanage controlvm VBOXNAME poweroff`

    + Go to the corresponding site directory (e.g. `sites/oracle6`)

    + Delete old `.vagrant` directory: `rm -rf .vagrant`

    + Start VBox: `vagrant up [--no-provision]`. The `--no-provision` flag has
    to be specified if the nightly  test inside the VBox should not be started
    automatically once  the machine has booted up.   During the initialization
    phase  the  `.box`  files  get   unpacked  to  a  temporary  directory  by
    Vagrant. Since  hdd space is low  on `/home`, this temporary  dir has been
    reset by  setting the `VAGRANT_HOME` environment variable  for testuser in
    `$HOME/.bashrc`.


Further Information
-------------------

A  flow  diagram  of  the  general  nightly test  workflow  is  documented  in
`nightly_docu.graphml`.   It   can  be  opened   with  the  free   [yEd  graph
editor](http://www.yworks.com).

For  more  information  concerning  the  setup  of  a  CDash  server  and  the
interaction please browse the references within the `ctest_*.cmake` files.

---

This       document       has       been      written       in       [Markdown
](http://daringfireball.net/projects/markdown/syntax) syntax  on 19/08/2013 by
Simon  Triebenbacher .   It can  be converted  e.g. to  a HTML  page  by using
[Pandoc](http://johnmacfarlane.net/pandoc):  `pandoc -f  markdown  -t html  -s
README.md > README.html`
