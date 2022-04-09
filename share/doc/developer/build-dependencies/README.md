Build Dependencies for CFS
==========================

The files in this directory contain instructions on how to setup your system for openCFS development.
Check the individaual markdown files for instructions for your distributaion, e.g. [`fedora_latest.md` for the latest stable Fedora](fedora_latest.md).
Just copy and paste the instuctionsthem into a (root) terminal one by one.
For the brave: use `sudo ./share/scripts/mdsh share/doc/developer/build-dependencies/<IMAGE>_<TAG>.md`.

The instructions are regularly tested in the build [pipeline](/.gitlab-ci.yml) in recent, official docker containes.

Maintainance
------------

Please follow the filename convention adapted from the naming of the official docker containers on [docker hub](https://hub.docker.com), i.e. `<IMAGE>:<tag>` e.g. fedora:latest for the latest stable fedora release.
We use `<IMAGE>_<tag>.md` since windows does not support `:` in file names.

Typically older instruction will keep on working for newer OS versions, thus we use symlinks to point to them.
Once an instruction breaks,
 1. replace the old-version symlink with the working instruction
 2. create a new instruction for the latest version

Instructions can be tested locally by using our [docker config](/share/docker/README.md). 
