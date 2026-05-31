Build Dependencies for CFS
==========================

The files in this directory contain instructions on how to set up your system for openCFS development.
Check the individual markdown files for instructions for your distribution, e.g. [`fedora.md` for the latest stable Fedora](fedora.md).
Just copy and paste the commands in them into a (root) terminal one by one.
For the brave: use `sudo ./share/scripts/mdsh share/doc/developer/build-dependencies/<IMAGE>_<TAG>.md`.

The instructions are regularly tested in the build [pipeline](/.gitlab-ci.yml) in recent, official docker containers.

Maintainance
------------

Please follow the filename convention adapted from the naming of the official docker containers on [docker hub](https://hub.docker.com), i.e. `<IMAGE>:<tag>` e.g. `fedora:latest` for the latest stable Fedora release.
We use `<IMAGE>_<tag>.md` since Windows does not support colons (`:`) in file names.
The instruction file for the current version has no version extension, e.g. is called `rhel.md` for Red Hat Enterprise Linux derivatives.
In order to indicate which instruction is usable for which docker-container we create symlinks, e.g. `rockylinux_8.md -> rhel.md`.

Typically, instructions will keep on working for newer OS versions.
Once an instruction breaks,
 1. replace the symlink with the last working version (and update symlinks for even older versions, if present)
 2. create a new instruction for the latest version

Instructions can be tested locally by using our [docker config](/share/docker/README.md).
Instructions are also tested in the _platform_ stage of the [CI pipeline](../../../../.gitlab-ci.yml).
When old instructions are removed from there, keep the last working copy as markdown file: In case the last tested version symlinks to the `*_latest.md`, one needs to create the file to ensure it is not changed when updating the latest version.
