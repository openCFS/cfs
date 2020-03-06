Build Dependencies for CFS
==========================

The files in this directory are tested in the build [pipeline](/.gitlab-ci.yml) in docker containes.
The filename convention follows form the naming of the official docker containers on [docker hub](https://hub.docker.com), i.e. `<IMAGE>:<tag>`, where the default tag `latest` can be omitted.

Typically older instruction will keep on working for newer OS versions, thus we use symlinks to point to them.
Once an instruction breaks,
 1. replace the old-version symlink with the working instruction
 2. create a new instruction for the latest version

Instructions can be tested locally by using our [docker config](/share/docker/README.md). 
