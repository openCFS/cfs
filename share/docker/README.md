Docker image for CFS++ development
==================================

First [install docker](https://docs.docker.com/install/).

Then you can use the prepared `Dockerfile`, which is also used to
check if the CFS build dependencies are up to date.
By default the latest version of the base image will be 
pulled from [docker hub](https://hub.docker.com/).

To select which image should be used set a environment variables, e.g.
`export IMAGE=fedora` and `export TAG=latest`.
Docker images are always built from the *context* i.e. the current working direcory.
In order to have files in the CFS source available in the expected locations, run docker
from the root CFS source directory. Here we have a [.dockerignore](../../.dockerignore) file,
which ignores all but the required sources.

To build the image go to the root directory of the CFS repo and run
```shell
docker build -t cfs-devel-$IMAGE-$TAG --build-arg BASE_IMAGE=$IMAGE:$TAG --build-arg IMAGE=$IMAGE --build-arg TAG=$TAG -f share/docker/Dockerfile .
```
This will execute our generic [Dockefile](Dockerfile) which, starting from a `BASE_IMAGE`
1. installs build dependencies
2. and creates a user

The `-t` option tags the image.

You can then list all images using
```shell
docker image ls
```

To use the image interactively run `docker run -it cfs-devel-$IMAGE bash`.
In order to have the source code availabe, mount a volume e.g. `docker run -v ~/cfs/CFS:/cfs -it cfs-devel-$IMAGE-$TAG bash` to mount the host dir `~/cfs/CFS` to `/cfs` in the container `cfs-devel-$IMAGE-$TAG`.

Updating Images
---------------

To update the base image use `docker pull $IMAGE:$TAG`.
You can then rebuild the image (see above).
If the instructions in `build-dependencies/$IMAGE_$TAG.md` fail,
you need to fix them.

1. First find a working image version by looking for availabe tags
   * try with `docker build -t cfs-devel-$IMAGE:<tag> --build-arg BASE_IMAGE=$IMAGE:<tag> --build-arg IMAGE=$IMAGE`
   * copy the original `build-dependencies/$IMAGE.md` to a fitting name `build-dependencies/$IMAGE_<tag>.md`
2. Fix the error in the latest version
3. Adapt the build pipeline in `.gitlab-ci.yml` to include both versions (or clean up if the old one is obsolete)

One can run the base image interactively to debug
1. start it by `docker run -it $IMAGE:$TAG bash`
2. run the instructions from `share/doc/developer/build-dependencies/$IMAGE_$TAG.md` in sequence

References and Notes
--------------------

* This instruction is [tested](/.gitlab-ci.yml) using `mdsh`.
* We use a [`.dockerignore` file](https://docs.docker.com/engine/reference/builder/#dockerignore-file)` as the docker *context* should be small.
* [docker-usage-doc](https://docs.docker.com/engine/reference/builder/#usage)
* One could also use [kaniko](https://docs.gitlab.com/ee/ci/docker/using_kaniko.html) to build the image.
