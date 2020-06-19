ARG BASE_IMAGE

# this is the base image to start from
FROM ${BASE_IMAGE}

# set the arguments from the command line
ARG IMAGE
ARG TAG

LABEL description="CFS++ devel image based on ${IMAGE}:${TAG}"

# install build dependencies
COPY share/scripts/mdsh .
COPY share/doc/developer/build-dependencies/${IMAGE}_${TAG}.md .
RUN ./mdsh --eval ${IMAGE}_${TAG}.md | bash -e -x

# add a user in the container
RUN useradd -ms /bin/bash developer
USER developer
