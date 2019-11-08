ARG BASE_IMAGE

# this is the base image to start from
FROM ${BASE_IMAGE}

# for some reason the BASE_IMAGE is not available any more
# we need to set the same thing again
ARG IMAGE

LABEL description="CFS++ devel image based on ${IMAGE}"

# install build dependencies
COPY share/scripts/mdsh .
COPY share/doc/developer/build-dependencies/${IMAGE}.md .
RUN ./mdsh --eval ${IMAGE}.md | bash -e

# add a user in the container
RUN useradd -ms /bin/bash developer
USER developer
