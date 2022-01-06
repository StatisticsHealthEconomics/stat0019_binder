FROM rocker/binder:4.1.2

## Declares build arguments
ARG NB_USER
ARG NB_UID

COPY --chown=${NB_USER} . ${HOME}

ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN echo "Checking for 'apt.txt'..." \
        ; if test -f "./binder/apt.txt" ; then \
        apt-get update --fix-missing > /dev/null\
        && xargs -a ./binder/apt.txt apt-get install --yes \
        && apt-get clean > /dev/null \
        && rm -rf /var/lib/apt/lists/* \
        && dpkg  -i ./binder/openbugs_3.2.2-1_amd64.deb && rm openbugs_3.2.2-1_amd64.deb \
        ; fi
USER ${NB_USER}

## Run an install.R script, if it exists.
RUN if [ -f ./binder/install.R ]; then R --quiet -f ./binder/install.R; fi
