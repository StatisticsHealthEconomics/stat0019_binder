FROM rocker/binder:latest 

## Declares build arguments
ARG NB_USER
ARG NB_UID
####ARG git_personal_token

COPY --chown=${NB_USER} . ${HOME}

ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN echo "Checking for 'apt.txt'..." \
        ; if test -f "apt.txt" ; then \
        rm -rf /var/lib/apt/lists && mkdir /var/lib/apt/lists && \ 
        apt-get update --fix-missing > /dev/null && \
        xargs -a apt.txt apt-get install --yes && \
        apt-get clean > /dev/null && \
        rm -rf /var/lib/apt/lists/* && \
        dpkg  -i openbugs_3.2.2-1_amd64.deb && rm openbugs_3.2.2-1_amd64.deb \
        ; fi
USER ${NB_USER}

## Copies the Rstudio configuration file
## Among other things, this tells the container to open Rstudio in ~/practical
COPY rstudio-prefs.json .config/rstudio/rstudio-prefs.json

## Run an install.R script, if it exists.
RUN if [ -f install.R ]; then R --quiet -f install.R; fi
