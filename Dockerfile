FROM rocker/binder:4.1.2 

## Declares build arguments
ARG NB_USER
ARG NB_UID

COPY --chown=${NB_USER} . ${HOME}

ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN echo "Checking for 'apt.txt'..." \
        ; if test -f "apt.txt" ; then \
        apt-get update --fix-missing > /dev/null\
        && xargs -a apt.txt apt-get install --yes \
        && apt-get clean > /dev/null \
        && rm -rf /var/lib/apt/lists/* \
        && dpkg  -i openbugs_3.2.2-1_amd64.deb && rm openbugs_3.2.2-1_amd64.deb \
#######################################################################################################
## This installs r2u for package managements through apt (source: https://github.com/eddelbuettel/r2u)
## NB: There are some issues with sudo --- need to finalise/check!
##        && apt-get update -qq && apt-get install --yes --no-install-recommends wget ca-certificates gnupg \
##        && wget -q -O- https://eddelbuettel.github.io/r2u/assets/dirk_eddelbuettel_key.asc | tee -a /etc/apt/trusted.gpg.d/cranapt_key.asc \
##        && echo "deb [arch=amd64] https://r2u.stat.illinois.edu/ubuntu jammy main" > /etc/apt/sources.list.d/cranapt.list \
##        && apt-get update -qq \
##        && apt-get install --yes --no-install-recommends python3-dbus python3-gi python3-apt \
##        && Rscript -e 'install.packages("bspm")' \
##        && RHOME=$(R RHOME) \
##        && echo "suppressMessages(bspm::enable())" >> ${RHOME}/etc/Rprofile.site \
##        && echo "options(bspm.version.check=FALSE)" >> ${RHOME}/etc/Rprofile.site \
#######################################################################################################
        ; fi
USER ${NB_USER}

## Run an install.R script, if it exists.
RUN if [ -f install.R ]; then R --quiet -f install.R; fi
