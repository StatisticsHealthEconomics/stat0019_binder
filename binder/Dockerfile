FROM rocker/binder:4.1.2

## Declares build arguments
ARG NB_USER
ARG NB_UID

COPY --chown=${NB_USER} . ${HOME}

ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN echo apt-get update --fix-missing > /dev/null \
         && apt-get install --yes \
         default-jdk \
         gfortran \
         libudunits2-dev \
         libgdal-dev \
         libv8-dev \
         libssl-dev \
         libsasl2-dev \
         autoconf \
         automake \
         libudunits2-dev \
         libicu-dev \
         libeigen3-dev \
         gcc-multilib \
         g++-multilib \
         jags \
         && apt-get clean > /dev/null \ 
         && rm -fr /var/lib/apt/lists/* \
         && dpkg -i openbugs_3.2.2-1_amd64.deb && rm openbugs_3.2.2-1_amd64.deb

#RUN echo "Checking for 'apt.txt'..." \
#        ; if test -f "apt.txt" ; then \
#        apt-get update --fix-missing > /dev/null\
#        && xargs -a apt.txt apt-get install --yes \
#        && apt-get clean > /dev/null \
#        && rm -rf /var/lib/apt/lists/* \
#        ; fi
USER ${NB_USER}

RUN R --quiet -f install.R

## Run an install.R script, if it exists.
#RUN if [ -f install.R ]; then R --quiet -f install.R; fi
