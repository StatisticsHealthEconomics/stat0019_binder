FROM rocker/verse:latest
MAINTAINER Damir Perisa 

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
# Installs JAGS
  jags \
  libc6-i386

# Installs OpenBUGS
RUN wget -nd -P /tmp http://pj.freefaculty.org/Debian/squeeze/amd64/openbugs_3.2.2-1_amd64.deb
RUN dpkg -i /tmp/openbugs_3.2.2-1_amd64.deb && rm /tmp/openbugs_3.2.2-1_amd64.deb 

# adding deps separately so it may build in dockerhub (works on my WS)
RUN apt-get install -y r-cran-rcpp r-cran-rcppeigen

#RUN install2.r --error \
#  --repos "https://stat.ethz.ch/CRAN/" \
#  rstan \
#  rstantools \ 
#  rstanarm \
#  bayesplot \
#  brms \
#  tidybayes
  
RUN install2.r --error \
  --repos "https://stat.ethz.ch/CRAN/" \
  rjags \
  R2jags \  
  R2OpenBUGS \
  BCEA
  
#RUN install2.r --error \
#  --repos "https://inla.r-inla-download.org/R/stable" \
#  INLA 

# Now runs the server
# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
USER ${USER}
