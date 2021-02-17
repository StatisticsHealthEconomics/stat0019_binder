# See https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile

# Must specigy a tag
FROM giabaio/stat0019:17022021

# Install recent version of Jupyter Notebook
RUN pip install --no-cache-dir notebook==5.*

# Install relevant packages not already available in the main image
RUN R -e "options(repos = \
   list(CRAN = 'http://cran.ma.imperial.ac.uk')); \
   install.packages('R2jags'); \
   install.packages('BCEA')"

# Creates user
ARG NB_USER=stat0019
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} ${NB_USER}

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
