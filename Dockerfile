# See https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile

# Must specigy a tag
FROM giabaio/stat0019:17022021

# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER rstudio
ENV HOME /home/rstudio

RUN echo $HOME
RUN echo $NB_USER

# Copy repo into ${HOME}, make user own $HOME
USER root
COPY . ${HOME}
RUN chown -R rstudio ${HOME}
USER rstudio

## run any install.R script we find
RUN if [ -f install.R ]; then R --quiet -f install.R; fi

