###############################################################################
# Modified from https://github.com/jupyter/docker-stacks/blob/master/base-notebook/Dockerfile
# Distributed under the terms of the Modified BSD License.
###############################################################################
FROM jupyter/scipy-notebook:c7fb6660d096
# FROM jupyter/minimal-notebook:033056e6d164
# FROM jupyter/all-spark-notebook:c7fb6660d096
# FROM rocker/tidyverse:3.4.2
# FROM jupyter/r-notebook

LABEL maintainer="Jedidiah Carlson <jed.e.carlson@gmail.com>"

# Configure environment
# ENV CONDA_DIR=/opt/conda \
#     NB_UID=1000 \
#     NB_GID=100 \
#     SHELL=/bin/bash \
#     HOME=/home/$NB_USER
# ENV PATH=$CONDA_DIR/bin:$PATH \

# create conda install directory
# RUN mkdir -p $CONDA_DIR && \
#   chown $NB_USER:$NB_GID $CONDA_DIR

# WORKDIR ${HOME}

# Force bash always
# RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# ADD fix-permissions /usr/local/bin/fix-permissions

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    fonts-dejavu \
    tzdata \
    gfortran \
    gcc && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER ${NB_USER}
COPY . ${HOME}

USER root
RUN chown -R ${NB_UID} ${HOME}

###############################################################################
# create environment from config file and activate
###############################################################################
RUN conda env create -n doomsayer -f env.yml && \
  conda clean -tipsy

# ENV CONDA_DS_ENV "doomsayer"
# ENV CONDA_ACTIVATE "source activate $CONDA_DS_ENV"
RUN ["/bin/bash", "-c", "source activate doomsayer"]
ENV PATH="/opt/conda/envs/doomsayer/bin:${PATH}"
# ENV CONDA_PREFIX /opt/conda/envs/doomsayer

###############################################################################
# Install R packages
# rmarkdown is installed separately via devtools to resolve issue
# with pandoc (https://github.com/rstudio/rmarkdown/issues/1120)
# --for some reason, this requires tar to be aliased as gtar
# per https://github.com/hadley/devtools/issues/379
###############################################################################
RUN ln -s /bin/tar /bin/gtar
RUN R --quiet -e "devtools::install_github('rstudio/rmarkdown')"
RUN if [ -f install.r ]; then R --quiet -f install.r; fi

USER ${NB_USER}
# FROM jupyter/scipy-notebook:c7fb6660d096
# ADD pip_reqs.txt ./
# ADD env.yml ./

# RUN conda create -n doomsayer-environment python=3.6 anaconda
# RUN source activate doomsayer-environment
# RUN conda env export > environment.yml

# create conda environment and install Python library dependencies


# activate the conda environment


# old version; installs via pip
# causes problems with cyvcf2--better to use bioconda
# RUN pip install -r pip_reqs.txt
