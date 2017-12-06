FROM rocker/tidyverse:3.4.2

LABEL maintainer="Jedidiah Carlson <jed.e.carlson@gmail.com>"

# Install all OS dependencies for notebook server that starts but lacks all
# features (e.g., download as all possible file formats)
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get -yq dist-upgrade \
 && apt-get install -yq --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    sudo \
    locales \
    fonts-liberation \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

# Install Tini
RUN wget --quiet https://github.com/krallin/tini/releases/download/v0.10.0/tini && \
    echo "1361527f39190a7338a0b434bd8c88ff7233ce7b9a4876f3315c22fce7eca1b0 *tini" | sha256sum -c - && \
    mv tini /usr/local/bin/tini && \
    chmod +x /usr/local/bin/tini

# Configure environment
ENV CONDA_DIR=/opt/conda \
    SHELL=/bin/bash \
    NB_USER=rstudio \
    NB_UID=1000 \
    NB_GID=100 \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8
ENV PATH=$CONDA_DIR/bin:$PATH \
    HOME=/home/$NB_USER

WORKDIR ${HOME}

# ADD fix-permissions /usr/local/bin/fix-permissions
RUN mkdir -p $CONDA_DIR && \
  chown $NB_USER:$NB_GID $CONDA_DIR

# ADD fix-permissions /usr/local/bin/fix-permissions
# # Create jovyan user with UID=1000 and in the 'users' group
# # and make sure these dirs are writable by the `users` group.
# RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \


# USER $NB_USER

# # Setup work directory for backward-compatibility
# RUN mkdir /home/$NB_USER/work && \
#     fix-permissions /home/$NB_USER

# Install conda as jovyan and check the md5 sum provided on the download site
ENV MINICONDA_VERSION 4.3.30
RUN cd /tmp && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh && \
    echo "0b80a152332a4ce5250f3c09589c7a81 *Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh" | md5sum -c - && \
    /bin/bash Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh && \
    $CONDA_DIR/bin/conda config --system --prepend channels conda-forge && \
    $CONDA_DIR/bin/conda config --system --set auto_update_conda false && \
    $CONDA_DIR/bin/conda config --system --set show_channel_urls true && \
    $CONDA_DIR/bin/conda update --all --quiet --yes && \
    conda clean -tipsy

# Install Jupyter Notebook and Hub
RUN conda install --quiet --yes \
    'notebook=5.2.*' \
    'jupyterhub=0.8.*' \
    'jupyterlab=0.29.*' \
    && conda clean -tipsy && \
    jupyter labextension install @jupyterlab/hub-extension@^0.6.0 && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging
# USER root

# EXPOSE 8888
# WORKDIR $HOME

# # Configure container startup
# ENTRYPOINT ["tini", "--"]
# CMD ["start-notebook.sh"]
#
# # Add local files as late as possible to avoid cache busting
# COPY start.sh /usr/local/bin/
# COPY start-notebook.sh /usr/local/bin/
# COPY start-singleuser.sh /usr/local/bin/
# COPY jupyter_notebook_config.py /etc/jupyter/
# RUN fix-permissions /etc/jupyter/

# add contents of repo to ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID}:${NB_UID} ${HOME}

USER ${NB_USER}

# run install.r script to load R package dependencies
# RUN R --quiet -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"
RUN if [ -f install.r ]; then R --quiet -f install.r; fi

# create conda environment and install Python library dependencies
RUN conda env create -n doomsayer -f env.yml

# activate the conda environment
ENV PATH /opt/conda/envs/doomsayer/bin:$PATH

# old version; installs via pip
# causes problems with cyvcf2--better to use bioconda
# RUN pip install -r pip_reqs.txt
