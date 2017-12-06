FROM jupyter/all-spark-notebook:c7fb6660d096

# FROM rocker/binder:3.4.2
#
# # Copy repo into ${HOME}, make user own $HOME
# USER root
# COPY . ${HOME}
# RUN chown -R ${NB_USER} ${HOME}
# USER ${NB_USER}
#
# ## run any install.R script we find
# RUN if [ -f install.R ]; then R --quiet -f install.R; fi
#
# FROM rocker/tidyverse:3.4.2
#
# RUN apt-get update && \
#     apt-get -y install python3-pip && \
#     pip3 install --no-cache-dir notebook==5.2 && \
#     apt-get purge && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/*
#
# ENV NB_USER rstudio
# ENV NB_UID 1000
# ENV HOME /home/rstudio
# WORKDIR ${HOME}
#
# USER ${NB_USER}
#
# # Set up R Kernel for Jupyter
# RUN R --quiet -e "install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest'))"
# RUN R --quiet -e "devtools::install_github('IRkernel/IRkernel')"
# RUN R --quiet -e "IRkernel::installspec()"
#
# # Make sure the contents of our repo are in ${HOME}
# COPY . ${HOME}
# USER root
# RUN chown -R ${NB_UID}:${NB_UID} ${HOME}
# USER ${NB_USER}

# Run install.r if it exists
RUN if [ -f install.r ]; then R --quiet -f install.r; fi

RUN pip install -r --no-cache-dir pip_reqs.txt

# RUN R --quiet -e "install.packages(c('heatmaply', 'viridis'), repos = 'http://cran.us.r-project.org')"
