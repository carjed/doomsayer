<!-- # DOOMSAYER -->

![](https://raw.githubusercontent.com/carjed/doomsayer/master/assets/doomsayer_logo2.png)

[![DOI](https://zenodo.org/badge/95728792.svg?style=flat)](https://zenodo.org/badge/latestdoi/95728792) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](http://www.jedidiahcarlson.com/docs/doomsayer) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://opensource.org/licenses/MIT)  [![Binder](https://img.shields.io/badge/launch-binder-d06681.svg?style=flat)](https://mybinder.org/v2/gh/carjed/doomsayer/master) [![Build Status](https://travis-ci.org/carjed/doomsayer.svg?style=flat?branch=master)](https://travis-ci.org/carjed/doomsayer)

# Introduction

_**DOOMSAYER** ( **D**etection **O**f **O**utliers using **M**utation **S**pectrum **A**nal**Y**sis in **E**xtremely **R**are variants) is a utility for analyzing patterns of rare, single-nucleotide variants (SNVs) in whole-genome (WGS) or whole-exome sequencing (WES) data._

_The basic intuition behind Doomsayer is that the non-somatic mutation spectra of rare SNVs should have little inter-individual heterogeneity. If an individual's mutation spectrum differs drastically from the expected distribution, it is likely due to cryptic error biases or batch effects rather than genuine biological variation. Doomsayer uses a series of statistical analyses to identify these outlier samples, and provides a diagnostic report summarizing the observed error signatures, helping to ensure rigor and reproducibility in the analysis of WGS/WES data._

_In addition to its purpose as a quality control program, Doomsayer can be applied more generally to study between-sample differences in somatic and germline mutation signatures._

<!-- ### Citation

If you use DOOMSAYER in your research, please cite our [paper](#) (citation pending). -->

------------------------------------

# Setup

## Using Conda (recommended)

The following commands will create a Conda environment named `doomsayer` and install all dependencies from `env.yml`. This environment will install the R binaries and most required R packages; however, some packages are not available via the Conda channels, so must be installed in the environment using the `install.r` script (note that this script will NOT install these packages outside of the `doomsayer` environment).

```{sh}
git clone https://github.com/carjed/doomsayer.git
cd doomsayer

conda env create -n doomsayer -f env.yml
source activate doomsayer

R --quiet -f install.r
```

## Local install

Prerequisites for Doomsayer can also be installed using `pip` and the included `install.r` script:

```{sh}
git clone https://github.com/carjed/doomsayer.git
cd doomsayer

pip install -r pip_reqs.txt

R --quiet -f install.r
```

Note that this method assumes you have `pip` and `R` (version 3.1 or higher) already installed. You will also need `pandoc` (version 1.19) installed in order to render RMarkdown reports. If [RStudio](https://www.rstudio.com/products/rstudio/download/#download) (or [RStudio Server]()) is installed on your system, the necessary pandoc binaries should already be available, and no additional acction is needed.

Debian/Ubuntu users may run the [`check_pandoc.sh`](check_pandoc.sh) script to confirm if the pandoc binaries are installed--if they are not found, the script will attempt to download the binaries to `doomsayer/pandoc/`.

```{sh}
bash check_pandoc.sh
```

Mac users will need to either install RStudio or manually install the binaries to `doomsayer/pandoc/` per the instructions described [here](https://github.com/jgm/pandoc/blob/master/INSTALL.md).

## Docker

For more flexible deployment options, Doomsayer is available as a Docker container. The following command will pull and run the preconfigured image from the [Docker Hub](https://hub.docker.com/):

```
docker run -d --name doomsayer \
  -v /path/to/local/data:/data \ # map directory containing input data
  -p 8888:8888 \ # expose jupyter notebook on port 8888
  start-notebook.sh --NotebookApp.token='' \ # start with token disabled
  carjed/doomsayer
```

You may also clone this repository and build the dockerfile locally, using the following commands:

```{sh}
git clone https://github.com/carjed/doomsayer.git
cd doomsayer

docker build -t latest --force-rm .

docker run -d --name doomsayer \
  -p 8888:8888 \
  start-notebook.sh --NotebookApp.token='' \
  doomsayer
```

In both cases, Doomsayer will be available as a Jupyter notebook server, accessible at http://[machine ip]:8888.

## Binder

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/carjed/doomsayer/develop)

A prebuilt Doomsayer docker image can be accessed via the cloud-based [Binder](https://mybinder.org/v2/gh/carjed/doomsayer/develop) platform.

When launched, this will spawn a Jupyter notebook with an interactive tutorial to guide new users through the various options and use cases for Doomsayer. 

Due to the resource constraints of the public BinderHub server, this should not be used to run Doomsayer on large datasets. However, if you have generated the subtype count matrix locally, you can easily upload this file into a Binder instance and run Doomsayer.

------------------------------------

