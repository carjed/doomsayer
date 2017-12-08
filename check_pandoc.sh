#!/bin/bash

os=${OSTYPE//[0-9.]/}

# install python dependencies
# may need to update path to pip
# by default installs to user library location (default ~/.local/bin)
# remove --user if installing system-wide
# echo "Installing python dependencies"
# /usr/bin/pip install --user cyvcf2 numpy pandas Biopython nimfa pyfaidx joblib scipy

# Check for pandoc libraries
if [ -d /usr/lib/rstudio ]; then
  echo "RStudio pandoc libraries found"
elif [ -d /usr/lib/rstudio-server ]; then
  echo "RStudio Server pandoc libraries found"
else
  if [ os == "linux-gnu"]; then
    echo "No pandoc binaries found on your system--attempting to download"
    wget https://github.com/jgm/pandoc/releases/download/1.19.2.1/pandoc-1.19.2.1-1-amd64.deb
    mkdir pandoc
    ar p pandoc-1.19.2.1-1-amd64.deb data.tar.gz | tar xvz --strip-components 2 -C ./pandoc/
  elif [ os == "darwin"]; then
    echo "Please install RStudio"
  else
    echo "Unrecognized operating system"
  fi
fi
