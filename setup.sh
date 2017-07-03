#!/bin/bash

# install python dependencies
echo "Installing python dependencies"
# /usr/local/bin/pip install --user cyvcf2 numpy Biopython scikit-learn

# Link rstudio pandoc libraries to path

if [ -d /usr/lib/rstudio ]; then
  echo "Linking RStudio pandoc libraries"
  ln -s /usr/lib/rstudio/bin/pandoc/pandoc /usr/local/bin
  ln -s /usr/lib/rstudio/bin/pandoc/pandoc-citeproc /usr/local/bin
elif [ -d /usr/lib/rstudio-server ]; then
  echo "Linking RStudio Server pandoc libraries"
  ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc /usr/local/bin
  ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc-citeproc /usr/local/bin
else
  echo "No RStudio pandoc binaries found--attempting to install standalone pandoc"
  sudo apt-get install pandoc pandoc-citeproc
fi
