#!/bin/bash

# install python dependencies
# may need to update path to pip
# by default installs to user library location (default ~/.local/bin)
# remove --user if installing system-wide
echo "Installing python dependencies"
/usr/bin/pip install --user cyvcf2 numpy Biopython nimfa pyfaidx

# # Link rstudio pandoc libraries to path
# if [ -d /usr/lib/rstudio ]; then
#   echo "Linking RStudio pandoc libraries"
#   ln -s /usr/lib/rstudio/bin/pandoc/pandoc /usr/local/bin
#   ln -s /usr/lib/rstudio/bin/pandoc/pandoc-citeproc /usr/local/bin
# elif [ -d /usr/lib/rstudio-server ]; then
#   echo "Linking RStudio Server pandoc libraries"
#   ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc /usr/local/bin
#   ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc-citeproc /usr/local/bin
# else
#   # may need to manually update package manager if not using apt-get
#   echo "No RStudio pandoc binaries found--attempting to install standalone pandoc"
#   sudo apt-get install pandoc pandoc-citeproc
# fi
