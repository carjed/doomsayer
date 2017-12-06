# install R dependencies from CRAN repository
require(devtools)
packages <- c("dplyr", "tidyr", "ggplot2", "yaml", "knitr", "corrr", "viridis", "RColorBrewer", "heatmaply", "dendextend", "e1071")

package.check <- suppressWarnings(suppressMessages(lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x,
          dependencies = TRUE,
          repos = "http://cran.us.r-project.org")
        library(x, character.only = TRUE, verbose=FALSE)
    }
})))
