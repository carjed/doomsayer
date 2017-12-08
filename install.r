# install R dependencies from CRAN repository
require(devtools)
packages <- c("dplyr", "tidyr", "ggplot2", "yaml", "knitr", "viridis", "RColorBrewer", "heatmaply", "dendextend")

package.check <- suppressWarnings(suppressMessages(lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x,
          dependencies = TRUE,
          repos = "http://cran.us.r-project.org")
    }
})))
