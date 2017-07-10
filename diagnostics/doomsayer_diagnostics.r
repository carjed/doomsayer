Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")

# check and load packages
packages <- c("rmarkdown", "knitr")

package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
    }
})

# receive YAML config file name as argument
cmd_args <- commandArgs(TRUE)
yaml_cfg <- cmd_args[1]
proj_dir <- dirname(yaml_cfg)

# Copy the RMarkdown template to the specified output folder
file.copy("./diagnostics/diagnostics.Rmd",
          paste0(proj_dir, "/diagnostics.Rmd"),
          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

# Render the RMarkdown document as both an HTML page and PDF,
# these will be located in the specified output folder
rmarkdown::render(paste0(proj_dir, "/diagnostics.Rmd"),
  c("html_document", "pdf_document"),
  params=list(yaml_cfg=yaml_cfg))
