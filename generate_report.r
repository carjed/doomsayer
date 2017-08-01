rstdir <- "/usr/lib/rstudio/bin/pandoc"
rstsrvdir <- "/usr/lib/rstudio-server/bin/pandoc"
pandocdir <- "./pandoc/bin"

if(dir.exists(rstdir)){
  Sys.setenv(RSTUDIO_PANDOC=rstdir)
} else if(dir.exists(rstsrvdir)){
  Sys.setenv(RSTUDIO_PANDOC=rstsrvdir)
} else {
  Sys.setenv(RSTUDIO_PANDOC=pandocdir)
}

# check and load packages
packages <- c("rmarkdown", "knitr")

package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE, silent=TRUE)
    }
})

# receive YAML config file name as argument
cmd_args <- commandArgs(TRUE)
yaml_cfg <- cmd_args[1]
template <- cmd_args[2]
proj_dir <- dirname(yaml_cfg)

templateRmd <- paste0("./report_templates/", template, ".Rmd")
outputRmd <- paste0(proj_dir, "/report.Rmd")

# Copy the RMarkdown template to the specified output folder
file.copy(templateRmd, outputRmd,
  overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

# Render the RMarkdown document as an HTML page, located in the specified
# output folder
# c("html_document", "pdf_document"),
setwd(proj_dir)
rmarkdown::render(outputRmd,
  c("html_document"),
  params=list(yaml_cfg=yaml_cfg),
  quiet=TRUE)
