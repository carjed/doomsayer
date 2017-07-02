## Diagnostic reports for Doomsayer

This directory contains R scripts for evaluating the results of a Doomsayer run where the `--diagnostics` option has been enabled. Assuming your output is contained in `~/doomsayer_output/` you can generate a diagnostic report with the following command:

```{sh id:"chj4m1wcp3"}
cd /path/to/doomsayer
Rscript diagnostics/doomsayer_diagnostics.r ~/doomsayer_output/config.yaml
```

Unless you are running this script from within the `doomsayer_output` directory, be sure to specify an absolute path of `config.yaml` (e.g., `/path/to/doomsayer_output/config.yaml`), not a relative path (e.g.,  `./doomsayer_output/config.yaml`); otherwise, the config file will not be found.

This script will copy an [RMarkdown template](diagnostics.Rmd) into `~/doomsayer_output/`, which reads the files listed in `~/doomsayer_output/diagnostics.yaml` and renders diagnostic reports (in Markdown, PDF, and HTML formats) detailing the results of your Doomsayer run.

An example of a final report is available [here](diagnostics.md).

All plots shown in this report will be available as standalone .png images, located in `~/doomsayer_output/diagnostics/diagnostics_files/figure-html/`

### Dependencies

These R scripts require the following packages:
- [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html)
- [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [devtools](https://cran.r-project.org/web/packages/devtools/index.html)

The scripts will check if they are installed, and attempt to install any that are not found. Packages will be installed to your default library location.
