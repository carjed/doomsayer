## Diagnostic reports for Doomsayer

This directory contains R scripts for evaluating the results of a Doomsayer run where the `--diagnostics` or `--autodiagnostics` option has been enabled.

The `doomsayer_diagnostics.r` script will copy an [RMarkdown template](diagnostics.Rmd) into `~/doomsayer_output/`, which reads the files listed in `~/doomsayer_output/diagnostics.yaml` and renders diagnostic reports (in Markdown and HTML formats) detailing the results of your Doomsayer run. If the `--autodiagnostics` option is enabled, this script should run automatically; otherwise, you can manually generate a diagnostic report with the following command:

```{sh id:"chj4m1wcp3"}
cd /path/to/doomsayer
Rscript diagnostics/doomsayer_diagnostics.r ~/doomsayer_output/config.yaml
```

Unless you are running this script from within the `doomsayer_output` directory, be sure to specify an absolute path of `config.yaml` (e.g., `/path/to/doomsayer_output/config.yaml`), not a relative path (e.g.,  `./doomsayer_output/config.yaml`); otherwise, the config file will not be found.

An example of a final report is available [here](diagnostics.md).

All plots shown in this report will be available as standalone .png images, located in `~/doomsayer_output/diagnostics/diagnostics_files/figure-html/`

### Dependencies

These R scripts require the following packages, which will be installed automatically to your default library location:
- [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html)
- [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
