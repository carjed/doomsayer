<!-- # DOOMSAYER -->

![](assets/doomsayer_logo.png)

## Introduction

_**DOOMSAYER** ( **D**etection **O**f **O**utliers using **M**utation **S**pectrum **A**nal**Y**sis in **E**xtremely **R**are variants) is a utility for analyzing patterns of genetic variation in next-generation sequencing (NGS) data._

_Doomsayer uses non-negative matrix factorization (NMF) to summarize the distribution of rare single-nucleotide variants within each sample. In studies of standing genetic variation, this information can be a valuable indicator of cryptic error biases or batch effects. By identifying and filtering these outlier samples, Doomsayer helps ensure rigor and reproducibility in the analysis of NGS data._

_In addition to its function as a quality control program, Doomsayer can be applied more generally to study between-sample differences in somatic and germline mutation signatures._

_Doomsayer is currently under active development, and should be considered as beta software._

<!-- ## Citation

If you use DOOMSAYER in your research, please cite our [paper](#) (citation pending). -->

## Contents
- [Download and installation](#download-and-installation)
- [Tutorial](#tutorial)
  - [Usage options](#usage-options)
  - [Input](#input)
    - [VCF input](#vcf-input)
    - [Fasta reference file](#fasta-reference-file)
  - [Output](#output)
  - [Optional output settings](#optional-output-settings)
    - [Custom output directory](#custom-output-directory)
    - [Generate auto-filtered VCF](#auto-filtered-vcf)
    - [Run NMF on all variants in the VCF](#run-nmf-on-all-variants-in-the-vcf)
    - [NMF output](#nmf-output)
    - [Diagnostic reports](#diagnostic-reports)
  - [Other options](#other-options)
- [Demonstration](#demonstration)
- [Compatibility and troubleshooting](#compatibility-and-troubleshooting)
- [Contact](#contact)

## Download and installation
To begin using Doomsayer, clone this repository with the following command:

```{sh id:"chj4lmgsky"}
git clone https://github.com/carjed/doomsayer.git
```

Enter the directory and run the `setup.sh` script, which will install dependencies and perform some additional setup tasks:

```{sh id:"chj4lmgbgi"}
cd doomsayer
bash setup.sh
```

`setup.sh` uses pip to check for and install the following python libraries and perform additional setup tasks:
- [CyVCF2](https://github.com/brentp/cyvcf2) (or [PyVCF](https://github.com/jamescasbon/PyVCF))
- [numpy](https://docs.scipy.org/doc/numpy/reference/)
- [Biopython](https://github.com/biopython/biopython)
- [scikit-learn](http://scikit-learn.org/stable/)

## Tutorial

Doomsayer requires two mandatory user-specified arguments: 1) a [Variant Call Format (VCF) input file](#vcf-input), and 2) a [fasta-format reference genome](#fasta-reference-file). Below is a basic, fully-functional Doomsayer command:

```{bash id:"chj4lkgfvz"}
python /path/to/doomsayer.py \
  --inputvcf /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta
```

### Usage options

```
usage: doomsayer.py [-h] -i [INPUTVCF] -f FASTAFILE [-p PROJECTDIR] [-o] [-n]
                    [-d] [-ad] [-r {2,3,4,5,6,7,8,9,10}] [-l {1,3,5,7}] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUTVCF], --inputvcf [INPUTVCF]
                        input vcf file (use - for STDIN)
  -f FASTAFILE, --fastafile FASTAFILE
                        fasta file name
  -p PROJECTDIR, --projectdir PROJECTDIR
                        directory to store output files (do NOT include a
                        trailing '/')
  -o, --outputtovcf     filter input VCF (writes to stdout--use standard
                        output redirection [ > out.vcf] to write the output
                        VCF to a file on disk)
  -n, --nofilter        disables generation of keep/drop lists, turns off
                        default filtering criteria, and evaluates all sites in
                        the input VCF. Useful if analyzing somatic data or
                        pre-filtering with another tool)
  -d, --diagnostics     write the NMF matrices to the output directory and
                        generate yaml config to be passed to diagnostic script
  -ad, --autodiagnostics
                        same as the --diagnostics option, but automatically
                        generates diagnostic report
  -r {2,3,4,5,6,7,8,9,10}, --rank {2,3,4,5,6,7,8,9,10}
                        rank for NMF decomposition
  -l {1,3,5,7}, --length {1,3,5,7}
                        motif length
  -v, --verbose         Enable verbose logging
```

### Input
-------------
#### VCF input
We assume the input VCF file is formatted with mandatory columns (#CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO), contains an allele count field in the INFO column ([AC=N]), and includes individual genotypes (i.e., cannot be a "sites-only" file). Both uncompressed and [bgzip](http://www.htslib.org/doc/tabix.html)-compressed VCF formats are supported.

The VCF input can either be read from a file on disk or piped from STDIN (using `--inputfile -`), enabling compatibility with a wide range of existing pipelines and workflows. In most cases, no extra preprocessing of the input VCF is necessary. The VCF can contain any combination of variant types (e.g., SNVs, indels, CNVs), allele counts, or filter flags. Doomsayer will only analyze singleton SNVs with a "PASS" value in the FILTER column of the VCF.

If your data are spread across multiple files (e.g., one VCF per chromosome), you will need to combine them into a single VCF file using [bcftools](https://samtools.github.io/bcftools/) or a similar utility. The following command demonstrates how to concatenate a set of VCFs with bcftools and pipe the output to Doomsayer:

```{sh id:"chj4lkh8gx"}
bcftools concatenate /path/to/input/vcfs/chr*.vcf.gz | \
  python /path/to/doomsayer.py \
    --inputvcf - \
    --fastafile /path/to/genome.fasta
```

#### Fasta reference file
Doomsayer also requires a fasta sequence file be specified with the `--fastafile` argument. This **must** corresponding to the same reference genome build used to call the variants contained in the input VCF file, with [consistently-formatted records](#reference-file-incompatibilities). The fasta file may be either uncompressed or bgzip-compressed.

### Output
-------------

The basic function of Doomsayer is to provide sample-level filtering by using NMF to identify individuals with abnormal distributions of singleton variants. Thus, the default output consists of two plain-text files: `keep_ids.txt`, containing the IDs of samples that passed the filter, and `drop_ids.txt`, containing the IDs of samples that failed the filter. These lists can be easily integrated into bcftools, vcftools, PLINK, and other downstream analysis tools.

By default, output files are stored in a subfolder named `doomsayer_output/` in your current directory.

### Optional output settings

#### Custom output directory
The name/location of the Doomsayer output directory can be customized with the `--projectdir` argument. The following example will send output to a folder named `my_doomsayer_output` in your home directory:

```{sh id:"chj4lkfg8u"}
python /path/to/doomsayer.py \
  --inputvcf /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta \
  --projectdir ~/my_doomsayer_output
```

Be sure to exclude the trailing "/" from any directory specified.

#### Auto-filtered VCF
With the `--outputtovcf` option, Doomsayer can automatically filter the input VCF and write a new VCF, removing all samples specified in the `drop_ids.txt` file.

```{sh id:"chj4lkcw3i"}
python /path/to/doomsayer.py \
  --inputvcf /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta \
  --projectdir ~/my_doomsayer_output \
  --outputtovcf
```

This command will print the VCF header and records to STDOUT, which can then be piped to compatible programs, or redirected to a file using standard UNIX redirection (adding `> /path/to/output.vcf` at the end of the above example), allowing Doomsayer to be used as an intermediate step of more complex VCF processing pipelines.

As far as possible, all records present in the original input will be preserved in the output, with three notable exceptions:
1. All samples in the drop list (and their genotypes) will be omitted
2.  **Variants that were unique to the set of dropped samples will be automatically removed.**
3. Four of the INFO fields will be updated to reflect the smaller set of retained samples:
    1. The total number of alleles (AN)
    2. allele count in genotypes (AC)
    3. number of samples (NS)
    4. combined depth (DP)

#### Run NMF on all variants in the VCF
The `--nofilter` option will suppress creation of the keep and drop lists and include **all** biallelic SNVs from the input VCF file in the NMF analysis (even those that fail a filter).

This option is useful in instances where interindividual variation is expected, and we wish to use NMF simply to characterize this heterogeneity (e.g., [comparing somatic mutation signatures in tumor samples](http://cancer.sanger.ac.uk/cosmic/signatures)). This option should usually be used in conjunction with the [`--diagnostics`](#nmf-output) flag.

#### NMF output
if you wish to analyze the results of the NMF decomposition, you must run Doomsayer with the `--diagnostics` option enabled. This will write three data files from the NMF decomposition to the output directory:
1. the M matrix of observed mutation spectra per sample
2. the H matrix containing loadings of the mutation subtypes in each signature
3. the W matrix containing the contributions of each signature within each sample's mutation spectrum

```{sh id:"chj4lkjmny"}
python doomsayer.py \
  --inputvcf /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta \
  --projectdir ~/my_doomsayer_output \
  --diagnostics
```

#### Diagnostic reports
The `--diagnostics` flag will also create a YAML configuration file that can optionally be passed to an R script for generating more detailed diagnostic reports from the NMF data files. See the [diagnostics](diagnostics) directory for further documentation and usage examples.

You can also change this to `--autodiagnostics` to have Doomsayer attempt to automatically generate the reports.

### Other options
-------------

- `-r {2,3,4,5,6,7,8,9}, --rank {2,3,4,5,6,7,8,9}` This option specifies the rank of the NMF decomposition (up to 10). The default rank (3) should be sufficient for most QC applications, but somatic signature analysis will likely require a higher rank.


- `-l {1,3,5,7}, --length {1,3,5,7}` This option specifies the motif length to be considered in determining the mutation subtypes. The default (3) produces 96 3-mer subtypes, which is standard for NMF-based mutation signature analysis. Higher values are currently experimental, and will likely lead to extremely noisy or unreliable results. Setting `-l 1` will improve speed, but reduces the resolution of the NMF input matrix, and may not differentiate outliers as accurately.


- `-v, --verbose` This enables verbose logging to the STDERR stream


- `-h --help` Access a brief description of all available options and exit:


## Demonstration

Sample VCF and reference files can be downloaded using the [download_demo_data.sh](download_demo_data.sh) script. Files will be placed in the `demo/input` subdirectory of the doomsayer install directory. You can also manually download the necessary files from [http://mutation.sph.umich.edu/share/doomsayer_demo/](http://mutation.sph.umich.edu/share/doomsayer_demo/).

The following command will run the demo, writing all output to `demo/output`:

```{sh id:"chj4ll4db1"}
cd /path/to/doomsayer
bash download_demo_data.sh

python doomsayer.py \
  --inputvcf demo/input/chr20.singletons.sample.vcf \
  --fastafile demo/input/chr20.fasta.gz \
  --projectdir demo/output \
  --diagnostics \
  --verbose \
  --outputtovcf > demo/output/chr20.filtered.vcf
```

Once complete, you can run the diagnostic script to generate diagnostic reports:
```{sh id:"chj4lm2fcv"}
Rscript diagnostics/doomsayer_diagnostics.r demo/output/config.yaml
```

An example of a final report is available [here](diagnostics/diagnostics.md).

## Compatibility and troubleshooting
Doomsayer was written for *nix-based systems.

The main Doomsayer program is cross-compatible with Python versions 2.7+ and 3.5+. Generation of diagnostic reports requires R version 3.1 or higher.

### What type of sequencing data can Doomsayer analyze?
Doomsayer's core functionality can be used to rapidly evaluate mutation signatures present in any collection of sequencing data with called variants, and offers an extremely efficient and almost entirely automated (albeit slightly less flexible) alternative to the functions provided by the [SomaticSignatures](http://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html) R package:

>Gehring JS, Fischer B, Lawrence M, Huber W. SomaticSignatures: inferring mutational signatures from single-nucleotide variants. Bioinformatics. 2015;31(22):3673-3675. doi:10.1093/bioinformatics/btv408.

Although we have only tested Doomsayer with human data, it should be fully capable of analyzing NGS data from any other organism.

Note that the sample-level filtering function of Doomsayer is exclusively designed for analyzing **population** variants in **unrelated** individuals. In its current state, this filtering is **not** appropriate for quality control in somatic sequencing data or *de novo* mutation analysis.

Also, Doomsayer does not currently include multiallelic variants in its analysis. This feature is under development.

### Doomsayer is not reading my input VCF
If you encounter any issues with reading in a VCF file, first confirm that it is properly formatted. [vcf-validator](https://github.com/EBIvariation/vcf-validator) is a useful utility for this. You can also try indexing the VCF with `tabix -p input.vcf.gz` and running a simple bcftools command, such as `bcftools view input.vcf.gz`. Also check that your VCF contains singleton SNVs, individual genotypes are included, and that the FILTER column contains "PASS" values.

### Reference file incompatibilities
Errors may occur if the chromosome records in your fasta reference file are formatted differently than the CHROM field of your VCF (e.g., if the fasta header is formatted as `>chrN`). You will need to modify the records in one or the other to be consistent.

### RMarkdown reports are not rendering
Rendering RMarkdown reports directly from the command line requires the pandoc binaries. If you have sudo access, running `sudo setup.sh` will look for an installation of RStudio or RStudio Server, and attempt to soft-link the included pandoc binaries from `/usr/lib/rstudio/bin/` or `/usr/lib/rstudio-server/bin/` to `/usr/local/bin`. Failing this, the setup will attempt to install compatible standalone pandoc binaries: `apt-get install pandoc pandoc-citeproc`. If you do not have sudo access and RStudio is installed on your server, you can try adding the RStudio paths to the $PATH variable in `~/.bashrc`.

Additional information for enabling command-line RMarkdown/pandoc rendering can be found [here](https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md).

## Contact
If you need any further help, come across a bug, or have a feature request, feel free to [open an issue](https://github.com/carjed/doomsayer/issues/new) or [contact me directly](mailto:jed.e.carlson@gmail.com). If you are reporting a bug or unexpected output, please run Doomsayer with the `--verbose` option and include this output in your inquiry.
