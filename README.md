<!-- # DOOMSAYER -->

![](assets/doomsayer_logo.png)

## Introduction

_**DOOMSAYER** ( **D**etection **O**f **O**utliers using **M**utation **S**pectrum **A**nal**Y**sis in **E**xtremely **R**are variants) is a utility for analyzing patterns of genetic variation in next-generation sequencing (NGS) data._

_Using non-negative matrix factorization (NMF), Doomsayer summarize the distribution of rare single-nucleotide variants within each sample. In studies of standing genetic variation, this information can be a valuable indicator of cryptic error biases or batch effects. By identifying and filtering these outlier samples, Doomsayer helps ensure rigor and reproducibility in the analysis of NGS data._

_In addition to its purpose as a quality control program, Doomsayer can be applied more generally to study between-sample differences in somatic and germline mutation signatures._

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
- [cyvcf2](https://github.com/brentp/cyvcf2) for VCF parsing
- [pyfaidx](https://github.com/mdshw5/pyfaidx) for reference genome parsing
- [numpy](https://docs.scipy.org/doc/numpy/reference/) for numeric arrays
- [Biopython](https://github.com/biopython/biopython) for various sequence parsing functions
- [nimfa](http://nimfa.biolab.si) for NMF algorithms

## Usage

Running `python doomsayer.py --help` will return the following usage information.

```
usage: doomsayer.py [-h] -i [INPUT] [-f FASTAFILE] [-p PROJECTDIR] [-o] [-n]
                    [-d] [-a] [-m [MMATRIXNAME]] [-s [SAMPLEFILE]]
                    [-g [GROUPFILE]] [-ns] [-r {2,3,4,5,6,7,8,9,10}]
                    [-l {1,3,5,7}] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        input file, usually a VCF. Can accept input from STDIN
                        with "--input -". If using in aggregation mode, input
                        should be a text file containing the file paths of the
                        M matrices to aggregate
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
                        the input VCF. (Useful if analyzing somatic data or
                        pre-filtering with another tool)
  -d, --diagnostics     write the NMF matrices to the output directory and
                        generate yaml config to be passed to diagnostic script
  -a, --autodiagnostics
                        same as the --diagnostics option, but automatically
                        generates diagnostic report
  -m [MMATRIXNAME], --mmatrixname [MMATRIXNAME]
                        custom filename for M matrix [without extension]
  -s [SAMPLEFILE], --samplefile [SAMPLEFILE]
                        file with sample IDs to include (one per line)
  -g [GROUPFILE], --groupfile [GROUPFILE]
                        two-column tab-delimited file containing sample IDs
                        (column 1) and group membership (column 2) for pooled
                        analysis
  -ns, --noscale        do not scale H and W matrices
  -r {2,3,4,5,6,7,8,9,10}, --rank {2,3,4,5,6,7,8,9,10}
                        rank for NMF decomposition
  -l {1,3,5,7}, --length {1,3,5,7}
                        motif length
  -v, --verbose         Enable verbose logging
```

## Tutorial

### Basic example
The basic purpose of Doomsayer is to provide sample-level filtering by using non-negative matrix factorization to identify individuals with abnormal distributions of singleton variants.

A basic, fully-functional Doomsayer command requires two mandatory arguments: 1) a [Variant Call Format (VCF) input file](#vcf-input), and 2) a [fasta-format reference genome](#fasta-reference-file):

```{bash id:"chj4lkgfvz"}
python /path/to/doomsayer.py \
  --input /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta
```

The default output consists of two plain-text files: `doomsayer_keep.txt`, containing the IDs of samples that passed the filter, and `doomsayer_drop.txt`, containing the IDs of samples that failed the filter. These lists can be easily integrated into bcftools, vcftools, PLINK, and other downstream analysis tools.

By default, output files are stored in a subfolder named `doomsayer_output/` in your current directory.

### Input options

#### Input file
`-i [INPUT], --input [INPUT]` **(MANDATORY)**

The basic input should be a VCF file formatted with mandatory columns (#CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO). This VCF must contain an allele count field in the INFO column ([AC=N]), and include individual genotypes (i.e., the VCF cannot be a "sites-only" file). Both uncompressed and [bgzip](http://www.htslib.org/doc/tabix.html)-compressed VCF formats are supported.

The VCF input can either be read from a file on disk or piped from STDIN (using `--input -`), enabling compatibility with a wide range of existing pipelines and workflows. In most cases, no extra preprocessing of the input VCF is necessary. The VCF can contain any combination of variant types (e.g., SNVs, indels, CNVs), allele counts, or filter flags. By default, Doomsayer will only analyze singleton SNVs with a "PASS" value in the FILTER column of the VCF.

If your data are spread across multiple files (e.g., one VCF per chromosome), you will need to combine them into a single VCF file using [bcftools](https://samtools.github.io/bcftools/) or a similar utility, or use Doomsayer in [aggregation mode](#). The following command demonstrates how to concatenate a set of VCFs with bcftools and pipe the output to Doomsayer:

```{sh id:"chj4lkh8gx"}
bcftools concat /path/to/input/vcfs/chr*.vcf.gz | \
  python /path/to/doomsayer.py \
    --input - \
    --fastafile /path/to/genome.fasta
```

#### Fasta reference file
`-f FASTAFILE, --fastafile FASTAFILE` **(MANDATORY)**

This parameter specifies the fasta-formatted reference genome used to look up the local sequence context for each SNV. This file **must** corresponding to the same reference genome build used to call the variants contained in the input VCF file, with [consistently-formatted records](#). The fasta file may be either uncompressed or bgzip-compressed.

#### Sample subsets
`-s [SAMPLEFILE], --samplefile [SAMPLEFILE]` **(optional)**

If you only wish to run Doomsayer on a subset of samples in the input VCF, this parameter will read a list of sample IDs to keep (one per line) and skip all other samples in the VCF. This will be much faster than pre-filtering the VCF with another tool (e.g., bcftools) and piping the input to Doomsayer.

#### Sample groupings
`-g [GROUPFILE], --groupfile [GROUPFILE]` **(optional)**

This parameter forces Doomsayer to evaluate mutation spectra across pooled groups of samples. The GROUPFILE should be a tab-delimited text file containing sample IDs in the first column and a grouping variable in the second column. This option is particularly useful if you wish to explicitly filter for batch effects in your data, e.g., by setting the grouping variable to be the sequencing date, study of origin, sequencing plate, etc.

Note that this option will assume the GROUPFILE contains all samples in the input VCF. If your GROUPFILE contains only a subset of samples, other samples will not be evaluated.

### Output options

#### Custom output directory
`-p PROJECTDIR, --projectdir PROJECTDIR` **(optional)**

This parameter will specify a custom name/location for the Doomsayer output directory. The following example will send output to a folder named `my_doomsayer_output` in your home directory:

```{sh id:"chj4lkfg8u"}
python /path/to/doomsayer.py \
  --input /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta \
  --projectdir ~/my_doomsayer_output
```

#### Auto-filtered VCF output
`-o, --outputtovcf` **(optional)**

This option will allow Doomsayer to automatically filter the input VCF and write a new VCF, removing any samples specified in the `doomsayer_drop.txt` file that has been generated.

```{sh id:"chj4lkcw3i"}
python /path/to/doomsayer.py \
  --input /path/to/input_data.vcf.gz \
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
`-n, --nofilter` **(optional)**

This option will force Doomsayer to include **all** biallelic SNVs from the input VCF file in the NMF analysis (even those that fail a filter).

This is useful in instances where interindividual variation is expected, and we wish to use NMF simply to characterize this heterogeneity (e.g., [comparing somatic mutation signatures in tumor samples](http://cancer.sanger.ac.uk/cosmic/signatures)), or where we have used another tool to pre-filter the set of variants to analyze. This option should usually be used in conjunction with the [`--diagnostics`](#) flag.

#### Detailed NMF output and diagnostic report
`-d, --diagnostics` **(optional)**

This option will write three data files from the NMF decomposition to the output directory:
1. the M matrix of observed mutation spectra per sample
2. the H matrix containing loadings of the mutation subtypes in each signature
3. the W matrix containing the contributions of each signature within each sample's mutation spectrum

```{sh id:"chj4lkjmny"}
python doomsayer.py \
  --input /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta \
  --projectdir ~/my_doomsayer_output \
  --diagnostics
```

The `--diagnostics` flag will also create a YAML configuration file that can optionally be passed to an R script for generating more detailed diagnostic reports from the NMF data files. See the [diagnostics](diagnostics) directory for further documentation and usage examples.

Running with the `--autodiagnostics` flag instead will force Doomsayer attempt to automatically generate the reports.

#### Write only M matrix with custom name
`-m [MMATRIXNAME], --mmatrixname [MMATRIXNAME]` **(optional)**

This parameter allows you to specify a custom name for the mutation spectra matrix file. This option was included for cases where large datasets are spread across several files and we wish to take advantage of the [aggregation](#) functionality of Doomsayer.  Note that using this option will force Doomsayer to skip the NMF decomposition and generation of keep/drop lists.


### Other optional settings
-------------


#### Specify NMF rank
`-r {2,3,4,5,6,7,8,9,10}, --rank {2,3,4,5,6,7,8,9,10}` **(optional)**

This parameter specifies the rank of the NMF decomposition (up to 10). The default rank (3) should be sufficient for most QC applications, but if you are using Doomsayer for somatic mutation signature analysis, you will likely require a higher rank.

#### Set motif length
`-l {1,3,5,7}, --length {1,3,5,7}` **(optional)**

This parameter specifies the motif length to be considered in determining the mutation subtypes. The default (3) produces 96 3-mer subtypes, which is likely sufficient for QC of whole-genome germline variants. If you are using Doomsayer for QC of whole-exome germline variants, you will likely need to specify `--length 1`, as most individuals will not have enough singleton variants to accurately infer higher-order mutation signatures. Values greater than 3 will likely lead to extremely noisy or unreliable QC results for germline variants, but may be useful for examining somatic sequencing data.

#### Verbose logging
`-v, --verbose` **(optional)**

This option enables detailed logging to the STDERR stream

### Aggregation mode
Doomsayer can be run independently and simultaneously on subsets of the same dataset and then act as an aggregator for these outputs. This is particularly necessary for very large datasets with millions of variants and thousands of samples where sequential processing of a single massive VCF file is not feasible.

The `doomsayer_by_chr.sh` script provides an example of how to implement this function. This shell script will simultaneously start 22 Doomsayer runs in the background (one per chromosome). Each run is specified with the `--mmatrixname chrN` argument, so only the mutation spectra matrices for each chromsomes are written, each with a unique file name.

This script also writes a file named `m_regions.txt` containing the file names/paths of the per-chromosome mutation spectra matrices.

Once the per-chromosome runs are complete, `m_regions.txt` can then be passed to Doomsayer with the `--input m_regions.txt` option--when Doomsayer detects that the input file is not a VCF, it will enter a routine to aggregate the mutation spectra matrices listed in `m_regions.txt` into a single M matrix, and then proceed with the NMF decomposition (and optional diagnostics) as normal.

<!-- If your data is split into N genomic regions (or if different files contain genotypes for N nonoverlapping subsets of individuals), you can run Doomsayer on each subset individually with a unique name for the M matrix, producing N unique matrix outputs.

 into Doomsayer may not always be practical, so we enabled Doomsayer to act as an aggregator for output from smaller datasets. -->

## Demonstration

Sample VCF and reference files can be downloaded using the [download_demo_data.sh](download_demo_data.sh) script. Files will be placed in the `demo/input` subdirectory of the doomsayer install directory. You can also manually download the necessary files from [http://mutation.sph.umich.edu/share/doomsayer_demo/](http://mutation.sph.umich.edu/share/doomsayer_demo/).

The following command will run the demo, writing all output to `demo/output`:

```{sh id:"chj4ll4db1"}
cd /path/to/doomsayer
bash download_demo_data.sh

python doomsayer.py \
  --input demo/input/chr20.singletons.sample.vcf \
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

## FAQ

### Is Doomsayer compatible with my system?
Doomsayer was written for *nix-based systems.

The main Doomsayer program is cross-compatible with Python versions 2.7+ and 3.5+. Generation of diagnostic reports requires R version 3.1 or higher.

### What are Doomsayer's hardware requirements?
Doomsayer is extremely lightweight, and should run with \<1GB of RAM

### What type of sequencing data can Doomsayer analyze?
Doomsayer's core functionality, the NMF-based mutation signature analysis, can be used to rapidly evaluate mutation patterns present in any collection of sequencing data, and offers an extremely efficient and almost entirely automated (albeit slightly less flexible) alternative to the functions provided by the [SomaticSignatures](http://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html) R package:

>Gehring JS, Fischer B, Lawrence M, Huber W. SomaticSignatures: inferring mutational signatures from single-nucleotide variants. Bioinformatics. 2015;31(22):3673-3675. doi:10.1093/bioinformatics/btv408.

Although we have only tested Doomsayer with human data, it should be fully capable of analyzing NGS data from any other organism.

Note that the sample-level filtering function of Doomsayer is exclusively designed for analyzing **population** variants in **unrelated** individuals. In its current state, this filtering is **not** appropriate for quality control in somatic sequencing data or *de novo* mutation analysis.

Also, Doomsayer does not currently include multiallelic variants in its analysis. This feature is under development.

### Doomsayer is not reading my input VCF
If you encounter any issues with reading in a VCF file, first confirm that it is properly formatted. [vcf-validator](https://github.com/EBIvariation/vcf-validator) is a useful utility for this. You can also try indexing the VCF with `tabix -p input.vcf.gz` and running a simple bcftools command, such as `bcftools view input.vcf.gz`. Also check that your VCF contains singleton SNVs, individual genotypes are included, and that the FILTER column contains "PASS" values.

### Reference file incompatibilities
Errors may occur if the chromosome records in your fasta reference file are formatted differently than the CHROM field of your VCF (e.g., if the fasta header is formatted as `>chrN`). You will need to modify the records in one or the other to be consistent.

### RMarkdown reports are not rendering
Rendering RMarkdown reports directly from the command line requires the pandoc binaries be accessible from your system path. If you have sudo access, running `sudo setup.sh` will look for an installation of RStudio or RStudio Server, and attempt to soft-link the included pandoc binaries from `/usr/lib/rstudio/bin/` or `/usr/lib/rstudio-server/bin/` to `/usr/local/bin`. Failing this, the setup will attempt to install compatible standalone pandoc binaries: `apt-get install pandoc pandoc-citeproc`. If you do not have sudo access and RStudio is installed on your server, you can try adding the RStudio paths to the $PATH variable in `~/.bashrc`.

Additional information for enabling command-line RMarkdown/pandoc rendering can be found [here](https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md).

## Contact
If you need any further help, come across a bug, or have a feature request, feel free to [open an issue](https://github.com/carjed/doomsayer/issues/new) or [contact me directly](mailto:jed.e.carlson@gmail.com). If you are reporting a bug or unexpected output, please run Doomsayer with the `--verbose` option and include this output in your inquiry.
