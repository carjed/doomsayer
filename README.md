<!-- # DOOMSAYER -->

![](assets/doomsayer_logo.png)

## Introduction

_**DOOMSAYER** ( **D**etection **O**f **O**utliers using **M**utation **S**pectrum **A**nal**Y**sis in **E**xtremely **R**are variants) is a quality control application for next-generation sequencing (NGS) data. Doomsayer uses non-negative matrix factorization (NMF) to identify samples with abnormal distributions of rare single-nucleotide variants. These outlying samples are likely a result of cryptic error biases or batch effects, which can negatively impact any number of downstream analyses. By identifying and filtering these problem samples, Doomsayer helps ensure rigor and reproducibility in the analysis of NGS data._

*This program is currently under active development*

<!-- ## Citation

If you use DOOMSAYER in your research, please cite our [paper](#) (citation pending). -->

## Contents

- [Download and Installation](#download-and-installation)
- [Tutorial](#tutorial)
  - [Input options](#input-options)
    - [VCF](#vcf)
    - [Fasta reference file](#fasta-reference-file)
  - [Output options](#output-options)
    - [Keep and drop files](#keep-and-drop-files)
    - [Auto-filtered VCF (optional)](#auto-filtered-vcf-optional)
    - [NMF output (optional)](#nmf-output-optional)
      - [Generating diagnostic reports (optional)](#generating-diagnostic-reports-optional)
  - [Other Options](#other-options)
- [Demonstration](#demonstration)
- [Contact](#contact)

## Download and Installation
To begin using Doomsayer, clone this repository with the following command:

```{sh id:"chj4lmgsky"}
git clone https://github.com/carjed/doomsayer.git
```

You will then need to enter the directory and run the `setup.py` script to install any missing dependencies:

```{sh id:"chj4lmgbgi"}
cd doomsayer
python setup.py install
```

<!-- (`pip install pyvcf cyvcf numpy Biopython`) -->
<!-- You can also install Doomsayer with `pip install doomsayer` -->

`setup.py` installs the following python libraries:
- [CyVCF2](https://github.com/brentp/cyvcf2) (or [PyVCF](https://github.com/jamescasbon/PyVCF))
- [numpy](https://docs.scipy.org/doc/numpy/reference/)
- [Biopython](https://github.com/biopython/biopython)
- [scikit-learn](http://scikit-learn.org/stable/)

Doomsayer is compatible with both Python 2.7+/Python 3.5+

## Tutorial

Doomsayer requires two mandatory user-specified arguments: 1) a [Variant Call Format](https://en.m.wikipedia.org/wiki/Variant_Call_Format) (VCF) input file, and 2) a fasta-format reference genome. A basic command line might look like this:

```{bash id:"chj4lkgfvz"}
python /path/to/doomsayer.py \
  --inputvcf /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta
```

### Input options
-------------
#### VCF
We assume the input VCF file is formatted with mandatory columns (CHR, POS, REF, ALT, QUAL, FILTER), and contains an allele count field in the INFO column ([AC=N]). The input VCF **must** contain individual genotypes (i.e., cannot be a "sites-only" file). Both uncompressed and [bgzip](http://www.htslib.org/doc/tabix.html)-compressed VCF formats are supported.

The VCF input can either be read from a file on disk or piped from STDIN (using `--inputfile -`), enabling compatibility with a wide range of existing pipelines and workflows. In most cases, no extra preprocessing of the input VCF is necessary. The VCF can contain any combination of variant types (e.g., SNVs, indels, CNVs), allele counts, or filter flags. Doomsayer will only analyze singleton SNVs with a "PASS" value in the FILTER column of the VCF.

If your data are spread across multiple files (e.g., one VCF per chromosome), you will need to combine them into a single VCF file using `bcftools` or a similar utility. The following command demonstrates how to concatenate a set of VCFs with bcftools and pipe the output to doomsayer:

```{sh id:"chj4lkh8gx"}
bcftools concatenate /path/to/input/vcfs/chr*.vcf.gz | \
  python /path/to/doomsayer.py \
    --inputvcf - \
    --fastafile /path/to/genome.fasta
```

#### Fasta reference file
Doomsayer also requires a fasta sequence file be specified with the `--fastafile` option. This **must** corresponding to the same reference genome build used in the input VCF. The fasta file may be either uncompressed or bgzip-compressed.

### Output options
-------------
By default, doomsayer creates a subfolder named `doomsayer_output/` in your current directory. To change the name/location of this output directory, use the `--projectdir` command line option. The following example will send output to a folder named `my_doomsayer_output` in your home directory:

```{sh id:"chj4lkfg8u"}
python /path/to/doomsayer.py \
  --inputvcf /path/to/input_data.vcf.gz \
  --fastafile /path/to/genome.fasta \
  --projectdir ~/my_doomsayer_output
```

Be sure to remove the trailing "/" from any directory specified.

#### Keep and drop files
The basic function of Doomsayer is to provide sample-level filtering by identifying individuals with abnormal distributions of singleton variants. The results of this analysis are parsed into two plain-text files: `doomsayer_output/keep_ids.txt`, containing the IDs of samples that passed the filter, and `doomsayer_output/drop_ids.txt`, containing the IDs of samples that failed the filter. These lists can be easily integrated into bcftools, vcftools, PLINK, and other downstream analysis tools.

#### Auto-filtered VCF (optional)
With the `--outputtovcf` flag, Doomsayer can automatically filter the input VCF and write a new VCF, removing all samples specified in the `doomsayer_output/drop_ids.txt` file.


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

#### NMF output (optional)
You may wish to analyze the intermediate data used in the NMF decomposition to understand why particular samples were determined to be outliers, or to derive more customized filtering rules.

Adding the `--diagnostics` flag to your command will write three data files to the output directory:
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

##### Generating diagnostic reports (optional)
The `--diagnostics` flag will also create a YAML configuration file that can optionally be passed to an R script for generating more detailed diagnostic reports from the NMF data files. See the [diagnostics](diagnostics) directory for further documentation and usage examples.


### Other Options
-------------
<!-- The PyVCF implementation is included for legacy and testing purposes. -->
- `-r {2,3,4,5,6,7,8,9}, --rank {2,3,4,5,6,7,8,9}` This option specifies the rank of the NMF decomposition (up to 10). The default rank (3) should be sufficient for most analyses


- `-l {1,3,5,7}, --length {1,3,5,7}` This option specifies the motif length to be considered in determining the mutation subtypes. The default (3) produces 96 3-mer subtypes, which is standard for NMF-based mutation signature analysis. Higher values are currently experimental, and will likely lead to extremely noisy or unreliable results. Setting `-l 1` will improve speed, but reduces the resolution of the NMF input matrix, and may not differentiate outliers as accurately.


- `-z`, `--nonoptimized` This option forces Doomsayer to use the PyVCF library for VCF parsing, instead of the default CyVCF2 library. CyVCF2 is a python wrapper for the [htslib](https://github.com/samtools/htslib) C++ libraries, and is recommended for most use cases due to the large speed improvements. Note that the PyVCF implementation currently only supports uncompressed VCF input files, and does not support writing to a filtered output VCF.


- `-v, --verbose` This enables verbose logging to the STDERR stream


- `-h --help` Access a brief description of all available options and exit:
```
usage: doomsayer.py [-h] -i [INPUTVCF] [-f FASTAFILE] [-p PROJECTDIR] [-o]
                    [-n] [-r {2,3,4,5,6,7,8,9,10}] [-l {1,3,5,7}] [-z] [-v]

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
  -n, --nmfoutput       write the NMF matrices to the output directory
  -r {2,3,4,5,6,7,8,9,10}, --rank {2,3,4,5,6,7,8,9,10}
                        rank for NMF decomposition
  -l {1,3,5,7}, --length {1,3,5,7}
                        motif length
  -z, --nonoptimized    Use non-optimized pyvcf library instead of cyvcf2
  -v, --verbose         Enable verbose logging
```

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
Rscript diagnostics/doomsayer_diagnostics.r demo/output/diagnostics.yaml
```

An example of a final report is available [here](diagnostics/diagnostics.md).

## Contact

jed.e.carlson@gmail.com, or submit an issue in this repository.
