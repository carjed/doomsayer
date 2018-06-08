# Usage

Running `python helmsman.py -h` will display a help message with detailed information for each option.

```
usage: doomsayer.py [-h] [-c [INT]] [-S [INT]] [-v] [-V] [-M [STR]] -i
                    [/path/to/input.vcf] [-f [/path/to/genome.fa]]
                    [-g [/path/to/sample_batches.txt]]
                    [-s [/path/to/kept_samples.txt]] [-C [INT]] [-X [INT]]
                    [-p [/path/to/project_directory]] [-m [STR]] [-o]
                    [-d [STR]] [-r [INT]] [-F [STR]] [-t [FLOAT]] [-l [INT]]
                    [-R] [-G] [-T [STR]]

optional arguments:
  -h, --help            show this help message and exit
  -c [INT], --cpus [INT]
                        number of CPUs. Must be integer value between 1 and 10
  -S [INT], --seed [INT]
                        random seed for NMF and outlier detection
  -v, --verbose         Enable verbose logging
  -V, --version         show program's version number and exit
  -M [STR], --mode [STR]
                        Mode for parsing input. Must be one of {vcf, agg,
                        txt}. Defaults to VCF mode.
  -i [/path/to/input.vcf], --input [/path/to/input.vcf]
                        In VCF mode (default) input file is a VCF or text file
                        containing paths of multiple VCFs. Defaults to accept
                        input from STDIN with "--input -". In aggregation
                        mode, input file is a text file containing mutation
                        subtype count matrices, or paths of multiple such
                        matrices. In plain text mode, input file is tab-
                        delimited text file containing 5 columns: CHR, POS,
                        REF, ALT, ID
  -f [/path/to/genome.fa], --fastafile [/path/to/genome.fa]
                        reference fasta file
  -g [/path/to/sample_batches.txt], --groupfile [/path/to/sample_batches.txt]
                        two-column tab-delimited file containing sample IDs
                        (column 1) and group membership (column 2) for pooled
                        analysis
  -s [/path/to/kept_samples.txt], --samplefile [/path/to/kept_samples.txt]
                        file with sample IDs to include (one per line)
  -C [INT], --minsnvs [INT]
                        minimum # of SNVs per individual to be included in
                        analysis. Default is 0.
  -X [INT], --maxac [INT]
                        maximum allele count for SNVs to keep in analysis.
                        Defaults to 1 (singletons) Set to 0 to include all
                        variants.
  -p [/path/to/project_directory], --projectdir [/path/to/project_directory]
                        directory to store output files (do NOT include a
                        trailing '/'). Defaults to ./doomsayer_output
  -m [STR], --matrixname [STR]
                        filename prefix for M matrix [without extension]
  -o, --filterout       in VCF or plain text modes, re-reads input file and
                        writes to STDOUT, omitting records that occur in the
                        detected outliers. To write to a new file, use
                        standard output redirection [ > out.vcf] at the end of
                        the doomsayer.py command
  -d [STR], --decomp [STR]
                        mode for matrix decomposition. Must be one of {nmf,
                        pca}. Defaults to pca.
  -r [INT], --rank [INT]
                        Rank for Matrix decomposition. If --decomp pca, will
                        select first r components. Default [0] will force
                        Doomsayer to iterate through multiple ranks to find an
                        optimal choice.
  -F [STR], --filtermode [STR]
                        Method for detecting outliers. Must be one of {ee,
                        lof, if, any, any2, all, none}. Defaults to ee.
  -t [FLOAT], --threshold [FLOAT]
                        threshold for fraction of potential outliers
  -l [INT], --length [INT]
                        motif length. Allowed values are 1,3,5,7
  -R, --report          automatically generates an HTML-formatted report in R.
  -G, --staticplots     use static ggplot figures instead of interactive
                        plotly figures
  -T [STR], --template [STR]
                        Template for diagnostic report. Must be one of
                        {diagnostics, msa}. Defaults to diagnostics.
```

## Program options

`--cpus [INT]`

Set the number of CPUs to use if input is a list of files.

`--seed [INT]`

Set the random seed to be used for NMF decomposition and outlier detection algorithms. If not specified, a random seed will be assigned and printed to the output if you need to reproduce results of a particular run.

`--verbose`

This flag enables debug logging for troubleshooting.

`--version`

Show program's version number and exit.

## Input options

These options specify the data to be processed by *Doomsayer*. By default, *Doomsayer* assumes input is a VCF file or text file containing file paths of multiple VCF files.

### Reference genome

`--fastafile [/path/to/genome.fa]`

In VCF mode, MAF mode, and plain text mode, you **must** specify the fasta-formatted reference genome to be used to look up the local sequence context for each SNV. This file must corresponding to the same reference genome build used to call the variants contained in the input file. The fasta file may be either uncompressed or bgzip-compressed.

### VCF mode

`--mode vcf`

*Doomsayer* can accept variant call format (VCF) or binary call format (BCF) files, either uncompressed or [bgzip](http://www.htslib.org/doc/tabix.html)-compressed. VCF/BCF files are assumed to be formatted with mandatory columns (`#CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO`). This VCF must contain an allele count field in the INFO column (`[AC=N]`), and include individual genotypes (i.e., the VCF cannot be a "sites-only" file).

#### Single VCF file

`--mode vcf --input /path/to/input.vcf`

The VCF input can be either a file on disk or piped from STDIN (using `--input -`), enabling compatibility with a wide range of existing pipelines and workflows. The VCF can contain any combination of variant types (e.g., SNVs, indels, CNVs), allele counts, or filter flags. By default, *Doomsayer* will only analyze singleton SNVs with a "PASS" value in the FILTER column of the VCF.

#### Multiple VCF files, same samples

`--mode vcf --input /path/to/vcfs.txt`

If your data are spread across multiple files, each containing the same samples (e.g., one VCF per chromosome), *Doomsayer* will accept a plain text file containing the absolute file paths (one per line) for the VCF files to process. This assumes all VCF files contain the exact same sample IDs and that the IDs in the VCF header are ordered identically across all subfiles. This option is best used in conjunction with the `--cpus {N}` argument to enable parallel processing of the files.

#### Multiple VCF files, different samples

`--mode vcf --input /path/to/vcfs.txt --rowwise`

If you instead have VCFs containing non-overlapping subsets of the sample (e.g., one per individual), use the `--rowwise` flag to properly generate the mutation spectra matrix.

#### Restrict to specific samples

`--samplefile /path/to/samples.txt`

If your input VCF file(s) contains samples that you do not wish to analyze, provide a tab-delimited text file containing the sample IDs you wish to **keep**. This file can contain any number of columns, as long as the sample IDs are contained in the first column, and the file has a header with column names (the first column must be named `ID`).

#### Pool samples together

`--samplefile /path/to/samples.txt --groupvar GROUP`

In some cases it may be necessary (or desired) to run *Doomsayer* with samples pooled together in some way. For example, if samples were sequenced on different dates and we want to explicitly look for batch effects, we can supply a file containing sample IDs in the `ID` column and the sequencing date in the `DATE`. If the `--groupvar DATE` option is used, *Doomsayer* will pool samples by the value present in the `DATE` column rather than the `ID` column. 

Note that the `--samplefile` option will operate as above--only samples present in this file will be considered when generating the mutation spectra matrix.

### MAF mode

`--mode maf --input /path/to/input.maf`

*Doomsayer* can also accept [Mutation Annotation Format](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) files. By default, SNVs will be counted per unique identifier in the `Tumor_Sample_Barcode` column, but any column can be specified using the `--groupvar` option.

### TXT mode

`--mode txt --input /path/to/input.txt`

*Doomsayer* can also accept a tab-delimited text file, containing a list of SNVs with the following 5 columns (order must be maintained):

- Chromosome
- Position
- Reference allele
- Alternative allele
- Sample ID (or other identifier variable)

### Aggregation mode

This mode can be used to run *Doomsayer* on one or more existing mutation spectra matrices, which can be obtained from the [output](#output-options) of previous *Doomsayer* runs.

#### Re-run on a single previously-generated matrix

`--mode agg --input /path/to/subtype_count_matrix.txt`

Suppose you have already run *Doomsayer* on a large dataset and generated the mutation spectra matrix, but wish to re-run it with different outlier detection parameters. Rather than re-running on the original data, you can use aggregation mode with the existing mutation spectra matrix as the input file.

#### Combine multiple matrices with different samples

`--mode agg --input /path/to/m_samples.txt`

If you have generated mutation spectra matrices from different **non-overlapping** subsamples (e.g., one VCF per individual), you can use aggregation mode to run a combined analysis on all samples together. Simply put the absolute file paths for each of the matrices (one per line) into a text file named `/path/to/m_samples.txt`

For example, if you have run *Doomsayer* on one sample of N individuals and another of M individuals, this will create an (N+M)x96 matrix and run the combined analysis.

We advise users to exercise caution with interpreting results of aggregated samples, as different datasets can have different SNV distributions.

#### Combine multiple matrices from different genomic regions

`--mode agg --input /path/to/m_regions.txt`

Similarly, if you have generated mutation spectra matrices from different **non-overlapping** genomic regions (e.g., one per chromsome), but in the same sample, put the absolute file paths for each of the matrices (one per line) into a text file named `/path/to/m_regions.txt`, and *Doomsayer* will sum the elements of each matrix. When combining by region, it is essential that the matrices have the **exact** same number of samples.

If you anticipate generating multiple mutation spectra matrices by sample or by region to be aggregated at a later point, we recommend setting a unique name for the output file using the [`--matrixname custom_output_name`](#customize-matrix-prefix) option and/or putting the results of each run in a different directory with the `--projectdir` option, to avoid mixing up which output files are being generated.

## Filtering options

*Doomsayer* is designed to integrate in existing VCF processing workflows, so if you wish to apply sample-level or variant-level filters to your input data prior to running *Doomsayer*, we advise pre-filtering your data (e.g., using bcftools) for greater flexibility and efficiency. However, we have included options that enable the user to perform a few common filtering tasks directly in *Doomsayer*.

### Exclude samples with too few SNVs

`--minsnvs [INT]`

In some cases, samples may appear as outliers simply because they have very few observed SNVs. The `--minsnvs` parameter forces *Doomsayer* to only evaluate individuals with at least X observed SNVs. By default, `minsnvs=0`, and all samples are retained. If this option is enabled, the low-SNV outliers will be written to a separate file named `doomsayer_snvs_lt{X}.txt` in the output directory, and the keep/drop lists will be derived from the remaining subset of samples.

Note that the `chisq` outlier detection mode is generally more robust to low-SNV outliers.

### Specify maximum allele count of SNVs to analyze

`--maxac [INT]`

By default, *Doomsayer* will evaluate only singleton SNVs (`--maxac 1`). To include more common SNVs, increase the value of this parameter. To include all SNVs, use `--maxac 0`. This is not recommended when using *Doomsayer* for QC functionality, but may be useful for more general mutation signature analysis.

## Output options

### Output directory

`--projectdir [/path/to/project_directory]`

By default, output files are stored in a subfolder named `/doomsayer_output/` in your current directory. The `--projectdir` parameter will specify a custom name/location for the *Doomsayer* output directory.

The output folder will contain the following data files (depending on which options are used):

1. the M matrix of observed mutation spectra per sample
2. the H matrix containing loadings of the mutation subtypes in each signature
3. the W matrix containing the contributions of each signature within each sample's mutation spectrum

The default output consists of two plain-text files: `doomsayer_keep.txt`, containing the IDs of samples that passed the filter, and `doomsayer_drop.txt`, containing the IDs of samples that failed the filter. These lists can be easily integrated into bcftools, vcftools, PLINK, and other utilities to exclude outlier samples from your downstream analyses.

### Customize matrix prefix

`--matrixname [STR]`

This parameter allows you to specify a custom name for the mutation spectra matrix file. The default is `--matrixname subtype_count_matrix`. This is particularly useful if you wish to generate an initial count matrix then play around with the outlier detection parameters.

### Auto-filter flagged outliers

`--filterout`

If the input is a single VCF or text file of sites, this flag will force *Doomsayer* to re-read the input file and write a new file of the same format, removing any samples specified in the `doomsayer_drop.txt` list that has been generated.

This command will print the VCF header and records to STDOUT, which can then be piped to compatible programs, or redirected to a file using standard UNIX redirection (i.e., adding `> /path/to/output.vcf` at the end of the command), allowing *Doomsayer* to be used as an intermediate step of more complex VCF processing pipelines.

As far as possible, all records present in the original input will be preserved in the output, with three notable exceptions:
1. All samples in the drop list (and their genotypes) will be omitted
2.  **Variants that were unique to the set of dropped samples will be automatically removed.**
3. Four of the INFO fields will be updated to reflect the smaller set of retained samples:
    1. The total number of alleles (AN)
    2. allele count in genotypes (AC)
    3. number of samples (NS)
    4. combined depth (DP)

## Matrix decomposition options

### NMF decomposition

`--decomp nmf`

*Doomsayer* will decompose the subtype count matrix using Nonnegative Matrix Factorization (NMF), as implemented in the [nimfa]() Python library.


### PCA decomposition

`--decomp pca`

*Doomsayer* will decompose the subtype count matrix using Principal Component Analysis (PCA).

### Motif length

`--length {1,3,5,7}`

This parameter specifies the (symmetric) motif length to be considered in determining the mutation subtypes. The default (3) produces 96 3-mer subtypes, which is likely sufficient for QC of whole-genome germline variants. If you are using *Doomsayer* for QC of whole-exome germline variants, you will likely need to specify `--length 1`, as most individuals will not have enough singleton variants to accurately infer higher-order mutation signatures. Values greater than 3 will likely lead to extremely noisy or unreliable QC results for germline variants, but may be useful for examining somatic sequencing data.

### Decomposition rank

`--rank [INT]`

Set the rank of the NMF decomposition to R components if using `--decomp nmf`, or select the first R principal components if using `--decomp pca`.

## Outlier detection options

`--filtermode [STR]`

Once we have reduced the input matrix to R components, we want to identify which individuals are outliers in that R-dimensional space. *Doomsayer* can apply three different outlier detection algorithms: elliptic envelopes, local outlier factors, and isolation forests. For a detailed overview of these methods, see the [scikit-learn outlier detection tutorial](http://scikit-learn.org/stable/modules/outlier_detection.html).

### Detection modes

#### Elliptic envelope
- `--filtermode ee`: Outliers are determined by assuming components follow a multivariate Gaussian distribution, and computing an elliptic envelope.
 

#### Local outlier factor

- `--filtermode lof`: Outliers are determined using local outlier factors. 

#### Isolation forests
- `--filtermode if`: Outliers are determined using isolation forests.


#### Omnibus filtering

In addition, *Doomsayer* can perform omnibus filtering using two or more of the above methods:
- `--filtermode any2`: Outliers must be called by at least two of the above methods
- `--filtermode all`: Outliers must be called by all three of the above methods (will call the fewest outliers)

#### Skip outlier detection

To disable outlier detection, use `--filtermode none` (useful if you plan on aggregating data from multiple runs)

### Detection threshold

`--threshold [FLOAT]`

This parameter specifies the fraction of the sample to flag as potential outliers by a given algorithm. If we set `--threshold 0.05` (the default) and use either the `ee`, `lof`, or `if` filter modes, *Doomsayer* will always flag ~5% of the sample as outliers.

When using the omnibus filtering modes, fewer outliers will be flagged, as they must meet more stringent criteria.

## Generate a diagnostic report
`--report` + `--template`

The `--report` option will tell *Doomsayer* to execute the `generate_report.r` script. This script will copy an [RMarkdown template](report_templates) (specified with the `-T, --template` parameter) into `~/doomsayer_output/` (or the user-specified project directory) and render an HTML-formatted report detailing the results of your *Doomsayer* run.

By default, plots are fully interactive using the Plotly engine. If you want to download static .png images, just mouse over the upper right corner of the figure and click the camera icon.

Alternatively, you can set the `--staticplots` option in your command, which will generate a report without any interactive elements. The standalone .png images will be saved to `~/doomsayer_output/report_files/figure-html/`. 

### Use static images in report

`--staticplots`

By default, figures in the *Doomsayer* report will be fully interactive. This option will force the report to include only static figures.

### Evaluate sample-level variables

`--svars [COLUMN1,COLUMN2,COLUMN3]`

If you have run *Doomsayer* with a sample file, you can include the `--svars` option, specifying a comma-separated list of variables to evaluate in the diagnostic report. This report will include boxplots showing the distribution of each variable for inliers compared to the different outlier clusters. Variables are assumed to be numeric.