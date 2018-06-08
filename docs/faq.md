## FAQ

### What are Doomsayer's hardware requirements?
Doomsayer has no specific hardware dependencies. We recommend ~2GB of free disk space if storing the demo data.

### What type of sequencing data can Doomsayer analyze?
Doomsayer's core functionality, the NMF-based or PCA-based mutation signature analysis, can be used to rapidly evaluate mutation patterns present in any collection of sequencing data, and offers an extremely efficient and almost entirely automated (albeit slightly less flexible) alternative to the functions provided by the [SomaticSignatures](http://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html) R package:

>Gehring JS, Fischer B, Lawrence M, Huber W. SomaticSignatures: inferring mutational signatures from single-nucleotide variants. Bioinformatics. 2015;31(22):3673-3675. doi:10.1093/bioinformatics/btv408.

Although we have only tested Doomsayer with human data, it should be fully capable of analyzing NGS data from any other organism.

Note that the outlier analysis function of Doomsayer is exclusively designed for analyzing **population** variants in **unrelated** individuals. In its current state, this filtering is **not** appropriate for quality control in somatic sequencing data or *de novo* mutation analysis, because genuine biological variability in these data will almost certainly be inferred as errors.

### Doomsayer is not reading my input VCF
Doomsayer uses the cyvcf2 wrapper for htslib libraries for parsing VCFs. If the input VCF is not formatted properly, htslib will often throw an error/warning or stop completely. If you encounter any issues with reading in a VCF file, first confirm that it is properly formatted. [vcf-validator](https://github.com/EBIvariation/vcf-validator) is a useful utility for this. You can also try indexing the VCF with `tabix -p input.vcf.gz` and running a simple bcftools command, such as `bcftools view input.vcf.gz`. Also check that your VCF contains singleton SNVs, individual genotypes are included, and that the FILTER column contains "PASS" values.

### Why am I getting errors about the reference genome?
To build the input matrix of subtype counts, Doomsayer must annotate each site with the surrounding sequence context, using a user-specified reference genome file.

The libraries for parsing fasta files require that the chromosome records in your fasta reference file are formatted identically to the CHROM field of your VCF (e.g., if the fasta header is formatted as `>chrN` but the VCF CHROM field is just `N`). If you encounter any errors, try modifying the the fasta file to either strip or add "chr" to each record.

### Why are my RMarkdown reports not rendering?
If any errors occur when trying to generate a report, it is likely because the pandoc binaries cannot be found on your system. Try installing RStudio; otherwise, you can install precompiled binaries into `doomsayer/pandoc/bin`

Alternatively, you can transfer the contents of the `doomsayer_output` directory to a computer that has RStudio installed and run `knit_with_parameters("report.Rmd")`. In the prompt window, enter `/path/to/doomsayer_output/config.yaml` in the `yaml_cfg` field.