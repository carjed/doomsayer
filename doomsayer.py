#!/usr/bin/python

from __future__ import print_function
import os
import sys
import argparse
import itertools
import timeit
import numpy as np
from subprocess import call
from pyfaidx import Fasta
# from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from sklearn.decomposition import NMF
from util import *

###############################################################################
# Parse arguments
###############################################################################
start = timeit.default_timer()

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputvcf",
                    help="input vcf file (use - for STDIN)",
                    required=True,
                    nargs='?',
                    type=str,
                    default=sys.stdin)

parser.add_argument("-f", "--fastafile",
                    help="fasta file name",
                    required=True,
                    type=str,
                    default="chr20.fasta.gz")

parser.add_argument("-p", "--projectdir",
                    help="directory to store output files \
                        (do NOT include a trailing '/')",
                    type=str,
                    default="doomsayer_output")

parser.add_argument("-o", "--outputtovcf",
                    help="filter input VCF (writes to stdout--use standard \
                        output redirection [ > out.vcf] to write the output \
                        VCF to a file on disk)",
                    action="store_true")

parser.add_argument("-n", "--nofilter",
                    help="disables generation of keep/drop lists, \
                        turns off default filtering criteria, \
                        and evaluates all sites in the input VCF. \
                        Useful if analyzing somatic data or pre-filtering \
                        with another tool)",
                    action="store_true")

parser.add_argument("-d", "--diagnostics",
                    help="write the NMF matrices to the output directory \
                        and generate yaml config to be passed to \
                        diagnostic script",
                    action="store_true")

parser.add_argument("-ad", "--autodiagnostics",
                    help="same as the --diagnostics option, but automatically \
                        generates diagnostic report",
                    action="store_true")

parser.add_argument("-r", "--rank",
                    help="rank for NMF decomposition",
                    type=int,
                    choices=range(2,11),
                    default=3)

parser.add_argument("-l", "--length",
                    help="motif length",
                    type=int,
                    choices=[1,3,5,7],
                    default=3)

parser.add_argument("-v", "--verbose",
                    help="Enable verbose logging",
                    action="store_true")

args = parser.parse_args()

###############################################################################
# Initialize project directory fasta and vcf files
###############################################################################
projdir = os.path.realpath(args.projectdir)
if args.verbose:
    eprint("checking if directory", projdir, "exists...")

if not os.path.exists(args.projectdir):
    if args.verbose:
        eprint("Creating output directory:", projdir)
    os.makedirs(args.projectdir)
else:
    eprint(projdir, "already exists") if args.verbose else None

eprint("Initializing reference genome...") if args.verbose else None
# fasta_reader = SeqIO.index(args.fastafile, 'fasta')
# fasta_reader = SeqIO.parse(args.fastafile, 'fasta')
fasta_reader = Fasta(args.fastafile, read_ahead=10000)

import cyvcf2 as vcf
from cyvcf2 import VCF
from cyvcf2 import Writer
vcf_reader = VCF(args.inputvcf, mode='rb', gts012=True)

###############################################################################
# index samples
###############################################################################
eprint("Indexing samples in", args.inputvcf, "...") if args.verbose else None
samples = vcf_reader.samples

samples_dict = {}
for i in range(len(samples)):
    samples_dict[samples[i]] = i

eprint(len(samples), "samples found") if args.verbose else None

###############################################################################
# index subtypes
###############################################################################
eprint("indexing subtypes...") if args.verbose else None

subtypes_dict = indexSubtypes(args)

M = np.zeros((len(samples), len(subtypes_dict)))

###############################################################################
# Query records in VCF and build matrix
###############################################################################
eprint("Parsing VCF records...") if args.verbose else None
numsites_keep = 0
numsites_skip = 0
chrseq = '1'

for record in vcf_reader:

    # debug--testing performance for triallelic sites
    # if(record.POS==91628): # triallelic site
    # if(record.POS==63549):
    #     eprint(acval)
    #     eprint(record.gt_types.tolist().index(1))

    # Filter by allele count, SNP status, and FILTER column
    if record.is_snp:
        acval = record.INFO['AC']
        if ((acval==1 and record.FILTER is None) or args.nofilter):

            # check and update chromosome sequence
            if record.CHROM != chrseq:
                if args.verbose:
                    eprint("Loading chromosome", record.CHROM, "reference...")
                # seq = fasta_reader[record.CHROM].seq
                # seq = SeqIO.to_dict(fasta_reader)[record.CHROM].seq
                sequence = fasta_reader[record.CHROM]
                chrseq = record.CHROM

            lseq = sequence[record.POS-2:record.POS+1].seq
            mu_type = record.REF + str(record.ALT[0])
            category = getCategory(mu_type)
            motif_a = getMotif(record.POS, lseq)
            subtype = str(category + "-" + motif_a)

            # use quick singleton lookup for default QC option
            if not args.nofilter:
                sample=samples[record.gt_types.tolist().index(1)]
                M[samples_dict[sample], subtypes_dict[subtype]] += 1

                # sample=np.where(record.gt_types == 1)[0]
                # M[sample, subtypes_dict[subtype]] += 1
            else:
                # sample=samples[record.gt_types.tolist().index(1)]
                samples_het = np.where(record.gt_types == 1)[0]
                M[samples_het, subtypes_dict[subtype]] += 1
                # for s1 in samples_het:
                #     M[s1, subtypes_dict[subtype]] += 1

                samples_hom = np.where(record.gt_types == 2)[0]
                M[samples_hom, subtypes_dict[subtype]] += 2
                # for s2 in samples_hom:
                #     M[s2, subtypes_dict[subtype]] += 2

                # eprint(record.POS, s2)
                # sample_gts=record.gt_types.tolist()
                # s = 0;
                # for gt in sample_gts:
                #     if gt == 1:
                #         M[s, subtypes_dict[subtype]] += 1
                #
                #         # if samples[s] == "1497-RMM-0968":
                #         #     print(record.CHROM, record.POS,
                #         #         record.REF, record.ALT[0], samples[s], subtype)
                #
                #     elif gt == 2:
                #         M[s, subtypes_dict[subtype]] += 2
                #     s += 1

            numsites_keep += 1
        else:
            numsites_skip += 1

        if args.verbose:
            if (numsites_keep > 0 and numsites_keep%100000==0):
                eprint("Processed", numsites_keep, "sites")

if numsites_keep == 0:
    eprint("No SNVs found. Please check your VCF file")
    sys.exit()

if args.verbose:
    eprint(numsites_keep, "sites kept")
    eprint(numsites_skip, "sites skipped")

vcf_reader.close()

###############################################################################
# Run NMF on final matrix
###############################################################################
# M /= np.max(np.abs(M),axis=1)
# drop rows with all zeroes
# M = M[~(M==0).all(1)]

# Convert count per entry to fraction of total per sample
M_f = M/(M.sum(axis=1)+1e-8)[:,None]

model = NMF(n_components=args.rank, init='random', random_state=0)
model.fit(M_f)

# Get loadings of subtypes per signature
H = model.components_

# Get signature contributions per sample
W = model.fit_transform(M)
W = W/W.sum(axis=1)[:,None]

# W= W[~np.isnan(W).any(axis=1)]

colmeans = np.mean(W, axis=0)
colstd = np.std(W, axis=0)
upper = colmeans+2*colstd
lower = colmeans-2*colstd

# output NMF results
if(args.diagnostics or args.autodiagnostics):
    if args.verbose:
        eprint("Mean signature contributions: ", colmeans)
        eprint("StdDev:", colstd)
        eprint("Upper:", upper)
        eprint("Lower:", lower)
        eprint("Writing NMF results")
    diagWrite(projdir, M, M_f, W, H, subtypes_dict, samples, args)

###############################################################################
# Write keep and drop lists
###############################################################################
if not args.nofilter:
    keep_samples = []
    drop_samples = []
    # Check for outliers
    i=0
    for n in W:
        if(np.greater(n, upper).any() or np.less(n, lower).any()):
            drop_samples.append(samples[i])
            # eprint("Adding", samples[i], "to drop list") if args.verbose else None
        else:
            keep_samples.append(samples[i])
        i += 1

    eprint("Printing keep and drop lists") if args.verbose else None
    keep_path = projdir + "/keep_samples.txt"
    np.savetxt(keep_path, keep_samples, delimiter='\t', fmt="%s")

    drop_path = projdir + "/drop_samples.txt"
    np.savetxt(drop_path, drop_samples, delimiter='\t', fmt="%s")
else:
    eprint("You are running with the --nofilter option. Keep and drop lists will not be generated")

###############################################################################
# write output vcf
###############################################################################
if args.outputtovcf:
    eprint("Filtering VCF by drop list...") if args.verbose else None
    keep_test = keep_samples[0:10]
    vcf = VCF(args.inputvcf, samples=keep_test, mode='rb')
    # vcf = VCF(args.inputvcf, samples=keep_samples, mode='rb')

    print(vcf.raw_header.rstrip())
    for v in vcf:
        v.INFO['AC'] = str(v.num_het + v.num_hom_alt*2)

        if int(v.INFO['AC']) > 0:
            v.INFO['NS'] = str(v.num_called)
            v.INFO['AN'] = str(2*v.num_called)
            v.INFO['DP'] = str(np.sum(v.format('DP')))
            print(str(v).rstrip())

    vcf.close()

if args.autodiagnostics:
    cmd = "diagnostics/doomsayer_diagnostics.r " + projdir + "/config.yaml"
    if args.verbose:
        eprint("Rscript will run the following command:")
        eprint("Rscript " + cmd)
        eprint("Auto-generating diagnostic report...")
    call("/usr/bin/Rscript --vanilla " + cmd, shell=True)


stop = timeit.default_timer()
tottime = round(stop - start, 2)
eprint("Total runtime:", tottime, "seconds") if args.verbose else None
