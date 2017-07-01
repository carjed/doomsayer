#!/usr/bin/python

from __future__ import print_function
import os
import sys
import argparse
import itertools
import timeit
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from sklearn.decomposition import NMF

###############################################################################
# print to stderr
###############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

###############################################################################
# collapse mutation types per strand symmetry
###############################################################################
def getCategory(mu_type):
    if (mu_type == "AC" or mu_type == "TG"):
        category = "A_C"
    if (mu_type == "AG" or mu_type == "TC"):
        category = "A_G"
    if (mu_type == "AT" or mu_type == "TA"):
        category = "A_T"
    if (mu_type == "CA" or mu_type == "GT"):
        category = "C_A"
    if (mu_type == "CG" or mu_type == "GC"):
        category = "C_G"
    if (mu_type == "CT" or mu_type == "GA"):
        category = "C_T"
    return category

###############################################################################
# query reference genome for local sequence motif
###############################################################################
def getMotif(pos):
    # get 3-mer motif
    motif = seq[pos-2:pos+1]
    altmotif = motif.reverse_complement()

    m1 = motif[1]
    m2 = altmotif[1]

    if m1 < m2:
        motif_a = motif
    else:
        motif_a = altmotif

    return motif_a

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

parser.add_argument("-d", "--diagnostics",
                    help="write the NMF matrices to the output directory \
                        and generate yaml config to be passed to \
                        diagnostic script",
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

parser.add_argument("-z", "--nonoptimized",
                    help="Use non-optimized pyvcf library instead of cyvcf2",
                    action="store_true")

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
fasta_reader = SeqIO.index(args.fastafile, 'fasta')

if args.nonoptimized:
    eprint("Using pyvcf library...") if args.verbose else None
    import vcf
    # vcf_reader = vcf.Reader(args.inputvcf, 'rb')
    if args.inputvcf=="-":
        vcf_reader = vcf.Reader(sys.stdin)
    else:
        vcf_reader = vcf.Reader(open(args.inputvcf), 'rb')
else:
    eprint("Using cyvcf library...") if args.verbose else None
    import cyvcf2 as vcf
    from cyvcf2 import VCF
    from cyvcf2 import Writer
    vcf_reader = VCF(args.inputvcf, mode='rb')

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
categories = ["A_C", "A_G", "A_T", "C_G", "C_T", "C_A"]
bases = ["A", "C", "G", "T"]

motiflength = args.length
flank = (motiflength-1)/2

kmers = itertools.product(bases, repeat=motiflength-1)

subtypes_list = []
# i = 0
for kmer in kmers:
    kmerstr = ''.join(kmer)
    for category in categories:
        ref = category[0]
        # subtype = category + "-" + kmer[0] + ref + kmer[1]
        if motiflength > 1:
            subtype = category + "-" \
                + kmerstr[0:flank] + ref + kmerstr[flank:(motiflength-1)]
        else:
            subtype = category

        subtypes_list.append(subtype)

i = 0
subtypes_dict = {}
for subtype in sorted(subtypes_list):
        subtypes_dict[subtype] = i
        i += 1

M = np.zeros((len(samples), len(subtypes_dict)))

###############################################################################
# Query records in VCF and build matrix
###############################################################################
eprint("Parsing VCF records...") if args.verbose else None
numsites_keep = 0
numsites_skip = 0
chrseq = '1'

for record in vcf_reader:

    # pyvcf2 sets the 'AC' field as an array
    if args.nonoptimized==True:
        acval = record.INFO['AC'][0]
    else:
        acval = record.INFO['AC']

    # Filter by allele count, SNP status, and FILTER column
    if (acval==1 and record.FILTER!='SVM' and record.is_snp):

        # use get_hets() function if using pyvcf, otherwise need to query
        # record.gt_types in cyvcf2
        if args.nonoptimized:
            sample = record.get_hets()[0].sample
        else:
            sample=samples[record.gt_types.tolist().index(1)]

        #
        mu_type = record.REF + str(record.ALT[0])
        category = getCategory(mu_type)

        # check and update chromosome sequence
        if record.CHROM != chrseq:
            if args.verbose:
                eprint("Loading chromosome", record.CHROM, "reference...")
            seq = fasta_reader[record.CHROM].seq
            chrseq = record.CHROM

        motif_a = getMotif(record.POS)
        subtype = str(category + "-" + motif_a)

        # if sample == "1497-RMM-0968":
        #     print(record.CHROM, record.POS, record.REF, record.ALT[0], sample, subtype)

        # query M matrix by sample and subtype; iterate by 1 if record matches
        M[samples_dict[sample], subtypes_dict[subtype]] += 1
        numsites_keep += 1
    else:
        numsites_skip += 1

if numsites_keep == 0:
    eprint("No singleton SNVs found. Please check your VCF file")
    sys.exit()

if args.verbose:
    eprint(numsites_keep, "sites kept")
    eprint(numsites_skip, "sites skipped")

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

if args.verbose:
    eprint("Mean signature contributions: ", colmeans)
    eprint("StdDev:", colstd)
    eprint("Upper:", upper)
    eprint("Lower:", lower)

keep_samples = []
drop_samples = []
# Check for outliers
i=0
for n in W:
    if(np.greater(n, upper).any() or np.less(n, lower).any()):
        drop_samples.append(samples[i])
    else:
        keep_samples.append(samples[i])
        # print(n)
    i += 1

vcf_reader.close()
# print(keep_samples[0:10])
# print(len(keep_samples))
# print(drop_samples[0:10])

###############################################################################
# Write output data
###############################################################################
keep_path = projdir + "/keep_samples.txt"
np.savetxt(keep_path, keep_samples, delimiter='\t', fmt="%s")

drop_path = projdir + "/drop_samples.txt"
np.savetxt(drop_path, drop_samples, delimiter='\t', fmt="%s")

if args.diagnostics:
    yaml = open(projdir + "/config.yaml","w+")
    print("# Config file for doomsayer_diagnostics.r", file=yaml)
    print("keep_path: " + keep_path, file=yaml)
    print("drop_path: " + drop_path, file=yaml)

    colnames = ["ID"]
    M_colnames = colnames + list(sorted(subtypes_dict.keys()))

    ###############################
    # M matrix (counts)
    ###############################
    if args.verbose:
        eprint("Saving M matrix (observed spectra counts)")

    # add ID as first column
    M_fmt = np.concatenate((np.array([samples]).T, M), axis=1)

    # add header
    M_fmt = np.concatenate((np.array([M_colnames]), M_fmt), axis=0)

    # write out
    M_path = projdir + "/NMF_M_spectra.txt"
    print("M_path: " + M_path, file=yaml)
    np.savetxt(M_path, M_fmt, delimiter='\t', fmt="%s")

    ###############################
    # M matrix (rates)
    ###############################
    if args.verbose:
        eprint("Saving M_f matrix (observed spectra rates)")

    # add ID as first column
    M_fmt = np.concatenate((np.array([samples]).T, M_f), axis=1)

    # add header
    M_fmt = np.concatenate((np.array([M_colnames]), M_fmt), axis=0)

    # write out
    M_path_rates = projdir + "/NMF_M_spectra_rates.txt"
    print("M_path_rates: " + M_path_rates, file=yaml)
    np.savetxt(M_path_rates, M_fmt, delimiter='\t', fmt="%s")

    ###############################
    # W matrix (contributions)
    ###############################
    if args.verbose:
        eprint("Saving W matrix (signature contributions per sample)")

    # add ID as first column
    W_fmt = np.concatenate((np.array([samples]).T, W), axis=1)

    # add header
    W_colnames = colnames + ["S" + str(i) for i in xrange(1,args.rank+1)]
    W_fmt = np.concatenate((np.array([W_colnames]), W_fmt), axis=0)

    # write out
    W_path = projdir + "/NMF_W_sig_contribs.txt"
    print("W_path: " + W_path, file=yaml)
    np.savetxt(W_path, W_fmt, delimiter='\t', fmt="%s")

    ###############################
    # H matrix (loadings)
    ###############################
    if args.verbose:
        eprint("Saving H matrix (feature loadings per signature)")

    # add signature ID as firsat column
    H_rownames = ["S" + str(i) for i in xrange(1,args.rank+1)]
    H_fmt = np.concatenate((np.array([H_rownames]).T, H), axis=1)

    # add header
    H_colnames = ["Sig"] + list(sorted(subtypes_dict.keys()))
    H_fmt = np.concatenate((np.array([H_colnames]), H_fmt), axis=0)

    # write out
    H_path = projdir + "/NMF_H_sig_loads.txt"
    print("H_path: " + H_path, file=yaml)
    np.savetxt(H_path, H_fmt, delimiter='\t', fmt="%s")

    yaml.close()

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

stop = timeit.default_timer()
tottime = round(stop - start, 2)
eprint("Total runtime:", tottime, "seconds") if args.verbose else None
