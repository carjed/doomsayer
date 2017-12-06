#!/usr/bin/python

from __future__ import print_function
import os
import sys

sys.path.append(os.getcwd())

import textwrap
import itertools
import timeit
import collections
import csv
import nimfa
import re
from pandas import *
import numpy as np
import cyvcf2 as vcf
from cyvcf2 import VCF
from scipy.stats import chisquare
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

###############################################################################
# print to stderr
###############################################################################
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

###############################################################################
# Custom class for args
###############################################################################
def restricted_float(x):
    x = float(x)
    if x < 1.0:
        raise argparse.ArgumentTypeError("%r must be greater than 1"%(x,))
    return x

###############################################################################
# collapse mutation types per strand symmetry
###############################################################################
def getCategory(mu_type):
    if re.match("^[ACGT]*$", mu_type):
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
    else:
        category = "unknown"
    return category

###############################################################################
# query reference genome for local sequence motif
###############################################################################
def getMotif(pos, sequence):
    motif = Seq(sequence, IUPAC.unambiguous_dna)
    altmotif = motif.reverse_complement()
    central_base = (len(motif)-1)//2

    m1 = motif[central_base]
    m2 = altmotif[central_base]

    if m1 < m2:
        motif_a = motif
    else:
        motif_a = altmotif

    return motif_a

###############################################################################
# define k-mer mutation subtypes
###############################################################################
def indexSubtypes(motiflength):
    categories = ["A_C", "A_G", "A_T", "C_G", "C_T", "C_A"]
    bases = ["A", "C", "G", "T"]
    flank = (motiflength-1)//2

    if motiflength > 1:
        kmers = itertools.product(bases, repeat=motiflength-1)

        subtypes_list = []

        for kmer in kmers:
            kmerstr = ''.join(kmer)

            for category in categories:
                ref = category[0]

                subtype = category + "." \
                    + kmerstr[0:flank] + ref + kmerstr[flank:(motiflength-1)]

                subtypes_list.append(subtype)
    else:
        ext = [".A", ".C"]
        extr = list(np.repeat(ext,3))
        subtypes_list = [m+n for m,n in zip(categories,extr)]

    i = 0
    subtypes_dict = {}
    for subtype in sorted(subtypes_list):
        subtypes_dict[subtype] = i
        i += 1

    return subtypes_dict

###############################################################################
# Build dictionary with sample ID as key, group ID as value
###############################################################################
def indexGroups(groupfile):
    sg_dict = {}
    with open(groupfile) as sg_file:
        for line in sg_file:
           (key, val) = line.split()
           sg_dict[key] = val

    samples = sorted(list(set(sg_dict.values())))
    return samples

###############################################################################
# get samples from VCF file
###############################################################################
def getSamplesVCF(args, inputvcf):
    if args.samplefile:
        with open(args.samplefile) as f:
            keep_samples = f.read().splitlines()

        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True, samples=keep_samples)
        # vcf_reader.set_samples(keep_samples) # <- set_samples() subsets VCF
    else:
        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True)

    if args.groupfile:
        all_samples = vcf_reader.samples
        samples = indexGroups(args.groupfile)
    else:
        samples = vcf_reader.samples

    vcf_reader.close()
    return samples

###############################################################################
# Main function for parsing VCF
###############################################################################
def processVCF(args, inputvcf, subtypes_dict, par):
    eprint("Initializing reference genome...") if args.verbose else None
    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)
    # record_dict = SeqIO.to_dict(SeqIO.parse(args.fastafile, "fasta"))

    # 'demo/input/keep.txt'
    if args.samplefile:
        with open(args.samplefile) as f:
            keep_samples = f.read().splitlines()
        eprint(len(keep_samples), "samples kept") if args.verbose else None

        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True, samples=keep_samples)
        # vcf_reader.set_samples(keep_samples) # <- set_samples() subsets VCF
    else:
        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True)

    nbp = (args.length-1)//2

    # index samples
    eprint("Indexing samples in", inputvcf, "...") if args.verbose else None

    if args.groupfile:
        all_samples = vcf_reader.samples
        samples = indexGroups(args.groupfile)
    else:
        samples = vcf_reader.samples

    samples_dict = {}
    for i in range(len(samples)):
        samples_dict[samples[i]] = i

    eprint(len(samples), "samples indexed") if args.verbose else None

    # Query records in VCF and build matrix
    eprint("Parsing VCF records...") if args.verbose else None
    M = np.zeros((len(samples), len(subtypes_dict)))
    numsites_keep = 0
    numsites_skip = 0
    chrseq = '0'

    batchit = 0
    sample_batch = []
    subtype_batch = []

    for record in vcf_reader:
        # debug--testing performance for triallelic sites
        # if(record.POS==91628): # triallelic site
        # if(record.POS==63549):
        #     eprint(acval)
        #     eprint(record.gt_types.tolist().index(1))

        # Filter by allele count, SNP status, and FILTER column
        # if len(record.ALT[0])==1:
        if record.is_snp:
            # eprint("SNP check: PASS")
            acval = record.INFO['AC']
            # eprint(record.POS, acval)

            if ((acval<=args.maxac or args.maxac==0) and record.FILTER is None):
                # eprint(record.CHROM, record.POS, record.REF, record.ALT[0],
                    # acval, record.FILTER)

                # check and update chromosome sequence
                if record.CHROM != chrseq:
                    sequence = fasta_reader[record.CHROM]
                    chrseq = record.CHROM

                if nbp > 0:
                    lseq = sequence[record.POS-(nbp+1):record.POS+nbp].seq
                else:
                    lseq = sequence[record.POS-1].seq

                mu_type = record.REF + str(record.ALT[0])
                category = getCategory(mu_type)
                motif_a = getMotif(record.POS, lseq)
                subtype = str(category + "." + motif_a)

                if subtype in subtypes_dict:
                    st = subtypes_dict[subtype]

                    if args.groupfile:
                        sample = all_samples[record.gt_types.tolist().index(1)]

                        if sample in sg_dict:
                            sample_gp = sg_dict[sample]
                            ind = samples.index(sample_gp)
                            M[ind,st] += 1
                    else:
                        M[:,st] = M[:,st]+record.gt_types

                    numsites_keep += 1

                else:
                    numsites_skip += 1

                if args.verbose:
                    if (numsites_keep%10000==0):
                        eprint("Processed", numsites_keep, "sites",
                            "(Skipped", numsites_skip, "sites)")
            else:
                numsites_skip += 1

    if args.verbose:
        eprint(numsites_keep, "sites kept")
        eprint(numsites_skip, "sites skipped")

    vcf_reader.close()
    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)

    if par:
        return M
    else:
        return out

###############################################################################
# process tab-delimited text file, containing the following columns:
# CHR    POS    REF    ALT    SAMPLE_ID
###############################################################################
def processTxt(args, subtypes_dict):
    eprint("Initializing reference genome...") if args.verbose else None
    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)

    nbp = (args.length-1)//2
    samples_dict = {}

    # M = np.zeros((len(samples), len(subtypes_dict)))
    numsites_keep = 0
    numsites_skip = 0
    chrseq = '0'

    with open(args.input, 'r') as f:
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            chrom = row[0]
            pos = int(row[1])
            ref = row[2]
            alt = row[3]
            sample = row[4]

            if chrom != chrseq:
                if args.verbose:
                    eprint("Loading chromosome", chrom, "reference...")

                sequence = fasta_reader[chrom]
                chrseq = chrom

            if(len(alt) == 1 and len(ref)==1):
                mu_type = ref + alt
                category = getCategory(mu_type)
                if nbp > 0:
                    lseq = sequence[pos-(nbp+1):pos+nbp].seq
                else:
                    lseq = sequence[pos-1].seq
                    # eprint("lseq:", lseq)
                motif_a = getMotif(pos, lseq)
                subtype = str(category + "-" + motif_a)
                st = subtypes_dict[subtype]

                if sample not in samples_dict:
                    samples_dict[sample] = {}

                if subtype not in samples_dict[sample]:
                    samples_dict[sample][subtype] = 1
                else:
                    samples_dict[sample][subtype] += 1

        M = DataFrame(samples_dict).T.fillna(0).values
        samples = sorted(samples_dict)

    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
    return out

###############################################################################
# get samples from input M matrix when using aggregation mode
###############################################################################
def getSamples(fh):
    samples = np.loadtxt(fh,
        dtype='S20',
        skiprows=1,
        delimiter='\t',
        usecols=(0,))

    return samples

###############################################################################
# aggregate M matrices from list of input files
###############################################################################
def aggregateM(inputM, subtypes_dict):
    colnames = ["ID"]
    M_colnames = colnames + list(sorted(subtypes_dict.keys()))
    colrange = range(1,len(M_colnames))

    if inputM.lower().endswith('nmf_m_spectra.txt'):
        samples = getSamples(inputM)
        M = np.loadtxt(inputM, skiprows=1, usecols=colrange)
        M = M.astype(np.float)
    else:
        with open(inputM) as f:
            file_list = f.read().splitlines()

        # M output by sample
        if inputM.lower().endswith('m_samples.txt'):

            M_out = np.array([M_colnames])

            for mfile in file_list:
                samples = getSamples(mfile)

                M_it = np.loadtxt(mfile, skiprows=1, usecols=colrange)
                M_it = np.concatenate((np.array([samples]).T, M_it), axis=1)
                M_out = np.concatenate((M_out, M_it), axis=0)

            M = np.delete(M_out, 0, 0)
            M = np.delete(M, 0, 1)
            M = M.astype(np.float)

        # M output by region
        elif inputM.lower().endswith('m_regions.txt'):
            samples = getSamples(file_list[0])

            M_out = np.zeros((len(samples), len(M_colnames)-1))
            for mfile in file_list:
                M_it = np.loadtxt(mfile, skiprows=1, usecols=colrange)
                M_out = np.add(M_out, M_it)

            M = M_out.astype(np.float)

    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
    return out

###############################################################################
# Generate keep/drop lists
###############################################################################
def detectOutliers(M, samples, filtermode, threshold):
    M_f = M/(M.sum(axis=1)+1e-8)[:,None]

    keep_samples = []
    drop_samples = []
    drop_indices = []
    # lowsnv_samples = []

    if filtermode == "fold":
        i=0
        M_err_d = np.divide(M_f, np.mean(M_f, axis=0))
        for row in M_err_d:
            if any(err > threshold for err in row):
                drop_samples.append(samples.flatten()[i])
                drop_indices.append(i)
            else:
                keep_samples.append(samples.flatten()[i])
            i += 1
    elif filtermode == "chisq":
        i=0
        mean_spectrum = np.mean(M, axis=0)
        # n_pass = sum(np.sum(M, axis=1) > minsnvs)
        for row in M:
            exp_spectrum = mean_spectrum*sum(row)/sum(mean_spectrum)
            pval = chisquare(row, f_exp=exp_spectrum)[1]
            if pval < 0.05/M.shape[0]:
                drop_samples.append(samples.flatten()[i])
                drop_indices.append(i)
            else:
                keep_samples.append(samples.flatten()[i])
            i += 1
    elif filtermode == "sd":
        i=0
        mean_spectrum = np.mean(M_f, axis=0)
        spec_std = np.std(M_f, axis=0, ddof=1)
        std_threshold = mean_spectrum + threshold*spec_std

        for row in M_f:
            if np.greater(row, std_threshold).any():
                drop_samples.append(samples.flatten()[i])
                drop_indices.append(i)
            else:
                keep_samples.append(samples.flatten()[i])
            i += 1

    out_handles = ['keep_samples',
        'drop_samples',
        'drop_indices']

    out = collections.namedtuple('Out', out_handles) \
        (keep_samples, drop_samples, drop_indices)
    return out

###############################################################################
# run NMF on input matrix
###############################################################################
def NMFRun(M_run, args, projdir, samples, subtypes_dict):
    if args.rank > 0:
        if args.verbose:
            eprint("Running NMF with rank =", args.rank)

        model = nimfa.Nmf(M_run,
            rank=args.rank,
            update="divergence",
            objective='div',
            n_run=1,
            max_iter=200)
        model_fit = model()
        evar = model_fit.fit.evar()
        maxind = args.rank

    elif args.rank == 0:
        if args.verbose:
            eprint("Finding optimal rank for NMF...")
        evarprev = 0
        for i in range(1,6):
            model = nimfa.Nmf(M_run,
                rank=i,
                update="divergence",
                objective='div',
                n_run=1,
                max_iter=200)
            model_fit = model()
            evar = model_fit.fit.evar()
            if args.verbose:
                eprint("Explained variance for rank " + str(i) + ":", evar)
            # if evar > 0.8:
            if(i > 2 and evar - evarprev < 0.001):
                if args.verbose:
                    eprint(textwrap.dedent("""\
                            Stopping condition met: <0.1 percent difference
                            in explained variation between ranks
                            """))
                    model = nimfa.Nmf(M_run,
                        rank=i-1,
                        update="divergence",
                        objective='div',
                        n_run=1,
                        max_iter=200)
                    model_fit = model()
                break
            evarprev = evar

    W = model_fit.basis()
    H = model_fit.coef()

    out = collections.namedtuple('Out', ['W', 'H'])(W, H)
    return out

###############################################################################
# write M matrix
###############################################################################
def writeM(M, M_path, subtypes_dict, samples):
    colnames = ["ID"]
    M_colnames = colnames + list(sorted(subtypes_dict.keys()))

    # add ID as first column
    #M_fmt = np.concatenate((np.array([samples]).T, M.astype('|S20')), axis=1)
    M_fmt = np.concatenate((samples.T, M.astype('|S20')), axis=1)

    # add header
    M_fmt = np.concatenate((np.array([M_colnames]), M_fmt), axis=0)

    # write out
    np.savetxt(M_path, M_fmt, delimiter='\t', fmt="%s")

###############################################################################
# write W matrix
###############################################################################
def writeW(W, W_path, samples):
    colnames = ["ID"]

    # add ID as first column
    W_fmt = np.concatenate((samples.T, W.astype('|S20')), axis=1)
    num_samples, num_sigs = W.shape

    # add header
    W_colnames = colnames + ["S" + str(i) for i in range(1,num_sigs+1)]
    W_fmt = np.concatenate((np.array([W_colnames]), W_fmt), axis=0)

    # write out
    np.savetxt(W_path, W_fmt, delimiter='\t', fmt="%s")

###############################################################################
# write H matrix
###############################################################################
def writeH(H, H_path, subtypes_dict):
    num_sigs, num_subtypes = H.shape

    # add signature ID as first column
    H_rownames = ["S" + str(i) for i in range(1,num_sigs+1)]
    H_fmt = np.concatenate((np.array([H_rownames]).T, H.astype('|S20')), axis=1)

    H_colnames = ["Sig"] + list(sorted(subtypes_dict.keys()))
    H_fmt = np.concatenate((np.array([H_colnames]), H_fmt), axis=0)

    # write out
    np.savetxt(H_path, H_fmt, delimiter='\t', fmt="%s")

###############################################################################
# write RMSE per sample
###############################################################################
def writeRMSE(M_rmse, rmse_path, samples):
	sample_col = samples.T
	rmse_col = np.array([M_rmse]).T
	rmse_arr = np.column_stack((sample_col, rmse_col))
	np.savetxt(rmse_path, rmse_arr, delimiter='\t', fmt="%s")

###############################################################################
# filter VCF input by kept samples
###############################################################################
def filterVCF(inputvcf, keep_samples):
    vcf = VCF(inputvcf, samples=keep_samples, mode='rb')

    print(vcf.raw_header.rstrip())
    for v in vcf:
        v.INFO['AC'] = str(v.num_het + v.num_hom_alt*2)

        if int(v.INFO['AC']) > 0:
            v.INFO['NS'] = str(v.num_called)
            v.INFO['AN'] = str(2*v.num_called)
            v.INFO['DP'] = str(np.sum(v.format('DP')))
            print(str(v).rstrip())

    vcf.close()

###############################################################################
# filter txt input by kept samples
###############################################################################
def filterTXT(inputtxt, keep_samples):
    with open(inputtxt, 'r') as f:
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            chrom = row[0]
            pos = row[1]
            ref = row[2]
            alt = row[3]
            sample = row[4]

            if sample in keep_samples:
                print("\t".join(row))
