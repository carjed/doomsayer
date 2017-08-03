#!/usr/bin/python

from __future__ import print_function
import os
import sys
import argparse
import itertools
import timeit
import collections
import csv
import nimfa
from pandas import *
import numpy as np
import cyvcf2 as vcf
from cyvcf2 import VCF
from cyvcf2 import Writer
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
def indexSubtypes(args):
    categories = ["A_C", "A_G", "A_T", "C_G", "C_T", "C_A"]
    bases = ["A", "C", "G", "T"]

    motiflength = args.length
    flank = (motiflength-1)//2

    if motiflength > 1:
        kmers = itertools.product(bases, repeat=motiflength-1)

        subtypes_list = []
        # i = 0
        for kmer in kmers:
            kmerstr = ''.join(kmer)
            # eprint(kmerstr)
            for category in categories:
                ref = category[0]

                subtype = category + "." \
                    + kmerstr[0:flank] + ref + kmerstr[flank:(motiflength-1)]

                subtypes_list.append(subtype)
    else:
        ext = [".A", ".C"]
        extr = list(np.repeat(ext,3))
        subtypes_list = [m+n for m,n in zip(categories,extr)]
        # eprint(subtypes_list)

    i = 0
    subtypes_dict = {}
    for subtype in sorted(subtypes_list):
        subtypes_dict[subtype] = i
        i += 1
    # eprint(subtypes_dict)
    return subtypes_dict

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
        sg_dict = {}
        with open(args.groupfile) as sg_file:
            for line in sg_file:
               (key, val) = line.split()
               sg_dict[key] = val

        all_samples = vcf_reader.samples
        samples = sorted(list(set(sg_dict.values())))
        # eprint(samples)
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
        sg_dict = {}
        with open(args.groupfile) as sg_file:
            for line in sg_file:
               (key, val) = line.split()
               sg_dict[key] = val

        all_samples = vcf_reader.samples
        samples = sorted(list(set(sg_dict.values())))
        # eprint(samples)
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

            if ((acval==1 and record.FILTER is None) or args.novarfilter):
                # eprint(record.CHROM, record.POS, record.REF, record.ALT[0],
                    # acval, record.FILTER)
                # check and update chromosome sequence
                if record.CHROM != chrseq:
                    if args.verbose:
                        eprint("Loading chromosome", record.CHROM,
                            "reference...")

                    sequence = fasta_reader[record.CHROM]
                    chrseq = record.CHROM

                mu_type = record.REF + str(record.ALT[0])
                category = getCategory(mu_type)
                if nbp > 0:
                    lseq = sequence[record.POS-(nbp+1):record.POS+nbp].seq
                else:
                    lseq = sequence[record.POS-1].seq
                    # eprint("lseq:", lseq)
                motif_a = getMotif(record.POS, lseq)
                subtype = str(category + "-" + motif_a)
                st = subtypes_dict[subtype]
                # eprint(record.CHROM, record.POS,
                    # record.REF, record.ALT[0], subtype)

                # testing for speed improvements by batching updates to M
                asbatch = False
                if(asbatch and not args.groupfile):
                    batchit += 1
                    sample = record.gt_types.tolist().index(1)
                    sample_batch.append(sample)
                    subtype_batch.append(st)

                    if batchit == 10000:
                        M[sample_batch, subtype_batch] += 1
                        batchit = 0

                elif args.groupfile:
                    sample = all_samples[record.gt_types.tolist().index(1)]

                    if sample in sg_dict:
                        sample_gp = sg_dict[sample]
                        ind = samples.index(sample_gp)
                        M[ind,st] += 1
                else:
                    M[:,st] = M[:,st]+record.gt_types

                numsites_keep += 1

                if args.verbose:
                    if (numsites_keep%10000==0):
                        eprint("Processed", numsites_keep, "sites",
                            "(Skipped", numsites_skip, "sites)")
            else:
                numsites_skip += 1

    if numsites_keep == 0:
        eprint("No SNVs found. Please check your VCF file")
        sys.exit()

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
# aggregate M matrices from list of input files
###############################################################################
def aggregateM(args, subtypes_dict):
    colnames = ["ID"]
    M_colnames = colnames + list(sorted(subtypes_dict.keys()))
    colrange = range(1,len(M_colnames))

    if args.input.lower().endswith('nmf_m_spectra.txt'):
        samples = getSamples(args.input)
        M = np.loadtxt(args.input, skiprows=1, usecols=colrange)
        # M_it = np.concatenate((np.array([samples]).T, M_it), axis=1)
        # M_out = np.concatenate((M_out, M_it), axis=0)

        # M = np.delete(M_out, 0, 0)
        # M = np.delete(M, 0, 1)
        M = M.astype(np.float)
    else:
        with open(args.input) as f:
            file_list = f.read().splitlines()

        # M output by sample
        if args.input.lower().endswith('m_samples.txt'):

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
        elif args.input.lower().endswith('m_regions.txt'):

            samples = getSamples(file_list[0])

            # eprint(len(samples))
            M_out = np.zeros((len(samples), len(M_colnames)-1))
            # eprint(M_out.shape)
            for mfile in file_list:
                M_it = np.loadtxt(mfile, skiprows=1, usecols=colrange)
                M_out = np.add(M_out, M_it)

            M = M_out.astype(np.float)

    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
    return out

###############################################################################
# process tab-delimited text file, containing the following columns:
# CHR    POS    REF    ALT    SAMPLE_ID
###############################################################################
def aggregateTxt(args, subtypes_dict):
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
            # seed="nndsvd",
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
    M_fmt = np.concatenate((np.array([samples]).T, M.astype('|S20')), axis=1)

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
    W_fmt = np.concatenate((np.array([samples]).T, W.astype('|S20')), axis=1)
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
    rmse = open(rmse_path, "w")
    i = 0
    for val in np.nditer(M_rmse):
        # if val > 0.002:
        line = str(samples[i]) + "\t" + str(val) + "\n"
        rmse.write("%s" % line)
        i += 1
    rmse.close()
