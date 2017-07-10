#!/usr/bin/python

from __future__ import print_function
import os
import sys
import argparse
import itertools
import timeit
import collections
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
    # get 3-mer motif
    # motif = Seq(sequence[pos-2:pos+1].seq, IUPAC.unambiguous_dna)
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

                subtype = category + "-" \
                    + kmerstr[0:flank] + ref + kmerstr[flank:(motiflength-1)]

                subtypes_list.append(subtype)
    else:
        ext = ["-A", "-C"]
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
# Main function for parsing VCF
###############################################################################
def processVCF(args, subtypes_dict):
    eprint("Initializing reference genome...") if args.verbose else None
    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)

    # 'demo/input/keep.txt'
    if args.samplefile:
        with open(args.samplefile) as f:
            keep_samples = f.read().splitlines()
        eprint(len(keep_samples), "samples kept") if args.verbose else None

        vcf_reader = VCF(args.input,
            mode='rb', gts012=True, lazy=True, samples=keep_samples)
        # vcf_reader.set_samples(keep_samples) # <- set_samples() subsets VCF
    else:
        vcf_reader = VCF(args.input,
            mode='rb', gts012=True, lazy=True)

    nbp = (args.length-1)//2

    ####################
    # index samples
    ####################
    eprint("Indexing samples in", args.input, "...") if args.verbose else None

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

    ############################################
    # Query records in VCF and build matrix
    ############################################
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

            if ((acval==1 and record.FILTER is None) or args.nofilter):
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

                # use quick singleton lookup for default QC option
                # if not args.nofilter:
                    # sample=samples[record.gt_types.tolist().index(1)]
                    # if sample == "1497-RMM-1269RD":
                    #     print(record.CHROM, record.POS,
                    #         record.REF, record.ALT[0], sample, subtype)
                    # eprint(record.gt_types.tolist())
                    # sample = record.gt_types.tolist().index(1)
                    # sample=np.where(record.gt_types == 1)[0]
                    # eprint(sample)

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
    return out

###############################################################################
# prepar data for diagnostics
###############################################################################
def diagWrite(projdir, M, M_f, W, H, subtypes_dict, samples, args):

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
    M_path = projdir + "/" + args.mmatrixname + ".txt"
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
    np.savetxt(M_path_rates, M_fmt, delimiter='\t', fmt="%s")

    ###############################
    # W matrix (contributions)
    ###############################
    if args.verbose:
        eprint("Saving W matrix (signature contributions per sample)")

    # add ID as first column
    W_fmt = np.concatenate((np.array([samples]).T, W), axis=1)

    # add header
    W_colnames = colnames + ["S" + str(i) for i in range(1,args.rank+1)]
    W_fmt = np.concatenate((np.array([W_colnames]), W_fmt), axis=0)

    # write out
    W_path = projdir + "/NMF_W_sig_contribs.txt"
    np.savetxt(W_path, W_fmt, delimiter='\t', fmt="%s")

    ###############################
    # H matrix (loadings)
    ###############################
    if args.verbose:
        eprint("Saving H matrix (feature loadings per signature)")

    # add signature ID as first column
    H_rownames = ["S" + str(i) for i in range(1,args.rank+1)]
    H_fmt = np.concatenate((np.array([H_rownames]).T,
        np.char.mod('%d', H)), axis=1)

    eprint(H_fmt)
    # add header
    H_colnames = ["Sig"] + list(sorted(subtypes_dict.keys()))
    H_fmt = np.concatenate((np.array([H_colnames]), H_fmt), axis=0)
    # eprint(H_fmt)
    # H_fmt = H
    # eprint(H_fmt)
    # write out
    H_path = projdir + "/NMF_H_sig_loads.txt"
    np.savetxt(H_path, H_fmt, delimiter='\t', fmt="%s")

    yaml = open(projdir + "/config.yaml","w+")
    print("# Config file for doomsayer_diagnostics.r", file=yaml)
    print("keep_path: " + projdir + "/doomsayer_keep.txt", file=yaml)
    print("drop_path: " + projdir + "/doomsayer_drop.txt", file=yaml)
    print("M_path: " + M_path, file=yaml)
    print("M_path_rates: " + M_path_rates, file=yaml)
    print("W_path: " + W_path, file=yaml)
    print("H_path: " + H_path, file=yaml)
    yaml.close()
