#!/usr/bin/python

from __future__ import print_function
import os
import sys
import argparse
import itertools
import timeit
import numpy as np
from pyfaidx import Fasta
# from Bio import SeqIO
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

    m1 = motif[1]
    m2 = altmotif[1]

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
    return subtypes_dict

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
    M_path = projdir + "/NMF_M_spectra.txt"
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

    # add signature ID as firsat column
    H_rownames = ["S" + str(i) for i in range(1,args.rank+1)]
    H_fmt = np.concatenate((np.array([H_rownames]).T, H), axis=1)

    # add header
    H_colnames = ["Sig"] + list(sorted(subtypes_dict.keys()))
    H_fmt = np.concatenate((np.array([H_colnames]), H_fmt), axis=0)

    # write out
    H_path = projdir + "/NMF_H_sig_loads.txt"
    np.savetxt(H_path, H_fmt, delimiter='\t', fmt="%s")

    yaml = open(projdir + "/config.yaml","w+")
    print("# Config file for doomsayer_diagnostics.r", file=yaml)
    print("keep_path: " + projdir + "/keep_samples.txt", file=yaml)
    print("drop_path: " + projdir + "/drop_samples.txt", file=yaml)
    print("M_path: " + M_path, file=yaml)
    print("M_path_rates: " + M_path_rates, file=yaml)
    print("W_path: " + W_path, file=yaml)
    print("H_path: " + H_path, file=yaml)
    yaml.close()
