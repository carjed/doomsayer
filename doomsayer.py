#!/usr/bin/python

from __future__ import print_function
import os
import sys
import textwrap
import argparse
import itertools
import timeit
import time
import numpy as np
from subprocess import call
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import nimfa
from util import *

###############################################################################
# Parse arguments
###############################################################################
start = timeit.default_timer()

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input",
                    help="input file, usually a VCF. Can accept input from \
                        STDIN  with \"--input -\". If using in aggregation \
                        mode, input should be a text file containing the file \
                        paths of the M matrices to aggregate",
                    required=True,
                    nargs='?',
                    type=str,
                    default=sys.stdin)

parser.add_argument("-f", "--fastafile",
                    help="fasta file name",
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
                        (Useful if analyzing somatic data or pre-filtering \
                        with another tool)",
                    action="store_true")

parser.add_argument("-d", "--diagnostics",
                    help="write the NMF matrices to the output directory \
                        and generate yaml config to be passed to \
                        diagnostic script",
                    action="store_true")

parser.add_argument("-a", "--autodiagnostics",
                    help="same as the --diagnostics option, but automatically \
                        generates diagnostic report",
                    action="store_true")

parser.add_argument("-m", "--mmatrixname",
                    help="custom filename for M matrix [without extension]",
                    nargs='?',
                    type=str,
                    default="NMF_M_spectra")

parser.add_argument("-s", "--samplefile",
                    help="file with sample IDs to include (one per line)",
                    nargs='?',
                    type=str)

parser.add_argument("-g", "--groupfile",
                    help="two-column tab-delimited file containing sample IDs \
                        (column 1) and group membership (column 2) for pooled \
                        analysis",
                    nargs='?',
                    type=str)

parser.add_argument("-b", "--subtypefile",
                    help="two-column tab-delimited file containing sample IDs \
                        (column 1) and group membership (column 2) for pooled \
                        analysis",
                    nargs='?',
                    type=str)

parser.add_argument("-ns", "--noscale",
                    help="do not scale H and W matrices",
                    action="store_true")

parser.add_argument("-t", "--threshold",
                    help="threshold for dropping samples, in standard \
                    deviations away from the mean signature contribution. \
                    The default is 2--lower values are more stringent",
                    type=int,
                    default=2)

parser.add_argument("-r", "--rank",
                    help="rank for NMF decomposition",
                    type=int,
                    choices=range(1,11),
                    default=0)

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
# Initialize project directory and index subtypes
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

eprint("indexing subtypes...") if args.verbose else None
subtypes_dict = indexSubtypes(args)

###############################################################################
# Check inputs
###############################################################################
# vcf
if(args.input.lower().endswith(('.vcf', '.vcf.gz')) or args.input == "-"):
    if args.verbose:
        eprint("Input detected as VCF file or VCF from STDIN")
    data = processVCF(args, subtypes_dict)
    M = data.M
    samples = data.samples
# testing--integrate per-chromosome parsing
# elif args.input.lower().endswith('inputs.txt'):
#     with open(args.input) as f:
#         file_list = f.read().splitlines()
#
#     i = 1
#     for vcf in file_list:
#         cmd = "python doomsayer.py" + \
#             " --input " + vcf + \
#             " --fastafile " + args.fastafile + \
#             " --projectdir " + args.projectdir + \
#             " --length " + str(args.length) + \
#             " --rank " + str(args.rank) + \
#             " --threshold " + str(args.threshold) + \
#             " --mmatrixname " + "NMF_" + str(i)
#         eprint("Running job:", cmd)
#         call(cmd + " &", shell=True)
#         i += 1
#
#     njobs = i
#     if args.verbose:
#         eprint("Waiting for " + str(njobs-1) + " subjobs to finish...")


elif(args.input.lower().endswith('m_samples.txt') or
        args.input.lower().endswith('m_regions.txt')):

    colnames = ["ID"]
    M_colnames = colnames + list(sorted(subtypes_dict.keys()))
    colrange = range(1,len(M_colnames))

    with open(args.input) as f:
        file_list = f.read().splitlines()

    # M output by sample
    if args.input.lower().endswith('m_samples.txt'):
        if args.verbose:
            eprint("Aggregating sample subset spectra matrices")

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
        if args.verbose:
            eprint("Aggregating regional subset spectra matrices")

        samples = getSamples(file_list[0])

        # eprint(len(samples))
        M_out = np.zeros((len(samples), len(M_colnames)-1))
        # eprint(M_out.shape)
        for mfile in file_list:
            M_it = np.loadtxt(mfile, skiprows=1, usecols=colrange)
            M_out = np.add(M_out, M_it)

        M = M_out.astype(np.float)
else:
    eprint("ERROR: invalid input detected. See documentation")
    sys.exit()

###############################################################################
# Run NMF on final matrix
###############################################################################
# M /= np.max(np.abs(M),axis=1)
# drop rows with all zeroes
# M = M[~(M==0).all(1)]

###############################
# M matrix (counts)
###############################
if args.mmatrixname != "NMF_M_spectra":
    if args.verbose:
        eprint("Saving M matrix (observed spectra counts)")

    colnames = ["ID"]
    M_colnames = colnames + list(sorted(subtypes_dict.keys()))

    # add ID as first column
    M_fmt = np.concatenate((np.array([samples]).T, M), axis=1)

    # add header
    M_fmt = np.concatenate((np.array([M_colnames]), M_fmt), axis=0)

    # write out
    M_path = projdir + "/" + args.mmatrixname + ".txt"
    np.savetxt(M_path, M_fmt, delimiter='\t', fmt="%s")
else:
    M_f = M/(M.sum(axis=1)+1e-8)[:,None]
    M_r = M

    # development option--pass file containing motif counts and
    # run NMF on relative rates
    if args.subtypefile:
        if args.verbose:
            eprint("Scaling M into relative rate matrix")
        st_dict = {}
        with open(args.subtypefile) as st_file:
            for line in st_file:
                (key, val) = line.split()
                st_dict[key] = int(val)
                M_r[:,subtypes_dict[key]] /= st_dict[key]
    # eprint(M)
    if args.noscale:
        M_run = M_r
    else:
        M_run = M_f

    if args.rank > 0:
        if args.verbose:
            eprint("Running NMF with specified rank=", args.rank)
        model = nimfa.Nmf(M_run,
            rank=args.rank,
            update="divergence",
            objective='div',
            n_run=1,
            max_iter=200)
        model_fit = model()
        evar = model_fit.fit.evar()
        maxind = args.rank
    else:
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
            maxind = i
            # if evar > 0.8:
            if(i > 2 and evar - evarprev < 0.001):
                if args.verbose:
                    eprint(textwrap.dedent("""\
                            Stopping condition met: <0.1 percent difference
                            in explained variation between ranks
                            """))
                maxind = i-1
                break
            # elif evar > 0.95:
            #     if args.verbose:
            #         eprint(textwrap.dedent("""\
            #                 Stopping condition met: rank explains >80 percent
            #                 of variation.
            #                 """))
            #     break
            evarprev = evar

    # if(maxind == 1 and evar > 0.95):
    #     stop = timeit.default_timer()
    #     tottime = round(stop - start, 2)
    #     if args.verbose:
    #         eprint(str(round(evar,2)*100) + \
    #             " percent of variance explained with 1 signature")
    #         eprint("Total runtime:", tottime, "seconds")
    #     sys.exit()
    # else:
    # maxind = evar_list.index(max(evar_list))+1
    # model = nimfa.Nmf(M_run, rank=maxind)
    # model_fit = model()
    W = model_fit.basis()
    H = model_fit.coef()
    # evar = model_fit.fit.evar()
    # eprint(evar)

    # eprint(H)
    # eprint(W)

    W_f = W
    W = W/np.sum(W, axis=1)
    # H = H/np.sum(H, axis=1)
    # W= W[~np.isnan(W).any(axis=1)]

    # output NMF results
    if(args.diagnostics or args.autodiagnostics):

        if args.verbose:
            colmeans = np.mean(W, axis=0)
            colstd = np.std(W, axis=0)
            upper = colmeans+args.threshold*colstd
            lower = colmeans-args.threshold*colstd

            eprint("Mean signature contributions: ", colmeans)
            eprint("StdDev:", colstd)
            eprint("Upper:", upper)
            eprint("Lower:", lower)
            eprint("Writing NMF results")

        diagWrite(projdir, M, M_run, W, H, subtypes_dict, samples, args)

###############################################################################
# Write keep and drop lists
###############################################################################
if(args.nofilter or args.mmatrixname != "NMF_M_spectra"):
    eprint(textwrap.dedent("""\
            You are running with the --nofilter or --mmatrixname option.
            Keep and drop lists will not be generated.
            """))
else:
    keep_samples = []
    drop_samples = []
    # Check for outliers
    i=0
    for n in W:
        if(np.greater(n, upper).any() or np.less(n, lower).any()):
            drop_samples.append(samples[i])
        else:
            keep_samples.append(samples[i])
        i += 1

    eprint("Printing keep and drop lists") if args.verbose else None
    keep_path = projdir + "/doomsayer_keep.txt"
    keeps = open(keep_path, "w")
    for sample in keep_samples:
        keeps.write("%s\n" % sample)
    # np.savetxt(keep_path, keep_samples, delimiter='\t', fmt="%s")

    drop_path = projdir + "/doomsayer_drop.txt"
    drops = open(drop_path, "w")
    for sample in drop_samples:
        drops.write("%s\n" % sample)
    # np.savetxt(drop_path, drop_samples, delimiter='\t', fmt="%s")

###############################################################################
# write output vcf
###############################################################################
if(args.outputtovcf and
        (args.input.lower().endswith(('.vcf', '.vcf.gz')) or
        args.input == "-")):
    eprint("Filtering VCF by drop list...") if args.verbose else None
    keep_test = keep_samples[0:10]
    vcf = VCF(args.input, samples=keep_test, mode='rb')
    # vcf = VCF(args.input, samples=keep_samples, mode='rb')

    print(vcf.raw_header.rstrip())
    for v in vcf:
        v.INFO['AC'] = str(v.num_het + v.num_hom_alt*2)

        if int(v.INFO['AC']) > 0:
            v.INFO['NS'] = str(v.num_called)
            v.INFO['AN'] = str(2*v.num_called)
            v.INFO['DP'] = str(np.sum(v.format('DP')))
            print(str(v).rstrip())

    vcf.close()

elif(args.outputtovcf and args.input.lower().endswith(('.txt'))):
    eprint(textwrap.dedent("""\
            WARNING: Doomsayer cannot write to a new VCF if running in
            aggregation mode. Please use the keep/drop lists to manually filter
            your VCF with bcftools or a similar utility
            """))

###############################################################################
# auto-generate diagnostic report in R
###############################################################################
if(args.autodiagnostics and args.mmatrixname == "NMF_M_spectra"):
    cmd = "diagnostics/doomsayer_diagnostics.r " + projdir + "/config.yaml"
    if args.verbose:
        eprint("Rscript will run the following command:")
        eprint("Rscript " + cmd)
        eprint("Auto-generating diagnostic report...")
    call("/usr/bin/Rscript --vanilla " + cmd, shell=True)

stop = timeit.default_timer()
tottime = round(stop - start, 2)
eprint("Total runtime:", tottime, "seconds") if args.verbose else None
