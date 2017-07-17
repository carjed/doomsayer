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

parser.add_argument("-b", "--baseline",
                    help="Transform individual mutation spectra relative to \
                    baseline average (mean across samples)",
                    action="store_true")

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
                    choices=range(2,11),
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
# build M matrix if vcf
if(args.input.lower().endswith(('.vcf', '.vcf.gz')) or args.input == "-"):
    if args.verbose:
        eprint("Input detected as VCF file or VCF from STDIN")
    data = processVCF(args, subtypes_dict)
    M = data.M
    samples = data.samples

# aggregate if M matrix file list
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
if args.mmatrixname != "NMF_M_spectra":
    if args.verbose:
        eprint("Saving M matrix (observed spectra counts)")

    colnames = ["ID"]
    M_colnames = colnames + list(sorted(subtypes_dict.keys()))

    # add ID as first column
    M_out = np.concatenate((np.array([samples]).T, M), axis=1)

    # add header
    M_out = np.concatenate((np.array([M_colnames]), M_out), axis=0)

    # write out
    M_path = projdir + "/" + args.mmatrixname + ".txt"
    np.savetxt(M_path, M_out, delimiter='\t', fmt="%s")
else:
    M_f = M/(M.sum(axis=1)+1e-8)[:,None]

    # eprint(M)
    if args.noscale:
        M_run = M
    else:
        M_run = M_f

    if args.verbose:
        eprint("Generating baseline signature")
    base_model = nimfa.Nmf(M_run,
        rank=1,
        update="divergence",
        objective='div',
        n_run=1,
        max_iter=200)

    base_model_fit = base_model()
    base_H = base_model_fit.coef()
    base_H = np.divide(base_H, np.sum(base_H))
    # base_H = abs(np.subtract(base_H, np.sum(base_H)))

    M_rmse = np.square(np.subtract(M_run, base_H))
    M_rmse = np.sqrt(M_rmse.sum(axis=1)/M.shape[1])

    if args.baseline:
        M_run = np.divide(M_run, base_H)

    if args.rank > 0:
        if args.verbose:
            eprint("Running NMF with specified rank =", args.rank)
        model = nimfa.Nmf(M_run,
            rank=args.rank,
            update="divergence",
            objective='div',
            n_run=1,
            max_iter=500)
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
                n_run=5,
                max_iter=500)
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
    W_f = W
    W = W/np.sum(W, axis=1)
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

        diagWrite(projdir, M, M_run, M_rmse, W, H, subtypes_dict, samples, args)

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
    keeps.close()
    # np.savetxt(keep_path, keep_samples, delimiter='\t', fmt="%s")

    drop_path = projdir + "/doomsayer_drop.txt"
    drops = open(drop_path, "w")
    for sample in drop_samples:
        drops.write("%s\n" % sample)
    drops.close()
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
