#!/usr/bin/python

from __future__ import print_function
import os
import sys

sys.path.append(os.getcwd())

import shutil
# from _version import __version__
import textwrap
import argparse
import itertools
import timeit
import time
import multiprocessing
import numpy as np
import cyvcf2 as vcf
from cyvcf2 import VCF
from cyvcf2 import Writer
from joblib import Parallel, delayed
from subprocess import call
from util import *

###############################################################################
# Parse arguments
###############################################################################
start = timeit.default_timer()

num_cores = multiprocessing.cpu_count()

parser = argparse.ArgumentParser()

mode_opts = ["vcf", "agg", "txt"]
parser.add_argument("-M", "--mode",
                    help="Mode for parsing input. Must be one of \
                        {"+", ".join(mode_opts)+ "}",
                    nargs='?',
                    type=str,
                    choices=mode_opts,
                    metavar='',
                    default="vcf")

parser.add_argument("-i", "--input",
                    help="In VCF mode (default) input file is a VCF \
                        or text file containing paths of multiple VCFs. \
                        Can accept input from STDIN  with \"--input -\". \
                        In aggregation mode, input file is a text file \
                        containing mutation subtype count matrices, \
                        or paths of multiple such matrices. \
                        In plain text mode, input file is tab-delimited text \
                        file containing 5 columns: CHR, POS, REF, ALT, ID",
                    required=True,
                    nargs='?',
                    type=str,
                    # metavar='',
                    default=sys.stdin)

parser.add_argument("-f", "--fastafile",
                    help="reference fasta file",
                    nargs='?',
                    type=str,
                    metavar='',
                    default="chr20.fasta.gz")

parser.add_argument("-s", "--samplefile",
                    help="file with sample IDs to include (one per line)",
                    nargs='?',
                    metavar='',
                    type=str)

parser.add_argument("-g", "--groupfile",
                    help="two-column tab-delimited file containing sample IDs \
                        (column 1) and group membership (column 2) for pooled \
                        analysis",
                    nargs='?',
                    metavar='',
                    type=str)

parser.add_argument("-p", "--projectdir",
                    help="directory to store output files \
                        (do NOT include a trailing '/')",
                    nargs='?',
                    type=str,
                    metavar='',
                    default="doomsayer_output")

parser.add_argument("-m", "--matrixname",
                    help="custom filename for M matrix [without extension]",
                    nargs='?',
                    type=str,
                    metavar='',
                    default="NMF_M_spectra")

parser.add_argument("-o", "--filterout",
                    help="in VCF or plain text modes, re-reads input \
                        file and writes to STDOUT, omitting records that occur \
                        in the detected outliers. To write to a new file, use \
                        standard output redirection [ > out.vcf] at the end of \
                        the doomsayer.py command",
                    action="store_true")

parser.add_argument("-a", "--allsamples",
                    help="disables generation of keep/drop lists. \
                        Forces NMF to run on the entire sample",
                    action="store_true")

parser.add_argument("-n", "--novarfilter",
                    help="turns off default variant filtering criteria \
                        and evaluates all sites in the input VCF. \
                        (Useful if analyzing somatic data or pre-filtering \
                        with another tool)",
                    action="store_true")

parser.add_argument("-R", "--report",
                    help="automatically generates an HTML-formatted report in \
                        R.",
                    action="store_true")

template_opts = ["diagnostics", "msa"]

parser.add_argument("-T", "--template",
                    help="Template for diagnostic report. Must be one of \
                    {"+", ".join(template_opts)+"}",
                    nargs='?',
                    type=str,
                    choices=template_opts,
                    metavar='',
                    default="diagnostics")

parser.add_argument("-t", "--threshold",
                    help="threshold for fold-difference RMSE cutoff, used to \
                        determine which samples are outliers. Must be a \
                        real-valued number > 1. The default is 2. \
                        higher values are more stringent",
                    nargs='?',
                    type=restricted_float,
                    metavar='',
                    default=2)

rank_opts = range(2,11)
ro_str = str(min(rank_opts)) + " and " + str(max(rank_opts))
parser.add_argument("-r", "--rank",
                    help="Rank for NMF decomposition. Must be an integer value \
                        between " + ro_str,
                    nargs='?',
                    type=int,
                    choices=rank_opts,
                    metavar='',
                    default=0)

motif_length_opts = [1,3,5,7]
mlo_str = ",".join(str(x) for x in motif_length_opts)

parser.add_argument("-l", "--length",
                    help="motif length. Allowed values are " + mlo_str,
                    nargs='?',
                    type=int,
                    choices=motif_length_opts,
                    metavar='',
                    default=3)

parser.add_argument("-c", "--cpus",
                    help="number of CPUs. Must be integer value between 1 \
                        and "+str(num_cores),
                    nargs='?',
                    type=int,
                    choices=range(1,num_cores+1),
                    metavar='',
                    default=1)

parser.add_argument("-v", "--verbose",
                    help="Enable verbose logging",
                    action="store_true")

args = parser.parse_args()

###############################################################################
# Initialize project directory
###############################################################################
projdir = os.path.realpath(args.projectdir)
if args.verbose:
    eprint("checking if directory", projdir, "exists...")

if not os.path.exists(args.projectdir):
    if args.verbose:
        eprint("Creating output directory:", projdir)
    os.makedirs(args.projectdir)
else:
    if args.verbose:
        eprint(projdir, "already exists")

if args.verbose:
    eprint("All output files will be located in ", projdir)

###############################################################################
# index subtypes
###############################################################################
eprint("indexing subtypes...") if args.verbose else None
subtypes_dict = indexSubtypes(args)

###############################################################################
# Build M matrix from inputs
###############################################################################
if args.mode == "vcf":
    eprint("Initializing reference genome...") if args.verbose else None
    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)

    if(args.input.lower().endswith(('.vcf', '.vcf.gz', '.bcf')) or
            args.input == "-"):
        par = False

        data = processVCF(args, args.input, fasta_reader, subtypes_dict, par)
        M = data.M
        samples = data.samples

    elif(args.input.lower().endswith(('.txt'))):
        par = True
        with open(args.input) as f:
            vcf_list = f.read().splitlines()
        results = Parallel(n_jobs=args.cpus) \
            (delayed(processVCF)(args, vcf, fasta_reader, subtypes_dict, par) for vcf in vcf_list)
        # eprint(results)
        # eprint(results[1].shape)
        nrow, ncol = results[1].shape
        M = np.zeros((nrow, ncol))

        for M_sub in results:
            M = np.add(M, M_sub)
        # data.M = M_comb
        eprint(M)
        samples = getSamplesVCF(args, vcf_list[1])

elif args.mode == "agg":
    data = aggregateM(args, subtypes_dict)
    M = data.M
    samples = data.samples

elif args.mode == "txt":
    data = aggregateTxt(args, subtypes_dict)
    M = data.M
    samples = data.samples

###############################################################################
# Write out M matrix if preparing for aggregation mode
###############################################################################
if args.matrixname != "NMF_M_spectra":
    eprint(textwrap.dedent("""\
            You are running with the --matrixname option. Keep and drop lists
            will not be generated.
            """))
    if args.verbose:
        eprint("Saving M matrix (spectra counts) to:", args.matrixname)

    M_path = projdir + "/" + args.matrixname + ".txt"
    writeM(M, M_path, subtypes_dict, samples)

###############################################################################
# Process M matrix
###############################################################################
else:
    # M_f is the relative contribution of each subtype per sample
    M_f = M/(M.sum(axis=1)+1e-8)[:,None]

    # M_err is N x K matrix of residual error profiles, used for RMSE calc
    M_err = np.subtract(M_f, np.mean(M_f, axis=0))
    M_rmse = np.sqrt(np.sum(np.square(M_err), axis=1)/M_err.shape[1])

    eprint("Writing M matrix and RMSE per sample") if args.verbose else None
    M_path = projdir + "/" + args.matrixname + ".txt"
    writeM(M, M_path, subtypes_dict, samples)

    M_path_rates = projdir + "/NMF_M_spectra_rates.txt"
    writeM(M_f, M_path_rates, subtypes_dict, samples)

    rmse_path = projdir + "/doomsayer_rmse.txt"
    writeRMSE(M_rmse, rmse_path, samples)

    if args.allsamples:
        if args.verbose:
            eprint("Using all samples--\
                keep and drop lists will not be generated")
        M_run = M_f
        samples = samples
    else:
        eprint("Printing keep and drop lists") if args.verbose else None
        M_err_d = np.divide(M_f, np.mean(M_f, axis=0))
        keep_samples = []
        drop_samples = []
        drop_indices = []
        i=0
        for row in M_err_d:
            if any(err > args.threshold for err in row):
            # if n > args.threshold:
            # if(np.greater(n, upper).any() or np.less(n, lower).any()):
                drop_samples.append(samples[i])
                drop_indices.append(i)
            else:
                keep_samples.append(samples[i])
            i += 1

        keep_path = projdir + "/doomsayer_keep.txt"
        keeps = open(keep_path, "w")
        for sample in keep_samples:
            keeps.write("%s\n" % sample)
        keeps.close()

        drop_path = projdir + "/doomsayer_drop.txt"
        drops = open(drop_path, "w")
        for sample in drop_samples:
            drops.write("%s\n" % sample)
        drops.close()

        if args.verbose:
            eprint(len(drop_samples), "potential outliers found.")

        M_run = M_f[np.array(drop_indices)]
        samples = drop_samples

    eprint("Running NMF model") if args.verbose else None
    NMFdata = NMFRun(M_run, args, projdir, samples, subtypes_dict)

    # W matrix (contributions)
    W_path = projdir + "/NMF_W_sig_contribs.txt"
    writeW(NMFdata.W, W_path, samples)

    # H matrix (loadings)
    H_path = projdir + "/NMF_H_sig_loads.txt"
    writeH(NMFdata.H, H_path, subtypes_dict)

    yaml = open(projdir + "/config.yaml","w+")
    print("# Config file for doomsayer_diagnostics.r", file=yaml)
    print("keep_path: " + projdir + "/doomsayer_keep.txt", file=yaml)
    print("drop_path: " + projdir + "/doomsayer_drop.txt", file=yaml)
    print("M_path: " + M_path, file=yaml)
    print("M_path_rates: " + M_path_rates, file=yaml)
    print("W_path: " + W_path, file=yaml)
    print("H_path: " + H_path, file=yaml)
    print("RMSE_path: " + rmse_path, file=yaml)
    yaml.close()

###############################################################################
# write output vcf
###############################################################################
if args.filterout:
    if args.mode == "vcf":
        eprint("Filtering input by drop list...") if args.verbose else None
        # keep_test = keep_samples[0:10]
        # vcf = VCF(args.input, samples=keep_test, mode='rb')
        vcf = VCF(args.input, samples=keep_samples, mode='rb')

        print(vcf.raw_header.rstrip())
        for v in vcf:
            v.INFO['AC'] = str(v.num_het + v.num_hom_alt*2)

            if int(v.INFO['AC']) > 0:
                v.INFO['NS'] = str(v.num_called)
                v.INFO['AN'] = str(2*v.num_called)
                v.INFO['DP'] = str(np.sum(v.format('DP')))
                print(str(v).rstrip())

        vcf.close()

    elif args.mode =="txt":
        eprint("Filtering input by drop list...") if args.verbose else None
        with open(args.input, 'r') as f:
            reader = csv.reader(f, delimiter='\t')

            for row in reader:
                chrom = row[0]
                pos = row[1]
                ref = row[2]
                alt = row[3]
                sample = row[4]

                if sample not in drop_samples:
                    print("\t".join(row))

    # elif(args.outputtovcf and args.input.lower().endswith(('.txt'))):
    elif args.mode == "agg":
        eprint(textwrap.dedent("""\
                WARNING: Doomsayer cannot write to a new VCF if running in
                aggregation mode. Please use the keep/drop lists to manually
                filter your VCF with bcftools or a similar utility
                """))

###############################################################################
# auto-generate diagnostic report in R
###############################################################################
if(args.report and args.matrixname == "NMF_M_spectra"):
    cmd_str = "Rscript --vanilla generate_report.r "
    param_str = projdir + "/config.yaml"

    shutil.copy("report_templates/" + args.template + ".Rmd",
        projdir + "/report.Rmd")

    cmd = cmd_str + param_str
    if args.verbose:
        eprint("Rscript will run the following command:")
        eprint(cmd)
    call(cmd, shell=True)

stop = timeit.default_timer()
tottime = round(stop - start, 2)
eprint("Total runtime:", tottime, "seconds") if args.verbose else None
