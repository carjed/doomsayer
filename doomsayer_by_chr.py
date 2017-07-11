#!/usr/bin/python

from __future__ import print_function
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input",
                    help="input file containing the file paths of VCFs to run",
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

projdir = os.path.realpath(args.projectdir)
regionfile = projdir + "/m_regions.txt"
mfile_list = [projdir + "/NMF_" + str(s) + ".txt" for s in range(1,njobs+1)]

myfile = open(regionfile, 'w')
for mfile in mfile_list:
    myfile.write("%s\n" % mfile)

myfile.close()

with open(args.input) as f:
    file_list = f.read().splitlines()

i = 1
for vcf in file_list:
    cmd = "python doomsayer.py" +
        " --input " + vcf +
        " --fastafile " + args.fastafile +
        " --projectdir " + args.projectdir +
        " --length " + args.length +
        " --rank " + args.length +
        " --threshold " + args.threshold +
        " --mmatrixname " + "NMF_" + i +
		" --autodiagnostics --verbose"
    print("Running job:", cmd)
    call(cmd + " &", shell=True)
    i += 1

print("Waiting for subjobs to finish")
njobs = i
mfile_list = [projdir + "/NMF_" + str(s) + ".txt" for s in range(1,njobs+1)]
mfile_check = [False for i in range(njobs)]
mfiles_exist = False
while not mfiles_exist:
    i = 0
    for mfile in mfile_list:
        mfile_check[i] = os.path.exists(mfile)
        i += 1
    if all(mfile_check):
        mfiles_exist = True
    else:
        time.sleep(5)

aggcmd = "python doomsayer.py" +
    " --input " + regionfile +
    " --projectdir " + args.projectdir +
	" --autodiagnostics --verbose"

print("Running aggregation script")
call(aggcmd, shell=True)
