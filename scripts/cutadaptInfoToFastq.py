#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description= 'Prepare fastq from cutadapt info file suitable for BLISS analysis')

parser.add_argument('--info', '-i', help= 'Input, tipycally the info file from cutadapt %(default)s', default= '-')
parser.add_argument('--umi-position', '-u', help= "Start (0-based) and end (1-based) position of UMI %(default)s", 
    default= [0, 8], type= int, nargs= 2, metavar= 'UMI')
parser.add_argument('--sep', '-s', help= 'Separator to use between read name and UMI [%(default)s]', default= '#')
parser.add_argument('--rname-index', help= "1-based column index of read name [%(default)s]", default= 1, type= int)
parser.add_argument('--umi-index', help= "1-based column index of UMI+barcode sequence [%(default)s]", default= 6, type= int)
parser.add_argument('--seq-index', help= "1-based column index of trimmed read sequence [%(default)s]", default= 7, type= int)
parser.add_argument('--qual-index', help= "1-based column index of trimmed read quality [%(default)s]", default= 11, type= int)

args = parser.parse_args()

if args.info == '-':
    fin = sys.stdin
else:
    fin= open(args.info)

for line in fin:
    line = line.strip().split('\t')
    if len(line) < (args.umi_index+1) or line[args.umi_index - 1] == '':
        # UMI-Barcode not found
        continue
    umi = line[args.umi_index- 1][args.umi_position[0] : args.umi_position[1]]
    rname = line[args.rname_index - 1].split()[0]
    assert args.sep not in rname
    rname = rname + args.sep + umi
    seq = line[args.seq_index - 1]
    qual = line[args.qual_index - 1]
    assert len(seq) == len(qual)
    out = '\n'.join(['@' + rname, seq, '+', qual])
    print(out)
