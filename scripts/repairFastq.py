#!/usr/bin/env python3

import sys
import argparse

def read_fastq_record(fin):
    rec = []
    while len(rec) < 4:
        line = fin.readline().strip()
        if line == '':
            return None
        if len(rec) == 0:
            assert line.startswith('@')
        rec.append(line)
    return(rec)

parser = argparse.ArgumentParser(description='Filter fastq by keeping only read names found in reference fastq.\
        Typically, the reference is a first-in-pair fastq of reads matching a known barcode (or with reads selected according to some criteria).\
        The fastq to be filtered is the second-in-pair with all reads.\
        Currently, only uncompressed files can be read directly. Use stdin to read gzip files.')
parser.add_argument('--input-fastq', '-i', help='Fastq file to be filtered [%(default)s]', default='-')
parser.add_argument('--reference-fastq', '-r', help='Keep read names found in this fastq [%(default)s]', default='-')
parser.add_argument('--read-name-sep', '-s', default='#', help='Match read names after stripping everything after and including this character [%(default)s]')
parser.add_argument('--keep-original-names', '-k', action='store_true', help='Output read names as found in input fastq. By default, replace with the corresponding names in reference')
args = parser.parse_args()

if args.reference_fastq == '-':
    r1 = sys.stdin
else:
    r1 = open(args.reference_fastq)

if args.input_fastq == '-':
    r2 = sys.stdin
else:
    r2 = open(args.input_fastq)

while True:
    read1 = read_fastq_record(r1)
    if read1 is None:
        break
    name_r1 = read1[0].split(' ')[0].split(args.read_name_sep)[0]
    mate_found = False
    while mate_found is False:
        read2 = read_fastq_record(r2)
        if read2 is None:
            sys.stderr.write('Read:\n%s\nnot found in fastq to be filtered or the order of reads in the two files is not the same\n' % '\n'.join(read1))
            sys.exit(1)
        name_r2 = read2[0].split(' ')[0]
        if name_r1 == name_r2:
            # If the name of this read from R2 matches the name from R1,
            # print it and move to the next read from R1.
            # We also replace the read name in R2 to be the same as R1,
            # i.e. including everything after the sepaprator (UMI etc.)
            if args.keep_original_names == False:
                read2[0] = read1[0]
            print('\n'.join(read2))
            break
