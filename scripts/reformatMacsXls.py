#!/usr/bin/env python3

import sys
import argparse
import pandas
import os

VERSION = '0.1.0'

parser = argparse.ArgumentParser(description= 'Reformat macs2 output format to BED-like')

parser.add_argument('input', help= 'Input *.xls file from macs2')


if __name__ == '__main__':
    args = parser.parse_args()

    with open(args.input) as fin:
        for line in fin:
            line = line.strip()
            if line.startswith('#'):
                line = line.lstrip('#')
                print('##' +  line)
            elif line == '':
                continue
            else:
                print('## Reformatted by %s v%s' % (os.path.basename(__file__), VERSION))
                break

    xls = pandas.read_csv(args.input, sep= '\t', comment= '#', dtype= {0: 'string'})
    xls.rename(columns= {'chr': '#chrom'}, inplace= True)
    xls['start'] = xls['start'] - 1
    xls['abs_summit'] = xls['abs_summit'] - 1
    
    csv = xls.to_csv(None, sep= '\t', index= False)
    if not csv.endswith('\n'):
        csv += '\n'
    sys.stdout.write(csv)

