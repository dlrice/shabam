#!/usr/bin/env python

import argparse

from shabam.plot import seqplot

def main():
    parser = argparse.ArgumentParser(description='A python/command tool to ' \
        'create sequence plots from bam/cram files.')
    parser.add_argument('--seqfiles', type=str,  nargs='+', required=True,
        help='BAM(s) or CRAM(s) to plot.')
    parser.add_argument('--fastafile', type=str, required=True,
        help='A reference FASTA file')
    parser.add_argument('--chrom', type=str, required=True,
        help='Chromosome')
    parser.add_argument('--start', type=int, required=True,
        help='Start base of plot')
    parser.add_argument('--end', type=int, required=True,
        help='End base of plot')
    parser.add_argument('--by-strand', default=False, action='store_true',
        help='whether to color reads by strand (default is not)')
    parser.add_argument('--out', type=str, required=True,
        help='Output file (extension determines type: png, pdf, jpg, etc.)')
    
    args = parser.parse_args()
    
    seqplot(args.seqfiles, args.fastafile, args.chrom, args.start, args.end,
        args.out, args.by_strand)

if __name__ == '__main__':
    main()
