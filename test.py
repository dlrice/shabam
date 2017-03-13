#!/usr/bin/env python

from shabam import plot

plot('test/example.bam',
    fastafile='test/human_g1k_v37.fasta',
    chrom=16,
    start=48000000,
    end=48000080,
    out='plot.png'
)