#!/usr/bin/env python
# import fire
import pysam
import numpy as np
import matplotlib.pyplot as plt
import os.path
import argparse


ENDCHAR = "-"
GAPCHAR = "."
DELCHAR = "*"

BAM_CMATCH = 0 # M
BAM_CINS = 1 # I
BAM_CDEL = 2 # D
BAM_CREF_SKIP = 3 # N
BAM_CSOFT_CLIP = 4 # S
BAM_CHARD_CLIP = 5 # H
BAM_CPAD = 6 # P
BAM_CEQUAL = 7 # =
BAM_CDIFF = 8 # X

# BASE2COLORS = {
#     'A': 'green', 
#     'C': 'blue', 
#     'G': 'orange', 
#     'T': 'red',
#     '-' : 'black',       # deletion
#     'I' : 'mediumpurple' # insertion
# }

BASE2COLORS = {
    'A' : np.array([0,110,0,255], dtype=np.uint8),
    'C' : np.array([0,0,255,255], dtype=np.uint8),
    'G' : np.array([255,150,50,255], dtype=np.uint8),
    'T' : np.array([255,0,0,255], dtype=np.uint8),
    'M' : np.array([232,232,232,255], dtype=np.uint8), # match
    '-' : 'pink',        # deletion
    'I' : 'mediumpurple' # insertion
}

def parseCigar(cigar, bases):
    """
    Return list of strings, each item corresponding to a single reference position
    """
    rep = []
    currentpos = 0
    wasinsert = False
    for operation, length in cigar:
        if operation in [BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF]: # match (M, X, =)
            if operation == BAM_CDIFF:
                print(operation)
            for i in range(length):
                if wasinsert:
                    rep[-1] = rep[-1] + bases[currentpos]
                else:
                    rep.append(bases[currentpos])
                wasinsert = False
                currentpos += 1
        elif operation == BAM_CINS: # put insertion in next base position (I)
            if wasinsert:
                rep[-1] = rep[-1] + nucs[currentpos:currentpos + length]
            else:
                rep.append(nucs[currentpos:currentpos + length])
            currentpos = currentpos + length
            wasinsert = True
        elif operation in [BAM_CDEL, BAM_CREF_SKIP]: # deletion (D) or skipped region from the reference (N)
            for i in range(length):
                if wasinsert:
                    rep[-1] = rep[-1] + DELCHAR
                else:
                    rep.append(DELCHAR)
                wasinsert = False
        elif operation == BAM_CPAD: # padding (silent deletion from padded reference) (P)
            if wasinsert:
                rep[-1] = rep[-1] + DELCHAR * length
            else:
                rep.append(DELCHAR * length)
            wasinsert = True
        elif operation not in [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]: # hard clipping or soft clipping
            sys.stderr.write("ERROR: Invalid CIGAR operation (%s) in read %s \n"%(operation, read.qname))
    return rep


def generateRGB(representations, plot_start, plot_end, reference):
    plot_length = plot_end - plot_start
    nrows = len(representations) + 1
    RGB = np.ones((nrows, plot_length, 4), dtype=np.uint8)*255
    
    # Plot reference
    for i, base in enumerate(reference):
            RGB[0, i, :] = BASE2COLORS[base]
    
    for row_index, representation in enumerate(representations):
        row_index += 1 # For reference
        bases = representation['bases']
        quals = representation['qualities']
        read_start = representation['position']
        read_end = read_start + len(bases) - 1
        start = max(plot_start, read_start)
        end = min(plot_end, read_end)
        for i in range(start, end):
            read_index = i - read_start
            plot_index = i - plot_start
            base = bases[read_index]
            qual = quals[read_index]
            if reference[plot_index] == base:
                base = 'M'
            color = BASE2COLORS[base]
            color[3] = to_alpha(qual)
            RGB[row_index, plot_index, :] = color

    return RGB

def to_alpha(qual):
    return np.floor(255 * (min(35, qual)/35))

def getRepresentations(reads):
    representations = []
    for read in reads:
        position = read.pos
        bases = read.query
        quals = read.query_qualities
        cigar = read.cigartuples
    #     if read.is_reverse:
    #         bases = bases.lower()
        representations.append({
            'position': position,
            'bases': parseCigar(cigar, bases),
            'qualities': quals,
        })
    return representations


def plot(seqfile, fastafile, chrom, start, end, out):
    seq = pysam.AlignmentFile(seqfile, 'rb')
    fasta = pysam.FastaFile(fastafile)
    reads = seq.fetch("16", start, end)
    representations = getRepresentations(reads)
    ref = fasta.fetch(start=start, end=end, region=str(chrom))
    RGB = generateRGB(representations, start, end, ref)

    fig, (ax) = plt.subplots(figsize=(RGB.shape[1]/10,RGB.shape[0]/10 + 5))
    ax.imshow(RGB)

    #Spacing between each line
    # plt.grid(b=True, which='minor', color='w',linestyle='-')

    xticks = np.arange(0.5,end - start + 1, 1)
    # ax.set_xticks(xticks)
    # major_labels = [str(start + d) for d in xticks]
    ax.set_xticks(xticks, minor=True)

    yticks = np.arange(0.5, len(representations))
    ax.set_yticks(yticks, minor=True)
  
    xticks=np.array(ax.get_xticks().tolist(), dtype=int) + start
    ax.set_xticklabels(xticks, rotation=90, ha='left')
    # plt.axhline(y=0.45, linewidth=1, color = 'k')
    plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off', labeltop='on')
    
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.tight_layout()
    _, ext = os.path.splitext(out)
    if not ext:
        ext = '.png'
    ext = ext[1:] # Cut off the period.
    fig.savefig(out, format=ext, dpi=200)

    return ax


# class Shabam(object):
#     """A python/command tool to create sequence plots from bam/cram files."""
#     def plot(self, seqfile, fastafile, chrom, start, end, out):
#         plot(seqfile, fastafile, chrom, start, end, out)


def main():
    parser = argparse.ArgumentParser(description='A python/command tool to create sequence plots from bam/cram files.')
    parser.add_argument('--seqfile', type=str, required=True,
        help='BAM or CRAM to plot')
    parser.add_argument('--fastafile', type=str, required=True,
        help='A reference FASTA file')
    parser.add_argument('--chrom', type=str, required=True,
        help='Chromosome')
    parser.add_argument('--start', type=int, required=True,
        help='Start base of plot')
    parser.add_argument('--end', type=int, required=True,
        help='End base of plot')
    parser.add_argument('--out', type=str, required=True,
        help='Output file (extension determines type: png, pdf, jpg, etc.)')

    args = parser.parse_args()
    seqfile = args.seqfile
    fastafile = args.fastafile
    chrom = args.chrom
    start = args.start
    end = args.end
    out = args.out

    plot(seqfile, fastafile, chrom, start, end, out)




if __name__ == '__main__':
    main()

    # fire.Fire(Shabam)
