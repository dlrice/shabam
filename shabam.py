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
    'M' : np.array([232,232,232,255], dtype=np.uint8), # match
    'M_f' : np.array([80,180,255,255], dtype=np.uint8), # match, forward color
    'M_r' : np.array([255,180,80,255], dtype=np.uint8), # match, reverse color
    '-' : 'pink',        # deletion
    'I' : 'mediumpurple' # insertion
}

def parseCigar(cigar, bases):
    """
    Return list of strings, each item corresponding to a single reference position

    Initial code lifted - with permission - from https://github.com/mgymrek/pybamview
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


def generateRGB(representations, plot_start, plot_end, reference, by_strand=False,
        use_ref=True):
    plot_length = plot_end - plot_start
    nrows = len(representations) + 1
    RGB = np.ones((nrows, plot_length, 4), dtype=np.uint8)*255

    # Plot reference
    if use_ref:
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
                if by_strand:
                    strand = {True: 'r', False: 'f'}[representation['is_reverse']]
                    base = 'M_{}'.format(strand)
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
        if not cigar:
            continue
        representations.append({
            'position': position,
            'bases': parseCigar(cigar, bases),
            'is_reverse': read.is_reverse,
            'qualities': quals,
        })
    return representations

def combine_bams(data):
    ''' combine data for one or more bams into a single array
    
    Args:
        data: list of numpy.ndarrays for bams
    
    Returns:
        numpy.ndarray, where the image data for the bams have been combined. We
        add gaps between successive bams by including white space, and a line.
    '''
    
    # define gap between BAMs as white space with a black line in the middle
    height, width = 2, data[0].shape[1]
    white = np.ones((height, width, 4), dtype=np.uint8) * 255
    line = np.ndarray((1, width, 4), np.uint8, np.array([[[0, 0, 0, 255]] * width]))
    spacer = np.concatenate((white, line, white))
    
    # insert spacer gaps between between successive bams
    data = [ item for x in zip(data, [spacer] * len(data)) for item in x ]
    data.pop()
    
    return np.concatenate(tuple(data))

def plot(seqfiles, fastafile, chrom, start, end, out=None, by_strand=False):
    
    if type(seqfiles) is not list:
        seqfiles = [seqfiles]
    
    chrom = str(chrom)
    fasta = pysam.FastaFile(fastafile)
    ref = fasta.fetch(start=start, end=end, region=chrom)
    
    use_ref = True
    RGBs = []
    for seqfile in seqfiles:
        seq = pysam.AlignmentFile(seqfile, 'rb')
        reads = seq.fetch(chrom, start, end)
        representations = getRepresentations(reads)
        RGB = generateRGB(representations, start, end, ref, by_strand, use_ref)
        
        use_ref = False
        RGBs.append(data)
    
    RGB = combine_bams(RGBs)

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

    if out:
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
    parser.add_argument('--by-strand', default=False, action='store_true',
        help='whether to color reads by strand (default is not)')
    parser.add_argument('--out', type=str, required=True,
        help='Output file (extension determines type: png, pdf, jpg, etc.)')

    args = parser.parse_args()
    seqfile = args.seqfile
    fastafile = args.fastafile
    chrom = args.chrom
    start = args.start
    end = args.end
    by_strand = args.by_strand
    out = args.out

    plot(seqfile, fastafile, chrom, start, end, out, by_strand)


if __name__ == '__main__':
    main()

    # fire.Fire(Shabam)
