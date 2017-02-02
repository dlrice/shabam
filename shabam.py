#!/usr/bin/env python
import pysam

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


def ParseCigar(cigar, nucs):
    """
    Return list of strings, each item corresponding to a single reference position
    """
    rep = []
    currentpos = 0
    wasinsert = False
    for operation, length in cigar:
#         print(operation, length)
        if operation in [BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF]: # match (M, X, =)
            for i in range(length):
                if wasinsert:
                    rep[-1] = rep[-1] + nucs[currentpos]
                else:
                    rep.append(nucs[currentpos])
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


def main():
    samfile = pysam.AlignmentFile("test/example.bam", "rb")

    chrom=16
    start=48000000
    end=48000080

    reads = samfile.fetch("16", start, end)
    representations = {}
    for read in reads:
        position = read.pos
        bases = read.query
        cigar = read.cigartuples
    #     if read.is_reverse:
    #         bases = bases.lower()
        representations[position] = ParseCigar(cigar, bases)
    print(representations)




if __name__ == '__main__':
    main()
