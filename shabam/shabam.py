
import os
from collections import OrderedDict

import pysam
import cairocffi as cairo

from shabam.parse_reads import parse_read
from shabam.plot_reads import plot_read
from shabam.plot_axis import plot_axis
from shabam.plot_gridlines import plot_grid
from shabam.utils import get_height, fileformat

def _plot(context, reads, start, end, axis_offset, height, ref_seq=None,
        by_strand=False):
    ''' plots reads to the Context
    
    Args:
        context: cairocffi.Context as a plotting device
        reads: iterator of parsed read Data
        start: nucleotide position at start of plotting window
        end: nucleotide position at end of plotting window
        axis_offset: height to account for in axis plotting.
        height: total height of the plotted image in pixels.
        ref_seq: reference sequence within plotting window (or None)
        by_strand: whether to shade reads by strand
    
    Returns:
        max end position at each row, indexed by row number
        e.g. {10: 100, 20: 150, 30: 50}
    '''
    
    if ref_seq is not None:
        plot_read(context, ref_seq, y_offset=axis_offset - 10)
    
    width = (end - start) * 10
    
    for read in reads:
        if read is None:
            continue
        
        plot_read(context, read['bases'], read['qualities'],
            read['position'] - start, read['offset'], width, read['is_reverse'],
            by_strand)
    
    plot_axis(context, start, end, axis_offset - 10)
    plot_grid(context, start, end, axis_offset, height)

def insert_spacer(context, coords, start, end):
    ''' combine data for one or more bams into a single array
    
    Args:
        context: cairocffi.Context as a plotting device
        coords: max end position at each row, indexed by row number
            e.g. {10: 100, 20: 150, 30: 50}
        start: initial nucleotide position of the region being plotted
        end: final nucleotide position of the region being plotted
        
    Returns:
        max end position at each row, indexed by row number
        e.g. {10: 100, 20: 150, 30: 50}
    '''
    
    # define gap between BAMs as white space with a black line in the middle
    lines = ['blank', 'blank', 'black', 'blank', 'blank']
    
    for i, line in enumerate(lines):
        coords[max(coords) + i * 10] = start + end
        
        if line != 'blank':
            width = (end - start) * 10
            context.rectangle(x=0, y=max(coords) + 10, width=width, height=10)
            context.set_source_rgba(0.4, 0.4, 0.4, 1.0)
            context.fill()

    return coords

def seqplot(seqfiles, chrom, start, end, fastafile, out=None, by_strand=False):
    ''' the plotting function
    
    Args:
        seqfiles: list of paths to sequence files
        chrom: chromosome to fetch reads from
        start: start nucleotide of plotting window.
        end: end nucleotide of plotting window.
        fastafile: path to reference FASTA file.
        out: path to write image file to, or None to return bytes-encoded png
        by_strand: whether to shade reads by strand
    
    Returns:
        None, or if out is None, returns image plot as bytes-encoded png
    '''
    
    if type(seqfiles) is not list:
        seqfiles = [seqfiles]
    
    chrom = str(chrom)
    with pysam.FastaFile(fastafile) as handle:
        reference = handle.fetch(start=start, end=end, region=chrom)
        ref = reference
    
    axis_offset = 75
    height = get_height(seqfiles, chrom, start, end, axis_offset)
    
    out_type, surface = fileformat(out, width=(end - start) * 10, height=height)
    context = cairo.Context(surface)
    
    depths = [axis_offset]
    for seqfile in seqfiles:
        seq = pysam.AlignmentFile(seqfile, 'rb')
        coords = OrderedDict({max(depths): -1e9})
        reps = ( parse_read(x, coords, ref, start) for x in seq.fetch(chrom, start, end) )
        
        _plot(context, reps, start, end, axis_offset, height, reference, by_strand)
        reference = None # don't plot the reference in subsequent BAMs
        
        if seqfiles.index(seqfile) < len(seqfiles) - 1:
            insert_spacer(context, coords, start, end)
        
        depths.append(max(coords))
    
    context.save()
    
    if out_type == 'png':
        surface.write_to_png(out)
    elif out_type is None:
        return surface.write_to_png()
