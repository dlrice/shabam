
import os

import pysam
import cairocffi as cairo

from shabam.parse_reads import parse_read
from shabam.plot_reads import plot_read
from shabam.plot_axis import plot_axis

def plot_grid(context, start, end, axis_offset, height):
    ''' plot vertical lines every ten bases, so we can orient ourselves
    
    Args:
        context: cairocffi.Context as a plotting device
        start: nucleotide position at start of plotting window
        end: nucleotide position at end of plotting window
        axis_offset: height to account for in axis plotting.
        height: total height of the plotted image in pixels.
    '''
    
    context.set_source_rgba(0.5, 0.5, 0.5, 0.5)
    context.set_line_width(1)
    
    pos = start + -start % 10
    while pos <= end:
        x_pos = (pos - start) * 10
        
        context.move_to(x_pos, axis_offset)
        context.line_to(x_pos, height)
        context.stroke()
        
        pos += 10

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

def get_height(seqfiles, chrom, start, end, axis_offset):
    ''' get the height of the output image in pixels
    
    This requires we loop through all of the reads in the various BAMs, before
    doing so a second time when we plot the reads, but I can't see an easier way.
    I can't get RecordingSurface to work, which would be one way to fix this,
    make a RecordingSurface, figure out the ink extents, then copy the surface
    to an output PDF or SVG surface.
    
    Args:
        seqfiles: list of paths to sequence files
        chrom: chromosome to fetch reads from
        start: start nucleotide of plotting window.
        end: end nucleotide of plotting window.
        axis_offset: height to account for in axis plotting.
    
    Returns:
        height of image plot in pixels.
    '''
    
    depths = [axis_offset]
    for seqfile in seqfiles:
        seq = pysam.AlignmentFile(seqfile, 'rb')
        
        coords = {}
        for read in seq.fetch(chrom, start, end):
            _ = parse_read(read, coords)
        
        depths.append(max(coords))
    
    spacer_depth = (len(depths) - 2) * 5 * 10
    
    return sum(depths) + spacer_depth + 10

def fileformat(filename, width, height):
    ''' figure out the output image format
    
    Args:
        filename: output filename, which will raise an error if it does not end
            in one of '.pdf', '.png', '.ps', or '.svg'
        width: width of the output image, in pixels
        height: height of the output image, in pixels
    
    Returns:
        tuple of cairo.Surface and filetype string e.g. 'pdf' or 'png'
    '''
    
    if filename is not None:
        _, ext = os.path.splitext(filename)
        if not ext:
            ext = '.png'
        ext = ext[1:].lower()
    else:
        ext = None
    
    assert ext in ['png', 'pdf', 'ps', 'svg', None], 'unknown format: .{}'.format(ext)
    
    if ext == 'pdf':
        surface = cairo.PDFSurface(filename, width, height)
    elif ext == 'svg':
        surface = cairo.SVGSurface(filename, width, height)
    elif ext == 'ps':
        surface = cairo.PSSurface(filename, width, height)
    else:
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width,  height)
    
    return ext, surface

def shabam(seqfiles, chrom, start, end, fastafile, out=None, by_strand=False):
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
        coords = {max(depths): -1e9}
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
