
import os
from collections import OrderedDict

import pysam
import cairocffi as cairo

from shabam.parse_reads import parse_read

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
        
        coords = OrderedDict()
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
