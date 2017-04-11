
import os
import math

import pysam
import cairocffi as cairo

COLORS = {
    'A':   [0.0, 0.4, 0.0], # green
    'C':   [0.0, 0.0, 1.0], # blue
    'G':   [1.0, 0.6, 0.2], # orange
    'T':   [1.0, 0.0, 0.0], # red
    'M':   [0.9, 0.9, 0.9], # match, gray
    'M_f': [0.3, 0.7, 1.0], # forward match, orange
    'M_r': [1.0, 0.7, 0.3], # reverse match, blue
    '-':   [0.9, 0.4, 0.8], # deletion, pink
    'I':   [0.4, 0.1, 0.8], # insertion, purple
    'N':   [0.0, 0.0, 1.0], # unknown, black
}

def parseCigar(cigar, bases):
    ''' get list of bases, each corresponding to a single reference position
    
    Initial code (with permission) from https://github.com/mgymrek/pybamview
    
    Args:
        cigar: list of cigar tuples in read.
        bases: nucleotide sequence of bases in read.
    
    Returns:
        list of bases, e.g ['A', 'T', 'M']. Matches to the reference are
        coded as 'M',  deletions as '-', and insertions as multinucleotide
        entries e.g. 'MATGC'.
    '''
    
    ENDCHAR = "-"
    GAPCHAR = "."
    DELCHAR = "-"
    
    BAM = {'CMATCH': 0, 'CINS': 1, 'CDEL': 2, 'CREF_SKIP': 3, 'CSOFT_CLIP': 4,
        'CHARD_CLIP': 5, 'CPAD': 6, 'CEQUAL': 7, 'CDIFF': 8,}
    
    rep = []
    currentpos = 0
    wasinsert = False
    for operation, length in cigar:
        if operation in [BAM['CMATCH'], BAM['CEQUAL'], BAM['CDIFF']]: # match (M, X, =)
            if operation == BAM['CDIFF']:
                print(operation)
            for i in range(length):
                if wasinsert:
                    rep[-1] = rep[-1] + bases[currentpos]
                else:
                    rep.append('M')
                wasinsert = False
                currentpos += 1
        elif operation == BAM['CINS']: # put insertion in next base position (I)
            if wasinsert:
                rep[-1] = rep[-1] + bases[currentpos:currentpos + length]
            else:
                rep.append(bases[currentpos:currentpos + length])
            currentpos = currentpos + length
            wasinsert = True
        elif operation in [BAM['CDEL'], BAM['CREF_SKIP']]: # deletion (D) or skipped region from the reference (N)
            for i in range(length):
                if wasinsert:
                    rep[-1] = rep[-1] + DELCHAR
                else:
                    rep.append(DELCHAR)
                wasinsert = False
        elif operation == BAM['CPAD']: # padding (silent deletion from padded reference) (P)
            if wasinsert:
                rep[-1] = rep[-1] + DELCHAR * length
            else:
                rep.append(DELCHAR * length)
            wasinsert = True
        elif operation not in [BAM['CSOFT_CLIP'], BAM['CHARD_CLIP']]: # hard clipping or soft clipping
            sys.stderr.write("ERROR: Invalid CIGAR operation (%s) in read %s \n"%(operation, read.qname))
    return rep

def plot_read(context, bases, quals=None, x_offset=0, y_offset=0, width=None,
        is_reverse=False, by_strand=False):
    ''' plots the bases in a read to a cairocffi.Context
    
    Args:
        context: cairocffi.Context as a plotting device
        bases: list of bases (per parseRead, so indels are odd)
        quals: list of quality scores for each base
        x_offset: x position to start plotting the read at
        y_offset: y position to plot the read at
        width: with of image in pixels
        is_reverse: whether the read is for the reverse strand
        by_strand: boolean for whether we want to shade reads by strand
    '''
    
    if quals is None:
        quals = [100] * len(bases)
    
    if width is None:
        width = len(bases) * 10
    
    for i, (base, qual) in enumerate(zip(bases, quals)):
        
        # TODO: This is an extremely crude adjustment to include insertions.
        # TODO: This should properly be handled by including the inserted
        # TODO: sequence at the insertion site.
        if len(base) > 1:
            base = 'I'
        
        if base == 'M' and by_strand:
            strand = {True: 'r', False: 'f'}[is_reverse]
            base = 'M_{}'.format(strand)
        
        x_pos = (x_offset + i) * 10
        if x_pos < 0 or x_pos > width - 1:
            # don't plot bases outside the required window. This is necessary
            # when plotting SVGs, otherwise the SVG includes the outside bases.
            continue
        
        context.rectangle(x=x_pos, y=y_offset, width=10, height=10)
        context.set_source_rgba(*COLORS[base], to_alpha(qual))
        context.fill()

def get_axis_increment(start, end):
    ''' figure out distance between axis ticks
    
    Args:
        start: nucleotide position at start of plotting window
        end: nucleotide position at end of plotting window
    
    Returns:
        distance in base-pairs between axis ticks, so as to get either three or
        four ticks along the axis
    '''
    
    assert start != end, 'the start position is the same as the end'
    
    # get plotting positions
    delta = abs(end - start)
    
    ratios = [1, 2, 5, 4]
    
    multiplier = 1
    while True:
        for ratio in ratios:
            if 2 < delta/(ratio * multiplier) < 5:
                return (ratio * multiplier)
        
        if ratio * multiplier > delta:
            multiplier /= 10
        else:
            multiplier *= 10

def plot_axis(context, start, end, axis_offset):
    ''' plot an axis, with nucldeotide positions for convenience
    
    Args:
        context: cairocffi.Context as a plotting device
        start: nucleotide position at start of plotting window
        end: nucleotide position at end of plotting window
        axis_offset: height to account for in axis plotting.
    '''
    
    increment = get_axis_increment(start, end)
    
    rotate = -90
    
    # select a font, and figure out the text sizes, so we can align text
    context.select_font_face('Arial')
    context.set_font_size(10)
    fascent, fdescent, fheight, fxadvance, fyadvance = context.font_extents()
    context.set_line_width(1)
    
    pos = start + -start % 10
    while pos <= end:
        xbearing, ybearing, width, height, xadvance, yadvance = \
            context.text_extents(str(pos))
        
        # center align the text
        x_pos = (pos - start) * 10 + 5
        x_text = x_pos + height / 2
        
        context.move_to(x_text, axis_offset - 10)
        context.rotate(math.radians(rotate))
        context.set_source_rgb(0, 0, 0)
        context.show_text(str(pos))
        context.rotate(-math.radians(rotate))
        
        # and plot a tick mark
        context.move_to(x_pos, axis_offset)
        context.line_to(x_pos, axis_offset - 5)
        context.stroke()
        
        pos += increment

def generateRGB(context, reads, start, end, axis_offset, height, ref_seq=None, by_strand=False):
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
        pattern = plot_read(context, ref_seq, y_offset=axis_offset - 10)
    
    width = (end - start) * 10
    
    for read in reads:
        if read is None:
            continue
        
        pattern = plot_read(context, read['bases'], read['qualities'],
            read['position'] - start, read['offset'], width, read['is_reverse'],
            by_strand)
    
    plot_axis(context, start, end, axis_offset - 10)

def to_alpha(qual, threshold=35):
    ''' convert base quality to an alpha transparency float
    '''
    return min(threshold, qual)/threshold

def parseRead(read, coords):
    ''' parse the read data into a useable structure
    
    This has an intended side-effect of modifying the coords dictionary.
    
    Args:
        read: pysam.AlignedSegment for sequence read
        coords: dictionary of end positions at each row, indexed by row number
            e.g. {10: 100, 20: 150, 30: 50}
    
    Returns:
        dictionary that includes start position of read, a list of bases
        (mapped to M if they match the reference), whether the read is on the
        reverse strand, a list of base quality scores, and a y-axis offset to
        avoid superimposing different reads.
    '''
    if not read.cigartuples:
        return None
    
    data = {'position': read.pos,
            'bases': parseCigar(read.cigartuples, read.query),
            'is_reverse': read.is_reverse,
            'qualities': read.query_qualities}
    
    y_pos = get_y_offset(data, coords)
    data['offset'] = y_pos
    
    if y_pos not in coords:
        coords[y_pos] = -1e9
    
    if data['position'] + len(data['bases']) > coords[y_pos]:
        coords[y_pos] = data['position'] + len(data['bases'])
    
    return data

def get_y_offset(read, coords):
    ''' get the y-position for the read, within the first row with open space
    
    This is decidedly less than optimal, but will produce a more compact plot.
    
    Args:
        read: dictionary with 'position' and 'bases' entries e.g.
            {'position': 15555, bases ['M', 'MTA', 'M']}
        positions: max end position at each row, indexed by row number
            e.g. {10: 100, 20: 150, 30: 50}
    
    Returns:
        y-axis offset.
    '''
    
    if len(coords) == 0:
        return 10
    
    for key in sorted(coords):
        if read['position'] > coords[key] + 1:
            return key
    
    return max(coords) + 10

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
            context.rectangle(x=0, y=(end - start) * 10, width=10, height=10)
            context.set_source_rgba(0.0, 0.0, 0.0, 1.0)
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
            unused = parseRead(read, coords)
        
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
    
    assert ext in ['png', 'pdf', 'ps', 'svg', None], 'unknown format: ".{}"'.format(ext)
    
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
    
    width = (end - start) * 10
    if type(seqfiles) is not list:
        seqfiles = [seqfiles]
    
    chrom = str(chrom)
    fasta = pysam.FastaFile(fastafile)
    reference = fasta.fetch(start=start, end=end, region=chrom)
    
    axis_offset = 70
    height = get_height(seqfiles, chrom, start, end, axis_offset)
    
    out_type, surface = fileformat(out, width, height)
    context = cairo.Context(surface)
    
    depths = [axis_offset]
    for seqfile in seqfiles:
        seq = pysam.AlignmentFile(seqfile, 'rb')
        coords = {max(depths): -1e9}
        reps = ( parseRead(x, coords) for x in seq.fetch(chrom, start, end) )
        
        generateRGB(context, reps, start, end, axis_offset, height, reference, by_strand)
        reference = None # don't plot the reference in subsequent BAMs
        
        if seqfiles.index(seqfile) < len(seqfiles) - 1:
            insert_spacer(context, coords, start, end)
        
        depths.append(max(coords))
    
    context.save()
    
    if out_type == 'png':
        surface.write_to_png(out)
    
    if out is None:
        return surface.write_to_png()
