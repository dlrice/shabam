
from shabam.colors import COLORS

def plot_read(context, bases, quals=None, x_offset=0, y_offset=0, width=None,
        is_reverse=False, by_strand=False):
    ''' plots the bases in a read to a cairocffi.Context
    
    Args:
        context: cairocffi.Context as a plotting device
        bases: list of bases (per parse_read, so indels are odd)
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
        x_pos = (x_offset + i) * 10
        if x_pos < 0 or x_pos > width - 1:
            # don't plot bases outside the required window. This is necessary
            # when plotting SVGs, otherwise the SVG includes the outside bases.
            continue
        
        if len(base) > 1:
            plot_insertion(context, base, x_pos, y_offset)
            base = 'M'
        
        if base == 'M' and by_strand:
            strand = {True: 'r', False: 'f'}[is_reverse]
            base = 'M_{}'.format(strand)
        
        context.rectangle(x=x_pos, y=y_offset, width=10, height=10)
        context.set_source_rgba(*(COLORS[base] + [to_alpha(qual)]))
        context.fill()

def to_alpha(qual, threshold=35):
    ''' convert base quality to an alpha transparency float
    '''
    return min(threshold, qual)/threshold

def plot_insertion(context, bases, x_pos, y_offset):
    ''' plot inserted bases at the insertion site
    
    Args:
        context: cairocffi.Context as a plotting device
        bases: string of inserted bases
        x_pos: position of insertion (in pixels)
        y_offset: y position to plot the read at
    '''
    
    # select a font, and figure out the text sizes, so we can align text
    context.select_font_face('Arial')
    context.set_font_size(7)
    _, _, width, _, _, _ = context.text_extents(bases)
    
    context.move_to(x_pos + 10 - width/2, y_offset - 3)
    context.set_source_rgb(0, 0, 0)
    context.show_text(bases)
    
    # plot an arrow to indicate the insertion point
    context.set_line_width(1)
    context.move_to(x_pos + 10 - 1.5, y_offset - 1.5)
    context.line_to(x_pos + 10, y_offset)
    context.line_to(x_pos + 10 + 1.5, y_offset - 1.5)
    context.close_path()
    context.stroke_preserve()
    context.set_source_rgb(0, 0, 0)
    context.fill()
