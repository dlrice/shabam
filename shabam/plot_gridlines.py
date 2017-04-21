
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
