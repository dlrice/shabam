
from __future__ import division

import math

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
    context.set_line_width(1)
    
    pos = start + -start % 10
    while pos <= end:
        _, _, _, height, _, _ = context.text_extents(str(pos))
        
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
