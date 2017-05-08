
import cairocffi as cairo
from shabam.colors import COLORS
from shabam.plot_reads import to_alpha

def plot_coverage(context, coverage, y_offset=0, x_offset=0, width=None, cov_height=50):
    ''' plots coverage along
    
    Args:
        context: cairo.context ploting device
        coverage: dict of dicts of read depths by base at each site in a region
            to be plotted
        y_offset: y-position to start plotting at
        x_offset: start base for the plotted region
        width: width of plotted region in pixels
        cov_height: height in pixels to plot coverage
    '''
    
    bottom = y_offset + cov_height
    max_cov = max(( sum(x.values()) for x in coverage.values() ))
    for pos in sorted(coverage):
        x_pos = (pos - x_offset) * 10
        if x_pos < 0 or x_pos > width - 1:
            continue
        
        cov = coverage[pos]
        height = sum(cov.values())/max_cov * cov_height
        
        # plot a gray box for the current base
        context.rectangle(x=x_pos, y=bottom-height, width=10, height=height)
        context.set_source_rgba(*(COLORS['M'] + [to_alpha(100)]))
        context.fill()
        
        # for each non-reference base, plot a shaded box
        y_pos = bottom - cov['M']/max_cov * cov_height
        for (depth, base) in sorted(zip(cov.values(), cov.keys()), reverse=True):
            if base == 'M':
                continue
            height = depth/max_cov * cov_height
            y_pos -= height
            
            context.rectangle(x=x_pos, y=y_pos, width=10, height=height)
            context.set_source_rgba(*(COLORS[base] + [to_alpha(100)]))
            context.fill()
