from pkg_resources import get_distribution

__version__ = get_distribution('shabam').version

from shabam.plot import seqplot
