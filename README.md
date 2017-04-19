# shabam
A python tool to plot BAM or CRAM sequence reads.

## Installation
install [cairo](https://www.cairographics.org/download/) if not already
installed:
```sh
sudo apt-get install libcairo2-dev
```

Install shabam:
```sh
pip install git+git://github.com/dlrice/shabam.git --user
```


## The dream
```py
from shabam import seqplot
seqplot('example.bam', chrom='1', start=30243, end=30321,
    fastafile='reference.fasta', out='plot.svg')
```

![Reality](/tests/data/reality.svg)

### Plotting options
- shade reads by strand with `by_strand=True`
- plot multiple sequence files together with a list of paths e.g.
  `['child.bam', 'mom.bam', 'dad.bam']`
- export PDF, PNG, SVG or PS formatted plots with matching filename extensions
- seqplot returns PNG data if you don't include an output filename

## Improvements to realize the dream
- [ ] Use consensus sequence when we don't provide a reference sequence
- [ ] Allow custom colors
- [ ] Compute proportion of variants at any site
- [ ] Plot read depth
- [ ] At positions with proportion of variants > threshold, reflect proportion
  with base colors

## Credit
Initial cigar parsing code lifted with permission from
[pybamview](https://github.com/mgymrek/pybamview).
