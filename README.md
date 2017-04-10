# shabam
A python tool to create sequence plots from bam/cram files.

The dream:

```
import shabam
shabam.plot('sequence.bam', chrom=2, start=48000000, end=48000080, file='plot.png')
```

![The dream](/dream.png?raw=true)

The reality:

![Reality](/reality.png?raw=true)

## Steps to realize the dream

#### Create text based view (similar to samtools tview)
- [x] Extract reads from a file
- [x] Read a reference file
- [x] Line up reads from a bam with a reference
- [x] Print reads in terminal


#### Sequence plot
- [x] Plot bases as a single matrix RGB triples
- [x] Save plot as png
- [ ] Deal with insertions/deletions
- [ ] Gridline every 10 bases
- [ ] Optimize vertical space
- [x] Scale colors based on qual values (Thanks @jeremymcrae)
- [x] Reflect positive/negative strand (Thanks @jeremymcrae)
- [x] Plot multiple sequence files at once (Thanks @jeremymcrae)
- [ ] Allow custom colors from JSON file


#### Historgram plot
- [ ] Compute proportion of variants at any site
- [ ] Plot read depth
- [ ] At positions with proportion of variants > threshold, reflect proportion with base colors


## Development setup
1. Install python ([try miniconda](http://conda.pydata.org/miniconda.html) if you haven't got it)
2. Install sampy, matplotlib, numpy
3. Ensure git is installed (try `git --version`)
4. Ensure you have a github account.
5. Fork this repository and clone onto your local machine.


## Credit
Initial cigar parsing code lifted with permission from [pybamview](https://github.com/mgymrek/pybamview).
