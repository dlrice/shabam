# shabam
A python tool to create sequence plots from bam/cram files.

The dream:

```
import shabam
shabam.plot('some_sequence.bam', chrom=2, start=2323234, end=4235235345, file='reads_for_sample_323.png')
```

![The dream](/dream.png?raw=true)

## The steps

#### Create text based view (similar to samtools tview)
- [ ] Extract reads from a file
- [ ] Read a reference file
- [ ] Line up reads from a bam with a reference
- [ ] Print reads in terminal
- [ ] Optimize vertical height

#### Sequence plot
- [ ] Plot a list of letters using matplotlib patches
- [ ] Scale colors based on qual values

#### Historgram plot
- [ ] Compute proportion of variants at any site
- [ ] Plot read depth
- [ ] At positions with proportion of variants > threshold, reflect proportion with base colors

## Development setup
1. Install python 2 [try miniconda](http://conda.pydata.org/miniconda.html)
2. Install [pysam](https://github.com/pysam-developers/pysam/)
3. Ensure git is installed (try `git --version`)
4. `git clone https://github.com/dlrice/shabam/`
