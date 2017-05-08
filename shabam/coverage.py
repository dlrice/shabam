
def get_coverage(reads):
    ''' get the coverage by base at each site in a set of reads
    
    Args:
        reads: iterable of reads. Each read is a dict containing 'position' and
            'bases' (among other fields).
    
    Returns:
        Returns a dict of read depths by base, indexed by nucleotide position.
        We previously checked for matches to the reference sequence, so we can
        include an additional 'base', 'M', for the base that matches the
        reference. For example:
            {
                1: {'A': 0, 'C': 0, 'G': 0, 'T': 20, 'M': 30},
                2: {'A': 30, 'C': 2, 'G': 0, 'T': 0, 'M': 20},
                ...
            }
    '''
    
    bases = set(['A', 'C', 'G', 'T', 'M'])
    coverage = {}
    
    for read in reads:
        if read is None:
            continue
        
        for i, base in enumerate(read['bases']):
            # only look at the first base of insertions, since this is the only
            # base that corresponds to a reference position
            if len(base) > 1:
                base = base[1]
            
            # ignore deletions, since those will correspond to coverage drop
            if base not in bases:
                continue
            
            pos = read['position'] + i
            if pos not in coverage:
                coverage[pos] = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'M': 0}
            
            coverage[pos][base] += 1
    
    return coverage
