
def parse_read(read, coords, ref=None, start=0):
    ''' parse the read data into a useable structure
    
    This has an intended side-effect of modifying the coords dictionary.
    
    Args:
        read: pysam.AlignedSegment for sequence read
        coords: ordered dictionary of end positions at each row, indexed by row
            number e.g. {10: 100, 20: 150, 30: 50}
        ref: reference sequence, or None
        start: position which the reference starts at
    
    Returns:
        dictionary that includes start position of read, a list of bases
        (mapped to M if they match the reference), whether the read is on the
        reverse strand, a list of base quality scores, and a y-axis offset to
        avoid superimposing different reads.
    '''
    if not read.cigartuples:
        return None
    
    data = {'position': read.pos,
            'bases': parse_cigar(read.cigartuples, read.query),
            'is_reverse': read.is_reverse,
            'qualities': parse_cigar(read.cigartuples, read.query_qualities)}
    
    # convert reference matches to 'M', so we can later color as reference bases
    if ref is not None:
        offset = read.pos - start
        for i, base in enumerate(data['bases']):
            if 0 <= offset + i < len(ref) and \
                    data['bases'][i] == ref[offset + i]:
                data['bases'][i] = 'M'
    
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
        positions: ordered dictionary of max end position at each row, indexed
            by row number e.g. {10: 100, 20: 150, 30: 50}
    
    Returns:
        y-axis offset.
    '''
    
    if len(coords) == 0:
        return 10
    
    for key in coords:
        if read['position'] > coords[key] + 1:
            return key
    
    return max(coords) + 10

def parse_cigar(cigar, bases):
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
    
    delchar = "-"
    bam = {'CMATCH': 0, 'CINS': 1, 'CDEL': 2, 'CREF_SKIP': 3, 'CSOFT_CLIP': 4,
        'CHARD_CLIP': 5, 'CPAD': 6, 'CEQUAL': 7, 'CDIFF': 8,}
    
    rep = []
    currentpos = 0
    wasinsert = False
    for operation, length in cigar:
        if operation in [bam['CMATCH'], bam['CEQUAL'], bam['CDIFF']]:
            # match (M, X, =)
            if operation == bam['CDIFF']:
                print(operation)
            for _ in range(length):
                if wasinsert:
                    rep[-1] += bases[currentpos]
                else:
                    rep.append(bases[currentpos])
                wasinsert = False
                currentpos += 1
        elif operation == bam['CINS']:
            # put insertion in next base position (I)
            if wasinsert:
                rep[-1] += bases[currentpos:currentpos + length]
            else:
                rep.append(bases[currentpos:currentpos + length])
            currentpos += length
            wasinsert = True
        elif operation in [bam['CDEL'], bam['CREF_SKIP']]:
            # deletion (D) or skipped region from the reference (N)
            for _ in range(length):
                if wasinsert:
                    rep[-1] += delchar
                else:
                    rep.append(delchar)
                wasinsert = False
        elif operation == bam['CPAD']:
            # padding (silent deletion from padded reference) (P)
            if wasinsert:
                rep[-1] += delchar * length
            else:
                rep.append(delchar * length)
            wasinsert = True
        elif operation not in [bam['CSOFT_CLIP'], bam['CHARD_CLIP']]:
            # hard clipping or soft clipping
            raise ValueError("Invalid CIGAR operation: {}".format(operation))
    return rep
