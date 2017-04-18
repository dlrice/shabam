
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
        for i in range(len(data['bases'])):
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

def parse_cigar(cigar, values):
    ''' get list of data, each corresponding to a single reference position
    
    Initial code (with permission) from https://github.com/mgymrek/pybamview
    
    Args:
        cigar: list of cigar tuples in read.
        values: list of data, either nucleotide sequence of values in read
            or list of quality scores.
    
    Returns:
        list of values, e.g ['A', 'T', 'C'] or [10, 5, 30]. Positions
        corresponding to deletions contain '-', and insertions have multiple
        data points e.g. ['A', 'TC'] or [10, [5, 30]].
    '''
    
    delchar = "-"
    bam = {'CMATCH': 0, 'CINS': 1, 'CDEL': 2, 'CREF_SKIP': 3, 'CSOFT_CLIP': 4,
        'CHARD_CLIP': 5, 'CPAD': 6, 'CEQUAL': 7, 'CDIFF': 8,}
    
    rep = []
    pos = 0
    wasinsert = False
    for code, length in cigar:
        if code in [bam['CMATCH'], bam['CEQUAL'], bam['CDIFF']]:
            # match (M, X, =)
            for _ in range(length):
                rep = cigar_join(rep, values[pos], wasinsert)
                wasinsert = False
                pos += 1
        elif code == bam['CINS']:
            # put insertion in next base position (I)
            rep = cigar_join(rep, values[pos:pos + length], wasinsert)
            pos += length
            wasinsert = True
        elif code in [bam['CDEL'], bam['CREF_SKIP']]:
            # deletion (D) or skipped region from the reference (N)
            for _ in range(length):
                rep = cigar_join(rep, delchar, wasinsert)
                wasinsert = False
        elif code == bam['CPAD']:
            # padding (silent deletion from padded reference) (P)
            rep = cigar_join(rep, delchar * length, wasinsert)
            wasinsert = True
        elif code not in [bam['CSOFT_CLIP'], bam['CHARD_CLIP']]:
            # hard clipping or soft clipping
            raise ValueError("Invalid CIGAR operation: {}".format(code))
    return rep

def cigar_join(initial, value, wasinsert):
    ''' adds a value to a list
    
    Args:
        initial: list to be added to
        value: data point to be added. This is either a string or an integer
        wasinsert: whether the previous position contained an insertion
    
    Returns:
        list with value added appropriately
    '''
    if not wasinsert:
        initial.append(value)
        return initial
    
    # now handle insertions. Concatenate strings (e.g. base-pairs), but
    # append integers (e.g. quality scores).
    try:
        initial[-1] += value
    except TypeError:
        initial[-1].append(value)
    
    return initial
