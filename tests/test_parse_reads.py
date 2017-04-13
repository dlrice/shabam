
import unittest

from collections import namedtuple, OrderedDict

from shabam.parse_reads import parse_read, get_y_offset, parse_cigar

class TestReadParsing(unittest.TestCase):
    
    Read = namedtuple('Read', ['pos', 'cigartuples', 'query', 'is_reverse',
        'query_qualities'])
    
    def test_parse_read(self):
        ''' check that a standard read parses cleanly
        '''
        coords = {}
        read = self.Read(0, [(0, 3)], 'AAA', False, [5, 5, 5])
        parsed = parse_read(read, coords)
        self.assertEqual(parsed, {'position': 0,
            'bases': ['A', 'A', 'A'], 'is_reverse': False,
            'qualities': [5, 5, 5], 'offset': 10})
        
        # ensure the coords dictionary has been updated in transit
        self.assertEqual(coords, {10: 3})
    
    def test_parse_read_missing_cigar(self):
        ''' check that a read without cigartuples returns None
        '''
        
        read = self.Read(0, None, 'AAA', False, [5, 5, 5])
        self.assertEqual(parse_read(read, {}), None)
    
    def test_parse_read_with_ref_sequence(self):
        ''' check we can specify a reference sequence to match against, and the
        relevant bases get converted to the match symbol.
        '''
        read = self.Read(0, [(0, 3)], 'AAA', False, [5, 5, 5])
        self.assertEqual(parse_read(read, {}, ref='AAG'), {'position': 0,
            'bases': ['M', 'M', 'A'], 'is_reverse': False,
            'qualities': [5, 5, 5], 'offset': 10})
    
    def test_parse_read_start_offset(self):
        ''' check we can use a reference sequence that extends around the
        current read, and providing a start position will let us access the
        correct sequence.
        '''
        read = self.Read(0, [(0, 3)], 'TGA', False, [5, 5, 5])
        self.assertEqual(parse_read(read, {}, ref='AATGGA', start=-2),
            {'position': 0, 'bases': ['M', 'M', 'A'], 'is_reverse': False,
            'qualities': [5, 5, 5], 'offset': 10})
    
    def test_get_y_offset_empty(self):
        ''' test that get_y_offset_works correctly
        '''
        
        # all we need from the parsed read for getting the y-offset is the start
        # position.
        read = {'position': 50}
        
        # test with empty coords dictionary
        self.assertEqual(get_y_offset(read, coords=OrderedDict({})), 10)
        
        # test with non-empty coords, but where the row ends before the read
        self.assertEqual(get_y_offset(read, OrderedDict({10: 48})), 10)
        
        # test when the row end is too close
        self.assertEqual(get_y_offset(read, OrderedDict({10: 49})), 20)
        
        # test when the row end is past the read
        self.assertEqual(get_y_offset(read, OrderedDict({10: 51})), 20)
        
        # test when multiple rows already exist, but the read starts beyond one
        self.assertEqual(get_y_offset(read, OrderedDict({10: 50, 20: 30})), 20)
        
        # test when the read starts beyond all existsing rows
        self.assertEqual(get_y_offset(read, OrderedDict({10: 50, 20: 50})), 30)
        
        # test that we can use extremely high depths if need be
        coords = OrderedDict(zip(range(100001), [100] * 100001))
        self.assertEqual(get_y_offset(read, coords), 100010)
    
    def test_parse_cigar(self):
        ''' test that cigar parsing works correctly
        '''
        
        bases = 'ATGC'
        self.assertEqual(parse_cigar([(0, 4)], bases), ['A', 'T', 'G', 'C'])
        
        # test when an insertion comes before matches
        self.assertEqual(parse_cigar([(1, 2), (0, 2)], bases), ['ATG', 'C'])
        
        # test when an insertion comes after matches
        self.assertEqual(parse_cigar([(0, 2), (1, 2)], bases), ['A', 'T', 'GC'])
        
        # test when a deletion is present
        self.assertEqual(parse_cigar([(0, 2), (2, 2), (0, 2)], bases),
            ['A', 'T', '-' , '-', 'G', 'C'])
        
        # tests for CPAD
        self.assertEqual(parse_cigar([(0, 2), (6, 1), (0, 2)], bases),
            ['A', 'T', '-G', 'C'])
        
        # test when we have some clipped sequence. Properly we should test when
        # the clipped region occurs at the end of reads.
        self.assertEqual(parse_cigar([(0, 2), (5, 1), (0, 1)], bases),
            ['A', 'T', 'G'])
        
        # check we raise an error with unknown cigar codes
        with self.assertRaises(ValueError):
            parse_cigar([(0, 2), (9, 2)], bases)
