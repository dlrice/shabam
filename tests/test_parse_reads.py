
import unittest

from collections import namedtuple, OrderedDict

from shabam.parse_reads import parse_read, get_y_offset, parse_cigar

class TestReadParsing(unittest.TestCase):
    
    Read = namedtuple('Read', ['pos', 'cigartuples', 'query', 'is_reverse',
        'query_qualities'])
    
    def setUp(self):
        
        pos = 50
        cigar = [(0, 3), (2, 2), (1, 2)]
        bases = 'ATGCC'
        quals = [30, 20, 40, 40, 10]
        reverse = False
        self.read = Read(pos, cigar, bases, reverse, quals)
    
    def test_parse_read(self):
        ''' test that we can parse a read
        '''
        pass
    
    def test_get_y_offset_empty(self):
        ''' test that get_y_offset_works correctly
        '''
        
        # all we need from the parsed read for getting the y-offset is the start
        # position.
        read = {'position': 50}
        
        # test with empty coords dictionary
        self.assertEqual(get_y_offset(read, coords=OrderedDict({}), 10)
        
        # test
        self.assertEqual(get_y_offset(read, OrderedDict({10: 48})), 10)
        self.assertEqual(get_y_offset(read, OrderedDict({10: 49})), 20)
        self.assertEqual(get_y_offset(read, OrderedDict({10: 50})), 20)
        self.assertEqual(get_y_offset(read, OrderedDict({10: 50, 20: 30})), 20)
        self.assertEqual(get_y_offset(read, OrderedDict({10: 50, 20: 50})), 30)
        
