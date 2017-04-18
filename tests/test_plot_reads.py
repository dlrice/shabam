
import unittest

from shabam.plot_reads import to_alpha

class TestPlotReads(unittest.TestCase):
    def test_to_alpha(self):
        ''' test that we can convert quality scores to alpha transparency
        '''
        
        self.assertEqual(to_alpha(35), 1.0)
        self.assertEqual(to_alpha(0), 0.0)
        self.assertEqual(to_alpha(17), 0.4857142857142857)
        self.assertEqual(to_alpha(40), 1.0)
        self.assertEqual(to_alpha('-'), 1.0)
        
