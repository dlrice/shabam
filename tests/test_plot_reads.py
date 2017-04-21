
import unittest

from shabam.plot_reads import to_alpha

class TestPlotReads(unittest.TestCase):
    def test_to_alpha(self):
        ''' test that we can convert quality scores to alpha transparency
        '''
        
        # check at top of range
        self.assertEqual(to_alpha(35), 1.0)
        
        # check at bottom of range
        self.assertEqual(to_alpha(0), 0.0)
        
        # check in middle of range
        self.assertEqual(to_alpha(17), 0.4857142857142857)
        
        # check above top of range gets set to top of range
        self.assertEqual(to_alpha(40), 1.0)
        
        # check if threshold is changed, the value scales appropriately
        self.assertEqual(to_alpha(40, 80), 0.5)
        
        # check with missing deletion-style quality score
        self.assertEqual(to_alpha('-'), 1.0)
        
