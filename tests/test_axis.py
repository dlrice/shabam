
import unittest

from shabam.plot_axis import get_axis_increment

class TestAxisPlotting(unittest.TestCase):
    def test_axis_increment(self):
        ''' test that we can get reasonable axis increments
        '''
        
        self.assertEqual(get_axis_increment(10, 30), 5)
        self.assertEqual(get_axis_increment(0, 100), 40)
        self.assertEqual(get_axis_increment(0, 1000), 400)
        self.assertEqual(get_axis_increment(0.1, 0.15), 0.01)
        
        with self.assertRaises(AssertionError):
            get_axis_increment(0, 0)
