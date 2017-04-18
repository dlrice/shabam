
import os
import unittest
import tempfile

import cairocffi as cairo

from shabam.utils import get_height, fileformat

class TestUtils(unittest.TestCase):
    
    def test_get_height(self):
        ''' check that we can get the image height in pixels
        '''
        
        path = os.path.join(os.path.dirname(__file__), 'example.bam')
        
        height = get_height([path], '1', 30000, 30400, 0)
        self.assertEqual(height, 540)
        
        # test that we can allow for axis spacing
        height = get_height([path], '1', 30000, 30400, 50)
        self.assertEqual(height, 590)
        
        # test if we use a region without reads, we get an error
        # TODO: this is possible not the best way to handle this. Could be
        # TODO: better to return depths of zero (plus the axis offset)
        with self.assertRaises(ValueError):
            get_height([path], '2', 30000, 30400, 50)
    
    def test_fileformat_png(self):
        ''' test we can pick png output format
        '''
        temp = tempfile.NamedTemporaryFile(suffix='.png')
        ext, surf = fileformat(temp.name, 50, 100)
        self.assertEqual(ext, 'png')
        self.assertEqual(type(surf), cairo.ImageSurface)
        
        # check that we've set the height and width correctly
        self.assertEqual(surf.get_height(), 100)
        self.assertEqual(surf.get_width(), 50)
        
    def test_fileformat_pdf(self):
        ''' test we can pick pdf output format
        '''
        temp = tempfile.NamedTemporaryFile(suffix='.pdf')
        ext, surf = fileformat(temp.name, 50, 100)
        self.assertEqual(ext, 'pdf')
        self.assertEqual(type(surf), cairo.PDFSurface)
        
    def test_fileformat_svg(self):
        ''' test we can pick svg output format
        '''
        
        temp = tempfile.NamedTemporaryFile(suffix='.svg')
        ext, surf = fileformat(temp.name, 50, 100)
        self.assertEqual(ext, 'svg')
        self.assertEqual(type(surf), cairo.SVGSurface)
        
    def test_fileformat_ps(self):
        ''' test we can pick ps output format
        '''
        temp = tempfile.NamedTemporaryFile(suffix='.ps')
        ext, surf = fileformat(temp.name, 50, 100)
        self.assertEqual(ext, 'ps')
        self.assertEqual(type(surf), cairo.PSSurface)
        
    def test_fileformat_unknown(self):
        ''' test we raise error for unknown output format
        '''
        with self.assertRaises(AssertionError):
            fileformat('a.zzz', 50, 100)
