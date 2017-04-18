
import tempfile
import unittest
import os
import hashlib

from shabam import seqplot

class TestPlot(unittest.TestCase):
    
    def test_seqplot(self):
        ''' test that we can make plots
        '''
        
        folder = os.path.dirname(__file__)
        path = os.path.join(folder, 'data', 'example.bam')
        fasta = os.path.join(folder, 'data', 'reference.fasta')
        
        # check that we make a plot with the correct image. check this by
        # returning png data, then getting the sha1 digest. This will fail if
        # a single bit chenges, so perhaps this is too stringent.
        data = seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta)
        checksum = hashlib.sha1(data).hexdigest()
        
        # now try writing a plot to PNG file.
        single = tempfile.NamedTemporaryFile(suffix='.png')
        seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta,
            out=single.name)
        
        self.assertTrue(os.path.getsize(single.name) > 0)
        # check we have written the same data as would be returned
        with open(single.name, 'rb') as handle:
            singlechecksum = hashlib.sha1(handle.read()).hexdigest()
        
        self.assertEqual(checksum, singlechecksum)
        
        # now check we can plot multiple sequence files in a single image
        double = tempfile.NamedTemporaryFile(suffix='.png')
        seqplot([path, path], chrom='1', start=30000, end=30400, fastafile=fasta,
            out=double.name)
        self.assertTrue(os.path.getsize(double.name) > os.path.getsize(single.name))
        
        # check that the sha1 hash does not match the earlier one. this may
        # be excessive.
        with open(double.name, 'rb') as handle:
            doublechecksum = hashlib.sha1(handle.read()).hexdigest()
        self.assertNotEqual(singlechecksum, doublechecksum)
    
    def test_seqplot_by_strand(self):
        ''' test that we can plot bam, and color reads by strand
        '''
        
        folder = os.path.dirname(__file__)
        path = os.path.join(folder, 'data', 'example.bam')
        fasta = os.path.join(folder, 'data', 'reference.fasta')
        
        # check that we can plot by strand
        data = seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta,
            by_strand=True)
        checksum = hashlib.sha1(data).hexdigest()
        self.assertEqual(checksum, '0b6cc0823c11da3dd56e49cd7052b4e3c8342836')
        
        data = seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta)
        checksum_unstranded = hashlib.sha1(data).hexdigest()
        self.assertNotEqual(checksum, checksum_unstranded)
    
    def test_seqplot_cram(self):
        ''' test that we can plot from CRAM files
        '''
        
        folder = os.path.dirname(__file__)
        path = os.path.join(folder, 'data', 'example.cram')
        fasta = os.path.join(folder, 'data', 'reference.fasta')
        
        single = tempfile.NamedTemporaryFile(suffix='.png')
        seqplot([path], chrom='1', start=30000, end=30400, fastafile=fasta,
            out=single.name)
        
        # check that we have written some data to the file
        self.assertTrue(os.path.getsize(single.name) > 0)
