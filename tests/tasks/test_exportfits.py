import os
import numpy
import sys
import shutil
import unittest

from casatools import image

from casatasks import exportfits

class exportfits_test(unittest.TestCase):
          
    def test_CAS3675(self):
        """ test fix for CAS 3675, outfile must be specified """
        name = "my.im"
        yy = image()
        yy.fromshape(name, [1,1,1,1])
        yy.done()
        self.assertRaises(Exception, exportfits, imagename=name, overwrite=True)
 
def suite():
    return [exportfits_test]        

if __name__ == '__main__':
    unittest.main()
