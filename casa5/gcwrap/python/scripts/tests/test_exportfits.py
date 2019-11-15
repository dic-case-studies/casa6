import os
import numpy
import sys
import shutil
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import image
    from casatasks import exportfits
else:
    from tasks import exportfits
    from taskinit import iatool as image

class exportfits_test(unittest.TestCase):
          
    def test_CAS3675(self):
        """ test fix for CAS 3675, outfile must be specified """
        name = "my.im"
        yy = image()
        yy.fromshape(name, [1,1,1,1])
        yy.done()
        if is_CASA6:
            self.assertRaises(Exception, exportfits, imagename=name, overwrite=True)
        else:
            # CASA5 returns False
            ret = exportfits(imagename=name,overwrite=True)
            self.assertFalse(ret)
            
def suite():
    return [exportfits_test]        

if __name__ == '__main__':
    unittest.main()
