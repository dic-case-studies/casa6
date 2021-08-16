from __future__ import absolute_import
from __future__ import print_function
import os
import shutil
import numpy as np
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys
    from casatasks import predictcomp, casalog

    ctsys_resolve = ctsys.resolve
else:
    from __main__ import default
    from tasks import predictcomp 
    from taskinit import *

    dataRoot = os.path.join(os.environ.get('CASAPATH').split()[0],'data/')
    def ctsys_resolve(apath):
        return os.path.join(dataRoot,apath)

''' Python unit tests for the predictcomp task

 - tests the following parameters:
 objname: error for unsupported object vs supported object
          non-visible case
 standard: wrong standard vs correct standard
 minfreq/maxfreq: wrong unit vs correct unit
 output: check for the cl file
 antennalist: use of the configuration file to plot 'observed' visibility
              amplitudes vs uvdist. GUI is turned off.
 

'''

datapath = ctsys_resolve('alma/simmos/')


class predictcomp_test(unittest.TestCase):

    def setUp(self):
        self.res=None
        if not is_CASA6:
            default(predictcomp) 

    def tearDown(self):
        #pass
        os.system('rm -rf *.cl')

        
    def test_default(self):
        '''predictcomp: test defaults'''
        with self.assertRaises(ValueError):
            predictcomp()
 
    def test_invalid_objname(self): 
        '''predictcomp: invalid objname'''
        with self.assertRaises(RuntimeError):
            predictcomp(objname='Moon', minfreq='100GHz',maxfreq='120GHz')

    def test_valid_objname(self):
        '''predictcomp: valid objname'''
        self.res=predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100GHz',maxfreq='120GHz',
                             standard='Butler-JPL-Horizons 2012') 
        print("type(self.res) = ",type(self.res))
        self.assertTrue(type(self.res)==dict)
        self.assertTrue(os.path.exists(self.res['clist']))
             
    def test_invalid_freqrange(self):
        '''predictcomp: invalid freqrange'''
        with self.assertRaises(RuntimeError):
            predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100', maxfreq='120' )
    
    def test_badmaxfreq(self):
        '''predictcomp: invalid maxfreq '''
        # ignore the maxfreq just calculate for minfreq
        self.res=predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100GHz', maxfreq='120BadUnit', standard='Butler-JPL-Horizons 2012') 
        self.assertTrue(type(self.res)==dict)
        self.assertTrue(os.path.exists(self.res['clist']))

    @unittest.skipIf(is_CASA6,"no plotting in casatasks")
    def test_predicted_visplot(self):
        '''predictcomp: generate visibility plot for a given array configuration''' 
        self.res=predictcomp( objname='Titan', epoch='2017/09/01/00:00', minfreq='100GHz',
                              maxfreq='120GHz', standard='Butler-JPL-Horizons 2012',
                              antennalist=ctsys_resolve(os.path.join(datapath,'alma.cycle5.1.cfg')),
                              showplot=False,savefig='visplot.png' )
        print("self.res : %s" % self.res)
        print("type : %s" % type(self.res))
        self.assertTrue(type(self.res)==dict)
        self.assertTrue(os.path.exists(self.res['clist']))
        self.assertTrue(os.path.exists('visplot.png'))
          
    def test_valid_but_not_visible_objname(self):
        '''predictcomp: valid but not visible objname'''
        with self.assertRaises(RuntimeError):
            predictcomp(objname='Mars', epoch='2018/09/01/00:00', minfreq='100GHz', maxfreq='120GHz',
                               antennalist=ctsys_resolve(os.path.join(datapath,'alma.cycle5.1.cfg')), showplot=False )

def suite():
    return [predictcomp_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()
