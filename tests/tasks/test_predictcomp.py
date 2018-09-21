import os
import shutil
import unittest
import numpy as np

from casatools import ctsys
from casatasks import predictcomp


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

datapath = 'alma/simmos'


class predictcomp_test(unittest.TestCase):

    def setUp(self):
        self.res=None

    def tearDown(self):
        #pass
        os.system('rm -rf *.cl')

        
    def test_default(self):
        '''predictcomp: test defaults'''
        self.assertRaises(Exception,predictcomp)
 
    def test_invalid_objname(self): 
        '''predictcomp: invalid objname'''
        self.assertRaises(Exception, predictcomp, objname='Moon', minfreq='100GHz',maxfreq='120GHz')
        
    def test_valid_objname(self):
        '''predictcomp: valid objname'''
        self.res=predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100GHz',maxfreq='120GHz',
                             standard='Butler-JPL-Horizons 2012') 
        print("type(self.res) = ",type(self.res))
        self.assertTrue(type(self.res)==dict)
        self.assertTrue(os.path.exists(self.res['clist']))
             
    def test_invalid_freqrange(self):
        '''predictcomp: invalid freqrange'''
        self.assertRaises(Exception, predictcomp, objname='Titan', epoch='2017/09/01/00:00', minfreq='100', maxfreq='120' )
    
    @unittest.skip('no plotting in this part of casa')
    def test_predicted_visplot(self):
        '''predictcomp: generate visibility plot for a given array configuration''' 
        self.res=predictcomp( objname='Titan', epoch='2017/09/01/00:00', minfreq='100GHz',
                              maxfreq='120GHz', standard='Butler-JPL-Horizons 2012',
                              antennalist=ctsys.resolve(os.path.join(datapath,'alma.cycle5.1.cfg')),
                              showplot=False,savefig='visplot.png' )
        self.assertTrue(type(self.res)==dict)
        self.assertTrue(os.path.exists(self.res['clist']))
        self.assertTrue(os.path.exists('visplot.png'))
          
    def test_valid_but_not_visible_objname(self):
        '''predictcomp: valid but not visible objname'''
        self.assertRaises( Exception, predictcomp, objname='Mars', epoch='2018/09/01/00:00', minfreq='100GHz', maxfreq='120GHz',
                           antennalist=ctsys.resolve(os.path.join(datapath,'alma.cycle5.1.cfg')), showplot=False )

def suite():
    return [predictcomp_test]

if __name__ == '__main__':
    unittest.main()
