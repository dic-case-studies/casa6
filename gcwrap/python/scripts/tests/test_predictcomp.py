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

    dataRoot = os.path.join(os.environ.get('CASAPATH'.split()[0],'data'))
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

datapath = 'alma/simmos'


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
        if is_CASA6:
            # casatasks throw exceptions
            self.assertRaises(Exception,predictcomp)
        else:
            # CASA 5 tasks do not
            self.res=predictcomp() 
            self.assertIsNone(self.res)
 
    def test_invalid_objname(self): 
        '''predictcomp: invalid objname'''
        if is_CASA6:
            # casatasks throw exceptions
            self.assertRaises(Exception, predictcomp, objname='Moon', minfreq='100GHz',maxfreq='120GHz')
        else:
            # CASA 5 tasks do not
            self.res=predictcomp(objname='Moon', minfreq='100GHz',maxfreq='120GHz') 
            self.assertIsNone(self.res)

    def test_valid_objname(self):
        '''predictcomp: valid objname'''
        self.res=predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100GHz',maxfreq='120GHz',
                             standard='Butler-JPL-Horizons 2012') 
        print("type(self.res) = ",type(self.res))
        self.assertTrue(type(self.res)==dict)
        self.assertTrue(os.path.exists(self.res['clist']))
             
    def test_invalid_freqrange(self):
        '''predictcomp: invalid freqrange'''
        if is_CASA6:
            # casatasks throw exceptions
            self.assertRaises(Exception, predictcomp, objname='Titan', epoch='2017/09/01/00:00', minfreq='100', maxfreq='120' )
        else:
            # CASA 5 tasks do not
            self.res=predictcomp(objname='Titan', epoch='2017/09/01/00:00', minfreq='100',maxfreq='120')
            self.assertIsNone(self.res)
    
    def test_predicted_visplot(self):
        # no plotting in this part of casa, skip this test
        if is_CASA6:
            casalog.post('test_predicted_visplot SKIPPED: no plotting in this part of casa','INFO')
        else:
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
        if is_CASA6:
            # casatasks throw exceptions
            self.assertRaises( Exception, predictcomp, objname='Mars', epoch='2018/09/01/00:00', minfreq='100GHz', maxfreq='120GHz',
                               antennalist=ctsys.resolve(os.path.join(datapath,'alma.cycle5.1.cfg')), showplot=False )
        else:
            # CASA 5 tasks do not
            self.res=predictcomp(objname='Mars', epoch='2018/09/01/00:00', minfreq='100GHz',maxfreq='120GHz',
                                 antennalist=datapath+'alma.cycle5.1.cfg',showplot=False) 
            self.assertIsNone(self.res)

def suite():
    return [predictcomp_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()
