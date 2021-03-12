#############################################################################
# $Id:$                                                                     #
# Test Name: test_regression_sim_components.py                              #
#    Regression Test Script for simobserve/simanalyze                       #
#                                                                           #
# Rationale for Inclusion:                                                  #
#  Original regression in casa 5 was testcompsim.py                         #
#                                                                           #
# Input data:                                                               #
#  sim_complist_6334.cl component list                                      #
#                                                                           #
# CAS-13086 JIRA                                                            #
#                                                                           #
#############################################################################

import shutil
import unittest
import os

CASA6 = False
try:
    from casatools import ctsys
    from casatasks import simobserve, simanalyze, casalog
    CASA6 = True

    def default(atask):
        pass
except ImportError:
    from tasks import simobserve, simanalyze
    from taskinit import casalog
    from __main__ import default

if CASA6:
    datadir = ctsys.resolve('regression/sim_components/')
    refdir = ctsys.resolve('regression/sim_components/sim_reference/')
else:
    repodir = os.path.join(os.environ['CASAPATH'].split()[0],'casatestdata/')
    datadir = repodir + 'regression/sim_components/'
    refdir = repodir + 'regression/sim_components/sim_reference/'

from casatestutils import testhelper as th

projname = "cltest"

def logprint(msg):
    print(msg)
    casalog.post(msg,origin='test_regression_sim_components')

class regression_sim_components_test(unittest.TestCase):

    def setUp(self):
        self.complist = 'sim_complist_6334.cl'
        if os.path.exists(self.complist):
            shutil.rmtree(self.complist)

        shutil.copytree(datadir+self.complist, self.complist)

    def tearDown(self):
        if os.path.exists(self.complist):
            shutil.rmtree(self.complist)

        if os.path.exists(projname):
            shutil.rmtree(projname)

    def test_regression(self):

        default("simobserve")
        simobserve(project=projname, complist=self.complist, compwidth="1.875GHz", setpointings=True, integration="10s",
                   direction="J2000 17h20m53.2s -35d47m00s", mapsize="13arcsec",maptype="ALMA", pointingspacing ="",
                   obsmode="int", refdate="2014/01/01", hourangle="transit",totaltime="7200s", antennalist="alma.cycle0.extended.cfg",
                   thermalnoise="",graphics='file', verbose=True, overwrite=True)

        default('simanalyze')
        simanalyze(project=projname, image=True, vis="default", imsize=300, imdirection="J2000 17h20m53.2s -35d47m00s",
                   cell="0.043294054arcsec", niter=500, threshold="0.1mJy", analyze=True,
                   graphics='file', verbose=True, overwrite=True)

        logprint('Regression results')
        
        regstate=True
        verbose=True
        
        img_suffix = projname+ '.alma.cycle0.extended.image'
        newImage=projname+"/" + img_suffix
        if verbose: logprint(newImage)
        templateImage = refdir + img_suffix
        regstate=regstate and th.compImages(newImage,templateImage,verbose=verbose)
        
        diff_suffix = projname+ '.alma.cycle0.extended.diff'        
        newImage=projname+"/" + diff_suffix
        if verbose: logprint(newImage)
        templateImage = refdir + diff_suffix
        regstate=regstate and th.compImages(newImage,templateImage,verbose=verbose)

        ms_suffix = projname+ '.alma.cycle0.extended.ms'
        newMS = projname+"/"+ ms_suffix
        if verbose: logprint(newMS)
        templateMS = refdir + ms_suffix
        regstate = regstate and th.compMS(newMS,templateMS,verbose=verbose)

        logprint('')
        if regstate:
            logprint('Regression PASSED')
        else:
            logprint('Regression FAILED')
            
        logprint('')
        
        self.assertTrue(regstate)

def suite():
    return[regression_sim_components_test]

if __name__ == '__main__':
    unittest.main()
            

