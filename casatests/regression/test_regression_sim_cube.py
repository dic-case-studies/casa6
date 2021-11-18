#############################################################################
# $Id:$                                                                     #
# Test Name: test_regression_sim_cube.py                                    #
#    Regression Test Script for simobserve                                  #
#                                                                           #
# Rationale for Inclusion:                                                  #
#  Original regression in casa 5 was testcube2_regression.py                #
#                                                                           #
# Input data:                                                               #
#  sim_alma_3dcube_128x128x10.fits : fits skymodel image                    #
#  sim_alma_3dcube_128x128x10.txt : ptg file                                #
#                                                                           #
# CAS-13086 JIRA                                                            #
#                                                                           #
#############################################################################

import os
import shutil
import numpy
import time
import unittest

from casatools import ctsys, ms
from casatasks import simobserve, importfits, casalog
_ms = ms()

datadir = ctsys.resolve("regression/sim_cube/")
cfgdir = ctsys.resolve("alma/simmos/")

projname = "tc2"


def logprint(msg):
    print(msg)
    casalog.post(msg,origin='test_regression_sim_cube')

class regression_sim_cube_test(unittest.TestCase):

    def setUp(self):
        if os.path.exists("testcube2"):
            shutil.rmtree("testcube2")
        importfits(fitsimage=datadir+'sim_alma_3dcube_128x128x10.fits',imagename="testcube2")

    def tearDown(self):
        if os.path.exists("testcube2"):
            shutil.rmtree("testcube2")

        if os.path.exists(projname):
            shutil.rmtree(projname)

    def test_regression(self):

        logprint("simobserve of test cube")

        startTime = time.time()
        startProc = time.perf_counter()

        simobserve(project=projname, skymodel="testcube2",inbright=".1",indirection="J2000 19h00m00s -40d00m00s",
                   incell="0.2arcsec",incenter="350GHz",inwidth="0.5MHz",setpointings=False,ptgfile=datadir+"sim_alma_3dcube_128x128x10.txt",
                   obsmode="int",antennalist=cfgdir+"alma.out01.cfg",refdate="2012/06/21/03:25:00",
                   totaltime="7200s",thermalnoise="",graphics="file",verbose=True,overwrite=True)

        endTime = time.time()
        endProc = time.perf_counter()

        # Regression
        logprint("********** Regression *****************")
        
        _ms.open(projname+"/"+projname+".alma.out01.ms")
        newdata = _ms.getdata(items="data")['data']
        _ms.close()

        refshape=[2,10,882000]
        
        refstats = { 'max': 2.05e-01 +  7.52e-03j,
                     'min':-1.90e-01 +  4.33e-02j,
                     'sum': 1.72e+04 + -1.53e+03j,
                     'std': 5.53e-02 }

        ### tight
        reftol   = {'max':  5e-3,
                    'min':  5e-3,
                    'sum':  5e-3,
                    'std':  5e-3}

        regstate = True
        logprint("Regression results")

        if max(abs(newdata.shape-numpy.array(refshape)))<=0:
            logprint("* Passed shape test with shape "+str(newdata.shape))
        else:
            logprint("* FAILED shape test, expecting %s, got %s" % (str(refshape),str(newdata.shape)))
            regstate = False

        cube_stats={'max':newdata.max(),
                    'min':newdata.min(),
                    'sum':newdata.sum(),
                    'std':newdata.std()}

        rskes = list(refstats.keys())
        rskes.sort()
        for ke in rskes:
            adiff=abs(cube_stats[ke] - refstats[ke])/abs(refstats[ke])
            if adiff < reftol[ke]:
                status="* Passed "
            else:
                status="* FAILED "
                regstate = False
            status=status+" %3s test, got " % ke
            if type(refstats[ke])==complex:
                status=status+"%9.2e + %9.2ej , expected %9.2e + %9.2ej." % (cube_stats[ke].real, cube_stats[ke].imag, refstats[ke].real, refstats[ke].imag)
            else:
                status=status+"%9.2e          , expected %9.2e." % (cube_stats[ke], refstats[ke])
            logprint(status)

        logprint("---")
        logprint('')
        if regstate:
            logprint('Regression PASSED')
        else:
            logprint('Regression FAILED')
        logprint('')
        logprint('---')

        logprint ('*********************************')
        logprint('')
        logprint('********** Benchmarking **************')
        logprint('')
        logprint('Total wall clock time was: %8.3f s.' % (endTime - startTime))
        logprint('Total CPU        time was: %8.3f s.' % (endProc - startProc))
        logprint('Wall processing  rate was: %8.3f MB/s.' % (17896.0 /
                                                             (endTime - startTime)))
        logprint('*************************************')
        
        logprint ('--Finished  simdata of test cube regression--')
        self.assertTrue(regstate)

def suite():
    return[regression_sim_cube_test]

if __name__ == '__main__':
    unittest.main()
