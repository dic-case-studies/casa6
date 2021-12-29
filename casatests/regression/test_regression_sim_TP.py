#############################################################################
# $Id:$                                                                     #
# Test Name: test_regression_TP.py                                          #
#    Regression Test Script for simobserve/simanalyze                       #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    Test the use of simobserve and simanalyze on simdata of a 2d image     #
#    Exercise simobserve of total power simdata                             #
#    Original regression in casa 5 was m51_tpsim_regression.py               #
#                                                                           #
#                                                                           #
# Input data:                                                               #
#    simdata of M51 (ALMA-12m INT + ACA-7m INT + 12m TP)                    #
#                                                                           #
# CAS-13086 JIRA                                                            #
#                                                                           #
#############################################################################

import os
import time
import shutil
import unittest


from casatools import ctsys, image, ms
from casatasks import simobserve, simanalyze, casalog

_ia = image()
_ms = ms()

datadir = ctsys.resolve('regression/sim_TP/')
cfgdir = ctsys.resolve('alma/simmos/')

projname = "m51sd_co32"

def logprint(msg):
    print(msg)
    casalog.post(msg,origin='test_regression_sim_TP')

class regression_sim_TP_test(unittest.TestCase):

    def setUp(self):
        self.modelname = "m51ha.model"
        if os.path.exists(self.modelname):
            shutil.rmtree(self.modelname)

        shutil.copytree(datadir+self.modelname, self.modelname)

    def tearDown(self):
        shutil.rmtree(projname)
        shutil.rmtree(self.modelname)

    def test_regression(self):
        '''test total power simulations '''

        logprint('sd total power simobserve of M51')

        startTime = time.time()
        startProc = time.perf_counter()
        
        simobserve(project = projname, skymodel = self.modelname, inbright = '0.004', indirection = 'B1950 23h59m59.96 -34d59m59.50',
                   incell = '0.5arcsec',incenter = '330.076GHz' , inwidth = '50MHz',
                   setpointings = True,integration = '10s', 
                   mapsize = '',maptype = "square", pointingspacing = '9arcsec',obsmode = "sd",
                   sdantlist = cfgdir+'aca.tp.cfg',sdant = 0,  refdate='2012/11/21/20:00:00',
                   totaltime = '31360s', thermalnoise = 'tsys-manual',
                   t_ground = 263.0, t_sky = 263.0,
                   graphics="file", verbose=True, overwrite = True)

        obsEndTime = time.time()
        obsEndProc = time.perf_counter()

        logprint('simanalyze of total power (M51)')
       
        simanalyze(project=projname, image=True, imsize=[512,512], cell='1.0arcsec',
                   imdirection = 'B1950 23h59m59.96 -34d59m59.50', analyze=True,
                   showpsf = False, showresidual = False, showconvolved = True, graphics='file',
                   verbose=True)

        endTime = time.time()
        endProc = time.perf_counter()

        # Regression

        m51sd_im = _ia.open(projname+"/"+projname +'.sd.image')
        m51sd_stats = _ia.statistics(verbose=False,list=False)
        _ia.close()

        # 2015/10/23: update after the fix to rasterutil._get_sampling (r35037)
        #refstats = {'max': 5.01661158,
        #             'min': -0.80671942,
        #             'rms': 0.5335443,
        #             'sigma': 0.48867723,
        #             'sum':  56116.85862354}
        #2016/06/15 after fixing B1950 component corruption
        # refstats = {'max': 5.2856,
        #              'min': -0.79999,
        #              'rms': 0.56074,
        #              'sigma': 0.51274,
        #              'sum':  59485}
        #2016/07/08 after fixing minweight bug in sdimaging
        refstats = {'max': 5.2856,
                    'min': -0.79999,
                    'rms': 0.56827,
                    'sigma': 0.5176,
                    'sum':  59391}

        m51sd_diff_im = _ia.open(projname+"/"+projname + '.sd.diff')
        m51sd_diffstats = _ia.statistics(verbose=False,list=False)
        _ia.close()
        
        # 2015/10/23: update after the fix to rasterutil._get_sampling (r35037)
        #diffstats = {'max': 0.91776484,
        #             'min': -0.89943612,
        #             'rms': 0.20261368,
        #             'sigma': 0.20169726,
        #             'sum': 5044.85529928}
        #2016/06/15 after fixing 1950 component corruption
        # diffstats = {'max': 0.91084,
        #              'min': -0.90267,
        #              'rms': 0.19858,
        #              'sigma': 0.19848,
        #              'sum': 1676.6}
        #2016/07/08 after fixing minweight bug in sdimaging
        #diffstats = {'max': 0.91084,
        #             'min': -0.90267,
        #             'rms': 0.19583,
        #             'sigma': 0.19572,
        #             'sum': 1640.8}
        # CAS-13086 : 5.8 values - results are identical in 6.2
        diffstats = {'max': 0.88841,
                     'min': -0.8805,
                     'rms': 0.191,
                     'sigma': 0.19088,
                     'sum': 1702.3}
        # relative tolerances to reference values
        reftol   = {'sum':  1e-2,
                    'max':  1e-2,
                    'min':  5e-2,
                    'rms':  1e-2,
                    'sigma': 1e-2}

        logprint('')
        logprint('********** Regression *****************')
        logprint('')
        logprint('Regression results')


        # more info
        _ms.open(projname+"/"+projname+".aca.tp.sd.ms")
        logprint("Noiseless MS, amp stats:")
        logprint("%s" % _ms.statistics('DATA','amp'))
        logprint("Noiseless MS, phase stats:")
        logprint("%s" % _ms.statistics('DATA','phase'))
        _ms.close()
       
        regstate = True
        rskes = list(refstats.keys())
        rskes.sort()
        for ke in rskes:
            adiff=abs(m51sd_stats[ke][0] - refstats[ke])/abs(refstats[ke])
            if adiff < reftol[ke]:
                logprint("* Passed %-5s image test, got % -11.5g expected % -11.5g." % (ke, m51sd_stats[ke][0], refstats[ke]))
            else:
                logprint("* FAILED %-5s image test, got % -11.5g instead of % -11.5g." % (ke, m51sd_stats[ke][0], refstats[ke]))
                regstate = False

        rskes = list(diffstats.keys())
        rskes.sort()
        for ke in rskes:
            adiff=abs(m51sd_diffstats[ke][0] - diffstats[ke])/abs(diffstats[ke])
            if adiff < reftol[ke]:
                logprint("* Passed %-5s  diff test, got % -11.5g expected % -11.5g." % (ke, m51sd_diffstats[ke][0], diffstats[ke]))
            else:
                logprint("* FAILED %-5s  diff test, got % -11.5g instead of % -11.5g." % (ke, m51sd_diffstats[ke][0], diffstats[ke]))
                regstate = False

        logprint ('---')
        if regstate:
            logprint ('Regression PASSED')
        else:
            logprint ('Regression FAILED')
            
        logprint ('---')
        logprint ('*********************************')

        logprint('')
        logprint('********** Benchmarking **************')
        logprint('')
        logprint('Total wall clock time was: %8.3f s.' % (endTime - startTime))
        logprint('Total CPU        time was: %8.3f s.' % (endProc - startProc))
        logprint('Wall processing  rate was: %8.3f MB/s.' % (17896.0 /
                                                             (endTime - startTime)))

        ### Get last modification time of .ms.
        msfstat = os.stat(projname+"/"+projname+'.aca.tp.sd.ms')
        #msfstat = os.stat(projname+"/"+projname+'.aca.tp.ms')
        logprint('* Breakdown:                           *')
        logprint('*  generating visibilities took %8.3fs,' % (msfstat[8] - startTime))
        logprint('*************************************')

        logprint('--Finished simdata of M51 (total power) regression--')
        self.assertTrue(regstate)

def suite():
    return[regression_sim_TP_test]

if __name__ == '__main__':
    unittest.main()
