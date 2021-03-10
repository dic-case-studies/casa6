#############################################################################
# $Id:$                                                                     #
# Test Name: test_regression_sim_components_and_skymodel.py                 #
#    Regression Test Script for simobserve/simanalyze                       #
#                                                                           #
# Rationale for Inclusion:                                                  #
#                                                                           #
#                                                                           #
# Input data:                                                               #
#    protoplanetary disks                                                   #
#                                                                           #
#############################################################################

import os
import time
import shutil
import unittest

CASA6 = False
try:
    from casatools import ctsys, image, componentlist
    from casatasks import simobserve, simanalyze, casalog
    CASA6 = True
    
    _cl = componentlist()
    _ia = image()
    def default(atask):
        pass
except ImportError:
    from tasks import simobserve, simanalyze
    from taskinit import cltool, iatool, casalog
    from __main__ import default
    
    _cl = cltool()
    _ia = iatool()

if CASA6:
    datadir = ctsys.resolve('regression/simdata/')
    cfgdir = ctsys.resolve('alma/simmos/')

else:
    repodir = os.path.join(os.environ['CASAPATH'].split()[0],'data/')
    datadir = repodir + 'regression/simdata/'
    cfgdir = repodir + 'alma/simmos/'

def logprint(msg):
    print(msg)
    casalog.post(msg,origin='test_regression_sim_components_and_skymodel')

logprint ('--Running simdata of input672GHz_50pc.image--')

#projname = "psim2"
#my_modelimage="diskmodel.im2"
class regression_components_skymodel_test(unittest.TestCase):

    def setUp(self):
        default("simanalyze")
        default("simobserve")
        self.projname = "psim2"
        self.modelimage = datadir + "/input50pc_672GHz.fits"
        self.antennalist = cfgdir + "/alma.out20.cfg"
        self.complist = "star672GHz.cl"
        if os.path.exists(self.complist):
            shutil.rmtree(self.complist)

        # NOTE: it will create the modelimage in the local directory, not in the project directory
        self.direction = "J2000 18h00m00.031s -22d59m59.6s"
        _cl.done()
        _cl.addcomponent(dir=self.direction,flux=0.0003,freq="672GHz")
        _cl.rename("star672GHz.cl")
        _cl.done()
                           
#        shutil.copytree(datadir+self.modelimage, self.modelimage)

 
    def tearDown(self):
        if os.path.exists(self.projname):
            shutil.rmtree(self.projname)
        if os.path.exists(self.complist):
            shutil.rmtree(self.complist)

    def test_regression(self):
        '''Test components and skymodel...'''

        startTime = time.time()
        startProc = time.clock()
        
        # Check if modelimage and antennalist are being read from original repository
        # ground must be set here to use appropriate defaults for "tsys-atm" when run without casalith sub-parameter defaults
        # tsys-atm -> t_ground = 269.0, usr_pwv=0.5, seed=11111 : seed already defaults to that value
        simobserve(project=self.projname, skymodel=self.modelimage, direction=self.direction, complist=self.complist,
                   setpointings=True, mapsize="0.76arcsec", pointingspacing="0.5arcsec", integration="10s", obsmode="int",
                   antennalist=self.antennalist, refdate="2012/06/21/03:25:00", totaltime="1200s",
                   thermalnoise="tsys-atm", user_pwv=0.5, t_ground=269.0, maptype="ALMA2012",  graphics="file", verbose=True)

        # Analyze and clean  results of project created above
        default(simanalyze)
        # skymodel is NOT set in log for ppdisk2_regression, which I don't understand, it doesn't matter as it picks up the correct one in psim2
        simanalyze(project=self.projname, skymodel="", image=True, cell="", niter=1000, 
                   threshold="1e-7Jy", imsize=[192, 192], stokes="I", weighting="natural", analyze=True, interactive=False,
                   graphics="file", verbose=True)

        endTime = time.time()
        endProc = time.clock()

        # Regression

        test_name_ppd = """simdata observation of Wolf & D'Angelo's protoplanetary disk"""

        ppdso_im= _ia.open(self.projname+"/"+self.projname + '.alma.out20.noisy.image')
        ppdso_stats = _ia.statistics(list=True, verbose=True)
        _ia.close()

        #refstats = { 'flux': 0.0359,
        #             'max': 6.1e-04,
        #             'min': -1.6e-04, 
        #             'rms': 1.9e-04,
        #             'sigma': 1.4e-04 }
        # r38847
        # refstats = { 'flux': 0.036714,
        #            'max': 0.00061958,
        #            'min': -0.0001894,
        #            'rms': 0.0001925,
        #            'sigma': 0.00013576 }
        # CAS-13086 : 5.8
        refstats = { 'flux': 0.037245,
                     'max': 0.00063979,
                     'min': -0.00019047,
                     'rms': 0.00019654,
                     'sigma': 0.00014007 }

        reftol   = {'flux':  0.01,
                    'max':   0.02,
                    'min':   0.02,
                    'rms':   0.01,
                    'sigma': 0.01}

        import datetime
        datestring = datetime.datetime.isoformat(datetime.datetime.today())        
        
        loghdr = """
        ********** Simulation Summary *********
        
        The disk input image is a simulation done by Wolf and D'Angelo, converted from
        900 GHz to 672 GHz
        
        ********** Regression *****************
        """
                
        regstate = True
        rskes = list(refstats.keys())
        rskes.sort()
        for ke in rskes:
            adiff=abs(ppdso_stats[ke][0] - refstats[ke])/abs(refstats[ke])
            if adiff < reftol[ke]:
                logprint("* Passed %-5s test, got % -11.5g , expected % -11.5g." % (ke, ppdso_stats[ke][0], refstats[ke]))
            else:
                logprint("* FAILED %-5s test, got % -11.5g instead of % -11.5g." % (ke, ppdso_stats[ke][0], refstats[ke]))
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
        
        logprint ('--Finished simdata of input672GHz_50pc.image regression--')
        self.assertTrue(regstate)
        
def suite():
    return[regression_components_skymodel_test]


if __name__ == '__main__':
    unittest.main()
       
        
        
