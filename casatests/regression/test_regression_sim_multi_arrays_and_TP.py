#############################################################################
# $Id:$                                                                     #
# Test Name: test_regression_sim_multi_arrays_and_TP.py                     #
#    Regression Test Script for simobserve/simanalyze                       #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    Test the use of simobserve and simanalyze on simdata of a 2d image     #
#    Exercise simobserve of IF, TP and ACA simdata                          #
#    Original regression in casa 5 was m51_3sim_regression.py               #
#                                                                           #
# Input data:                                                               #
#    m51ha.model : model for simdata of                                     #
#                  M51 (ALMA-12m INT + ACA-7m INT + 12m TP)                 #
#                                                                           #
# CAS-13086 JIRA                                                            #
#                                                                           #
#                                                                           #
#############################################################################
 
import os
import shutil
import unittest


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
    datadir = ctsys.resolve('regression/sim_multi_arrays_and_TP/')
    cfgdir = ctsys.resolve('alma/simmos/')
    refdir = ctsys.resolve('regression/sim_multi_arrays_and_TP/m51c_reference/')

else:
    repodir = os.path.join(os.environ['CASAPATH'].split()[0],'casatestdata/')
    datadir = repodir + 'regression/sim_multi_arrays_and_TP/'
    cfgdir = repodir + 'alma/simmos/'
    refdir = repodir + 'regression/sim_multi_arrays_and_TP/m51c_reference/'

from casatestutils import testhelper as th

noise = False # add noise
my_verbose = True

projname = "m51c"

def logprint(msg):
    print(msg)
    casalog.post(msg,origin='test_regression_sim_multi_arrays_and_TP')

class regression_sim_multiarrays_test(unittest.TestCase):

    def setUp(self):
        self.modelname="m51ha.model"
        if os.path.exists(self.modelname):
            shutil.rmtree(self.modelname)
                    
        shutil.copytree(datadir+self.modelname, self.modelname)

 
    def tearDown(self):
        shutil.rmtree(projname)
        shutil.rmtree(self.modelname)

    def test_regression(self):
        '''test 12m IF, TP'''
        
        ############################  12m IF  ############################
        logprint('12m - Interferometry simobserve')

        if noise:
            thermalnoise = 'tsys-atm'  #w/ noise 
#            user_pwv=3.0
        else:
            thermalnoise=""


        simobserve(project = projname, skymodel = self.modelname, inbright = '0.004', indirection = 'B1950 23h59m59.96 -34d59m59.50',
                   incell = '0.1arcsec',incenter = '330.076GHz' , inwidth = '50MHz', setpointings = True,integration = '10s',
                    mapsize = '1arcmin',maptype = "hex", pointingspacing = '9arcsec',obsmode = "int", refdate='2012/11/21/20:00:00',
                    totaltime = '3600s', antennalist="alma;0.5arcsec", thermalnoise = thermalnoise, 
                    graphics="file", verbose=my_verbose, overwrite = True)


        ############################ 12m TP  ############################
        logprint('12m - Total Power simobserve')
        default(simobserve)
        project = projname
        
        if noise:
            thermalnoise = 'tsys-atm'  #w/ noise 
#            user_pwv=3.0
        else:
            thermalnoise=""

        simobserve(project = projname, skymodel = self.modelname, inbright = '0.004', indirection = 'B1950 23h59m59.96 -34d59m59.50',
                   incell = '0.1arcsec',incenter = '330.076GHz' , inwidth = '50MHz', setpointings = True,integration = '10s',
                    mapsize = '1arcmin',maptype = "square", pointingspacing = '9arcsec',obsmode = "sd", refdate='2012/11/21/20:00:00',
                    totaltime = '2h', sdantlist = cfgdir+'aca.tp.cfg', sdant = 0, thermalnoise = thermalnoise, 
                    graphics="file", verbose=my_verbose, overwrite = True)
                


        ############################ ACA  ############################
        logprint('ACA - simobserve')
        default(simobserve)
        if noise:
            thermalnoise = 'tsys-atm'  #w/ noise 
#            user_pwv=3.0
        else:
            thermalnoise=""

        simobserve(project = projname, skymodel = self.modelname, inbright = '0.004', indirection = 'B1950 23h59m59.96 -34d59m59.50',
                   incell = '0.1arcsec',incenter = '330.076GHz' , inwidth = '50MHz', setpointings = True, integration = '10s',
                    mapsize = '1arcmin',maptype = "hex", pointingspacing = '15arcsec',obsmode = "int", refdate='2012/11/21/20:00:00',
                    totaltime = '3', antennalist="aca.i.cfg", thermalnoise = thermalnoise, graphics="file", 
                    verbose=my_verbose, overwrite = True)


        ############################ clean ACA with SD model  ############################
        logprint('clean ACA with SD model')
        default(simanalyze)
        if noise:
            myvis = '$project.aca.i.noisy.ms,$project.aca.tp.sd.noisy.ms'  #w/ noise
        else:
            myvis = '$project.aca.i.ms,$project.aca.tp.sd.ms'  #w/ noise

        simanalyze(project=projname, vis=myvis, image=True, imsize = [512,512], cell = '0.2arcsec', modelimage='$project.sd.image',
                   analyze = True, showpsf = False, showresidual = False, showconvolved = True, graphics="file")


        ############################ clean ALMA with ACA+SD model  ############################
        logprint('clean ACA with SD model')
        default(simanalyze)
        if noise:
            myvis = '$project.alma_0.5arcsec.noisy.ms'
        else:
            myvis = '$project.alma_0.5arcsec.ms'

        simanalyze(project=projname, vis=myvis, image=True, imsize = [512,512], cell = '0.2arcsec', modelimage="$project.aca.i.image",
                   analyze = True, showpsf = False, showresidual = False, showconvolved = True, graphics="file")

        # Regression
        logprint('')
        logprint('********************')
        logprint('')
        logprint('Regression Results')

        # It seems that it will compare what is created by the regression inside m51c with
        # what is available in casa-data/regression/sim_m51c/reference/m51c
        
        regstate=True
        
        # test SD first so that we can check that indep. from INT
        # Compare  newMS with templateMS
        newMS=project+"/"+project+".aca.tp.sd.ms"
        if my_verbose: logprint (newMS)
        templateMS = refdir + project + ".aca.tp.sd.ms"
        regstate=regstate and th.compMS(newMS,templateMS,verbose=my_verbose)
        
        # Compare newImage with templateImage
        newImage=project+"/"+project + '.sd.image'
        if my_verbose: logprint (newImage)
        templateImage = refdir + project + ".sd.image"
        regstate=regstate and th.compImages(newImage,templateImage,verbose=my_verbose)
        
        # then INT
        # Compare  newMS with templateMS
        newMS=project+"/"+project+".alma_0.5arcsec.ms"
        if my_verbose: logprint (newMS)
        templateMS = refdir + project + ".alma_0.5arcsec.ms"
        regstate=regstate and th.compMS(newMS,templateMS,verbose=my_verbose)
        
        # Compare newImage with tempateImage
        newImage=project+"/"+project + '.alma_0.5arcsec.image'
        if my_verbose: logprint (newImage)
        templateImage = refdir + project + ".alma_0.5arcsec.image"
        regstate=regstate and th.compImages(newImage,templateImage,verbose=my_verbose)
        
        # Compare newImage with templateImage
        newImage=project+"/"+project + '.alma_0.5arcsec.diff'
        if my_verbose: logprint (newImage)
        templateImage = refdir + project + '.alma_0.5arcsec.diff'
        regstate=regstate and th.compImages(newImage,templateImage,verbose=my_verbose)


        logprint ('')
        if regstate:
            logprint ('Regression PASSED')
        else:
            logprint ('Regression FAILED')
            
        logprint('')
        logprint('********************')
        logprint('')

        self.assertTrue(regstate)

def suite():
    return[regression_sim_multiarrays_test]


if __name__ == '__main__':
    unittest.main()


