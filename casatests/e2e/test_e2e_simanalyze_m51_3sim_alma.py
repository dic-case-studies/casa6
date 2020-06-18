#############################################################################
# $Id:$                                                                     #
# Test Name: test_e2e_simanalyze_m51_3sim_alma.py                           #
#    Regression Test Script for simobserve/simanalyze                       #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    Test the use of simobserve and simanalyze on simdata of a 2d image     #
#    Exercise simobserve of IF, TP and ACA simdata                          #
#                                                                           #
#                                                                           #
#                                                                           #
# Input data:                                                               #
#    simdata of M51 (ALMA-12m INT + ACA-7m INT + 12m TP)                    #
#                                                                           #
#############################################################################
 
import os
import shutil
import unittest


CASA6 = False
try:
    from casatools import ctsys
    from casatasks import simobserve, simanalyze
    CASA6 = True
    
    def default(atask):
        pass
except ImportError:
    from tasks import simobserve, simanalyze
    from __main__ import default

if CASA6:
    ctsys_resolve = ctsys.resolve
    datadir = ctsys.resolve('regression/simdata/')
    cfgdir = ctsys.resolve('alma/simmos/')
    regdir = ctsys.resolve('regression/sim_m51c/reference/m51c/')

else:
    repodir = os.path.join(os.environ['CASAPATH'].split()[0],'data/')
    datadir = repodir + 'regression/simdata/'
    cfgdir = repodir + 'alma/simmos/'
    regdir = repodir + 'regression/sim_m51c/reference/m51c/'

from casatestutils import testhelper as th

noise = False # add noise
my_verbose = True

projname = "m51c"

class regression_m51_3sim_test(unittest.TestCase):

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
        print('12m - Interferometry simobserve')

        if noise:
            thermalnoise = 'tsys-atm'  #w/ noise 
#            user_pwv=3.0
        else:
            thermalnoise=""

        simobserve(project = projname, skymodel = self.modelname, inbright = '0.004', indirection = 'B1950 23h59m59.96 -34d59m59.50',
                   incell = '0.1arcsec',incenter = '330.076GHz' , inwidth = '50MHz', setpointings = True,integration = '10s',
                    mapsize = '1arcmin',maptype = "hex", pointingspacing = '9arcsec',obsmode = "int", refdate='2012/11/21/20:00:00',
                    totaltime = '3600s', antennalist="alma;0.5arcsec", thermalnoise = thermalnoise, 
                    graphics="file", verbose=False, overwrite = True)


        ############################ 12m TP  ############################
        print('12m - Total Power simobserve')
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
                    graphics="file", verbose=False, overwrite = True)
                


        ############################ ACA  ############################
        print('ACA - simobserve')
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
                    verbose=False, overwrite = True)


        ############################ clean ACA with SD model  ############################
        print('clean ACA with SD model')
        default(simanalyze)
        if noise:
            myvis = '$project.aca.i.noisy.ms,$project.aca.tp.sd.noisy.ms'  #w/ noise
        else:
            myvis = '$project.aca.i.ms,$project.aca.tp.sd.ms'  #w/ noise

        simanalyze(project=projname, vis=myvis, image=True, imsize = [512,512], cell = '0.2arcsec', modelimage='$project.sd.image',
                   analyze = True, showpsf = False, showresidual = False, showconvolved = True, graphics="file")


        ############################ clean ALMA with ACA+SD model  ############################
        print('clean ACA with SD model')
        default(simanalyze)
        if noise:
            myvis = '$project.alma_0.5arcsec.noisy.ms'
        else:
            myvis = '$project.alma_0.5arcsec.ms'

        simanalyze(project=projname, vis=myvis, image=True, imsize = [512,512], cell = '0.2arcsec', modelimage="$project.aca.i.image",
                   analyze = True, showpsf = False, showresidual = False, showconvolved = True, graphics="file")


        # Regression
        print('Regression Results')

        # It seems that it will compare what is created by the regression inside m51c with
        # what is available in casa-data/regression/sim_m51c/reference/m51c
        
        regstate=True
        verbose=my_verbose
        
        # test SD first so that we can check that indep. from INT
        # Compare  newMS with templateMS
        newMS=project+"/"+project+".aca.tp.sd.ms"
        if verbose: print (newMS)
        # templateMS is /casa-data/regression/sim_m51c/reference/m51c/m51c.aca.tp.sd.ms
        templateMS = regdir + project + ".aca.tp.sd.ms"        
        regstate=regstate and th.compMS(newMS,templateMS,verbose=verbose)
        
        # Compare newImage with templateImage
        newImage=project+"/"+project + '.sd.image'
        if verbose: print (newImage)
        # templateImage is /casa-data/regression/sim_m51c/reference/m51c/m51c.sd.image
        templateImage = regdir + project + ".sd.image"
        regstate=regstate and th.compImages(newImage,templateImage,verbose=verbose)
        
        # then INT
        # Compare  newMS with templateMS
        newMS=project+"/"+project+".alma_0.5arcsec.ms"
        if verbose: print (newMS)
        templateMS = regdir + project + ".alma_0.5arcsec.ms"
        regstate=regstate and th.compMS(newMS,templateMS,verbose=verbose)
        
        # Compare newImage with tempateImage
        newImage=project+"/"+project + '.alma_0.5arcsec.image'
        if verbose: print (newImage)
        templateImage = regdir + project + ".alma_0.5arcsec.image"
        regstate=regstate and th.compImages(newImage,templateImage,verbose=verbose)
        
        # Compare newImage with templateImage
        newImage=project+"/"+project + '.alma_0.5arcsec.diff'
        if verbose: print (newImage)
        templateImage = regdir + project + '.alma_0.5arcsec.diff'
        regstate=regstate and th.compImages(newImage,templateImage,verbose=verbose)


        if regstate:
        #    print >> logfile, 'Passed',
            print ('')
            print ('Regression PASSED')
            print ('')
        else:
        #    print >> logfile, 'FAILED',
            print ('')
            print ('Regression FAILED')
            print ('')

def suite():
    return[regression_m51_3sim_test]


if __name__ == '__main__':
    unittest.main()


