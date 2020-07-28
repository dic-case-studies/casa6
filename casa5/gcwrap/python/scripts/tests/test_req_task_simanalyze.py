#test_req_task_simanalyze.py
#
# Copyright (C) 2020
# Associated Universities, Inc. Washington DC, USA.
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# https://open-jira.nrao.edu/browse/CAS-3669
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_simanalyze/about
#
# Test case: requirement
#
##########################################################################
 
 
####    Imports     ####
import os
import unittest
import shutil

try:
    from casatasks.private.casa_transition import is_CASA6
except ImportError:
    is_CASA6 = False

if is_CASA6:
    import casatools # not a good idea inside the casashell...perhaps os.path.exists(casatools.__file__) instead
    from casatasks import tclean, simobserve, simanalyze 
    from casatasks import casalog
else:
    from __main__ import default
    from tasks import *
    from taskinit import *
 
# DATA #
if is_CASA6:
    dataroot = casatools.ctsys.resolve()
    configpath_int  = casatools.ctsys.resolve((os.path.join(dataroot, 'alma/simmos/vla.a.cfg')))
    imagepath_int   = casatools.ctsys.resolve((os.path.join(dataroot, 'nrao/VLA/CalModels/3C286_Q.im/')))
    configpath_sd   = casatools.ctsys.resolve((os.path.join(dataroot, 'alma/simmos/aca.tp.cfg')))
    mspath_sd       = casatools.ctsys.resolve((os.path.join(dataroot, 'regression/unittest/clean/refimager/refim_twopoints_twochan.ms')))
else:
    dataroot = os.environ.get('CASAPATH').split()[0]
    configpath_int  = os.path.join(dataroot, 'data/alma/simmos/vla.a.cfg')
    imagepath_int   = os.path.join(dataroot, 'data/nrao/VLA/CalModels/3C286_Q.im/')
    configpath_sd   = os.path.join(dataroot, 'data/alma/simmos/aca.tp.cfg')
    mspath_sd       = os.path.join(dataroot, 'data/regression/unittest/clean/refimager/refim_twopoints_twochan.ms')

logpath = casalog.logfile()

int_project = 'sim_interferometric'
sd_project = 'sim_single_dish'
imagepath_sd = 'sd_model.image'
#both_project = 'sim_both'

####    interferometric    ####
class simanalyze_main_usage_modes_test_int(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not is_CASA6: # needs similar branch condition for casalith
            default(simobserve)
            default(tclean)

        simobserve(project=int_project, skymodel=imagepath_int, complist='', setpointings=True, direction=[], 
                   mapsize='', maptype='square', pointingspacing='',caldirection='',calflux='1Jy',obsmode='int', 
                   refdate='2020/02/13', hourangle='transit', totaltime='100s', antennalist=configpath_int,
                   outframe='LSRK', thermalnoise='', leakage=0.0, graphics='none',verbose=False, overwrite=False)

    @classmethod
    def tearDownClass(cls):
        os.system("rm -rf {}".format(int_project))

    def setUp(self):
        if not is_CASA6: # needs similar branch condition for casalith
            default(simanalyze)

    def tearDown(self):
        os.system("rm -rf *.last")

    def test_imaging_False_analysis_False(self):
        '''test_imaging_False_analysis_False: 
        -------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        The imagename parameter takes the name of an already synthesized image
        '''

        visname = str(int_project +'/'+ int_project + '.' + configpath_int.split('/')[-1][:-3] +'ms')

        try:
            simanalyze(project=int_project, image=False, vis=visname,
                       imagename=imagepath_int, analyze=False,
                       graphics='none', verbose=False, overwrite=True, dryrun=False,
                       logfile=logpath)
        except Exception:
            self.fail()
        
    @unittest.skip("ModuleNotFoundError: No module named clean")
    def test_imaging_True_interferometric_analysis_False(self):
        '''test_imaging_True_interferometric_analysis_False:
        ----------------------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        '''

        visname_int = str(int_project +'/'+ int_project + '.' + configpath_int.split('/')[-1][:-3] +'ms')

        simanalyze(project=int_project, image=True, vis=visname_int, modelimage='',imsize = [],
                   imdirection ='',cell = '',interactive = False, niter = 0,threshold = '0.1mJy',
                   weighting = 'natural',mask = [],outertaper = [],pbcor = False,stokes = 'I', 
                   featherimage = '',analyze=False, graphics='none', verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        self.assertTrue(os.path.isdir(visname_int[:-2]+'image'))

    @unittest.skipUnless(os.path.isdir(int_project+'/'+int_project+'.'+configpath_int.split('/')[-1][:-3]+'image'),"Analysis-only mode presumes the existence of some image data, such as that generated by test_imaging_True_interferometric*")
    def test_imaging_False_interferometric_analysis_True_showfidelity_True(self):
        '''test_imaging_False_interferometric_analysis_True_showfidelity_True:
        ----------------------------------------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        The showfidelity parameter displays the fidelity image
        '''

        visname_int = str(int_project +'/'+ int_project + '.' + configpath_int.split('/')[-1][:-3] +'ms')

        simanalyze(project=int_project, image=False, analyze=True, showfidelity=True,graphics='none', 
                   verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        # confirm that fidelity image was generated
        self.assertTrue(os.path.isdir(visname_int[:-2]+'fidelity')) 

####    Single Dish     ####
class simanalyze_main_usage_modes_test_sd(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        '''Called before tests in this class of test cases are run. Executes simobserve to create expected directory structure. Minimally, f"{project}/{project}.{suffix}" where suffix in ['skymodel','newmodel','compskymodel'].'''
        if not is_CASA6: # needs similar branch condition for casalith
            default(simobserve)
            default(tclean)

        # create reference image > 2.5*PB to use for SD sim
        tclean(vis=mspath_sd, imagename=imagepath_sd.split('.')[0],
               specmode='cube',nchan=1,imsize=750,cell='20arcsec',niter=100,
               gridder='wproject', wprojplanes=5,pblimit=-0.1, 
               phasecenter='J2000 19h59m28.449 40d44m01.199')

        # now that a model image exists, we can call simobserve to set up the directory structure
        simobserve(project=sd_project, skymodel=imagepath_sd, complist='', setpointings=True, direction=[], 
                   mapsize='', maptype='square', pointingspacing='',caldirection='',calflux='1Jy',obsmode='sd', 
                   refdate='2020/02/13', hourangle='transit', totaltime='100s', 
                   antennalist=configpath_sd, sdantlist=configpath_sd, sdant=0, 
                   outframe='LSRK', thermalnoise='', leakage=0.0, graphics='none',verbose=False, overwrite=False)

    @classmethod
    def tearDownClass(cls):
        os.system("rm -rf {}*".format(imagepath_sd[:-5]))
        os.system("rm -rf {}".format(sd_project))

    def setUp(self):
        if not is_CASA6: # needs similar branch condition for casalith
            default(simanalyze)

    def tearDown(self):
        os.system("rm -rf *.last")

    def test_imaging_True_single_dish_analysis_False(self):
        '''test_imaging_True_single_dish_analysis_False:
        ------------------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        The vis parameter can accept one or more MSs that can be interferometric or single dish data
        Modelimage will not be used if the MS is in total power
        '''

        visname_sd = str(sd_project +'/'+ sd_project + '.' + configpath_sd.split('/')[-1][:-3] +'sd.ms')

        simanalyze(project=sd_project, image=True, vis=visname_sd, modelimage='', imsize = [],
                   imdirection ='',cell = '',interactive = False, niter = 0,threshold = '0.01mJy',
                   weighting = 'natural',mask = [],outertaper = [],pbcor = False,stokes = 'I', 
                   featherimage = '',analyze=False, graphics='none', verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        # Observed value
        a = (os.path.isdir(sd_project+'/'+sd_project+'.sd.image') and 
             os.path.isdir(sd_project+'/'+sd_project+'.sd.image0') and 
             os.path.isdir(sd_project+'/'+sd_project+'.sd.image0.scaled') and 
             os.path.isdir(sd_project+'/'+sd_project+'.sd.image0.weight'))
        b = True # Expected value
        self.assertEqual(a,b) 

    @unittest.skipUnless(os.path.isdir(sd_project+'/'+sd_project+'.sd.image'), "Analysis-only mode presumes the existence of some image data, such as that generated by test_imaging_True_single_dish*")
    def test_imaging_False_single_dish_analysis_True_showfidelity_True(self):
        '''test_imaging_False_single_dish_analysis_True_showfidelity_True:
        ------------------------------------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        The showfidelity parameter displays the fidelity image
        '''

        simanalyze(project=sd_project, image=False, imagename='', analyze=True, showfidelity=True,
                   graphics='none', verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        # confirm that fidelity image was generated
        self.assertTrue(os.path.isdir(sd_project +'/'+ sd_project + '.' +'sd.fidelity')) 

####    both Interferometric and Single Dish    ####
class simanalyze_main_usage_modes_test_both(unittest.TestCase):

    @classmethod
    def setUpClass(cls): # TODO
        '''Called before tests in this class of test cases are run. Executes simobserve to create expected directory structure. Minimally, f"{project}/{project}.{suffix}" where suffix in ['skymodel','newmodel','compskymodel'].'''
        if not is_CASA6: # needs similar branch condition for tests to run in casalith
            default(simobserve)
            default(tclean)

        ## create images to serve as skymodel inputs to reference simulations?
        # gaussian = np.random.multivariate_normal(mean=[0,0], cov=[[1, 0], [0, 100]], size=[128,128])
        # gaus_img = gaussian[:, :, np.newaxis, :] # add degenerate axis to represent stokes
        # ia.newimagefromarray('demo.im', gaus_imag, csys, linear=False, overwrite=True, log=True

        # # create reference simulation with both int and sd
        # simobserve(project=both_project, 
        #            skymodel=imagepath_sd, 
        #            complist='', 
        #            setpointings=True, 
        #            direction=[], 
        #            mapsize=['1deg'], 
        #            maptype='square', 
        #            pointingspacing='',
        #            caldirection='',
        #            calflux='1Jy',
        #            obsmode='sd', 
        #            refdate='2020/02/13', 
        #            hourangle='transit', 
        #            totaltime='100s', 
        #            antennalist=configpath_int,
        #            sdantlist=configpath_sd, 
        #            sdant=0,
        #            outframe='LSRK', 
        #            thermalnoise='', 
        #            leakage=0.0, 
        #            graphics='none',
        #            verbose=False, 
        #            overwrite=False)
        # simobserve(project=both_project, 
        #            skymodel=imagepath_sd, 
        #            complist='', 
        #            setpointings=True, 
        #            direction=[], 
        #            mapsize=['1deg'], 
        #            maptype='square', 
        #            pointingspacing='',
        #            caldirection='',
        #            calflux='1Jy',
        #            obsmode='int', 
        #            refdate='2020/02/13', 
        #            hourangle='transit', 
        #            totaltime='100s', 
        #            antennalist=configpath_int,
        #            sdantlist=configpath_sd, 
        #            sdant=0,
        #            outframe='LSRK', 
        #            thermalnoise='', 
        #            leakage=0.0, 
        #            graphics='none',
        #            verbose=False, 
        #            overwrite=False)

    def setUp(self):
        if not is_CASA6:
            default(simobserve)
            default(tclean)

    @unittest.skip("Still need data to generate reference simulation for single dish + interferometry case")
    def test_imaging_True_interferometric_and_single_dish_analysis_False(self):
        '''test_imaging_True_interferometric_and_single_dish_analysis_False:
        --------------------------------------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        The vis parameter can accept one or more MSs that can be interferometric or single dish data        
        Modelimage will not be used if the MS is in total power
        '''

        simanalyze(project=both_project, image=True, vis= 'default',modelimage='', imsize = [1024],
                   imdirection ='',cell = '',interactive = False, niter = 0,threshold = '0.01mJy',weighting = 'natural',
                   mask = [],outertaper = [],pbcor = False,stokes = 'I', featherimage = '',analyze=False, 
                   graphics='none', verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        visname_both_int = str(both_project +'/'+ both_project + '.' + configpath_int.split('/')[-1][:-3] +'ms')
        visname_both_sd = str(both_project +'/'+ both_project + '.' + configpath_sd.split('/')[-1][:-3] +'sd.ms')

        # confirm that both the IF part and SD part generated output
        a = (os.path.isdir(visname_int[:-2]+'image') and
             os.path.isdir(both_project+'/'+both_project+'.sd.image') and 
             os.path.isdir(both_project+'/'+both_project+'.image0') and 
             os.path.isdir(both_project+'/'+both_project+'.sd.image0.scaled') and 
             os.path.isdir(both_project+'/'+both_project+'.sd.image0.weight'))
        b = True # Expected value
        self.assertEqual(a,b) 

####    Suite: Required for CASA5     ####

def suite():
    return[simanalyze_main_usage_modes_test_int, simanalyze_main_usage_modes_test_sd]
    #return[simanalyze_main_usage_modes_test_int, simanalyze_main_usage_modes_test_sd, simanalyze_main_usage_modes_test_both]
  
####    Main     ####
if __name__ == '__main__':
    unittest.main()
