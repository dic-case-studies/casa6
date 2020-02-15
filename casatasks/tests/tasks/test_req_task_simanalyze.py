##########################################################################
# test_req_task_simanalyze.py
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
 
CASA6 = False
 
try:
    import casatools # not a good idea inside the casashell...
    import casatasks # perhaps os.path.exists(casatools.__file__) instead
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
 
# DATA #
if CASA6:
    dataroot = casatools.ctsys.resolve()
    configpath = casatools.ctsys.resolve((os.path.join
                                          (dataroot, 
                                           'alma/simmos/vla.a.cfg')))
    imagepath_int = casatools.ctsys.resolve((os.path.join
                                             (dataroot,
                                              'nrao/VLA/CalModels/3C286_Q.im/')))
    mspath_sd = casatools.ctsys.resolve((os.path.join
                                         (dataroot,
                                          'regression/unittest/clean/refimager/'\
                                          'refim_twopoints_twochan.ms')))

else:
    dataroot = os.environ.get('CASAPATH').split()[0]
    configpath = (dataroot + 
                  os.path.join(dataroot, 'alma/simmos/vla.a.cfg'))
    imagepath_int = (dataroot + 
                     os.path.join(dataroot,'nrao/VLA/CalModels/3C286_Q.im/'))
    mspath_sd = casatools.ctsys.resolve((os.path.join
                                         (dataroot,
                                          'regression/unittest/clean/refimager/'\
                                          'refim_twopoints_twochan.ms')))

logpath = casalog.logfile()

####    Tests     ####
class simanalyze_main_usage_modes_test(unittest.TestCase):

    ### Set Up
    @classmethod
    def setUpClass(cls):
        '''Called before tests in this class of test cases are run. Executes simobserve to create expected directory structure. Minimally, f"{project}/{project}.{suffix}" where suffix in ['skymodel','newmodel','compskymodel'].'''

        # create first reference simulation
        int_project = 'sim_interferometric'
        simobserve(project=int_project, 
                   skymodel=imagepath_int, 
                   complist='', 
                   setpointings=True, 
                   direction=[], 
                   mapsize='', 
                   maptype='square', 
                   pointingspacing='',
                   caldirection='',
                   calflux='1Jy',
                   obsmode='int', 
                   refdate='2020/02/13', 
                   hourangle='transit', 
                   totaltime='100s', 
                   antennalist=configpath,
                   outframe='LSRK', 
                   thermalnoise='', 
                   leakage=0.0, 
                   graphics='none',
                   verbose=False, 
                   overwrite=False)
        
        # create reference image > 2.5*PB to use for SD sim
        tclean(vis=mspath_sd, 
               imagename='sim_sd_model',
               imsize=512,
               cell='20arcsec',
               niter=0,
               gridder='wproject', 
               wprojplanes=5,pblimit=-0.1, 
               phasecenter='J2000 19h59m28.449 40d44m01.199')

        imagepath_sd = 'sim_sd_model.image'

        # create second reference simulation
        sd_project = 'sim_single_dish'
        simobserve(project=sd_project, 
                   skymodel=imagepath_sd, 
                   complist='', 
                   setpointings=True, 
                   direction=[], 
                   mapsize='', 
                   maptype='square', 
                   pointingspacing='',
                   caldirection='',
                   calflux='1Jy',
                   obsmode='sd', 
                   refdate='2020/02/13', 
                   hourangle='transit', 
                   totaltime='100s', 
                   sdantlist=configpath, 
                   sdant=0,
                   outframe='LSRK', 
                   thermalnoise='', 
                   leakage=0.0, 
                   graphics='none',
                   verbose=False, 
                   overwrite=False)
  
    def setUp(self):
        '''Method called to prepare the test fixture.  This is called immediately before calling the test method'''
        if not CASA6:
            default(simanalyze)
        pass
 
    ### Teardown
    def tearDown(self):
        '''A method called after an individual test in a class has run'''
        try:
            os.remove('simanalyze.last')
        except FileNotFoundError:
            pass
 
    @classmethod
    def tearDownClass(cls):
        '''Called after tests in this class of test cases are run. Undoes the execution of simobserve outputs determined by the setUpClass method run.'''

        # Remove the interferometric reference simulation
        shutil.rmtree(int_project)

        # Remove the products of the tclean call required 
        to create the reference image for SD simulation
        shutil.rmtree(imagepath_sd[:-5]+'image')
        shutil.rmtree(imagepath_sd[:-5]+'model')
        shutil.rmtree(imagepath_sd[:-5]+'pb')
        shutil.rmtree(imagepath_sd[:-5]+'psf')
        shutil.rmtree(imagepath_sd[:-5]+'residual')
        shutil.rmtree(imagepath_sd[:-5]+'sumwt')

        # Remove the single dish reference simulation
        shutil.rmtree(sd_project)

        # Remove (CASA<6 or casashell) parameter storage files
        try:
            os.remove('simobserve.last')
            os.remove('tclean.last')
        except FileNotFoundError:
            pass
 
    ### Test Cases
    def test_imaging_False_analysis_False(self):
        '''
        test_simanalyze: Not all of the output files will be generated depending on parameter selections
        test_simanalyze: The image parameter determines if an image is produced (False for no .image file)
        test_simanalyze: The analyze parameter turns on or off the creation of analytical images
        test_simanalyze: The imagename parameter takes the name of an already synthesized image
        '''
        simanalyze(project=int_project, 
                   image=False, 
                   imagename='', 
                   skymodel='', 
                   analyze=False, 
                   graphics='none', 
                   verbose=False, 
                   overwrite=False, 
                   dryrun=False, 
                   logfile=logpath)
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_interferometric_only_analysis_False(self):
        '''
        test_simanalyze: Not all of the output files will be generated depending on parameter selections
        test_simanalyze: The image parameter determines if an image is produced (False for no .image file)
        test_simanalyze: The analyze parameter turns on or off the creation of analytical images
        '''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_total_power_only_analysis_False(self):
        '''test_simanalyze: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_interferometric_and_total_power_analysis_False(self):
        '''test_simanalyze: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_False_analysis_True_showfidelity_True(self):
        '''
        test_simanalyze: Not all of the output files will be generated depending on parameter selections
        test_simanalyze: The image parameter determines if an image is produced (False for no .image file)
        test_simanalyze: The analyze parameter turns on or off the creation of analytical images
        test_simanalyze:  The showfidelity parameter displays the fidelity image
        '''
        simanalyze(project=int_project, 
                   image=False, 
                   imagename=imagepath_int, 
                   skymodel=int_project+'/'+int_project+".skymodel", 
                   analyze=True, 
                   showuv=False,
                   showpsf=False,
                   showmodel=False,
                   showconvolved=False,
                   showclean=False,
                   showresidual=False,
                   showdifference=False,
                   showfidelity=True,
                   graphics='none', 
                   verbose=False, 
                   overwrite=False, 
                   dryrun=False, 
                   logfile=logpath)

        # A fidelity image should have been produced
        a = os.path.exists() # Observed value
        b = True # Expected value
        self.assertEqual(a,b) 
 
    def test_imaging_True_interferometric_only_analysis_True(self):
        '''test_simanalyze: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_total_power_only_analysis_True(self):
        '''test_simanalyze: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 

    def test_imaging_True_interferometric_and_total_power_analysis_True(self):
        '''test_simanalyze: Requirement'''
        # Some Single Test
        a = 1 # Observed value
        b = 1 # Expected value
        self.assertEqual(a,b) 


####    Suite: Required for CASA5     ####
def suite():
    return[simanalyze_main_usage_modes_test]
  
####    Imports     ####
if __name__ == '__main__':
    unittest.main()
