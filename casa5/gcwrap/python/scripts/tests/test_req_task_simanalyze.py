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
    from casatools import ctsys, image, coordsys, componentlist, quanta
    ia = image()
    cs = coordsys()
    cl = componentlist()
    qa = quanta()
    from casatasks import tclean, simobserve, simanalyze, casalog
else:
    from __main__ import default
    from tasks import *
    from taskinit import *
 
# DATA #
if is_CASA6:
    ctsys_resolve = ctsys.resolve
    configpath_int      = ctsys_resolve('alma/simmos/vla.a.cfg')
    imagepath_int       = ctsys_resolve('nrao/VLA/CalModels/3C286_Q.im/')
    configpath_both_int = ctsys_resolve('alma/simmos/aca.cycle7.cfg')
    configpath_sd       = ctsys_resolve('alma/simmos/aca.tp.cfg')
else:
    dataroot = os.environ.get('CASAPATH').split()[0]
    configpath_int      = os.path.join(dataroot, 'data/alma/simmos/vla.a.cfg')
    imagepath_int       = os.path.join(dataroot, 'data/nrao/VLA/CalModels/3C286_Q.im/')
    configpath_both_int = os.path.join(dataroot, 'data/alma/simmos/aca.cycle7.cfg')
    configpath_sd       = os.path.join(dataroot, 'data/alma/simmos/aca.tp.cfg')
    def ctsys_resolve(apath):
        dataPath = os.path.join(os.environ['CASAPATH'].split()[0],'casatestdata/')
        return os.path.join(dataPath,apath)
    
# Path to MS data
mspath_sd = ctsys_resolve('unittest/simanalyze/refim_twopoints_twochan.ms')

logpath = casalog.logfile()

int_project = 'sim_interferometric'
sd_project = 'sim_single_dish'
imagepath_sd = 'sd_model.image'
both_project = 'sim_both'
both_component_list = 'J0319+4130.cl'
both_model_image_int = both_component_list[:-3]+'.im'
both_model_image_sd = both_component_list[:-3]+'_rebin.im'

####    interferometric    ####
class simanalyze_main_usage_modes_test_int(unittest.TestCase):

    def setUp(self):
        if not is_CASA6: # needs similar branch condition for casalith
            default(simanalyze)
            default(simobserve)
            default(tclean)

        simobserve(project=int_project, skymodel=imagepath_int, complist='', setpointings=True, direction=[], 
                   mapsize='', maptype='square', pointingspacing='',caldirection='',calflux='1Jy',obsmode='int', 
                   refdate='2020/02/13', hourangle='transit', totaltime='100s', antennalist=configpath_int,
                   outframe='LSRK', thermalnoise='', leakage=0.0, graphics='none',verbose=False, overwrite=False)

    def tearDown(self):
        os.system("rm -rf {}".format(int_project))
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

    def test_imaging_False_interferometric_analysis_True_showfidelity_True(self):
        '''test_imaging_False_interferometric_analysis_True_showfidelity_True:
        ----------------------------------------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        The showfidelity parameter displays the fidelity image
        '''

        visname_int = str(int_project +'/'+ int_project + '.' + configpath_int.split('/')[-1][:-3] +'ms')

        # this is not the subject of this test
        simanalyze(project=int_project, image=True, vis=visname_int, modelimage='',imsize = [],
                   imdirection ='',cell = '',interactive = False, niter = 0,threshold = '0.1mJy',
                   weighting = 'natural',mask = [],outertaper = [],pbcor = False,stokes = 'I', 
                   featherimage = '',analyze=False, graphics='none', verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        # this is the subject of this test
        simanalyze(project=int_project, image=False, analyze=True, showfidelity=True,graphics='none', 
                   verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        # confirm that fidelity image was generated
        self.assertTrue(os.path.isdir(visname_int[:-2]+'fidelity')) 

####    Single Dish     ####
class simanalyze_main_usage_modes_test_sd(unittest.TestCase):

    def setUp(self):
        """Executes simobserve to create expected directory structure. Minimally, f"{project}/{project}.{suffix}" where suffix in ['skymodel','newmodel','compskymodel']."""
        if not is_CASA6: # needs similar branch condition for casalith
            default(simobserve)
            default(tclean)
            default(simanalyze)

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

    def tearDown(self):
        os.system("rm -rf {}*".format(imagepath_sd[:-5]))
        os.system("rm -rf {}".format(sd_project))
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

    def test_imaging_False_single_dish_analysis_True_showfidelity_True(self):
        '''test_imaging_False_single_dish_analysis_True_showfidelity_True:
        ------------------------------------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        The showfidelity parameter displays the fidelity image
        '''

        visname_sd = str(sd_project +'/'+ sd_project + '.' + configpath_sd.split('/')[-1][:-3] +'sd.ms')
        print(visname_sd)

        # this is not the subject of this test
        simanalyze(project=sd_project, image=True, vis=visname_sd, modelimage='', imsize = [],
                   imdirection ='',cell = '',interactive = False, niter = 0,threshold = '0.01mJy',
                   weighting = 'natural',mask = [],outertaper = [],pbcor = False,stokes = 'I', 
                   featherimage = '',analyze=False, graphics='none', verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        # this is the subject of this test
        simanalyze(project=sd_project, image=False, imagename='', analyze=True, showfidelity=True,
                   graphics='none', verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        # confirm that fidelity image was generated
        self.assertTrue(os.path.isdir(sd_project +'/'+ sd_project + '.' +'sd.fidelity')) 

####    both Interferometric and Single Dish    ####
class simanalyze_main_usage_modes_test_both(unittest.TestCase):

    def setUp(self):
        """Executes simobserve to create expected directory structure. Minimally, f"{project}/{project}.{suffix}" where suffix in ['skymodel','newmodel','compskymodel']."""
        if not is_CASA6:
            default(simobserve)
            default(tclean)
            default(simanalyze)

        ## create image to serve as skymodel inputs to reference simulations
        # build a point source component and convert to image
        cl.done()
        # J0319+4130 at band 3
        cl.addcomponent(dir="J2000 03h19m48.160s +41d30m42.11s",
                        flux=14.38, fluxunit='Jy', freq='115.271GHz',
                        shape="point", spectrumtype="spectral index", index=-1.0)
        cl.rename(both_component_list)
    
        ia.fromshape(both_model_image_int,[256,256,1,128],overwrite=True)
        cs=ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])
        cell_rad=qa.convert(qa.quantity("0.1arcsec"),"rad")['value']
        cs.setincrement([-cell_rad,cell_rad],'direction')
        cs.setreferencevalue([qa.convert("03h19m48.160s",'rad')['value'],
                              qa.convert("41d30m42.11s",'rad')['value']],type="direction")
        cs.setreferencevalue("115.271GHz",'spectral')
        cs.setincrement('15.1368MHz','spectral')
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(),subtract=False)
        ia.done()
        cs.done()

        # make another one with larger pixel size and FOV
        ia.fromshape(both_model_image_sd,[256,256,1,128],overwrite=True)
        cs=ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])
        cell_rad=qa.convert(qa.quantity("1.0arcsec"),"rad")['value']
        cs.setincrement([-cell_rad,cell_rad],'direction')
        cs.setreferencevalue([qa.convert("03h19m48.160s",'rad')['value'],
                              qa.convert("41d30m42.11s",'rad')['value']],type="direction")
        cs.setreferencevalue("115.271GHz",'spectral')
        cs.setincrement('15.1368MHz','spectral')
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(),subtract=False)
        cl.done()
        ia.done()
        cs.done()

        # create reference simulation with both int and sd
        simobserve(project=both_project, 
                   skymodel=both_model_image_int,
                   incell='0.1arcsec',
                   complist='', 
                   setpointings=True, 
                   direction=[], 
                   mapsize=['20arcsec'], 
                   maptype='square', 
                   pointingspacing='',
                   caldirection='',
                   calflux='1Jy',
                   obsmode='int', 
                   refdate='2020/02/13', 
                   hourangle='transit', 
                   totaltime='720s', 
                   antennalist=configpath_both_int,
                   sdantlist=configpath_sd, 
                   sdant=0,
                   outframe='LSRK', 
                   thermalnoise='', 
                   leakage=0.0, 
                   graphics='none',
                   verbose=False, 
                   overwrite=True)
        simobserve(project=both_project, 
                    skymodel=both_model_image_sd, 
                    complist='', 
                    setpointings=True, 
                    direction=[], 
                    mapsize=['200arcsec'], 
                    maptype='square', 
                    pointingspacing='',
                    caldirection='',
                    calflux='1Jy',
                    obsmode='sd', 
                    refdate='2020/02/13', 
                    hourangle='transit', 
                    totaltime='7200s', 
                    antennalist=configpath_sd,
                    sdantlist=configpath_sd, 
                    sdant=0,
                    outframe='LSRK', 
                    thermalnoise='', 
                    leakage=0.0, 
                    graphics='none',
                    verbose=False, 
                    overwrite=True)

    def tearDown(self):
        os.system("rm -rf {}".format(both_project))
        os.system("rm -rf *.last")
        os.system("rm -rf *.cl")
        os.system("rm -rf *.im")

    @unittest.skip("Task call exits with SEVERE error from feather task: Failed AlwaysAssert chan >=0 && chan < Int(nchan()) && stokes >= 0 && stokes < Int(nstokes())")
    def test_imaging_True_interferometric_and_single_dish_analysis_False(self):
        '''test_imaging_True_interferometric_and_single_dish_analysis_False:
        --------------------------------------------------------------------
        Not all of the output files will be generated depending on parameter selections
        The image parameter determines if an image is produced (False for no .image file)
        The analyze parameter turns on or off the creation of analytical images
        The vis parameter can accept one or more MSs that can be interferometric or single dish data        
        Modelimage will not be used if the MS is in total power
        '''

        simanalyze(project=both_project, image=True, vis='default', modelimage='', imsize=[1458],
                   imdirection='',cell='', interactive=False, niter=0, threshold='0.01mJy', weighting='natural',
                   mask=[], outertaper=[], pbcor=False, stokes='I', featherimage='',
                   analyze=False, graphics='none', verbose=False, overwrite=True, dryrun=False, logfile=logpath)

        visname_both_int = str(both_project +'/'+ both_project + '.' + configpath_both_int.split('/')[-1][:-3] +'ms')
        visname_both_sd = str(both_project +'/'+ both_project + '.' + configpath_sd.split('/')[-1][:-3] +'sd.ms')

        # confirm that both the IF and SD imaging generated output
        a = (os.path.isdir(visname_both_int[:-2]+'image') and
             os.path.isdir(both_project+'/'+both_project+'.sd.image') and 
             os.path.isdir(both_project+'/'+both_project+'.image0') and 
             os.path.isdir(both_project+'/'+both_project+'.sd.image0.scaled') and 
             os.path.isdir(both_project+'/'+both_project+'.sd.image0.weight'))
        # perhaps we should check for feather as well...?
        b = True # Expected value
        self.assertEqual(a,b) 

####    Suite: Required for CASA5     ####

def suite():
    return[simanalyze_main_usage_modes_test_int, simanalyze_main_usage_modes_test_sd, simanalyze_main_usage_modes_test_both]
  
####    Main     ####
if __name__ == '__main__':
    unittest.main()
