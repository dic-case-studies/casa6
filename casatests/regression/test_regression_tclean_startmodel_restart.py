##########################################################################
#
##########################################################################
#
#  Datasets
#
##########################################################################

from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import unittest
import inspect

from casatestutils.imagerhelpers.imagetesthelpers import TestHelpers

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
     from casatools import ctsys
     from casatasks import casalog, tclean
     from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
     from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

     refdatapath = ctsys.resolve('regression/unittest/clean/refimager/')
else:
     from __main__ import default
     from tasks import *
     from taskinit import *
     from parallel.parallel_task_helper import ParallelTaskHelper
     from imagerhelpers.parallel_imager_helper import PyParallelImagerHelper

     refdatapath = os.environ.get('CASAPATH').split()[0] + '/data/regression/unittest/clean/refimager/'
     
## List to be run
def suite():
     return [test_csys_startmodel]
 
class test_csys_startmodel(unittest.TestCase):

     def setUp(self):
          self.msfile = ""
          self.parallel = False
          self.nnode = 0
          if ParallelTaskHelper.isMPIEnabled():
              self.parallel = True
              self.PH = PyParallelImagerHelper()
              self.nnode = len(self.PH.getNodeList())

          self.th = TestHelpers()

     def tearDown(self):
          """ don't delete it all """
          self.delData()

     # Separate functions here, for special-case tests that need their own MS.
     def prepData(self,msname=""):
          if msname != "":
               self.msfile=msname
          if (os.path.exists(self.msfile)):
               os.system('rm -rf ' + self.msfile)
          shutil.copytree(os.path.join(refdatapath,self.msfile), self.msfile)
          
     def delData(self,msname=""):
          if msname != "":
               self.msfile=msname
          os.system('rm -rf ' + self.msfile)
          os.system('rm -rf savemod.*')

     def checkfinal(self,pstr=""):
#          pstr += "["+inspect.stack()[1][3]+"] : To re-run this test :  runUnitTest.main(['test_regression_tclean_startmodel_restart["+ inspect.stack()[1][3] +"]'])"
          casalog.post(pstr,'INFO')
          print("posted %s" % pstr)
          if( pstr.count("(Fail") > 0 ):
               self.fail("\n"+pstr)

     ## Run this test only in parallel mode as it tests combinations of serial and parallel runs. 
     @unittest.skipIf(not ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_csys_startmodel_restart_mfs(self):
          """ [startmodel] test_csys_startmodel_restart_cube

          Run a sequence of tclean runs to trigger a complicated situation of restarts, mixing serial/parallel and model writes. 
          This sequence, coupled with the algorithm options listed below in tests #1 through #6 trigger three different errors that
          have been fixed in this branch, and one that will be addressed via CAS-9386 (cube refactor). 

          tclean call sequence : 
          --- (a) Parallel run for niter=0
          --- (b) Serial/Parallel run for niter=10, with calcres=F, calcpsf=F  : to reuse images from prev.
          --- (c) Serial model-predict run (without/with startmodel) : in one case it reuses prev image-set. in other case it reuses only 'model'. 
          --- (d) Impbcor on the output of (b)

          Note that this is not a full fix of the various instances of the 'latpole' inconsistency, but only a workaround. 
          Hence it needs a test to ensure this keeps working. 
          """
          self.prepData('refim_point.ms')

          specmode='mfs'
          runpar=self.parallel

          # Test 1 :  mfs + hogbom/mtmfs + usestartmod=False : Triggers both the latpole errors 
          #                                                                      (a) when the serial modelpredict runs on the output of the prev parallel run.
          #                                                                      (b) when the ouput of a parallel run is used by impbcor.
          # Test 2 : mfs + hogbom/mtmfs + usestartmod=True : Triggers the problem of 'model' being insufficient to create a valid ImStore. 
          # Run for deconvolver='hogbom' and 'mtmfs' because it exercises both SIImageStore and SIImageStoreMultiTerm 
          #      where some infrastructure code is re-implemented for multi-term vs single-term. 

          report = ''

          for runtype in ['serial','parallel']:
               for deconvolver in ['hogbom','mtmfs']:
                    for usestartmod in [False,True]:
                         infostr = 'Run with '+specmode +' - ' + deconvolver + ' - usestartmodel = ' + str(usestartmod) + ' - imaging in ' + runtype 
                         print(infostr)
                         report = report + infostr+'\n'

                         # Clear the model column.
                         self.th.delmodels(msname=self.msfile,modcol='reset0')

                         os.system('rm -rf savemod*')

                         ## (a) Always parallel run for niter=0
                         tclean( vis = self.msfile, imagename='savemod.par', spw='0:0~5,0:15~19',   niter=0,  imsize=100, cell='10.0arcsec', deconvolver=deconvolver,  specmode=specmode,  parallel=True)
                         
                         ## (b) Serial/Parallel run for niter=10, with calcres=F, calcpsf=F  : to reuse images from prev.
                         if runtype=='serial':
                              runpar=False
                         else:
                              runpar=True
                         tclean( vis = self.msfile, imagename='savemod.par', spw='0:0~5,0:15~19', niter=10, imsize=100, cell='10.0arcsec',deconvolver=deconvolver, specmode=specmode,  calcpsf=False, calcres=False,  parallel=runpar)
                         
                         # (c) Serial model-predict run (without/with startmodel) : 
                         if usestartmod==False:    # Do the restart by re-using the same image name as before.
                              imname2 = 'savemod.par'
                              startmodel=[]
                         else:
                              imname2 = 'savemod.ser'  # Do the restart with a new imagename, and using 'startmodel'
                              if deconvolver=='mtmfs':
                                   startmodel=['savemod.par.model.tt0', 'savemod.par.model.tt1']
                              else:
                                   startmodel= ['savemod.par.model']

                         tclean( vis=self.msfile, imagename=imname2, spw='0:0~5,0:15~19', niter=0,  imsize=100, cell='10.0arcsec',deconvolver=deconvolver, specmode=specmode, savemodel='modelcolumn', startmodel=startmodel, calcres=False,  calcpsf= False, restoration=False, parallel=False)
                         
                         # Check the values of the model column
                         report = report + self.th.check_chanvals(self.msfile,[(0,">",0.5),(19,"<",0.9)])
                         
                         ## (d) Impbcor on the output of (b) + check the output of pbcor
                         if deconvolver=='hogbom':
                              impbcor(imagename='savemod.par'+'.image', pbimage='savemod.par'+'.pb', outfile='savemod.par'+'.impbcor',overwrite=True)
                              report=report +self.th.checkall(imgexist=['savemod.par.impbcor'], imgval=[  ('savemod.par.impbcor',1.1,[50,50,0,0]) ])
                         else:
                              impbcor(imagename='savemod.par'+'.image.tt0', pbimage='savemod.par'+'.pb.tt0', outfile='savemod.par'+'.impbcor.tt0',overwrite=True)
                              report=report +self.th.checkall(imgexist=['savemod.par.impbcor.tt0'], imgval=[  ('savemod.par.impbcor.tt0',1.1,[50,50,0,0]) ])

          self.checkfinal(report)

     @unittest.skipIf(not ParallelTaskHelper.isMPIEnabled(), "Skip the test temporarily")
     def test_csys_startmodel_restart_cube(self):
          """ [startmodel] test_csys_startmodel_restart_cube : Check that csys differences w.r.to latpoles for parallel vs serial runs are appropriately squashed. 

          Run a sequence of tclean runs to trigger a complicated situation of restarts, mixing serial/parallel and model writes. 
          This sequence, coupled with the algorithm options listed below in tests #1 through #6 trigger three different errors that
          have been fixed in this branch, and one that will be addressed via CAS-9386 (cube refactor). 

          tclean call sequence : 
          --- (a) Parallel run for niter=0
          --- (b) Serial/Parallel run for niter=10, with calcres=F, calcpsf=F  : to reuse images from prev.
          --- (c) Serial model-predict run (without/with startmodel) : in one case it reuses prev image-set. in other case it reuses only 'model'. 
          --- (d) Impbcor on the output of (b)

          Note that this is not a full fix of the various instances of the 'latpole' inconsistency, but only a workaround. 
          Hence it needs a test to ensure this keeps working. 
          """
          self.prepData('refim_point.ms')

          specmode='cube'
          deconvolver='hogbom'
          runpar=self.parallel

          report=''

          # Test  : cube + hogbom + usestartmod=True/False :  Triggers the problem of refconcat image outputs being incompatible with restarts in serial later.
          ####for runtype in ['serial','parallel']:     ## This will fail in 'serial' because we cannot mix and match refconcat images with regular ones. 
          ####                                                                 This must be revisited after CAS-9386. 
          for runtype in ['parallel']:
               for usestartmod in [False,True]:
                    print('Run with '+specmode +' - ' + deconvolver + ' - ' + str(usestartmod))
                    
                    infostr = 'Run with '+specmode +' - ' + deconvolver + ' - usestartmodel = ' + str(usestartmod) + ' - imaging in ' + runtype 
                    print(infostr)
                    report = report + infostr+'\n'
                    
                    # Clear the model column.
                    self.th.delmodels(msname=self.msfile,modcol='reset0')
                    
                    os.system('rm -rf savemod*')
                    
                    ## (a) Always parallel run for niter=0
                    tclean( vis = self.msfile, imagename='savemod.par', spw='0:0~5,0:15~19',   niter=0,  imsize=100, cell='10.0arcsec', deconvolver=deconvolver,  specmode=specmode,  parallel=True)
                    
                    ## (b) Serial/Parallel run for niter=10, with calcres=F, calcpsf=F  : to reuse images from prev.
                    if runtype=='serial':
                         runpar=False
                    else:
                         runpar=True
                    tclean( vis = self.msfile, imagename='savemod.par', spw='0:0~5,0:15~19', niter=10, imsize=100, cell='10.0arcsec',deconvolver=deconvolver, specmode=specmode,  calcpsf=False, calcres=False,  parallel=runpar)
                         
                    # (c) Serial model-predict run (without/with startmodel) : 
                    if usestartmod==False:    # Do the restart by re-using the same image name as before.
                         imname2 = 'savemod.par'
                         startmodel=[]
                    else:
                         imname2 = 'savemod.ser'  # Do the restart with a new imagename, and using 'startmodel'
                         if deconvolver=='mtmfs':
                              startmodel=['savemod.par.model.tt0', 'savemod.par.model.tt1']
                         else:
                              startmodel= ['savemod.par.model']
                              
                    tclean( vis=self.msfile, imagename=imname2, spw='0:0~5,0:15~19', niter=0,  imsize=100, cell='10.0arcsec',deconvolver=deconvolver, specmode=specmode, savemodel='modelcolumn', startmodel=startmodel, calcres=False,  calcpsf= False, restoration=False, parallel=False)
                              
                    # Check the values of the model column
                    report = report + self.th.check_chanvals(self.msfile,[(0,">",0.5),(19,"<",0.9)])
                    
                    ## (d) Impbcor on the output of (b) + check the output of pbcor (first channel)
                    impbcor(imagename='savemod.par'+'.image', pbimage='savemod.par'+'.pb', outfile='savemod.par'+'.impbcor',overwrite=True)
                    report=report +self.th.checkall(imgexist=['savemod.par.impbcor'], imgval=[  ('savemod.par.impbcor',1.5,[50,50,0,0]) ])


          self.checkfinal(report)

if is_CASA6:
     if __name__ == '__main__':
          unittest.main()
