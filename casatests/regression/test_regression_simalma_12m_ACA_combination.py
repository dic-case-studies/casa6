############################################################################
# $Id:$                                                                     #
# Test Name: test_regression_simalma_12m_ACA_combination                    #
#    Regression Test Script for                                             #
#       https://casaguides.nrao.edu/index.php?title=Simalma_(CASA_5.4)      #
#                                                                           #
# Rationale for Inclusion:                                                  #
#    Simple check that simalma works in CASA                                #
#                                                                           #
# Data used in test:                                                        #
#    The script can retrieve the data using curl. If it fails it            #
#    will get the data from the casatestdata repository                     #
#                                                                           #
#  The regression will work in 9 steps as follows                           #
#    Step 1: simulating 12m ALMA array (run simobserve)                     #
#    Step 2: simulating 7m ACA array (run simobserve)                       #
#    Step 3: simulating Total Power (run simobserve)                        #
#    Step 4: generating a total power image (run sdimaging)                 #
#    Step 5: imaging and analyzing alma.cycle6.3.cfg (run simanalyze/tclean)#
#    Step 6: imaging and analyzing aca.cycle6.cfg (run simanalyze/tclean)   #
#    Step 7: concatenating interferometric visibilities (run concat)        #
#    Step 8: imaging and analyzing m51/m51.concat.ms(run simanalyze/tclean) #
#    Step 9: combining a total power and synthesis image (run imregrid)     #
# Input data:                                                               #
#     M51ha.fits                                                            #
#                                                                           #
#############################################################################

# 
import os
import shutil
import unittest

CASA6 = False
try:
    from casatools import ctsys
    from casatasks import simalma
    CASA6 = True
    datapath = ctsys.resolve('regression/simalma_12m_ACA_combination/')
except ImportError:
    from __main__ import *
    from tasks import *
    from taskinit import *
    datapath = os.environ['CASAPATH'].split()[0] + '/casatestdata/regression/simalma_12m_ACA_combination/'

#os.system('curl https://casaguides.nrao.edu/images/3/3f/M51ha.fits.txt -f -o M51ha.fits')

class regression_simalma_12m_ACA_test(unittest.TestCase):

    def setUp(self):
        if not CASA6:
            default(simalma)
        
        self.skymodel = 'M51ha.fits'
        shutil.copyfile(os.path.join(datapath,self.skymodel), self.skymodel)
    
    def tearDown(self):
        os.system('rm -rf m51')
        os.system('rm -rf '+self.skymodel)


    def test_regression(self):
        '''Test regression 'simalma' with Main 12m Array and the ACA: manual combination of the data'''
        simalma(
            project="m51", 
            overwrite=True, 
            skymodel=self.skymodel, 
            indirection="J2000 23h59m59.96s -34d59m59.50s",
            incell="0.1arcsec",
            inbright="0.004",
            incenter="330.076GHz",
            inwidth="50MHz",
            antennalist=["alma.cycle6.3.cfg","aca.cycle6.cfg"],
            totaltime="1800s",
            tpnant = 2,
            tptime="7200s",
            pwv=0.6,
            imsize=[128,128],
            mapsize="1arcmin",
            dryrun = False )

        images = [
            './m51/m51.aca.cycle6.noisy.analysis.png','./m51/m51.aca.cycle6.noisy.image.png',
            './m51/m51.aca.cycle6.observe.png','./m51/m51.aca.cycle6.skymodel.png',
            './m51/m51.aca.tp.observe.png','./m51/m51.aca.tp.skymodel.png',
            './m51/m51.alma.cycle6.3.noisy.analysis.png','./m51/m51.alma.cycle6.3.noisy.image.png',
            './m51/m51.alma.cycle6.3.observe.png','./m51/m51.alma.cycle6.3.skymodel.png',
            './m51/m51.analysis.png','./m51/m51.combine.png',
            './m51/m51.concat.analysis.png','./m51/m51.concat.image.png'
                ]

        for image in images: self.assertTrue(os.path.isfile(image))

        M51MSs = [
            './m51/m51.aca.cycle6.ms','./m51/m51.alma.cycle6.3.ms',
            './m51/m51.aca.cycle6.noisy.ms','./m51/m51.concat.ms',
            './m51/m51.alma.cycle6.3.noisy.ms'
                ]
        for M51MS in M51MSs: self.assertTrue(os.path.isdir(M51MS))



def suite():
    return[regression_simalma_12m_ACA_test]


if __name__ == '__main__':
    unittest.main()

