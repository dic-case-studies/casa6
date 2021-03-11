##########################################################################
# test_req_task_phaseshift.py
#
# Copyright (C) 2018
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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_phaseshift/about
#
#
##########################################################################
import sys
import os
import unittest
import shutil
import numpy

CASA6 = False
try:
    import casatools
    from casatasks import phaseshift, tclean
    CASA6 = True

    tb = casatools.table()
    ia = casatools.image()
    ms = casatools.ms()

except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
    ia = iatool()


# Data ###
if CASA6:
    d = os.path.join('unittest', 'phaseshift')
    datapath = casatools.ctsys.resolve(
        os.path.join(d, 'refim_twopoints_twochan.ms')
    )
    datapath_Itziar = casatools.ctsys.resolve(
        os.path.join(d, 'Itziar.ms')
    )
    datapath_ngc = casatools.ctsys.resolve(
        os.path.join(d, 'ngc7538_ut.ms')
    )
    datapath_nep = casatools.ctsys.resolve(
        os.path.join(d, 'nep2-shrunk.ms')
    )
    datapath_mms = casatools.ctsys.resolve(
        os.path.join(d, 'uid___X02_X3d737_X1_01_small.mms')
    )
else:
    d = os.path.join(
        os.environ.get('CASAPATH').split()[0], 'casatestdata',
        'unittest', 'phaseshift'
    )
    datapath = os.path.join(d, 'refim_twopoints_twochan.ms')
    datapath_Itziar = os.path.join(d, 'Itziar.ms')
    datapath_ngc = os.path.join(d, 'ngc7538_ut.ms')
    datapath_nep = os.path.join(d, '/nep2-shrunk.ms')
    datapath_mms = os.path.join(d, '/uid___X02_X3d737_X1_01_small.mms')


def change_perms(path):
    os.chmod(path, 0o777)
    for root, dirs, files in os.walk(path):
        for d in dirs:
            os.chmod(os.path.join(root,d), 0o777)
        for f in files:
            os.chmod(os.path.join(root,f), 0o777)


datacopy = 'datacopy.ms'
datacopy_Itziar = 'Itziar_copy.ms'
datacopy_ngc = 'ngc_copy.ms'
datacopy_nep = 'nep_copy.ms'
datacopy_mms = 'mms_copy.mms'
output = 'phaseshiftout.ms'

class phaseshift_test(unittest.TestCase):

    def setUp(self):
        if not CASA6:
            default(phaseshift)
        shutil.copytree(datapath, datacopy)
        shutil.copytree(datapath_Itziar, datacopy_Itziar)
        shutil.copytree(datapath_ngc, datacopy_ngc)
        shutil.copytree(datapath_nep, datacopy_nep)
        shutil.copytree(datapath_mms, datacopy_mms)
        
        change_perms(datacopy)
        change_perms(datacopy_Itziar)
        change_perms(datacopy_ngc)
        change_perms(datacopy_nep)
        change_perms(datacopy_mms)
        
    def tearDown(self):
        shutil.rmtree(datacopy)
        shutil.rmtree(datacopy_Itziar)
        shutil.rmtree(datacopy_ngc)
        shutil.rmtree(datacopy_nep)
        shutil.rmtree(datacopy_mms)
        
        if os.path.exists('post_phaseshift.ms'):
            shutil.rmtree('post_phaseshift.ms')
        
        if os.path.exists(output):
            shutil.rmtree(output)
            
    def test_takesVis(self):
        ''' Check that the task requires a valid input MS '''
        result = phaseshift(datacopy, outputvis=output, phasecenter='J2000 19h53m50 40d06m00')
        # Completion without throwing an exception indicates success in CASA 6
        if not CASA6:
            self.assertTrue(result)
        
    def test_outvis(self):
        ''' Check that the outvis parameter specifies the name of the output '''
        phaseshift(datacopy, outputvis=output, phasecenter='J2000 19h53m50 40d06m00')
        
        self.assertTrue(os.path.exists(output))
        
    def test_fieldSelect(self):
        ''' Check the field selection parameter '''
        phaseshift(datacopy_Itziar, outputvis=output, phasecenter='J2000 00h00m01 -29d55m40', field='2')
        tb.open(output)
        data_selected = len(tb.getcol('FIELD_ID'))
        tb.close()
        
        self.assertTrue(data_selected == 6125)
        
    def test_spwSelect(self):
        ''' Check the spw selection parameter '''
        phaseshift(datacopy_ngc, outputvis=output, phasecenter='B1950_VLA 23h11m54 61d10m54', spw='1')
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()
        
        self.assertTrue(data_selected == 13338, msg=data_selected)
        
    def test_intentSelect(self):
        ''' Check the intent selection parameter '''
        phaseshift(datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35', intent='*FLUX*')
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()
        
        self.assertTrue(data_selected == 570)
        
    def test_arraySelect(self):
        ''' Check the array selection parameter '''
        msg = "specified array incorrectly found"
        if CASA6:
            with self.assertRaises(RuntimeError, msg=msg):
                phaseshift(
                    datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35',
                    array='1'
                )
        else:
            self.assertFalse(
                phaseshift(
                    datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35',
                    array='1'
                ), msg=msg
            )
        phaseshift(datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35', array='0')
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()
        
        self.assertTrue(data_selected == 6270, "Incorrect number of rows found")
        
    def test_observationSelect(self):
        ''' Check the observation selection parameter '''
        msg = "Observation not out of range"
        if CASA6:
            with self.assertRaises(RuntimeError, msg=msg):
                phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35', observation='1'
                )
        else:
            self.assertFalse(
                phaseshift(
                    datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35',
                    observation='1'
                ), msg=msg
            )
        phaseshift(datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35', observation='0')
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()
        
        self.assertTrue(data_selected == 6270, "Incorrect number of rows found")
        
    def test_keepsMMS(self):
        ''' Test the keepmms paramter creates the output as an MMS if the input is one as well '''
        phaseshift(datacopy_mms, outputvis=output, phasecenter='J2000 05h30m48 13d31m48', keepmms=False)
        ms.open(output)
        is_mms = ms.ismultims()
        ms.close()
        
        self.assertFalse(is_mms)
        
    def test_datacolumn(self):
        ''' Check that this parameter selects which datacolumns to write to the output MS '''
        msg = "Data column incorrectly present"
        if CASA6:
            with self.assertRaises(RuntimeError, msg=msg):
                phaseshift(
                    datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35',
                    datacolumn='MODEL'
                )
            # running to completion indicates success in CASA 6
            phaseshift(
                datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35',
                datacolumn='DATA'
            )
        else:
            self.assertFalse(
                phaseshift(
                    datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35',
                    datacolumn='MODEL'
                ), msg=msg
            )
            self.assertTrue(
                phaseshift(
                    datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -06d23m35',
                    datacolumn='DATA'
                ), msg="phaseshift unexpectedly failed"
            )        
        
    def test_phasecenter(self):
        ''' Check that this parameter sets the sky coordinates of the new phasecenter '''
        phaseshift(datacopy_nep, outputvis=output, phasecenter='ICRS 00h06m14 -08d23m35')
        tb.open(output)
        data_mean = numpy.mean(tb.getcol('DATA'))
        tb.close()

        self.assertTrue(numpy.isclose(data_mean, -0.00968202886279957-0.004072808512879953j))

    def test_shiftAndCompare(self):
        '''
            Check that changing the phasecenter with phaseshift and
            reverting with tclean results in the correct flux values at
            selected pixel locations
        '''
        # Run phaseshift to shift the MS phasecenter to a new location.
        post_vis = 'post_phaseshift.ms'
        os.system('rm -rf ' + post_vis)
        phaseshift(
            vis=datacopy, outputvis=post_vis,
            phasecenter='J2000 19h53m50 40d06m00'
        )

        # (1) Imaging on the original dataset
        os.system('rm -rf im2_pre*')
        tclean(vis=datacopy, imagename='im2_pre',
               imsize=2048,cell='5arcsec',niter=0,
               gridder='wproject', wprojplanes=128,pblimit=-0.1, 
               phasecenter='J2000 19h59m28.449 40d44m01.199')

        # (2) Image the phaseshifted dataset at it's new phasecenter as
        # image center.
        post_image = ''
        post_image = 'im2_post_phaseshift'
        os.system('rm -rf im2_post_phaseshift*')

        tclean(
            vis=post_vis, imagename='im2_post_phaseshift',
            imsize=2048, cell='5arcsec', niter=0,
            gridder='wproject', wprojplanes=128, pblimit=-0.1
        )

        # (3) Imaging on phaseshifted dataset, but with the imaging phasecenter
        # set. If this is working correctly, it should shift back to the same
        # source positions as the previous tclean result.
        post_image = 'im2_post_phaseshift_tclean_phasecenter'
        os.system('rm -rf im2_post_phaseshift_tclean_phasecenter*')

        tclean(
            vis=post_vis, imagename=post_image,
            imsize=2048, cell='5arcsec', niter=0, gridder='wproject',
            wprojplanes=128, pblimit=-0.1,
            phasecenter='J2000 19h59m28.449 40d44m01.199' 
        )

        # In the above 3 images, (1) has the correct locations.
        # Both (2) and (3) show the offset error when viewed in world
        # coordinates. Open in the viewer as an image stack, and step
        # through.
        # For comparisons with (1), pick the result from (3) because when
        # this works correctly the sources should appear at the same
        # pixel location as in (1).  This test is encoded below.

        ia.open('im2_pre.image')
        src1_pre = ia.pixelvalue([1024, 1024, 0, 0] )['value']
        src2_pre = ia.pixelvalue([1132, 1168, 0, 0] )['value']
        ia.close()
        ia.open(post_image+'.image')
        src1_post = ia.pixelvalue([1024, 1024, 0, 0] )['value']
        src2_post = ia.pixelvalue([1132, 1168, 0, 0] )['value']
        ia.close()

        print("Image value at source locations")
        print("Original MS : "+str(src1_pre) + " and " + str(src2_pre))
        print("Fixvis'd MS : "+str(src1_post) + " and " + str(src2_post))

        os.system('rm -rf im2_pre*')
        os.system('rm -rf im2_post_phaseshift*')
        os.system('rm -rf im2_post_phaseshift_tclean_phasecenter*')

        self.assertTrue(numpy.isclose(src1_pre['value'], src1_post['value'], rtol=0.01))
        self.assertTrue(numpy.isclose(src2_pre['value'], src2_post['value'], rtol=0.01))


def suite():
    return[phaseshift_test]


if __name__ == '__main__':
    unittest.main()
