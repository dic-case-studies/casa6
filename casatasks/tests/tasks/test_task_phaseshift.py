##########################################################################
# test_task_phaseshift.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.phaseshift.html
#
#
##########################################################################
import glob
import numpy as np
import os
import shutil
import unittest

from casatools import (
        componentlist, ctsys, image, measures, ms, msmetadata, quanta,
        regionmanager, simulator, table
    )
from casatasks import flagdata, imstat, phaseshift, tclean
from casatasks.private import simutil
cl = componentlist()
ia = image()
md = msmetadata()
me = measures()
ms = ms()
qa = quanta()
rg = regionmanager()
sm = simulator()
tb = table()
datadir = os.path.join('unittest', 'phaseshift')
ctsys_resolve = ctsys.resolve

datadir = os.path.join('unittest', 'phaseshift')
datapath = ctsys_resolve(os.path.join(datadir, 'refim_twopoints_twochan.ms'))
datapath_Itziar = ctsys_resolve(os.path.join(datadir, 'Itziar.ms'))
datapath_ngc = ctsys_resolve(os.path.join(datadir, 'ngc7538_ut.ms'))
datapath_nep = ctsys_resolve(os.path.join(datadir, 'nep2-shrunk.ms'))
datapath_mms = ctsys_resolve(
    os.path.join(datadir, 'uid___X02_X3d737_X1_01_small.mms')
)


def change_perms(path):
    os.chmod(path, 0o777)
    for root, dirs, files in os.walk(path):
        for d in dirs:
            os.chmod(os.path.join(root, d), 0o777)
        for f in files:
            os.chmod(os.path.join(root, f), 0o777)


datacopy = 'datacopy.ms'
datacopy_Itziar = 'Itziar_copy.ms'
datacopy_ngc = 'ngc_copy.ms'
datacopy_nep = 'nep_copy.ms'
datacopy_mms = 'mms_copy.mms'
output = 'phaseshiftout.ms'


class phaseshift_test(unittest.TestCase):

    def setUp(self):
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
        result = phaseshift(
            datacopy, outputvis=output,
            phasecenter='J2000 19h53m50 40d06m00'
        )

    def test_outvis(self):
        '''
        Check that the outvis parameter specifies the name of the output
        '''
        phaseshift(
            datacopy, outputvis=output,
            phasecenter='J2000 19h53m50 40d06m00'
        )

        self.assertTrue(os.path.exists(output))

    def test_fieldSelect(self):
        ''' Check the field selection parameter '''
        phaseshift(
            datacopy_Itziar, outputvis=output,
            phasecenter='J2000 00h00m01 -29d55m40', field='2'
        )
        tb.open(output)
        data_selected = len(tb.getcol('FIELD_ID'))
        tb.close()

        self.assertTrue(data_selected == 6125)

    def test_spwSelect(self):
        ''' Check the spw selection parameter '''
        phaseshift(
            datacopy_ngc, outputvis=output,
            phasecenter='B1950_VLA 23h11m54 61d10m54', spw='1'
        )
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()

        self.assertTrue(data_selected == 13338, msg=data_selected)

    def test_intentSelect(self):
        ''' Check the intent selection parameter '''
        phaseshift(
            datacopy_nep, outputvis=output,
            phasecenter='ICRS 00h06m14 -06d23m35', intent='*FLUX*'
        )
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()

        self.assertTrue(data_selected == 570)

    def test_arraySelect(self):
        ''' Check the array selection parameter '''
        msg = "specified array incorrectly found"
        with self.assertRaises(RuntimeError, msg=msg):
            phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35',
                    array='1'
            )
        phaseshift(
            datacopy_nep, outputvis=output,
            phasecenter='ICRS 00h06m14 -06d23m35', array='0'
        )
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()

        self.assertTrue(
            data_selected == 6270,
            "Incorrect number of rows found"
        )

    def test_observationSelect(self):
        ''' Check the observation selection parameter '''
        msg = "Observation not out of range"
        with self.assertRaises(RuntimeError, msg=msg):
            phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35', observation='1'
            )

        phaseshift(
            datacopy_nep, outputvis=output,
            phasecenter='ICRS 00h06m14 -06d23m35', observation='0'
        )
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()

        self.assertTrue(
            data_selected == 6270, "Incorrect number of rows found"
        )

    def test_keepsMMS(self):
        '''
        Test the keepmms paramter creates the output as an MMS
        if the input is one as well
        '''
        phaseshift(
            datacopy_mms, outputvis=output,
            phasecenter='J2000 05h30m48 13d31m48', keepmms=False
        )
        ms.open(output)
        is_mms = ms.ismultims()
        ms.close()

        self.assertFalse(is_mms)

    def test_datacolumn(self):
        '''
        Check that this parameter selects which datacolumns to write
        to the output MS
        '''
        msg = "Data column incorrectly present"
        with self.assertRaises(RuntimeError, msg=msg):
            phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35',
                    datacolumn='MODEL'
            )
            # running to completion indicates success in CASA 6
            phaseshift(
                datacopy_nep, outputvis=output,
                phasecenter='ICRS 00h06m14 -06d23m35', datacolumn='DATA'
            )

    def test_phasecenter(self):
        '''
        Check that this parameter sets the sky coordinates of the new
        phasecenter
        '''
        msg = 'Empty phasecenter param incorrectly runs'
        with self.assertRaises(ValueError, msg=msg):
            phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter=''
            )

        phaseshift(
            datacopy_nep, outputvis=output,
            phasecenter='ICRS 00h06m14 -08d23m35'
        )
        tb.open(output)
        data_mean = np.mean(tb.getcol('DATA'))
        tb.close()

        self.assertTrue(np.isclose(
            data_mean, -0.00968202886279957-0.004072808512879953j)
        )

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
        tclean(
            vis=datacopy, imagename='im2_pre',
            imsize=2048, cell='5arcsec', niter=0,
            gridder='wproject', wprojplanes=128, pblimit=-0.1,
            phasecenter='J2000 19h59m28.449 40d44m01.199'
        )

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
        src1_pre = ia.pixelvalue([1024, 1024, 0, 0])['value']
        src2_pre = ia.pixelvalue([1132, 1168, 0, 0])['value']
        ia.close()
        ia.open(post_image+'.image')
        src1_post = ia.pixelvalue([1024, 1024, 0, 0])['value']
        src2_post = ia.pixelvalue([1132, 1168, 0, 0])['value']
        ia.close()

        print("Image value at source locations")
        print("Original MS : "+str(src1_pre) + " and " + str(src2_pre))
        print("Phase shifted MS : "+str(src1_post) + " and " + str(src2_post))

        os.system('rm -rf im2_pre*')
        os.system('rm -rf im2_post_phaseshift*')
        os.system('rm -rf im2_post_phaseshift_tclean_phasecenter*')

        self.assertTrue(
            np.isclose(src1_pre['value'], src1_post['value'], rtol=0.01)
        )
        self.assertTrue(
            np.isclose(src2_pre['value'], src2_post['value'], rtol=0.01)
        )


class reference_frame_tests(unittest.TestCase):
    # much of this code is adapted from that of rurvashi
    # https://gitlab.nrao.edu/rurvashi/simulation-in-casa-6/-/blob/master/Simulation_Script_Demo.ipynb

    comp_list = 'sim_onepoint.cl'
    orig_ms = 'sim_data.ms'
    orig_im = 'im_orig'
    pshift_ms = 'sim_data_pshift.ms'
    pshift_im = 'im_post_phaseshift'
    pshift_shiftback_im = 'im_post_phaseshift_tclean_shiftback'
    exp_flux = 5

    @classmethod
    def __delete_intermediate_products(cls):
        if os.path.exists(cls.pshift_ms):
            shutil.rmtree(cls.pshift_ms)
        for im in (cls.pshift_im, cls.pshift_shiftback_im):
            for path in glob.glob(im + '*'):
                shutil.rmtree(path)

    @classmethod
    def __delete_data(cls):
        cls.__delete_intermediate_products()
        for x in (cls.orig_ms, cls.comp_list):
            if os.path.exists(x):
                shutil.rmtree(x)
        for path in glob.glob(cls.orig_im + '*'):
            shutil.rmtree(path)

    def setUp(self):
        self.__delete_data()

    def tearDown(self):
        self.__delete_data()

    @classmethod
    def tearDownClass(cls):
        cls.__delete_data()

    @classmethod
    def __phase_center_string(cls, ra, dec, frame):
        return ' '.join([frame, ra, dec])

    @classmethod
    def __makeMSFrame(cls, radir, decdir, dirframe):
        """
        Construct an empty Measurement Set that has the desired
        observation setup.
        """
        # Open the simulator
        sm.open(ms=cls.orig_ms)

        # Read/create an antenna configuration.
        # Canned antenna config text files are located at
        # /home/casa/data/trunk/alma/simmos/*cfg
        antennalist = os.path.join(
            ctsys.resolve("alma/simmos"), "vla.d.cfg"
        )


        # Fictitious telescopes can be simulated by specifying x, y, z, d,
        # an, telname, antpos.
        # x,y,z are locations in meters in ITRF (Earth centered)
        # coordinates.
        # d, an are lists of antenna diameter and name.
        # telname and obspos are the name and coordinates of the
        # observatory.
        (x, y, z, d, an, an2, telname, obspos) = (
            simutil.simutil().readantenna(antennalist)
        )

        # Set the antenna configuration
        sm.setconfig(
            telescopename=telname, x=x, y=y, z=z, dishdiameter=d,
            mount=['alt-az'], antname=an, coordsystem='global',
            referencelocation=me.observatory(telname)
        )

        # Set the polarization mode (this goes to the FEED subtable)
        sm.setfeed(mode='perfect R L', pol=[''])

        # Set the spectral window and polarization (one
        # data-description-id).
        # Call multiple times with different names for multiple SPWs or
        # pol setups.
        sm.setspwindow(
            spwname="LBand", freq='1.0GHz', deltafreq='0.1GHz',
            freqresolution='0.2GHz', nchannels=1, stokes='RR LL'
        )

        # Setup source/field information (i.e. where the observation phase
        # center is) Call multiple times for different pointings or source
        # locations.
        sm.setfield(
            sourcename="fake", sourcedirection=me.direction(
                rf=dirframe, v0=radir, v1=decdir
            )
        )

        # Set shadow/elevation limits (if you care). These set flags.
        sm.setlimits(shadowlimit=0.01, elevationlimit='1deg')

        # Leave autocorrelations out of the MS.
        sm.setauto(autocorrwt=0.0)

        # Set the integration time, and the convention to use for timerange
        # specification
        # Note : It is convenient to pick the hourangle mode as all times
        #   specified in sm.observe() will be relative to when the source
        #   transits.
        sm.settimes(
            integrationtime='2000s', usehourangle=True,
            referencetime=me.epoch('UTC', '2019/10/4/00:00:00')
        )

        # Construct MS metadata and UVW values for one scan and ddid
        # Call multiple times for multiple scans.
        # Call this with different sourcenames (fields) and spw/pol
        # settings as defined above.
        # Timesteps will be defined in intervals of 'integrationtime',
        # between starttime and stoptime.
        sm.observe(
            sourcename="fake", spwname='LBand', starttime='-5.0h',
            stoptime='+5.0h'
        )
        # Close the simulator
        sm.close()
        # Unflag everything (unless you care about elevation/shadow flags)
        flagdata(vis=cls.orig_ms, mode='unflag')

    @classmethod
    def __makeCompList(cls, ra, dec, frame):
        # Add sources, one at a time.
        # Call multiple times to add multiple sources.
        # ( Change the 'dir', obviously )
        cl.addcomponent(
            dir=cls.__phase_center_string(ra, dec, frame),
            flux=cls.exp_flux,      # For a gaussian, this is the
                                    # integrated area.
            fluxunit='Jy', freq='1.5GHz', shape='point',
            spectrumtype="constant"
        )
        # Save the file
        cl.rename(filename=cls.comp_list)
        cl.done()

    @classmethod
    def __sim2fields(cls, radir, decdir, dirframe, offset):
        """
        Construct an empty Measurement Set with two fieldsthat has
        the desired observation setup.
        """
        # Open the simulator
        sm.open(ms=cls.orig_ms)

        # Read/create an antenna configuration.
        # Canned antenna config text files are located at
        # /home/casa/data/trunk/alma/simmos/*cfg
        antennalist = os.path.join(
            ctsys.resolve("alma/simmos"), "vla.d.cfg"
        )


        # Fictitious telescopes can be simulated by specifying x, y, z, d,
        # an, telname, antpos.
        # x,y,z are locations in meters in ITRF (Earth centered)
        # coordinates.
        # d, an are lists of antenna diameter and name.
        # telname and obspos are the name and coordinates of the
        # observatory.
        (x, y, z, d, an, an2, telname, obspos) = (
            simutil.simutil().readantenna(antennalist)
        )

        # Set the antenna configuration
        sm.setconfig(
            telescopename=telname, x=x, y=y, z=z, dishdiameter=d,
            mount=['alt-az'], antname=an, coordsystem='global',
            referencelocation=me.observatory(telname)
        )

        # Set the polarization mode (this goes to the FEED subtable)
        sm.setfeed(mode='perfect R L', pol=[''])

        # Set the spectral window and polarization (one
        # data-description-id).
        # Call multiple times with different names for multiple SPWs or
        # pol setups.
        sm.setspwindow(
            spwname="LBand", freq='1.0GHz', deltafreq='0.1GHz',
            freqresolution='0.2GHz', nchannels=1, stokes='RR LL'
        )

        # Setup source/field information (i.e. where the observation phase
        # center is) Call multiple times for different pointings or source
        # locations.
        sm.setfield(
            sourcename="fake", sourcedirection=me.direction(
                rf=dirframe, v0=radir, v1=decdir
            )
        )
        # the second field is 10deg north of the first, so should be
        # emission free
        sm.setfield(
            sourcename="pretend", sourcedirection=me.direction(
                rf=dirframe, v0=radir, v1=qa.tos(
                    qa.add(qa.quantity(decdir), qa.quantity(offset))
                )
            )
        )
        # Set shadow/elevation limits (if you care). These set flags.
        sm.setlimits(shadowlimit=0.01, elevationlimit='1deg')

        # Leave autocorrelations out of the MS.
        sm.setauto(autocorrwt=0.0)

        # Set the integration time, and the convention to use for timerange
        # specification
        # Note : It is convenient to pick the hourangle mode as all times
        #   specified in sm.observe() will be relative to when the source
        #   transits.
        sm.settimes(
            integrationtime='2000s', usehourangle=True,
            referencetime=me.epoch('UTC', '2019/10/4/00:00:00')
        )

        # Construct MS metadata and UVW values for one scan and ddid
        # Call multiple times for multiple scans.
        # Call this with different sourcenames (fields) and spw/pol
        # settings as defined above.
        # Timesteps will be defined in intervals of 'integrationtime',
        # between starttime and stoptime.
        sm.observe(
            sourcename="fake", spwname='LBand', starttime='-5.0h',
            stoptime='-2.5h'
        )
        sm.observe(
            sourcename="pretend", spwname='LBand', starttime='-2.5h',
            stoptime='0h'
        )
        sm.observe(
            sourcename="fake", spwname='LBand', starttime='0h',
            stoptime='2.5h'
        )
        sm.observe(
            sourcename="pretend", spwname='LBand', starttime='2.5h',
            stoptime='5h'
        )

        # Close the simulator
        sm.close()
        # Unflag everything (unless you care about elevation/shadow flags)
        flagdata(vis=cls.orig_ms, mode='unflag')

    @classmethod
    def __predictSimFromComplist(cls):
        sm.openfromms(cls.orig_ms)
        # Predict from a component list
        sm.predict(complist=cls.comp_list, incremental=False)
        sm.close()

    @classmethod
    def __createImage(cls, msname, imagename, phasecenter):
        for path in glob.glob(imagename + '*'):
            shutil.rmtree(path)
        tclean(
            vis=msname, imagename=imagename, datacolumn='data',
            imsize=256, cell='8.0arcsec', gridder='standard',
            niter=20, gain=0.3, pblimit=-0.1,
            phasecenter=phasecenter
        )

    def __compare(self, imagename, radir, decdir, dirframe):
        ia.open('.'.join([imagename, 'image']))
        stats = ia.statistics()
        maxpos = stats['maxpos']
        (xc, yc) = maxpos[0:2]
        blc = [xc-10, yc-10, 0, 0]
        trc = [xc+10, yc+10, 0, 0]
        fit = ia.fitcomponents(region=rg.box(blc=blc, trc=trc))
        ia.done()
        cl.fromrecord(fit['deconvolved'])
        pos = cl.getrefdir(0)
        flux = cl.getfluxvalue(0)
        cl.done()
        expec = me.direction(dirframe, radir, decdir)
        diff = me.separation(pos, expec)
        self.assertTrue(
            qa.lt(diff, qa.quantity('0.15arcsec')),
            'position difference is too large for ' + str(pos)
            + ': ' + qa.tos(qa.convert(diff, 'arcsec'))
        )
        self.assertAlmostEqual(
            flux[0], self.exp_flux,
            msg='flux differs by too much: got: ' + str(flux[0])
            + ' expected: ' + str(self.exp_flux), delta=0.01
        )

    def __compare_ms(self, dirframe):
        tb.open(
            ctsys_resolve(
                os.path.join(
                    datadir, 'phaseshift_test_frames_expected.ms'
                )
            )
        )
        expuvw = tb.getcol('UVW')
        tb.done()
        tb.open(self.pshift_ms)
        gotuvw = tb.getcol('UVW')
        gotdata = tb.getcol('DATA')
        tb.done()
        # absolute tolerance in meters
        self.assertTrue(
            np.allclose(gotuvw, expuvw, atol=3e-6),
            'UVW do not match for ' + dirframe
        )
        # point source at phase center, so data should have
        # const amplitude and zero phase
        self.assertTrue(
            np.allclose(np.abs(gotdata), self.exp_flux, rtol=1e-6),
            'Amplitudes do not match for ' + dirframe
        )
        # phases in radians
        self.assertTrue(
            np.allclose(np.angle(gotdata), 0, atol=1e-4),
            'Phases do not match for ' + dirframe
        )

    def __run_direction_test(self, p, radir, decdir, dirframe):
        pr = [p['lon'], qa.time(p['lon'], 10)[0]]
        pd = [p['lat'], qa.time(p['lat'], 10)[0]]
        for unit in ['deg', 'rad']:
            pr.append(qa.tos(qa.convert(qa.toangle(p['lon']), unit)))
            pd.append(qa.tos(qa.convert(qa.toangle(p['lat']), unit)))
        for lon in pr:
            for lat in pd:
                shifted_pcenter = self.__phase_center_string(
                    lon, lat, p['frame']
                )
                # run phaseshift
                phaseshift(
                    vis=self.orig_ms, outputvis=self.pshift_ms,
                    phasecenter=shifted_pcenter
                )
                self.__compare_ms(p['frame'])
                # create image from phaseshifted MS
                self.__createImage(self.pshift_ms, self.pshift_im, "")
                self.__compare(self.pshift_im, radir, decdir, dirframe)
                x = imstat(self.pshift_im + '.image')
                # source should be at image center after phase shift
                self.assertTrue(
                    ((x['maxpos'] == [128, 128, 0, 0]).all()),
                    msg='Source stats' + str(x) + ' for ' + str(p)
                )
                self.__delete_intermediate_products()

    def test_frames(self):
        # This is the source position
        radir = '19h53m50'
        decdir = '40d06m00'
        dirframe = 'J2000'
        # this is the field center
        fra = '19h59m28.5'
        fdec = '+40.40.01.5'
        fframe = 'J2000'

        def create_input_ms():
            # this is how self.orig_ms was created. Do not delete
            # this or related code even if the data set is now
            # stored in the data repos. The code is useful to have
            # as a record.
            # make the MS
            self.__makeMSFrame(fra, fdec, fframe)
            # Make the component list
            self.__makeCompList(radir, decdir, dirframe)
            # Predict Visibilities
            self.__predictSimFromComplist()

        # create_input_ms()
        shutil.copytree(
            ctsys_resolve(
                os.path.join(datadir, 'phaseshift_test_frames_input.ms')
            ), self.orig_ms
        )
        # image simulated MS, the source is offset from the phase center
        # of the image
        tclean(
            vis=self.orig_ms, imagename=self.orig_im, datacolumn='data',
            imsize=2048, cell='8.0arcsec', gridder='wproject',
            niter=20, gain=0.3, wprojplanes=128, pblimit=-0.1
        )
        x = imstat(self.orig_im + '.image')
        self.assertTrue((x['maxpos'] == [1509, 773, 0, 0]).all())
        # self.__createImage(self.orig_ms, self.orig_im, orig_pcenter)
        self.__compare(self.orig_im, radir, decdir, dirframe)
        # J2000
        j2000 = {'lon': '19h53m50', 'lat': '40d06m00', 'frame': 'J2000'}
        # ICRS coordinates of the above
        icrs = {'lon': '19h53m49.9980', 'lat': '40d06m0.0019', 'frame': 'ICRS'}
        # GALACTIC coordinates of the above
        galactic = {
            'lon': '05h00m21.5326', 'lat': '+006d21m09.7433',
            'frame': 'GALACTIC'
        }
        # B1950_VLA coordinates of the above
        b1950_vla = {
            'lon': '19h52m05.65239', 'lat': '39d58m05.8512',
            'frame': 'B1950_VLA'
        }
        for p in (j2000, icrs, galactic, b1950_vla):
            self.__run_direction_test(p, radir, decdir, dirframe)
        self.__delete_data()

    def test_field(self):
        """Test that a field is correctly chosen in a multi-field MS"""
        def shift_and_clean(myfield, expdir):
            phaseshift(
                vis=self.orig_ms, outputvis=self.pshift_ms,
                phasecenter=pcenter, field=myfield
            )
            md.open(self.pshift_ms)
            exp_nfields = 2 if len(myfield) == 0 else 1
            self.assertEqual(
                md.nfields(), exp_nfields,
                msg='Wrong number of fields for field ' + myfield
            )
            sep = me.separation(md.refdir(field=0), expdir)
            md.done()
            self.assertEqual(
                qa.getvalue(sep), 0,
                msg='Ref direction is wrong for field ' + myfield
                + ' separation is ' + qa.tos(qa.convert(sep, 'arcsec'))
            )
            # check times and baselines
            field_id = ''
            if len(myfield) > 0:
                try:
                    field_id = int(myfield)
                except ValueError:
                    md.open(self.orig_ms)
                    field_id = md.fieldsforname(myfield)[0]
                    md.done()
            if field_id == 0:
                exp_ms = ctsys_resolve(
                    os.path.join(
                        datadir, 'phaseshift_test_field_0_expected.ms'
                    )
                )
            elif field_id == 1:
                exp_ms = ctsys_resolve(
                    os.path.join(
                        datadir, 'phaseshift_test_field_1_expected.ms'
                    )
                )
            elif len(field_id) == 0:
                exp_ms = ctsys_resolve(
                    os.path.join(
                        datadir, 'phaseshift_test_field_0_1_expected.ms'
                    )
                )
            tb.open(exp_ms)
            # myfilter = '' if len(myfield) == 0 else 'FIELD_ID='
            # + str(field_id)
            # x = tb.query(myfilter, columns='TIME, ANTENNA1, ANTENNA2')
            exptime = tb.getcol('TIME')
            expant1 = tb.getcol('ANTENNA1')
            expant2 = tb.getcol('ANTENNA2')
            expuvw = tb.getcol('UVW')
            expdata = tb.getcol('DATA')
            # x.done()
            tb.done()
            tb.open(self.pshift_ms)
            gottime = tb.getcol('TIME')
            gotant1 = tb.getcol('ANTENNA1')
            gotant2 = tb.getcol('ANTENNA2')
            gotuvw = tb.getcol('UVW')
            gotdata = tb.getcol('DATA')
            tb.done()
            self.assertTrue(
                (gottime == exptime).all(),
                msg='Failed TIME column test for "' + myfield + '"'
            )
            self.assertTrue(
                (gotant1 == expant1).all(),
                msg='Failed ANTENNA1 column test for "' + myfield + '"'
            )
            self.assertTrue(
                (gotant2 == expant2).all(),
                msg='Failed ANTENNA2 column test for "' + myfield + '"'
            )
            self.assertTrue(
                np.allclose(gotuvw, expuvw),
                msg='Failed UVW column test for "' + myfield + '"'
            )
            self.assertTrue(
                (gotdata == expdata).all(),
                msg='Failed DATA column test for "' + myfield + '"'
            )
            tclean(
                vis=self.pshift_ms, imagename=self.pshift_im,
                datacolumn='data', imsize=256, cell='8.0arcsec',
                gridder='standard', niter=20, gain=0.3, pblimit=-0.1
            )

        # This is the source position
        radir = '19h53m50'
        decdir = '40d06m00'
        dirframe = 'J2000'
        offset = '10deg'

        def create_ms():
            # do not delete this code, even if the MS is now in the data
            # repos; the code is useful to have as a record and a guide.
            # make the MS
            self.__sim2fields(radir, decdir, dirframe, offset)
            # Make the component list
            self.__makeCompList(radir, decdir, dirframe)
            # Predict Visibilities
            self.__predictSimFromComplist()

        # create_ms()
        shutil.copytree(
            ctsys_resolve(
                os.path.join(datadir, 'phaseshift_test_field_input.ms')
            ), self.orig_ms
        )
        # shift first field by 4 pixels north of source
        pcenter = self.__phase_center_string(
            radir,
            qa.tos(qa.add(qa.quantity(decdir), qa.quantity("32arcsec"))),
            dirframe
        )
        expdir = me.direction(dirframe, radir, decdir)
        # test both incarnations of the first field as well as both
        # fields together ('')
        for myfield in ('0', 'fake', ''):
            shift_and_clean(myfield, expdir)
            x = imstat(self.pshift_im + '.image')
            self.assertTrue(
                (x['maxpos'] == [128, 124, 0, 0]).all(),
                msg='maxpos is incorrect'
            )
            self.assertTrue(
                np.isclose(x['max'][0], self.exp_flux, 1e-6),
                msg='max is incorrect, expected ' + str(self.exp_flux)
                + ' got ' + str(x['max'][0])
            )
            self.__delete_intermediate_products()
        # 4 pixel shift of second field, which contains no signal, but
        # just sidelobes of source 10 degrees away
        decref = qa.add(qa.quantity(decdir), qa.quantity(offset))
        decdir = qa.tos(qa.add(decref, qa.quantity("32arcsec")))
        pcenter = self.__phase_center_string(radir, decdir, dirframe)
        expdir = me.direction(dirframe, radir, decref)
        for myfield in ('1', 'pretend'):
            shift_and_clean(myfield, expdir)
            x = imstat(self.pshift_im + '.image')
            self.assertTrue(
                x['max'][0]/x['rms'][0] < 5,
                msg='Incorrectly found signal in empty field, got S/N of '
                + str(x['max'][0]/x['rms'][0])
            )
            self.__delete_intermediate_products()

if __name__ == '__main__':
    unittest.main()
