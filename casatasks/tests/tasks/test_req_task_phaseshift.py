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
import glob
import numpy as np
import os
import shutil
import unittest

CASA6 = False
try:
    from casatools import (
        componentlist, ctsys, image, measures, ms, quanta, regionmanager,
        simulator, table
    )
    from casatasks import flagdata, phaseshift, tclean
    from casatasks.private import simutil
    CASA6 = True
    cl = componentlist()
    ia = image()
    me = measures()
    ms = ms()
    qa = quanta()
    rg = regionmanager()
    sm = simulator()
    tb = table()
    d = os.path.join('unittest', 'phaseshift')
    datapath = ctsys.resolve(
        os.path.join(d, 'refim_twopoints_twochan.ms')
    )
    datapath_Itziar = ctsys.resolve(
        os.path.join(d, 'Itziar.ms')
    )
    datapath_ngc = ctsys.resolve(
        os.path.join(d, 'ngc7538_ut.ms')
    )
    datapath_nep = ctsys.resolve(
        os.path.join(d, 'nep2-shrunk.ms')
    )
    datapath_mms = ctsys.resolve(
        os.path.join(d, 'uid___X02_X3d737_X1_01_small.mms')
    )
except ImportError:
    from __main__ import default
    from tasks import flagdata, phaseshift, tclean
    from taskinit import (
        cltool, iatool, metool, mstool, qatool, rgtool, smtool, tbtool
    )
    import simutil
    cl = cltool()
    ia = iatool()
    me = metool()
    ms = mstool()
    qa = qatool()
    rg = rgtool()
    sm = smtool()
    tb = tbtool()
    d = os.path.join(
        os.environ.get('CASAPATH').split()[0], 'casatestdata',
        'unittest', 'phaseshift'
    )
    datapath = os.path.join(d, 'refim_twopoints_twochan.ms')
    datapath_Itziar = os.path.join(d, 'Itziar.ms')
    datapath_ngc = os.path.join(d, 'ngc7538_ut.ms')
    datapath_nep = os.path.join(d, 'nep2-shrunk.ms')
    datapath_mms = os.path.join(d, 'uid___X02_X3d737_X1_01_small.mms')


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
        result = phaseshift(
            datacopy, outputvis=output,
            phasecenter='J2000 19h53m50 40d06m00'
        )
        # Completion without throwing an exception indicates success in CASA 6
        if not CASA6:
            self.assertTrue(result)

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
        if CASA6:
            with self.assertRaises(RuntimeError, msg=msg):
                phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35',
                    array='1'
                )
        else:
            self.assertFalse(
                phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35',
                    array='1'
                ), msg=msg
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
        if CASA6:
            with self.assertRaises(RuntimeError, msg=msg):
                phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35', observation='1'
                )
        else:
            self.assertFalse(
                phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35',
                    observation='1'
                ), msg=msg
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
        if CASA6:
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
        else:
            self.assertFalse(
                phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35', datacolumn='MODEL'
                ), msg=msg
            )
            self.assertTrue(
                phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35',
                    datacolumn='DATA'
                ), msg="phaseshift unexpectedly failed"
            )

    def test_phasecenter(self):
        '''
        Check that this parameter sets the sky coordinates of the new
        phasecenter
        '''
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
        if CASA6:
            antennalist = os.path.join(
                ctsys.resolve("alma/simmos"), "vla.d.cfg"
            )
        else:
            antennalist = os.path.join(
                '/home', 'casa', 'data', 'trunk', 'alma', 'simmos', "vla.d.cfg"
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

    def __run_direction_test(self, p, radir, decdir, dirframe):
        pr = [p['lon'], qa.time(p['lon'])[0]]
        pd = [p['lat'], qa.time(p['lat'])[0]]
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
                # create image from phaseshifted MS
                self.__createImage(self.pshift_ms, self.pshift_im, "")
                # create image from phaseshifted MS, using tclean to shift
                # phase center back to original position
                # self.__createImage(
                #    self.pshift_ms, self.pshift_shiftback_im, orig_pcenter
                # )
                self.__compare(self.pshift_im, radir, decdir, dirframe)
                # self.__compare(
                #    self.pshift_shiftback_im, radir, decdir, dirframe
                # )
                self.__delete_intermediate_products()

    def test_frames(self):
        # This is the source position
        radir = '19h53m50'
        decdir = '40d06m00'
        dirframe = 'J2000'
        # make the MS
        self.__makeMSFrame(radir, decdir, dirframe)
        # Make the component list
        self.__makeCompList(radir, decdir, dirframe)
        # Predict Visibilities
        self.__predictSimFromComplist()
        # image simulated MS, the source is offset from the phase center
        # of the image
        orig_pcenter = self.__phase_center_string(
            '19:59:28.5', '+40.40.01.5', 'J2000'
        )
        tclean(
            vis=self.orig_ms, imagename=self.orig_im, datacolumn='data',
            imsize=2048, cell='8.0arcsec', gridder='wproject',
            niter=20, gain=0.3, wprojplanes=128, pblimit=-0.1,
            phasecenter=orig_pcenter
        )
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
            'lon': '19h52m05.65239', 'lat': '40d06m00', 'frame': 'B1950_VLA'
        }
        for p in (j2000, icrs, galactic, b1950_vla):
            self.__run_direction_test(p, radir, decdir, dirframe)
        self.__delete_data()


def suite():
    return[phaseshift_test, reference_frame_tests]


if __name__ == '__main__':
    unittest.main()
