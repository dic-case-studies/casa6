##########################################################################
# test_simulator
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
# [CAS-13439]
#
#
##########################################################################
 
 
####    Imports     ####
import glob
import os
import sys
import shutil
import json
import numpy
import unittest
  
is_CASA6 = False
 
try:
    from casatools import ctsys, simulator, componentlist, table, agentflagger
    from casatasks.private.simutil import *

    is_CASA6 = True
    _sm = simulator()
    _cl = componentlist()
    _tb = table()
    _af = agentflagger()
 
    # Location of input data
    datapath = ctsys.resolve('unittest/simulator/')
    refpath = ctsys.resolve('unittest/simulator/smtool_reference/')
except ImportError:
    from __main__ import default
    from taskinit import * 
    _sm = smtool()
    _cl = cltool()
    _tb = tbtool()
    _af = agentflagger()
    
    # Location of input data
    datapath = os.path.join(os.environ['CASAPATH'].split()[0], 'casatestdata/unittest/simulator/')
    refpath = datapath+'smtool_reference/'

class sm_predict_test(unittest.TestCase):
    """
    """   
    comp_list = 'mycomplist.cl'
    orig_ms = 'myms.ms'
    ref_flux = 5
    # Antenna list created from CASA's distro data in alma/simmos/vla.d.cfg
    # after executing simutil.readantenna()
    antennalist = 'antlist_simutil.json'

    @classmethod
    def __delete(cls):
        m = glob.glob(cls.orig_ms + '*')
        m.append(cls.comp_list)
        for x in m:
            if os.path.exists(x):
                shutil.rmtree(x)
        

    def setUp(self):
        self.__delete()
        # Copy the antennalist json file to local directory
        shutil.copy(os.path.join(refpath,self.antennalist), self.antennalist)


    def tearDown(self):
        self.__delete()
        
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.antennalist)

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
        _sm.open(ms=cls.orig_ms)


        # Load the json file for a fictitious telescope that  can be simulated 
        # by specifying x, y, z, d, an, an2, telname, antpos.
        # x,y,z are locations in meters in ITRF (Earth centered)
        # coordinates.
        # d, an are lists of antenna diameter and name.
        # telname and obspos are the name and coordinates of the
        # observatory.
        with open(cls.antennalist, 'r') as readfile:
            jobj = json.load(readfile)
        
        x = jobj['x']
        y = jobj['y']
        z = jobj['z']
        d = jobj['diameter']
        an = jobj['ant_names']
        an2 = jobj['ant_names2']
        telname = jobj['telname']
        obspos = jobj['obspos']
        
        # Set the antenna configuration
        _sm.setconfig(
            telescopename=telname, x=x, y=y, z=z, dishdiameter=d,
            mount=['alt-az'], antname=an, coordsystem='global',
            referencelocation=me.observatory(telname)
        )

        # Set the polarization mode (this goes to the FEED subtable)
        _sm.setfeed(mode='perfect R L', pol=[''])

        # Set the spectral window and polarization (one
        # data-description-id).
        # Call multiple times with different names for multiple SPWs or
        # pol setups.
        _sm.setspwindow(
            spwname="LBand", freq='1.0GHz', deltafreq='0.1GHz',
            freqresolution='0.1GHz', nchannels=5, stokes='RR LL'
        )

        # Setup source/field information (i.e. where the observation phase
        # center is) Call multiple times for different pointings or source
        # locations.
        _sm.setfield(
            sourcename="fake", sourcedirection=me.direction(
                rf=dirframe, v0=radir, v1=decdir
            )
        )

        # Set shadow/elevation limits (if you care). These set flags.
        _sm.setlimits(shadowlimit=0.01, elevationlimit='1deg')

        # Leave autocorrelations out of the MS.
        _sm.setauto(autocorrwt=0.0)

        # Set the integration time, and the convention to use for timerange
        # specification
        # Note : It is convenient to pick the hourangle mode as all times
        #   specified in sm.observe() will be relative to when the source
        #   transits.
        _sm.settimes(
            integrationtime='2000s', usehourangle=True,
            referencetime=me.epoch('UTC', '2019/10/4/00:00:00')
        )

        # Construct MS metadata and UVW values for one scan and ddid
        # Call multiple times for multiple scans.
        # Call this with different sourcenames (fields) and spw/pol
        # settings as defined above.
        # Timesteps will be defined in intervals of 'integrationtime',
        # between starttime and stoptime.
        _sm.observe(
            sourcename="fake", spwname='LBand', starttime='-5.0h',
            stoptime='+5.0h'
        )
        # Close the simulator
        _sm.close()
        
        # Unflag everything (unless you care about elevation/shadow flags)
        _af.open(cls.orig_ms)
        _af.selectdata()
        agentUnflag={'apply':True,'mode':'unflag'}
        _af.parseagentparameters(agentUnflag)
        _af.init()
        _af.run(writeflags=True)
        _af.done()

    @classmethod
    def __makeCompList(cls, ra, dec, frame):
        _cl.addcomponent(
            dir=cls.__phase_center_string(ra, dec, frame),
            flux=cls.ref_flux,      # For a gaussian, this is the
                                    # integrated area.
            fluxunit='Jy', freq='1.5GHz', shape='point',
            spectrumtype="plp",  index=[3,2]
        )
        # Save the file
        _cl.rename(filename=cls.comp_list)
        _cl.done()

    def test_plp(self):
        """CAS-13439 verify support for plp, spectral curvature model"""
        # This is the source position
        radir = '19h53m50'
        decdir = '40d06m00'
        dirframe = 'J2000'
        # this is the field center
        fra = radir
        fdec = decdir
        fframe = dirframe

        self.__makeMSFrame(fra, fdec, fframe)
        # Make the component list
        self.__makeCompList(radir, decdir, dirframe)
        # Predict Visibilities
        _sm.openfromms(self.orig_ms)
        # Predict from a component list
        _sm.predict(complist=self.comp_list, incremental=False)
        _sm.close()
        _tb.open(self.orig_ms)
        x = _tb.getcol('DATA')[0,:,50]
        _tb.done()
        r = numpy.array([1, 1.1, 1.2, 1.3, 1.4])/1.5
        expec = 5*r**(3 + 2*numpy.log(r))
        self.assertTrue(
            numpy.allclose(numpy.real(x), expec, rtol=5e-8),
            'Incorrect visibility amplitiudes'
        )
        
def suite():
    return [sm_predict_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()

