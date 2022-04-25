##########################################################################
# test_tool_coordsys.py
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
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.coordsys.html
#
#
##########################################################################
import shutil
import unittest
import os

from casatools import image as iatool
from casatools import coordsys
from casatools import ctsys
from casatools import table
from casatools import measures
from casatools import quanta
import numpy as np
import re
import math

me = measures()
qa = quanta()
cs = coordsys()
ia = iatool()

myim = "center_0.fits"

datapath='unittest/coordsys/'

class coordsys_test(unittest.TestCase):
    
    def setUp(self):
        shutil.copy(ctsys.resolve(datapath + myim), myim)
        
    def tearDown(self):
        cs.done()
        os.remove(myim)
        tb = table( )
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done( )
        
    def test_CAS_2724(self):
        myia = iatool()
        myia.open(myim)
        mycsys = myia.coordsys()
        self.assertRaises(Exception, mycsys.toworldmany, [1])

    def test_constructor(self):
        """Test constructors"""
        x = coordsys( )
        
    def test_findaxisbyname(self):
        myia = iatool()
        myia.fromshape("", [4, 4, 4, 4])
        csys = myia.coordsys( )
        myia.done()

        self.assertTrue(csys.findaxisbyname("ri") == 0)
        self.assertRaises(Exception, csys.findaxisbyname, "ra", False)
        self.assertTrue(csys.findaxisbyname("ra", True) == 0)
        
    def test_transpose(self):
        """Test the transpose() method"""
        myia = iatool()
        myia.fromshape("", [4, 4, 4, 4])
        csys = myia.coordsys()
        myia.done()
        
        orig = csys.names()
        self.assertTrue(csys.transpose([3,2,1,0]))
        new = csys.names()
        self.assertTrue(orig[0] == new[3])
        self.assertTrue(orig[1] == new[2])
        self.assertTrue(orig[2] == new[1])
        self.assertTrue(orig[3] == new[0])
        
    def test_findaxis(self):
        """Test the findaxis() method"""
        myia = iatool()
        myia.fromshape("", [20, 20, 4, 20])
        csys = myia.coordsys()
        self.assertRaises(Exception, csys.findaxis, True, -1)
        self.assertRaises(Exception, csys.findaxis, True, 4)
        got = csys.findaxis(True, 0)
        self.assertTrue(got == {'axisincoordinate': 0, 'coordinate': 0})
        got = csys.findaxis(True, 1)
        self.assertTrue(got == {'axisincoordinate': 1, 'coordinate': 0})
        got = csys.findaxis(True, 2)
        self.assertTrue(got == {'axisincoordinate': 0, 'coordinate': 1})
        got = csys.findaxis(True, 3)
        self.assertTrue(got == {'axisincoordinate': 0, 'coordinate': 2})
        myia.done()
        
    def test_findcoordinate(self):
        """Test the findcoordinate() method"""
        myia = iatool()
        myia.fromshape("", [20, 20, 4, 20])
        csys = myia.coordsys()
        myia.done()
        self.assertFalse(csys.findcoordinate('linear')['return'])
        self.assertFalse(
            csys.findcoordinate('direction', 2)['return']
        )
        got = csys.findcoordinate("direction")
        self.assertTrue(
            got['return'] and (got['pixel'] == [0,1]).all()
            and (got['world'] == [0,1]).all()
        )
        got = csys.findcoordinate("spectral")
        self.assertTrue(
            got['return'] and (got['pixel'] == [3]).all()
            and (got['world'] == [3]).all()
        )
        got = csys.findcoordinate("stokes")
        self.assertTrue(
            got['return'] and (got['pixel'] == [2]).all()
            and (got['world'] == [2]).all()
        )

    # ===== Test cases from casa5 coordsys regression test =====

    def test_constructedSetUp(self):
        """Test that the coordsys was set up with the right number of coords and type"""

        csys = coordsys()
        mycs = csys.newcoordsys()

        self.assertTrue(mycs.ncoodrdinate() == 0)
        self.assertTrue(mycs.type() == 'coordsys')
        mycs.done()

    def test_createDirectionAxes(self):
        """Test that a coordsys with direction axes can be created properly"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction = True)
        t1 = mycs.axiscoordinatetypes(True)
        t2 = mycs.axiscoordinatetypes(False)

        self.assertTrue(mycs.ncoodrdinates() == 1)
        self.assertTrue(mycs.coordinatetype(0) == ['Direction'])
        self.assertTrue(t1[1] == 'Direction' and t1[1] == 'Direction')
        self.assertTrue(t2[0] == 'Direction' and t2[1] == 'Direction')
        mycs.done()


    def test_createSpectralAxes(self):
        """Test that coordsys with spectral axes can be created properly"""

        csys = coordsys()
        mycs = csys.newcoordsys(spectral=True)
        t1 = mycs.axiscoordinatetypes(True)
        t2 = mycs.axiscoordinatetypes(False)

        self.assertTrue(mycs.ncoodrdinates() == 1)
        self.assertTrue(mycs.coordinatetype(0) == ['Spectral'])
        self.assertTrue(t1 == 'Spectral' and t2 == 'Spectral')
        mycs.done()

    def test_createStokesAxes(self):
        """Test that coordsys with stokes axes can be created properly"""

        csys = coordsys()
        mycs = csys.newcoordsys(stokes="I Q U V")
        t1 = mycs.axiscoordinatetypes(True)
        t2 = mycs.axiscoordinatetypes(False)

        self.assertTrue(mycs.ncoodrdinates() == 1)
        self.assertTrue(mycs.coordinatetype(0) == ['Stokes'])
        self.assertTrue(t1 == 'Stokes' and t2 == 'Stokes')
        mycs.done()

    def tet_createLinearAxes(self):
        """"Test that coordsys with linear axes can be created properly"""

        csys = coordsys()
        mycs = csys.newcoordsys(linear=3)
        t1 = mycs.axiscoordinatetypes(True)
        t2 = mycs.axiscoordinatetypes(False)

        self.assertTrue(mycs.ncoodrdinates() == 1)
        self.assertTrue(mycs.coordinatetype(0) == ['Linear'])
        self.assertTrue(t1[0] == 'Linear' and t2[0] == 'Linear')
        mycs.done()

    def test_createTabularAxes(self):
        """Test that coordsys with tabular axes can be created properly"""

        csys = coordsys()
        mycs = csys.newcoordsys(tabular=True)
        t1 = mycs.axiscoordinatetypes(True)
        t2 = mycs.axiscoordinatetypes(False)

        self.assertTrue(mycs.ncoodrdinates() == 1)
        self.assertTrue(mycs.coordinatetype(0) == ['Tabular'])
        self.assertTrue(t1 == 'Tabular' and t2 == 'Tabular')
        mycs.done()

    def test_createMixedAxes(self):
        """Test that a coordsys with mixed axes types can be created properly"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True, spectral=True, stokes="I Q U V",
                              linear=1, tabular=True)
        types = mycs.coordinatetype()
        typesRes = ['Direction', 'Stokes', 'Spectral', 'Linear', 'Tabular']

        self.assertTrue(mycs.ncoodrdinates() == 5)
        self.assertTrue(np.all(types == typesRes))

        t1 = mycs.axiscoordinatetypes(True)
        t2 = mycs.axiscoordinatetypes(False)

        tRes = ['Direction', 'Direction', 'Stokes', 'Spectral', 'Linear', 'Tabular']

        self.assertTrue(np.all(t1 == tRes))
        self.assertTrue(np.all(t2 == tRes))
        mycs.done()

    def test_refCode(self):
        """Test that a reference code can be set"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True, spectral=True)
        d = me.frequency("LSRK")
        list = me.listcode(d)

        for i in list['normal']:
            if (i!='REST'):
                mycs.setreferencecode(type='spectral', value=i, adjust=True)
                if (mycs.referencecode(type='spectral')!=[i]):
                    self.fail()
        self.assertTrue(True)

    def test_setDirection(self):
        """Test that direction can be set"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True, spectral=True)
        d = me.direction('J2000')
        list = me.listcodes(d)
        for i in list['normal']:
            if not re.search('AZEL', i, re.IGNORECASE):
                mycs.setreferencecode(type='direction', value=i, adjust=False)
                if (mycs.referencecode(type='direction')!=[i]):
                    self.fail()
        mycs = csys.newcoordsys(direction=True, spectral=True, linear=1)
        mycs.setreferencecode(type='direction', value='B1950')
        mycs.setreferencecode(type='spectral', value='BARY')
        c = mycs.referencecode()

        self.assertTrue(len(c) == 3)
        self.assertTrue(c[0] == 'B1950')
        self.assertTrue(c[1] == 'BARY')
        self.assertTrue(c[2] == '')

    def test_projection(self):
        """Test that the projection can be set"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True, spectral=True, linear=1)
        mycs.setprojection('SIN', [1.0, 2.0])
        p = mycs.projection()

        self.assertTrue(p['type'] == 'SIN' and len(p['parameters']) == 2)
        self.assertTrue(p['parameters'][0] == 1.0)
        self.assertTrue(p['parameters'][1] == 2.0)

        p = mycs.projection('all')['types']
        self.assertTrue(len(p) == 27)
        for i in p:
            n = mycs.projectio(i)
            if not n: self.fail()

    def test_restFreq(self):
        """Test that the rest frequency can be set"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True, spectral=True, linear=1)
        rf1 = qa.quantity('1.2GHz')
        rf2 = mycs.restfrequency()
        rf2 = qa.convert(rf2, rf1['unit'])

        self.assertTrue(abs(qa.getvalue(rf1) - qa.getvalue(rf2)) < 1.0e-6)
        self.assertTrue(qa.getunit(rf1) == qa.getunit(rf2))

        unit = qa.getunit(mycs.restfrequency())
        rf1 = 2.0
        mycs.setrestfrequency(rf1)
        rf2 = mycs.restfrequency()
        rf1 = qa.unit(rf1, unit)

        self.assertTrue(abs(qa.getvalue(rf1)-qa.getvalue(rf2))<1.0e-6)
        self.assertTrue(qa.getunit(rf1)==qa.getunit(rf2))

        rf1 = qa.quantity([1e9, 2e9], 'Hz')
        # Select second freq
        mycs.setrestfrequency(value=rf1, which=1, append=False)
        rf2 = qa.convert(mycs.restfrequency(), qa.getunit(rf1))
        v1 = qa.getvalue(rf1)
        v2 = qa.getvalue(rf2)

        self.assertTrue(abs(v1[0] - v2[1]) < 1e-6)
        self.assertTrue(abs(v1[1] - v2[0]) < 1e-6)
        self.assertTrue(qa.getunit(rf1) == qa.getunit(rf2))

    def test_toRecord(self):
        """Test that the coordsys can be converted to a record"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True, spectral=True, stokes="I Q U V", linear=3)

        rec = mycs.torecord()

        self.assertTrue(rec.has_key('direction0'))
        self.assertTrue(rec.has_key('stokes1'))
        self.assertTrue(rec.has_key('spectral2'))
        self.assertTrue(rec.has_key('linear3'))

    def test_fromRecord(self):
        """Test that a record can be converted to a coordsys"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True, spectral=True, stokes="I Q U V", linear=3)

        rec = mycs.torecord()
        cs2 = csys.newcoordsys(direction=False, spectral=False, stokes="", linear=0)
        self.assertTrue(cs2.fromrecord(rec))
        mycs.done()
        cs2.done()

    def test_copy(self):
        """Test that the coordsys can be copied"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True)
        cs2 = mycs.copy()

        self.assertTrue(cs2)
        closed = mycs.done()
        self.assertTrue(closed and cs2)
        cs2.done()

    def test_epoch(self):
        """Test that you can set and get the epoch"""

        csys = coordsys()
        mycs = csys.newcoordsys(direction=True)
        epoch1 = me.epoch('UTC', 'today')
        mycs.setepoch(epoch1)
        epoch2 = mycs.epoch()

        self.assertTrue(abs(me.getvalue(epoch1)['m0']['value']-me.getvalue(epoch2)['m0']['value'])<1.0e-6)
        self.assertTrue(me.getvalue(epoch1)['m0']['unit'] == me.getvalue(epoch2)['m0']['unit'])
        self.assertTrue(me.gettype(epoch1) == me.gettype(epoch2))
        self.assertTrue(me.getref(epoch1) == me.getref(epoch2))
        mycs.done()

    def test_observer(self):
        """Test that you can set and get the observer"""
        csys = coordsys()
        mycs = csys.newcoordsys(direction=True)
        obs1 = 'Biggles'
        mycs.setobserver(obs1)
        obs2 = mycs.observer()

        self.assertTrue(obs1 == obs2)
        mycs.done()

    def test_telescope(self):
        """Test that you can set and get the telescope"""
        csys = coordsys()
        mycs = csys.newcoordsys(direction=True)
        tel1 = 'VLA'
        mycs.settelescope(tel1)
        tel2 = mycs.telescope()

        self.assertTrue(tel1 == tel2)
        mycs.done()

    def test_parentName(self):
        """Test that you can set and get the parent name"""
        csys = coordsys()
        mycs = csys.newcoordsys(direction=True)
        pn1 = 'Biggles.image'
        mycs.setparentname(pn1)
        pn2 = mycs.parentname()

        self.assertTrue(pn1 == pn2)
        mycs.done()

    def test_names(self):
        """Test that you can set and get names"""

        mycs = cs.newcoordsys(direction=True, spectral=True)
        val1 = "a b c"
        mycs.setnames(value=val1)
        val2 = mycs.names()
        self.assertTrue(val1.split() == val2)

        val2 = mycs.names('spec')
        self.assertTrue(val2[0] == val1.split()[2])

        val1 = 'fish'
        mycs.setnames(type='spec', value=val1)
        val2 = mycs.names('spec')
        self.assertTrue(val2==[val1])
        mycs.done()

    def test_units(self):
        """Test that you can set and get unit names"""

        mycs = cs.newcoordsys(direction=True, spectral=True)
        val1 = "deg rad GHz"
        val2 = mycs.units()

        self.assertTrue(np.all(val1.split() == val2))
        # self.assertFalse(mycs.setunits(value="Hz Hz Hz"))
        # self.assertFalse(mycs.setunits(value="m"))

        mycs.setunits(value=val1)
        val2 = mycs.units('spec')
        self.assertTrue(val2[0] == val2.split()[2])

        val1 = "kHz"
        mycs.setunits(type='spec', value=val1)
        val2 = mycs.units('spec')
        self.assertTrue(val2[0] == val1)

        mycs = cs.newcoordsys(direction=True, linear=2)
        val1 = "Hz kHz"
        mycs.setunits(type='linear', value=val1, overwrite=True)
        val2 = mycs.units()

        self.assertTrue(val1.split()[0] == val2[2] and val1.split()[1] == val2[3])
        mycs.done()

    def test_refPixel(self):
        """Test that you can get and set the reference pixel"""

        mycs = cs.newcoordsys(direction=True, spectral=True)
        val1 = [0, 1, 2]
        mycs.setreferencepixel(value=val1)
        val2 = mycs.referencepixel()['numeric']
        self.assertTrue(abs(val1[0]-val2[0]) < 1.0e-6)
        self.assertTrue(abs(val1[1]-val2[1]) < 1.0e-6)
        self.assertTrue(abs(val1[2]-val2[2]) < 1.0e-6)

        val2 = mycs.referencepixel('dir')['numeric']
        self.assertTrue(abs(val2[0] - val1[0]) < 1.0e-6)
        self.assertTrue(abs(val2[1] - val1[1]) < 1.0e-6)

        val2 = mycs.referencepixel('spec')['numeric']
        self.assertTrue(abs(val2 - val1[2]) < 1.0e-6)

        val1 = [0, 0]
        mycs.setreferencepixel(type='dir', value=val1)
        val2 = mycs.referencepixel('dir')['numeric']
        self.assertTrue(abs(val2[0] - val1[0]) < 1e-6)
        self.assertTrue(abs(val2[1] - val1[1]) < 1e-6)

        mycs.done()

    def test_linearTransformDirection(self):
        """Test set linear transform type to direction"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes='IQ', tabular=True, linear=3)
        val1 = ia.makearray(0, [2, 2])
        val1[0, 0] = math.cos(0.5)
        val1[0, 1] = -math.sin(0.5)
        val1[1, 0] = math.sin(0.5)
        val1[1, 1] = math.cos(0.5)
        type = 'direction'
        mycs.setlineartransform(value=val1, type=type)
        val2 = mycs.lineartransform(type=type)

        self.assertTrue(np.all(val1 == val2 == 1.0e-6))

    def test_linearTransformSpectral(self):
        """Test set linear transform type to spectral """
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes='IQ', tabular=True, linear=3)
        val1 = ia.makearray(2, [1, 1])
        type = 'spectral'
        mycs.setlineartransform(value=val1, type=type)
        val2 = mycs.lineartransform(type=type)
        self.assertTrue(np.all(val1 == val2 == 1.0e-6))

    def test_linearTransformStokes(self):
        """Test set linear transform type to stokes"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes='IQ', tabular=True, linear=3)
        val1 = ia.makearray(2, [1, 1])
        type = 'stokes'
        mycs.setlineartransform(value=val1, type=type)
        val2 = mycs.lineartransform(type=type)
        self.assertTrue(np.all(val1 == val2 == 1.0e-6))

    def test_linearTransformTabular(self):
        """Test set linear transform type to tabular"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes='IQ', tabular=True, linear=3)
        val1 = ia.makearray(4, [1, 1])
        type = 'tabular'
        mycs.setlineartransform(value=val1, type=type)
        val2 = mycs.lineartransform(type=type)
        self.assertTrue(np.all(val1 == val2 == 1.0e-6))

    def test_linearTransformLinear(self):
        """Test set linear transform type to linear"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes='IQ', tabular=True, linear=3)
        val1 = ia.makearray(0, [3, 3])
        val1[0, 0] = 2.0
        val1[1, 1] = 3.0
        val1[2, 2] = 4.0
        type = 'linear'
        mycs.setlineartransform(value=val1, type=type)
        val2 = mycs.lineartransform(type=type)
        self.assertTrue(np.all(val1 == val2 == 1.0e-6))

    def test_referenceValueSetValRad(self):
        """Test set reference value in radians"""
        mycs = cs.newcoordsys(direction=True)
        mycs.setunits(value="rad rad")
        val1 = mycs.referencevalue(format='q')
        val1['quantity']['*1'] = qa.quantity('0.01rad')
        val1['quantity']['*2'] = qa.quantity('-0.01rad')
        mycs.setreferencevalue(value=val1)
        val2 = mycs.referencevalue(format='q')

        self.assertTrue(abs(val1['quantity']['*1']['value'] - val2['quantity']['*1']['value']) < 1e-6)
        self.assertTrue(abs(val1['quantity']['*2']['value'] - val2['quantity']['*2']['value']) < 1e-6)
        self.assertTrue(val1['quantity']['*1']['unit'] == val2['quantity']['*1']['unit'])
        self.assertTrue(val1['quantity']['*2']['unit'] == val2['quantity']['*2']['unit'])
        mycs.done()

    def test_referenceValueSetValDefaultUnits(self):
        """Test set reference value using default units"""
        mycs = cs.newcoordsys(direction=True)
        units = mycs.units()
        val1 = [1.0, 2.0]
        mycs.setreferencevalue(value=val1)
        val2 = mycs.referencevalue(format='q')
        self.assertTrue(abs(val1[0] - val2['quantity']['*1']['value']) < 1e-6)
        self.assertTrue(abs(val1[1] - val2['quantity']['*2']['value']) < 1e-6)
        self.assertTrue(units[0] == val2['quantity']['*1']['unit'])
        self.assertTrue(units[1] == val2['quantity']['*2']['unit'])
        mycs.done()

    def test_referenceValueSetIncorrect(self):
        """Test setting the value incorrectly"""
        mycs = cs.newcoordsys(spectral=True)
        try:
            mycs.setreferencevalue(value='i like doggies')
        except Exception, e:
            pass
        self.fail()
        mycs.done()

    def test_referenceValueSetSpectral(self):
        """Test set reference value with spectral type"""
        mycs = cs.newcoordsys(direction=True, spectral=True)
        val1 = mycs.referencevalue(format='q')
        val2 = mycs.referencevalue(type='spec', format='q')
        self.assertTrue(abs(val1['quantity']['*3']['value'] - val2['quantity']['*1']['value']) < 1e-6 and val1['quantity']['*3'][
            'unit'] == val2['quantity']['*1']['unit'])

        val1 = [-10]
        mycs.setreferencevalue(type='spec', value=val1)
        val2 = mycs.referencevalue(type='spec', format='n')
        self.assertTrue(abs(val1[0] - val2['numeric'][0]) < 1e-6)

        val1 = mycs.referencevalue(format='n')
        val2 = mycs.referencevalue(type='spec', format='n')
        self.assertTrue(abs(val1['numeric'][2] - val2['numeric'][0]) < 1e-6)
        mycs.done()

    def test_incrementRad(self):
        """Test increment with radian units"""
        mycs = cs.newcoordsys(direction=True)
        mycs.setunits(value="rad rad")
        val1 = mycs.increment(format='q')
        val1['quantity']['*1'] = qa.quantity('0.01rad')
        val1['quantity']['*2'] = qa.quantity('-0.01rad')
        mycs.setincrement(value=val1)
        val2 = mycs.increment(format='q')
        #
        self.assertTrue(abs(val1['quantity']['*1']['value'] - val2['quantity']['*1']['value']) < 1e-6)
        self.assertTrue(abs(val1['quantity']['*2']['value'] - val2['quantity']['*2']['value']) < 1e-6)
        self.assertTrue(val1['quantity']['*1']['unit'] == val2['quantity']['*1']['unit'])
        self.assertTrue(val1['quantity']['*2']['unit'] == val2['quantity']['*2']['unit'])
        mycs.done()

    def test_incrementDefaultUnits(self):
        """Test increment with default coordsys units"""
        mycs = cs.newcoordsys(direction=True)
        units = mycs.units()
        val1 = [1.0, 2.0]
        mycs.setincrement(value=val1)
        val2 = mycs.increment(format='q')

        self.assertTrue(abs(val1[0] - val2['quantity']['*1']['value']) < 1e-6)
        self.assertTrue(abs(val1[1] - val2['quantity']['*2']['value']) < 1e-6)
        self.assertTrue(units[0] == val2['quantity']['*1']['unit'])
        self.assertTrue(units[1] == val2['quantity']['*2']['unit'])
        mycs.done()
        #

    def test_incrementSetIncorrect(self):
        """Test increment with incorrect values"""
        mycs = cs.newcoordsys(spectral=True)
        try:
            mycs.setincrement(value='i like doggies')
        except Exception, e:
            pass
        mycs.close()
        self.fail()

    def test_incrementSpectral(self):
        """Test increment using spectral type"""
        mycs = cs.newcoordsys(direction=True, spectral=True)
        val1 = [1.0, 2.0, 3.0]
        mycs.setincrement(value=val1)
        val2 = mycs.increment(format='n', type='dir')
        self.assertTrue(abs(val2['numeric'][0] - val1[0]) < 1e-6 and abs(val2['numeric'][1] - val1[1]) < 1e-6)
        val2 = mycs.increment(type='spe', format='n')
        self.assertTrue(abs(val2['numeric'][0] - val1[2]) < 1e-6)
        try:
            val2 = mycs.increment(type='lin', format='q')
        except Exception, e:
            val2 = False
            self.assertTrue(True)
        if val2: self.fail()
        #
        val1 = [-10]
        mycs.setincrement(type='spec', value=val1)
        val2 = mycs.increment(type='spec', format='n')
        self.assertTrue(abs(val1[0] - val2['numeric'][0]) < 1e-6)
        mycs.done()

    def test_stokes(self):
        """Test that you can create a coordsys with the given stokes values"""
        mycs = cs.newcoordsys(stokes="I RL")
        stokes = mycs.stokes()

        self.assertTrue(stokes[0] == 'I' and stokes[1] == 'RL', msg='Stokes recovered wrong values')
        mycs.done()

    def test_stokesSetValue(self):
        """Test that you can set the stokes values after creation"""
        mycs = cs.newcoordsys(stokes="I RL")
        mycs.setstokes("XX V")
        stokes = mycs.stokes()

        self.assertTrue(stokes[0] == 'XX' and stokes[1] == 'V', msg='Stokes recovered wrong values')
        mycs.done()

    def test_stokesWrongValues(self):
        """Test that incorrect values raise exceptions"""
        mycs = cs.newcoordsys(direction=True)
        try:
            stokes = True
            stokes = mycs.stokes()
        except Exception, e:
            stokes = False
        self.assertFalse(stokes)

        try:
            ok = True
            ok = mycs.setstokes("I V")
        except Exception, e:
            ok = False
        self.assertFalse(ok)
        mycs.done()

    def test_toWorld(self):
        """Test the toworld conversion"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes="I V", linear=2)
        rp = mycs.referencepixel()['numeric']
        rv = mycs.referencevalue(format='n')

        self.assertTrue(len(rp) == 6)
        self.assertTrue(np.isclose(mycs.toworld(value=rp, format='n')['numeric'], rv['numeric'], atol=1e-6))
        ###
        d = mycs.toworld(value=list(rp), format='q')
        u = mycs.units()

        self.assertTrue(len(d)['quantity'] == len(rv['numeric']))

        for i in range(len(d['quantity'])):
            self.assertTrue(abs(d['quantity']['*' + str(i + 1)]['value'] - rv['numeric'][i]) > 1e-6)
            self.assertTrue(d['quantity']['*' + str(i + 1)]['unit'] == u[i])
        ###
        q = mycs.toworld(value=rp, format='q')
        m = mycs.toworld(value=rp, format='m')
        m = m['measure']

        self.assertTrue(m.has_key('direction') and m.has_key('spectral'))
        self.assertTrue(m['spectral'].has_key('frequency'))
        self.assertTrue(m['spectral'].has_key('opticalvelocity'))
        self.assertTrue(m['spectral'].has_key('radiovelocity'))
        self.assertTrue(m['spectral'].has_key('betavelocity'))
        self.assertTrue(m.has_key('stokes'))
        self.assertTrue(m.has_key('linear'))

        d = m['direction']
        f = m['spectral']['frequency']
        l = m['linear']
        s = m['stokes']

        v = me.getvalue(d)
        q['quantity']['*1'] = qa.convert(q['quantity']['*1'], v['m0']['unit'])
        q['quantity']['*2'] = qa.convert(q['quantity']['*2'], v['m1']['unit'])

        self.assertTrue(abs(v['m0']['value'] - q['quantity']['*1']['value']) < 1e-6)
        self.assertTrue(abs(v['m1']['value'] - q['quantity']['*2']['value']) < 1e-6)
        self.assertTrue(v['m0']['unit'] == q['quantity']['*1']['unit'])
        self.assertTrue(v['m1']['unit'] == q['quantity']['*2']['unit'])

        v = me.getvalue(f)
        q['quantity']['*4'] = qa.convert(q['quantity']['*4'], v['m0']['unit'])

        self.assertTrue(abs(v['m0']['value'] - q['quantity']['*4']['value']) < 1e-6)
        self.assertTrue(v['m0']['unit'] == q['quantity']['*4']['unit'])

        q['quantity']['*5'] = qa.convert(q['quantity']['*5'], l['*1']['unit'])
        q['quantity']['*6'] = qa.convert(q['quantity']['*6'], l['*2']['unit'])

        self.assertTrue(abs(l['*1']['value'] - q['quantity']['*5']['value']) < 1e-6)
        self.assertTrue(abs(l['*2']['value'] - q['quantity']['*6']['value']) < 1e-6)
        self.assertTrue(l['*1']['unit'] == q['quantity']['*5']['unit'])
        self.assertTrue(l['*2']['unit'] == q['quantity']['*6']['unit'])
        mycs.done()


    def test_toWorldMany(self):
        """Test the toworldmany parameter"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes="I V", linear=2)

        p = mycs.referencepixel()['numeric']
        w = mycs.referencevalue()
        rIn = ia.makearray(0, [len(p), 10])
        for i in range(10):
            for j in range(len(p)):
                rIn[j, i] = p[j]
        rOut = mycs.toworldmany(rIn)

        self.assertTrue(len(rOut['numeric']) == len(rIn))
        for i in range(10):
            for j in range(len(p)):
                self.assertTrue((rOut['numeric'][j, i] - w['numeric'][j]) < 1e-6)
        mycs.done()

    def test_toPixel(self):
        """Test topixel parameter"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes="I V", linear=2)
        tol = 1.0e-6
        rp = mycs.referencepixel()['numeric']

        self.assertTrue(len(rp) == 6)

        rv = mycs.referencevalue(format='n')
        p = mycs.topixel(value=rv)['numeric']

        self.assertTrue(len(p) == 6)
        self.assertTrue(np.all(np.isclose(p, rp, atol=tol)))

        for format in ["n", "q", "m", "s", "nqms"]:
            for i in range(len(rp)): p[i] = rp[i] + 1
            w = mycs.toworld(value=p, format=format)
            p2 = mycs.topixel(value=w)['numeric']

            self.assertTrue(len(p2) == 6)
            self.assertTrue(np.all(np.isclose(p, rp, atol=tol)),
                            msg='toworld/topixel reflection failed for format "' + format + '"')
        mycs.done()

    def test_toPixelMany(self):
        """Test topixelmany parameter"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes="I V", linear=2)
        n = 10
        p = mycs.referencepixel()['numeric']
        w = mycs.toworld(p, 'n')
        w = w['numeric']
        rIn = ia.makearray(0, [len(w), n])
        for i in range(n):
            for j in range(len(w)):
                rIn[j, i] = w[j]
        r2 = mycs.topixelmany(rIn)

        self.assertTrue(len(r2['numeric']) == len(rIn))
        for i in range(n):
            for j in range(len(w)):
                self.assertTrue(abs(p[j] - r2['numeric'][j, i]) < 1e-6)
        mycs.done()

    def test_naxes(self):
        """Test the naxes param with 0 or multiple axes"""
        mycs = cs.newcoordsys()
        n = mycs.naxes()

        self.assertTrue(n == 0)
        self.assertTrue(mycs.naxes() == 0)

        mycs = cs.newcoordsys(direction=True, spectral=True, stokes="I V", linear=2)

        self.assertTrue(mycs.naxes == 6)

    def test_axesmap(self):
        """Test that maps are the same"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes="I V", linear=2)
        toworld = mycs.axesmap(toworld=True)
        topixel = mycs.axesmap(toworld=False)

        idx = range(0, len(mycs.referencepixel()['numeric']))

        self.assertTrue(np.all(toworld == idx))
        self.assertTrue(np.all(topixel == idx))
        mycs.done()

    def test_reorder(self):
        """Test that the coords are reordered"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes='I V', linear=1)
        order = [3, 2, 1, 0]
        mycs.reorder(order)

        self.assertTrue(mycs.coordinatetype(0) == ['Linear'])
        self.assertTrue(mycs.coordinatetype(1) == ['Spectral'])
        self.assertTrue(mycs.coordinatetype(2) == ['Stokes'])
        self.assertTrue(mycs.coordinatetype(3) == ['Direction'])

    def test_reorderInvalid(self):
        """Test that invalid inputs for reorder raises exceptions"""
        mycs = cs.newcoordsys(direction=True, spectral=True, stokes='I V', linear=1)
        try:
            ok = mycs.reorder([1, 2])
        except Exception, e:
            ok = False
        self.assertFalse(ok)

        try:
            ok = mycs.reorder([1, 2, 3, 10])
        except Exception, e:
            ok = False
        self.assertFalse(ok)
        mycs.done()


if __name__ == '__main__':
    unittest.main()
