#########################################################################
# test_task_fixplanets.py
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
#
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.fixplanets.html
#
##########################################################################
import os
import sys
import shutil
import unittest

from casatools import ctsys
from casatools import ms as mstool
from casatools import msmetadata as msmdtool
from casatasks import fixplanets

'''
Unit tests for task fixplanets.

Features tested:                                                       
  1. Does a standard fixplanets work on an MS imported from an ASDM from April 2011
  2. Does the setting of a given direction work on an MS imported from an ASDM from April 2011
  3. Does the setting of a given direction with ref !=J2000 and != sol.sys. object give the expected error?
  4. Does the setting of a given direction work with a sol system ref frame
  5. Does the use of an ephemeris via the direction parameter work

'''
datapath = ctsys.resolve('unittest/fixplanets/')
outms = 'uid___A002_X1c6e54_X223-thinned.ms'
inpms = os.path.join(datapath, outms)
outms2 = 'uid___A002_X1c6e54_X223-thinned.mms/'
inpms2 = os.path.join(datapath,outms2)

mymst = mstool()
mymsmdt = msmdtool()


class fixplanets_test1(unittest.TestCase):
    def setUp(self):
        res = None
        shutil.rmtree(outms, ignore_errors=True)
        shutil.copytree(inpms, outms)
        shutil.rmtree(outms2, ignore_errors=True)
        os.system('cp -R '+ inpms2 + ' ' + outms2)
        
    def tearDown(self):
        shutil.rmtree(outms, ignore_errors=True)
        shutil.rmtree(outms2, ignore_errors=True)

    def verify(self, thems, thefield, theref):
        therval = True
        mymsmdt.open(thems)
        thefieldids = mymsmdt.fieldsforname(thefield)
        mymsmdt.close()
        mymst.open(thems)
        thedir = mymst.getfielddirmeas(fieldid=thefieldids[0])
        mymst.close()
        print("Read direction result %s" % thedir)
        if not (thedir['refer']==theref):
            print("ERROR: reference not as expected: expected %s, got %s" % (theref,thedir['refer']))
            therval = False
        return therval
    
    def test1(self):
        '''test1: Does a standard fixplanets work on an MS imported from an ASDM from April 2011'''
        for myms in [outms,outms2]:
            fixplanets(myms, 'Titan', True)

    def test2(self):
        '''test2: Does the setting of a given direction work on an MS imported from an ASDM from April 2011'''
        for myms in [outms,outms2]:
            fixplanets(myms, 'Titan', False, 'J2000 0h0m0s 0d0m0s')
            self.assertTrue(self.verify(myms, 'Titan', 'J2000'))

    def test3(self):
        '''test3: Does the setting of a given direction with ref !=J2000 and != sol.sys. object give the expected error?'''
        for myms in [outms,outms2]:
            with self.assertRaises(RuntimeError):
                fixplanets(myms, 'Titan', False, 'B1950 0h0m0s 0d0m0s')

    def test4(self):
        '''test4: Does the setting of a given direction work with a sol system ref frame?'''
        for myms in [outms,outms2]:
            fixplanets(myms, 'Titan', False, 'SATURN 0h0m0s 0d0m0s')
            self.assertTrue(self.verify(myms, 'Titan', 'SATURN'))

    def test5(self):
        '''test5: Does a standard fixplanets work on an MS imported from an ASDM from April 2011 with parameter reftime'''
        for myms in [outms,outms2]:
            fixplanets(vis=myms, field='Titan', fixuvw=True, reftime='median')

    def test6(self):
        '''test6: Does a standard fixplanets with put of bounds parameter reftime give the expected error'''
        for myms in [outms,outms2]:
            with self.assertRaises(TypeError):
                fixplanets(vis=myms, field='Titan', fixuvw=True, reftime='2012/07/11/08:41:32')

    def test7(self):
        '''test7: Does a standard fixplanets with wrong parameter reftime give the expected error'''
        for myms in [outms,outms2]:
            with self.assertRaises(TypeError):
                fixplanets(vis=myms, field='Titan', fixuvw=True, reftime='MUDIAN')

    def test8(self):
        '''test8: Does a fixplanets with an ephemeris work'''
        for myms in [outms,outms2]:
            fixplanets(vis=myms, field='Titan', fixuvw=True,
                       direction=os.path.join(datapath,'Titan_55437-56293dUTC.tab') )
                
            self.assertTrue(os.path.exists(myms+'/FIELD/EPHEM0_Titan.tab'))
            self.assertTrue(self.verify(myms, 'Titan', 'APP'))

    def test9(self):
        '''test9: Does a fixplanets with an ephemeris in mime format work'''
        os.system('cp '+ os.path.join(datapath,'titan.eml')+' .')
        for myms in [outms,outms2]:
            os.system("rm -rf titan.eml.tab")
            fixplanets( vis=myms, field='Titan', fixuvw=True, direction='titan.eml' )

            self.assertTrue(os.path.exists(myms+'/FIELD/EPHEM0_Titan.tab'))
            self.assertTrue(self.verify(myms, 'Titan', 'J2000'))

if __name__ == '__main__':
    unittest.main()
