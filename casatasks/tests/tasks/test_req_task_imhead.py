##########################################################################
# test_req_task_imhead.py
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
# https://casa.nrao.edu/casadocs/casa-5.4.1/global-task-list/task_imhead/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import casalog, imhead, rmtables
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import numpy as np

    #DATA#
if CASA6:
    datapath = casatools.ctsys.resolve('image/ngc5921.clean.image/')
else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/image/ngc5921.clean.image/'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/image/ngc5921.clean.image/'

logfile = casalog.logfile()
datacopy = 'clean.image'
testlog = 'testlog.log'

expectedKeys = [
    'beammajor', 'beamminor', 'beampa', 'bunit', 'cdelt1', 'cdelt2', 'cdelt3', 'cdelt4',
    'crpix1', 'crpix2', 'crpix3', 'crpix4', 'crval1', 'crval2', 'crval3', 'crval4',
    'ctype1', 'ctype2', 'ctype3', 'ctype4', 'cunit1', 'cunit2', 'cunit3', 'cunit4',
    'datamax', 'datamin', 'date-obs', 'equinox', 'imtype', 'masks', 'maxpixpos',
    'maxpos', 'minpixpos', 'minpos', 'object', 'observer', 'projection', 'reffreqtype',
    'restfreq', 'shape', 'telescope'
]

class imhead_test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass
 
    def setUp(self):
        if not CASA6:
            default(imhead)
        shutil.copytree(datapath, datacopy)
        os.chmod(datacopy, 493)
        for root, dirs, files in os.walk(datacopy):
            for d in dirs:
                os.chmod(os.path.join(root, d), 493)
            for f in files:
                os.chmod(os.path.join(root, f), 493)
 
    def tearDown(self):
        rmtables(datacopy)
        if os.path.exists(testlog):
            os.remove(testlog)
 
    @classmethod
    def tearDownClass(cls):
        pass
    
    
    def test_listkeys(self):
        '''
            test_listkeys
            ------------------------
            This test checks that the list mode displays all the expected keys from the imhead function
            A list is generated and the keys of that list are compared to what is expected.
            All items in the expected keys must be in the listed keys to pass.
            The length of the expected and recieved keys must also be the same.
        '''
        listed = imhead(datacopy, mode='list')
        self.assertTrue(all([True for item in expectedKeys if item in listed.keys()]))
        self.assertTrue(len(expectedKeys) == len(listed.keys()))
        
    def test_listkeysHdkeyVal(self):
        '''
            test_listkeysHdkeyVal
            ------------------------------
            This test checks to make sure that hdkey and hdvalue have no effect on the output of mode list
            A key and value is provided and then the keys are checked to make sure they are the expected values.
            The lists of expected values and recived ones must also be of equal length.
        '''
        listed = imhead(datacopy, mode='list', hdkey='beammajor', hdvalue=1)
        self.assertTrue(all([True for item in expectedKeys if item in listed.keys()]))
        self.assertTrue(len(expectedKeys) == len(listed.keys()))
    
    def test_history(self):
        '''
            test_history
            ---------------------
            This checks to make sure that when the mode is history a log is populated
            When the function is ran with mode history a logfile is populated with information
            Check that info is written to the file.
        '''
        casalog.setlogfile(testlog)
        imhead(datacopy, mode='history')
        # this is a history entry in the image
        msg = 'Tue Aug 7 20:48:18 2012  HISTORY imager::clean() [] Using Clark clean'
        self.assertTrue(msg in open(testlog).read())
        
    def test_historyHdkeyVal(self):
        '''
            test_historyHdkeyVal
            ------------------------
            This checks to make sure this task still functions when provided with Hdkey and value.
            These parameters should have no effect, so it should have the identical output to test_history.
        '''
        casalog.setlogfile(testlog)
        imhead(datacopy, mode='history', hdkey='beammajor', hdvalue=1)
        self.assertTrue('INFO' in open(testlog).read())
        
    def test_del(self):
        '''
            test_del
            ---------------
            This checks that the result of running the delete mode on all the hd keys has the
            effect that we expected
            Certain values should be removed where others will remain unchanged.
            The assertion checks that the image header after removal looks like the expected dictionary
        ''' 
        endDict = {
            'bunit': '', 'cdelt1': -7.27220521664304e-05, 'cdelt2': 7.27220521664304e-05, 'cdelt3': 1.0,
            'cdelt4': 24414.0625, 'crpix1': 128.0, 'crpix2': 128.0, 'crpix3': 0.0, 'crpix4': 0.0,
            'crval1': 4.022983925846928, 'crval2': 0.08843001543437938, 'crval3': 1.0,
            'crval4': 1412787144.0812755, 'ctype1': 'Right Ascension', 'ctype2': 'Declination',
            'ctype3': 'Stokes', 'ctype4': 'Frequency', 'cunit1': 'rad', 'cunit2': 'rad', 'cunit3': '',
            'cunit4': 'Hz', 'datamax': 0.054880399256944656, 'datamin': -0.011138656176626682,
            'date-obs': '1995/04/13/09:33:00', 'equinox': 'J2000', 'imtype': 'Intensity',
            'masks': np.array([], dtype='<U16'),
            'maxpixpos': np.array([134, 134,   0,  38], dtype='int32'),
            'maxpos': '15:21:53.976 +05.05.29.998 I 1.41371e+09Hz',
            'minpixpos': np.array([230,   0,   0,  15], dtype='int32'),
            'minpos': '15:20:17.679 +04.31.59.470 I 1.41315e+09Hz', 'object': '', 'observer': '',
            'projection': 'SIN', 'reffreqtype': 'LSRK', 'restfreq': np.array([1420405752.0]),
            'shape': np.array([256, 256,   1,  46], dtype='int32'), 'telescope': ''
        }
        for key in expectedKeys:
            imhead(datacopy, mode='del', hdkey=key)
        self.assertTrue(len(endDict.keys()) == len(imhead(datacopy, mode='list').keys()))
        self.assertTrue(
            np.all([
                np.all(imhead(datacopy, mode='list')[key]==endDict[key]) 
                for key in imhead(datacopy, mode='list').keys()
            ])
        )
        
    def test_add(self):
        '''
            test_add
            ---------------
            This test makes sure that add can add to all the keys specified in the documentation
            The endDict is the dictionary we expect out of the list function.
            we first delete all the values we can, and then try to add them back.
            Upon failing to add a section back the dictionaries will become differnt, so my
            asserting the dictionary is the same as it started assures add worked.
            This relies on del working as well, however if that fails so will the test_del function.
            
            imtype and restfreq cannot be tested, because there is no way to remove the values for them.
            Also adding a beam auto matically sets beamminor and beampa, those values cannot be all chosen using add
        '''
        # There are fields you can add to but can't delete from?
        # How are you supposed to use add when add requires th field to be empty
    
        endDict = {
            'beammajor': {'unit': 'arcsec', 'value': 51.7},
            'beamminor': {'unit': 'arcsec', 'value': 51.7},
            'beampa': {'unit': 'deg', 'value': 0.0},
            'bunit': 'Jy/beam', 'cdelt1': -7.27220521664304e-05, 'cdelt2': 7.27220521664304e-05,
            'cdelt3': 1.0, 'cdelt4': 24414.0625, 'crpix1': 128.0, 'crpix2': 128.0,
            'crpix3': 0.0, 'crpix4': 0.0, 'crval1': 4.022983925846928,
            'crval2': 0.08843001543437938, 'crval3': 1.0, 'crval4': 1412787144.0812755,
            'ctype1': 'Right Ascension', 'ctype2': 'Declination', 'ctype3': 'Stokes',
            'ctype4': 'Frequency', 'cunit1': 'rad', 'cunit2': 'rad', 'cunit3': '',
            'cunit4': 'Hz', 'datamax': 0.054880399256944656, 'datamin': -0.011138656176626682,
            'date-obs': '1995/04/13/09:33:00', 'equinox': 'J2000', 'imtype': 'Intensity',
            'masks': np.array([], dtype='<U16'),
            'maxpixpos': np.array([134, 134,   0,  38], dtype='int32'),
            'maxpos': '15:21:53.976 +05.05.29.998 I 1.41371e+09Hz',
            'minpixpos': np.array([230,   0,   0,  15], dtype='int32'),
            'minpos': '15:20:17.679 +04.31.59.470 I 1.41315e+09Hz', 'object': 'N5921_2',
            'observer': 'TEST', 'projection': 'SIN', 'reffreqtype': 'LSRK',
            'restfreq': np.array([1420405752.0]),
            'shape': np.array([256, 256,   1,  46], dtype='int32'), 'telescope': 'VLA'
        }
        for key in expectedKeys:
            imhead(datacopy, mode='del', hdkey=key)
            
        for key in expectedKeys:
            imhead(datacopy, mode='add', hdkey=key, hdvalue=endDict[key])
        
        self.assertTrue(
            len(endDict.keys()) == len(imhead(datacopy, mode='list').keys())
        )
        self.assertTrue(
            np.all([
                np.all(imhead(datacopy, mode='list')[key]==endDict[key])
                for key in imhead(datacopy, mode='list').keys()
            ])
        )
        
    def test_put(self):
        '''
            test_put
            ----------------
            This test checks to see that the put mode will replace the specified key with the specified
            value. Not all keys are replacable, and the test checks for that. The put values are all different
            from the initial values to enable more robust testing, except for non-mutable keys in which case
            they are the same values
        '''
        # axis 3 is the stokes axis and those values cannot be changed because its a Stokes axis
        immutable = [
            'maxpixpos', 'minpos', 'minpixpos', 'maxpos', 'masks', 'datamax',
            'shape', 'datamin', 'cdelt3', 'crpix3', 'cunit3'
        ] 
        myDict = {
            'crpix4': 1.0, 'equinox': 'B1950', 'crpix1': 129.0, 'crpix2': 129.0, 'crpix3': 0.0,
            'maxpixpos': np.array([134, 134,   0,  38], dtype=np.int32),
            # minpos will change not because we change its value directly but because we muck with
            # axes metadata
            'minpos': '00:11:59.709 +00.19.25.824 U 1.41414e+09kHz',
            'minpixpos': np.array([230,   0,   0,  15], dtype=np.int32),
            # you can put the Stokes axis name, but it remains Stokes
            'ctype4': 'Freq', 'ctype3': 'Stokes', 'ctype2': 'Dec',
            'ctype1': 'RA', 'beamminor': {'value': 51.5, 'unit': 'arcsec'},
            # maxpos will change not because we change its value directly but because we muck with
            # axes metadata
            'maxpos': '00:12:00.683 +00.19.46.197 U 1.41473e+09kHz',
            'projection': 'TAN', 'observer': 'FRED', 'telescope': 'ALMA',
            'date-obs': '1995/04/15/09:33:00', 'restfreq': {'value': 1520405752.0, 'unit': 'kHz'},
            'imtype': 'Optical Depth',
            'crval4': {'value': 1413787144.0812755, 'unit': 'kHz'},
            'crval2': {'value': 0.08943001543437938, 'unit': 'crad'},
            'crval3': np.array(['U']),
            'crval1': {'value': 4.122983925846928, 'unit': 'crad'},
            'beampa': {'value': 11.0, 'unit': 'deg'}, 'object': 'N5921_2_bogus',
            'masks': np.array([], dtype='|S1'),
            'cunit4': 'kHz', 'cunit1': 'crad', 'datamax': 0.054880399256944656,
            'cunit3': '', 'cunit2': 'crad',
            'beammajor': {'value': 55.7, 'unit': 'arcsec'}, 'reffreqtype': 'BARY',
            'cdelt1': {'value': -7.372e-05, 'unit': 'crad'},
            'cdelt2': {'value': 7.3722e-05, 'unit': 'crad'},
            'cdelt3': {'value': 1.0, 'unit': ''}, 'cdelt4': {'value': 25414.0625, 'unit': 'kHz'},
            'shape': np.array([256, 256,   1,  46], dtype=np.int32),
            'bunit': 'mJy/beam', 'datamin': -0.011138656176626682
        }        
        for k in myDict.keys():
            if k in immutable:
                self.assertFalse(
                    imhead(
                        datacopy, mode='put', hdkey=k, hdvalue=myDict[k]
                    ), "Incorrectly put hdkey " + k
                )
            else:
                self.assertTrue(
                    imhead(datacopy, mode='put', hdkey=k, hdvalue=myDict[k]),
                    "Did not put hdkey " + k
                )
        # do after first loop because we are not gauranteed the order things
        # are put, so everything must be put before we can check values
        for k in myDict.keys():
            # check the value is correct
            x = imhead(datacopy, mode='get', hdkey=k)
            msg = "Failed get for hdkey " + k + " got " + str(x) + " expected " + str(myDict[k])
            if type(x) == np.ndarray:
                try:
                    self.assertEqual(np.all(x, myDict[k]), msg)
                except:
                    # string arrays can throw "Cannot reduce flexible type" exceptions
                    self.assertEqual(list(x), list(myDict[k]), msg)
            else:
                if k == 'cdelt1':
                    # because its annoying
                    self.assertTrue(
                        np.isclose(x['value'], myDict[k]['value']),
                        "got wrong value for key " + k + " expected "
                        + str(myDict[k]) + " got " + str(x)
                    )
                    self.assertEqual(x['unit'], myDict[k]['unit'], 'wrong unit cdelt1')
                else:
                    self.assertEqual(x, myDict[k], msg)
       
    def test_get(self):
        '''
            test_get
            --------------
            Test test_get to verify mode='get' works and returns correct values
        '''
        exp_dict = {
            'crpix4': 0.0, 'equinox': 'J2000', 'crpix1': 128.0, 'crpix2': 128.0, 'crpix3': 0.0,
            'maxpixpos': np.array([134, 134,   0,  38], dtype=np.int32),
            'minpos': '15:20:17.679 +04.31.59.470 I 1.41315e+09Hz',
            'minpixpos': np.array([230,   0,   0,  15], dtype=np.int32),
            'ctype4': 'Frequency', 'ctype3': 'Stokes', 'ctype2': 'Declination',
            'ctype1': 'Right Ascension', 'beamminor': {'value': 47.2, 'unit': 'arcsec'},
            'maxpos': '15:21:53.976 +05.05.29.998 I 1.41371e+09Hz', 'projection': 'SIN',
            'observer': 'TEST', 'telescope': 'VLA', 'date-obs': '1995/04/13/09:33:00',
            'restfreq': {'value': 1420405752.0, 'unit': 'Hz'}, 'imtype': 'Intensity',
            'crval4': {'value': 1412787144.0812755, 'unit': 'Hz'},
            'crval2': {'value': 0.08843001543437938, 'unit': 'rad'},
            'crval3': np.array(['I']),
            'crval1': {'value': 4.022983925846928, 'unit': 'rad'},
            'beampa': {'value': 9.0, 'unit': 'deg'}, 'object': 'N5921_2',
            'masks': np.array([], dtype='|S1'), 'cunit4': 'Hz', 'cunit1': 'rad',
            'datamax': 0.054880399256944656, 'cunit3': '', 'cunit2': 'rad',
            'beammajor': {'value': 51.7, 'unit': 'arcsec'}, 'reffreqtype': 'LSRK',
            'cdelt1': {'value': -7.27220521664304e-05, 'unit': 'rad'},
            'cdelt2': {'value': 7.27220521664304e-05, 'unit': 'rad'},
            'cdelt3': {'value': 1.0, 'unit': ''},
            'cdelt4': {'value': 24414.0625, 'unit': 'Hz'},
            'shape': np.array([256, 256,   1,  46], dtype=np.int32), 'bunit': 'Jy/beam',
            'datamin': -0.011138656176626682
        }
        for k in expectedKeys:
            x = imhead(datacopy, mode='get', hdkey=k)
            print(k, x, exp_dict[k])
            msg = "Failed get for hdkey " + k + " got " + str(x) + " expected " + str(exp_dict[k])
            if type(x) == np.ndarray:
                try:
                    self.assertEqual(np.all(x, exp_dict[k]), msg)
                except:
                    # string arrays can throw "Cannot reduce flexible type" exceptions
                    self.assertEqual(list(x), list(exp_dict[k]), msg)
            else:
                self.assertEqual(x, exp_dict[k], msg)
        # test bogus key fails
        k = 'boguskey'
        self.assertFalse(imhead(datacopy, mode='get', hdkey=k), 'Incorrectly found bogus key')

    def test_summary(self):
        '''
            test_summary
            -----------------
        '''
        ###NOTE: to test the function of the verbose parameter an image with multiple beams is required
        summaryDictnoV = {
            'axisnames': np.array(['Right Ascension', 'Declination', 'Stokes', 'Frequency'], dtype='<U16'),
            'axisunits': np.array(['rad', 'rad', '', 'Hz'], dtype='<U16'), 'defaultmask': '',
            'hasmask': False, 'imagetype': 'Intensity',
            'incr': np.array([
                -7.27220521664304e-05,  7.27220521664304e-05,  1.00000000e+00,  24414.0625
            ]), 'masks': np.array([], dtype='<U16'), 'messages': np.array([], dtype='<U16'),
            'ndim': 4, 'refpix': np.array([128., 128.,   0.,   0.]),
            'refval': np.array([
                4.022983925846928, 0.08843001543437938, 1.00000000e+00, 1412787144.0812755
            ]), 'restoringbeam': {
                'major': {'unit': 'arcsec', 'value': 51.7}, 'minor': {'unit': 'arcsec', 'value': 47.2},
                'positionangle': {'unit': 'deg', 'value': -171.0}
            }, 'shape': np.array([256, 256,   1,  46], dtype='int32'),
            'tileshape': np.array([64, 64,  1,  8], dtype='int32'), 'unit': 'Jy/beam'
        }
        compared = imhead(datacopy, mode='summary')
        self.assertTrue(len(compared.keys()) == len(summaryDictnoV.keys()))
        self.assertTrue(np.all([np.all(compared[key]==summaryDictnoV[key]) for key in compared.keys()]))

    def test_bad_mode(self):
        """Test unupported mode fails"""
        if CASA6:
            self.assertRaises(Exception, imhead, datacopy, mode='bogus')
        else:
            self.assertFalse(
                imhead(datacopy, mode='bogus'),
                'application did not fail as expected for bad mode test'
            )
        
def suite():
    return[imhead_test]
 
# Main #
if __name__ == '__main__':
    unittest.main()
