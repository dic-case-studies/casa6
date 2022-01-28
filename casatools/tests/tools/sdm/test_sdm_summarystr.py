
# Trivial tests of asdmsummary.
# Tests that asdmsummary works without dying on a few ASDMs from 
# different telescopes.  Also checks that the number of new log rows written
# by asdmsummary matches the expected number for that ASDM.


import os
from casatools import ctsys, sdm
import unittest

datapath = ctsys.resolve('unittest/sdmtool/')
class sdm_summarystr_test(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_alma_asdm(self):
        ''' ALMA M51 data'''
        # used in test_importasdm, test_importasdm_mms.
        mysdm = sdm(datapath + 'uid___X5f_X18951_X1')
        #self.doASDMSummary(asdmpath,174)
        summary = mysdm.summarystr( ).splitlines( )
        self.assertTrue(len(summary) == 166, summary[0] if len(summary) == 1 else None)
        
    def test_vla_asdm(self):
        '''VLA data'''
        # used in test_importevla, test_importasdm_mms, test_importasdm
        mysdm = sdm(datapath + 'X_osro_013.55979.93803716435')
        summary = mysdm.summarystr( ).splitlines( )
        self.assertTrue(len(summary) == 246, summary[0] if len(summary) == 1 else None)

    def test_aca_asdm(self):
        '''ACA with mixed pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        mysdm = sdm(datapath + 'uid___A002_X72bc38_X000')
        summary = mysdm.summarystr( ).splitlines( )
        self.assertTrue(len(summary) == 2513, summary[0] if len(summary) == 1 else None)

    def test_12m_asdm(self):
        ''' 12m with mixedl pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        mysdm = sdm(datapath + 'uid___A002_X71e4ae_X317_short')
        summary = mysdm.summarystr( ).splitlines( )
        self.assertTrue(len(summary) == 1017, summary[0] if len(summary) == 1 else None)

    def test_bogus_file(self):
        ''' 12m with mixedl pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        passes = False
        try:
            mysdm = sdm('/tmp/some_nonexistent_directory_for_our_test')
            summary = mysdm.summarystr( )
        except:
            passes = True
        self.assertTrue(passes,"non-existent file fails to throw an exception")

if __name__ == '__main__':
    unittest.main()
