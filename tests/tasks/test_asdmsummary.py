
# Trivial tests of asdmsummary.
# Tests that asdmsummary works without dying on a few ASDMs from 
# different telescopes.  Also checks that the number of new log rows written
# by asdmsummary matches the expected number for that ASDM.


import os
from casatasks import asdmsummary, casalog
import unittest

def logfileLen():
    # count the lines in the current log file
    result = 0
    logfile = casalog.logfile()
    if os.path.isfile(logfile):
        with open(logfile) as f:
            for result, l in enumerate(f,1):
                pass
    return result

class asdmsummary_test(unittest.TestCase):

    def doASDMSummary(self, asdmpath, expectedLogLines):
        # run asddumsummary, expepctedLogLines is the expected number of new log lines
        logLength = logfileLen()
        asdmsummary(asdmpath)
        newLines = logfileLen()-logLength
        self.assertEqual(newLines,expectedLogLines)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_alma_asdm(self):
        ''' ALMA M51 data'''
        # used in test_importasdm, test_importasdm_mms.
        self.doASDMSummary('regression/asdm-import/input/uid___X5f_X18951_X1',166)
        
    def test_vla_asdm(self):
        '''VLA data'''
        # used in test_importevla, test_importasdm_mms, test_importasdm
        self.doASDMSummary('regression/unittest/importevla/X_osro_013.55979.93803716435',246)

    def test_aca_asdm(self):
        '''ACA with mixed pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        self.doASDMSummary('regression/asdm-import/input/uid___A002_X72bc38_X000',2513)

    def test_12m_asdm(self):
        ''' 12m with mixedl pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        self.doASDMSummary('regression/asdm-import/input/uid___A002_X71e4ae_X317_short',1017)

def suite():
    return [asdmsummary_test]

if __name__ == '__main__':
    unittest.main()
