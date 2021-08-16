
# Trivial tests of asdmsummary.
# Tests that asdmsummary works without dying on a few ASDMs from 
# different telescopes.  Also checks that the number of new log rows written
# by asdmsummary matches the expected number for that ASDM.

from __future__ import absolute_import
import os
import sys

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys
    from casatasks import asdmsummary, casalog

else:
    from __main__ import default
    from tasks import asdmsummary
    from taskinit import casalog

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

def taskLogRange(startat):
    # return the start and end line numbers in the current logfile for "asdmsummary" task start and end after startat line
    # task start is a line with "Begin Task" and "asdmsummary" in it
    # task end is a line with "End Task" and "asdmsummary" in it
    # the first line matching each is returned
    # each value is -1 if no matching line is found
    # the end line isn't checked unless a start line has already been found
    firstLine = -1
    lastLine = -1
    logfile = casalog.logfile()
    if os.path.isfile(logfile):
        with open(logfile) as f:
            # lineNum will be the number of lines read, l is the most recent line read
            for lineNum, l in enumerate(f,1):
                # nothing happens until at startat line
                if lineNum >= startat:
                    if firstLine < 0:
                        if ((l.find('Begin Task') >= 0) and (l.find('asdmsummary')>=0)):
                            firstLine = lineNum
                    else:
                        if ((l.find('End Task') >= 0) and (l.find('asdmsummary')>=0)):
                            lastLine = lineNum
                            break
    return (firstLine,lastLine)

class asdmsummary_test(unittest.TestCase):

    # trivial tests that just demonstrate it doesn't fail completely
    # also now counts the number of lines written to the log file against expected count

    # CASA5 needs to know where the data is
    # CASA5 spits out 8 additional lines to the log that CASA6 does not produce
    dataPath = ""
    extraLines = 0
    if is_CASA6:
        dataPath = ctsys.resolve('unittest/asdmsummary/')
    else:
        dataPath = os.path.join(os.environ.get('CASAPATH').split()[0],'/casatestdata/unittest/asdmsummary/')
        extraLines = 8

    def doASDMSummary(self, asdmpath, expectedLogLines):
        # run asdmsummary, expepctedLogLines is the expected number of new log lines
        logLength = logfileLen()
        asdmsummary(os.path.join(self.dataPath,asdmpath))
        (firstLine,lastLine) = taskLogRange(logLength)
        newLines = 0
        if (firstLine >= 0) and (lastLine >= 0):
            newLines = lastLine-firstLine+1
        self.assertEqual(newLines,expectedLogLines+self.extraLines)

    def setUp(self):
        if is_CASA6:
            pass
        else:
            default(asdmsummary)

    def tearDown(self):
        pass

    def test_alma_asdm(self):
        ''' ALMA M51 data'''
        # used in test_importasdm, test_importasdm_mms.
        nlines = 170
        self.doASDMSummary('uid___X5f_X18951_X1', nlines)
        
    def test_vla_asdm(self):
        '''VLA data'''
        # used in test_importevla, test_importasdm_mms, test_importasdm
        nlines = 250
        self.doASDMSummary('X_osro_013.55979.93803716435',nlines)

    def test_aca_asdm(self):
        '''ACA with mixed pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        nlines = 2517
        self.doASDMSummary('uid___A002_X72bc38_X000',nlines)

    def test_12m_asdm(self):
        ''' 12m with mixedl pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        nlines = 1021
        self.doASDMSummary('uid___A002_X71e4ae_X317_short',nlines)

def suite():
    return [asdmsummary_test]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()
