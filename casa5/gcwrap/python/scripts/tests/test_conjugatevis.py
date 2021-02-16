#############################################################################
# $Id:$
# Test Name:                                                                #
#    Regression Test Script for the conjugatevis task
#
#                                                                           #
#############################################################################
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import shutil
import glob
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys, table, ms
    from casatasks import conjugatevis

    _tb = table( )
    _ms = ms( )
    datapath=ctsys.resolve('unittest/conjugatevis/')
else:
    from __main__ import default
    from tasks import *
    from taskinit import *
    _tb = tb
    _ms = ms
    datapath=os.environ.get('CASAPATH').split()[0]+'/casatestdata/unittest/conjugatevis/'

myname = 'test_conjugatevis'

# name of the resulting MS
msname = 'conjugated.ms'

# Pick up alternative data directory to run tests on MMSs
testmms = False
if 'TEST_DATADIR' in os.environ:
    testmms = True
    DATADIR = str(os.environ.get('TEST_DATADIR'))+'/concat/input/'
    if os.path.isdir(DATADIR):
        datapath = DATADIR
    print('conjugatevis tests will use data from %s' % datapath)


def checktable(thename, theexpectation):
    global msname, myname
    _tb.open(msname+"/"+thename)
    if thename == "":
        thename = "MAIN"
    for mycell in theexpectation:
        print("%s: comparing %s"% (myname,mycell))
        value = _tb.getcell(mycell[0], mycell[1])
        # see if value is array
        try:
            isarray = value.__len__
        except:
            # it's not an array
            # zero tolerance?
            if mycell[3] == 0:
                in_agreement = (value == mycell[2])
            else:
                in_agreement = ( abs(value - mycell[2]) < mycell[3])
        else:
            # it's an array
            # zero tolerance?
            if mycell[3] == 0:
                in_agreement =  (value == mycell[2]).all()
            else:
                try:
                    in_agreement = (abs(value - mycell[2]) < mycell[3]).all()
                except:
                    in_agreement = False
        if not in_agreement:
            print("%s:  Error in MS subtable %s:" % (myname,thename))
            print ("     column %s row %s contains %s" % (mycell[0],mycell[1],value))
            print("     expected value is %s" % mycell[2])
            _tb.close()
            return False
    _tb.close()
    print("%s: table %s as expected." %  (myname, thename))
    return True


###########################
# beginning of actual test

class test_conjugatevis(unittest.TestCase):

    def setUp(self):
        res = None

        cpath = os.path.abspath(os.curdir)
        filespresent = sorted(glob.glob("*.ms"))
        os.chdir(datapath)
        mymsname = 'shortpart1.ms'
        if not mymsname in filespresent:
            print("Copying %s" % mymsname)
            shutil.copytree(mymsname, cpath+'/'+mymsname)
        os.chdir(cpath)

        if not is_CASA6:
            default(conjugatevis)
        
    def tearDown(self):
        shutil.rmtree(msname,ignore_errors=True)

    def test1(self):
        '''Conjugatevis 1: '''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }

        self.res = conjugatevis(vis='shortpart1.ms', spwlist=[5,7], outputvis=msname)
        self.assertEqual(self.res,None)

        print("%s: Success! Now checking output ..." % myname)
        mscomponents = set(["table.dat",
                            "table.f0"
                            ])
        for name in mscomponents:
            if not os.access(msname+"/"+name, os.F_OK):
                print("%s: Error  %s/%s doesn't exist ..." % (myname,msname,name))
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+msname+'/'+name+' does not exist'
            else:
                print("%s: %s present." % (myname,name))
        print("%s: MS exists. Try opening as MS ..." %  myname)
        try:
            _ms.open(msname)
        except:
            print("%s: Error  Cannot open MS table %s" % (myname,tablename))
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']+'Cannot open MS table '+tablename
        else:
            _ms.close()
            print("%s: OK. Checking tables in detail ..." % myname)
            retValue['success']=True

            # check main table
            name = ""
            #             col name, row number, expected value, tolerance
            expected = [
                ['DATA',           4000,
                 [[-0.00426177-0.00387163j],
                  [ 0.00058119+0.00283016j]],
                 0.00000001]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'
            expected = [
                ['DATA',           3000,
                 [[ 0.00347826-0.00406267j],
                  [ 0.00458098-0.00508398j]],
                 0.00000001]
                ]
            results = checktable(name, expected)
            if not results:
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']+'Check of table '+name+' failed'


class conjugatevis_cleanup(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        os.system('rm -rf *.ms')

    def testrun(self):
        '''Conjugatevis: Cleanup'''
        pass

def suite():
    return [test_conjugatevis,conjugatevis_cleanup]

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()
