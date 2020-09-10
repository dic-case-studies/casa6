##########################################################################
# test_req_task_accor.py
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
# https://casa.nrao.edu/casadocs-devel/stable/global-task-list/task_accor/about
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import accor, rmtables
    tb = casatools.table()
    CASA6 = True
except ImportError:
    from __main__ import default
    from tasks import *
    from taskinit import *
import sys
import os
import unittest
import shutil
import filecmp
import numpy

### DATA ###

if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/alma/uid___X02_X3d737_X1_01_small.ms/')
    libpath = casatools.ctsys.resolve('text/testcallib.txt')
    vladata = casatools.ctsys.resolve('visibilities/vla/ngc5921.ms/')
    cdfdata = casatools.ctsys.resolve('regression/evn/n08c1.ms/')

else:
    if os.path.exists(os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req'):
        datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
        libpath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/text/testcallib.txt'
        vladata = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/ngc5921.ms/'
    else:
        datapath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
        libpath = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/text/testcallib.txt'
        vladata = os.environ.get('CASAPATH').split()[0] + '/casa-data-req/visibilities/vla/ngc5921.ms/'
    cdfdata = os.environ.get('CASAPATH').split()[0] + '/data/regression/evn/n08c1.ms/'


caltab = 'cal.A'
cal_default = 'cal.default'
datacopy = 'uid_copy.ms'
vlacopy = 'vla_copy.ms'
cdfcopy = 'cdf_copy.ms'

def cal_size(cal):
    '''
        This returns the size of a dir
    '''
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(cal):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size

def getmean(data):
    
    tb.open(data)
    datamean = numpy.mean(tb.getcol('CPARAM'))
    tb.close()
    return datamean
    
    
class accor_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        shutil.copytree(datapath, datacopy)
        shutil.copytree(vladata, vlacopy)
    
    def setUp(self):
        if not CASA6:
            default(accor)
            
    def tearDown(self):
        rmtables(caltab)
        if os.path.exists('cal.B'):
            rmtables('cal.B')
        if os.path.exists(cdfcopy):
            shutil.rmtree(cdfcopy)
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(datacopy)
        shutil.rmtree(vlacopy)
        rmtables(cal_default)
    
    
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_makesTable(self):
        ''' Test that when accor is run it creates a calibration table '''
        
        accor(vis=datacopy, caltable=caltab)
        
        self.assertTrue(os.path.exists(caltab))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_fieldSelect(self):
        ''' Test that a calibration table generated with a field selection is different than one generated with no selection parameters '''
        
        accor(vis=datacopy, caltable=caltab, field='1')
        tb.open(caltab)
        fields = tb.getcol('FIELD_ID')
        tb.close()
        
        self.assertTrue(numpy.all(fields == 1))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_spwSelect(self):
        ''' Test that a calibration table generated with a spectral window selection is different than one generated with no selection parameters '''
        
        accor(vis=datacopy, caltable=caltab, spw='1')
        tb.open(caltab)
        spws = tb.getcol('SPECTRAL_WINDOW_ID')
        tb.close()
        
        self.assertTrue(numpy.all(spws == 1))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_intentSelect(self):
        ''' Test that a calibration table generated with an intent selection is different than one generated with no selection parameters '''
        
        accor(vis=datacopy, caltable=caltab, intent='*AMPLI*')
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean,1.1453650693098703+0j))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_selectData(self):
        ''' Test when selectdata = True selections may be used, and while selectdata = False they may not be used '''
        
        accor(vis=datacopy, caltable=caltab, selectdata=False, scan='2')
        tb.open(caltab)
        data1 = tb.getcol('CPARAM')
        tb.close()
        rmtables(caltab)
        
        accor(vis=datacopy, caltable=caltab, selectdata=True, scan='2')
        tb.open(caltab)
        data2 = tb.getcol('CPARAM')
        tb.close()
        rmtables(caltab)
        
        self.assertFalse(numpy.shape(data1) == numpy.shape(data2))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_timeRangeSelect(self):
        ''' Test that a calibration table generated with a timerange selection is different than one generated with no selection parameters '''
        
        accor(vis=datacopy, caltable=caltab, timerange='03:01:32.7~03:04:59.5')
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 1.1656893293062844+0j))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_antennaSelect(self):
        ''' Test that a calibration table generated with an antenna selection is different than one generated with no selection parameters '''
        
        accor(vis=datacopy, caltable=caltab, antenna='1&&')
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 1.075232790576087+0j))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_scanSelect(self):
        ''' Test that a calibration table generated with a scan selection is different than one generated with no selection parameters '''
        
        accor(vis=datacopy, caltable=caltab, scan='2')
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 1.1453650693098703+0j))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_obsSelect(self):
        ''' Test that a calibration table generated with an observation selection is different than one generated with no seletion parameters '''
        
        accor(vis=datacopy, caltable=caltab, observation='0')
        self.assertTrue(os.path.exists(caltab))
        
    #CAS-12736   
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_solint(self):
        ''' Test that the solint parameter changes he solution interval (?) '''
        
        accor(vis=datacopy, caltable=caltab, solint='10s')
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 1.1667321394651364+0j))
        
    #CAS-12736   
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_combineSelect(self):
        ''' Test that a calibration table generated with a combine selection is different than one generated with no selection parameteres '''
        
        accor(vis=datacopy, caltable=caltab, combine='scan, spw')
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 1.136721501747767+0j))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_append(self):
        ''' Generates a cal table and then attempts to append to it '''
        
        accor(datacopy, caltab)
        before_append = cal_size(caltab)
        
        accor(datacopy, caltab, append=True)
        after_append = cal_size(caltab)
        
        self.assertTrue(after_append > before_append)
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_docallib(self):
        ''' Test that the do callib parameter allows for the selection of caltables '''

        accor(datacopy, caltable='cal.B')
        accor(vis=datacopy, caltable=caltab, docallib=True, callib=libpath)
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 1.0017881503811588+0j))
        
        
    #CAS-12736   
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_gaintable(self):
        ''' Test that providing the gaintable will yeild a different final cal table than the default '''
        
        accor(datacopy, caltable='cal.B')
        accor(vis=datacopy, caltable=caltab, gaintable='cal.B')
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 1.0017881503811588+0j))
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_gainfield(self):
        ''' Test that adding a field selection to the gaintable will yeild a different cal table than gaintable with no field selection '''
        
        accor(datacopy, caltable='cal.B')
        accor(vis=datacopy, caltable=caltab, gaintable='cal.B', gainfield='1')
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 0.9921940355389206+0j))

    #CAS-12736   
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_interp(self):
        ''' Test that adding an interp selection to the gaintable will yeild a different cal table than gaintable with standard interp (linear, linear) '''
        accor(vlacopy, caltable='cal.B')
        
        accor(vlacopy, caltable='cal.A', gaintable=['cal.B'])
        accor(vlacopy, caltable='cal.C', gaintable=['cal.B'],  interp='nearest')
        
        print(len(filecmp.dircmp('cal.A', 'cal.C').diff_files),(filecmp.dircmp('cal.A', 'cal.C').diff_files))
        
        self.assertTrue(len(filecmp.dircmp('cal.A', 'cal.C').diff_files) > 1)
        
        rmtables('cal.C')
        
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_spwmap(self):
        ''' Test that adding a spwmap selection to the gaintable will yeild a different cal table than gaintable with no spwmap '''
        
        accor(datacopy, caltable='cal.B')
        accor(vis=datacopy, caltable=caltab, gaintable='cal.B', spwmap=[0,0])
        datamean = getmean(caltab)
        
        self.assertTrue(numpy.isclose(datamean, 1.282986655279442+0j))
        
    #CAS-13184
    @unittest.skipIf(sys.platform == "darwin", "Disabled for OSX")
    def test_corrdepflags(self):
        ''' Test that adding corrdepflags=True finds more solutions '''

        shutil.copytree(cdfdata, cdfcopy)

        # Cross autocorrelations are flagged; no solutions should be found
        accor(vis=cdfcopy, caltable=caltab)
        self.assertFalse(os.path.exists(caltab))

        # Cross autocorrelations are flagged; solutions should be found
        accor(vis=cdfcopy, caltable=caltab, corrdepflags=True)
        self.assertTrue(os.path.exists(caltab))
        
                
def suite():
    return[accor_test]

if __name__ == '__main__':
    unittest.main()
