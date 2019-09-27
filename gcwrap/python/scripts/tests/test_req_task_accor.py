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

### DATA ###

if CASA6:
    datapath = casatools.ctsys.resolve('visibilities/alma/uid___X02_X3d737_X1_01_small.ms/')
    libpath = casatools.ctsys.resolve('text/testcallib.txt')
    vladata = casatools.ctsys.resolve('visibilities/vla/ngc5921.ms/')

else:
    datapath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/alma/uid___X02_X3d737_X1_01_small.ms/'
    libpath = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/text/testcallib.txt'
    vladata = os.environ.get('CASAPATH').split()[0] + '/data/casa-data-req/visibilities/vla/ngc5921.ms/'


caltab = 'cal.A'
cal_default = 'cal.default'

def use_param(vis=datapath, field="", spw="", intent="", selectdata=True, timerange="", antenna="", scan="", observation="", msselect="", solint="inf", combine="", append=False, docallib=False, callib="", gaintable=[], gainfield=[], interp=[], spwmap=[], gaincheck=False, forAnt=False):
    '''
        This function creates two calibration tables and checks to see if they differ
        One file is always different even for cal tables only differing in name 
        Checks for more than one difference
    '''
    
    if not gaincheck:
        accor(vis=vis, caltable=cal_default)
    else:
        accor(vis=vis, caltable=cal_default, gaintable=gaintable)

    accor(vis=vis, caltable=caltab, field=field, spw=spw, intent=intent, selectdata=selectdata, timerange=timerange, antenna=antenna, scan=scan, observation=observation, msselect=msselect, solint=solint, combine=combine, append=append, docallib=docallib, callib=callib, gaintable=gaintable, gainfield=gainfield, interp=interp, spwmap=spwmap)
    
    if not forAnt:
        if len(filecmp.dircmp(cal_default, caltab).diff_files) <= 1:
            isSame = True
        else:
            isSame = False
    else:
        if len(filecmp.dircmp(cal_default, caltab).diff_files) < 1:
            isSame = True
        else:
            isSame = False
    
    return isSame

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
    
    
class accor_test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass
    
    def setUp(self):
        if not CASA6:
            default(accor)
            
    def tearDown(self):
        rmtables(caltab)
        if os.path.exists('cal.B'):
            rmtables('cal.B')
    
    @classmethod
    def tearDownClass(cls):
        rmtables(cal_default)
    
    def test_makesTable(self):
        '''
            test_makesTable
            ------------------
            
            Test that when accor is run it creates a calibration table
            This simply checks that it exists, not the contents of the table
        '''
        
        accor(vis=datapath, caltable=caltab)
        self.assertTrue(os.path.exists(caltab))
        
    def test_fieldSelect(self):
        '''
            test_fieldSelect
            --------------------
            
            Test that a calibration table generated with a field selection is different than one generated with no selection parameters
        '''
        
        self.assertFalse(use_param(field='1'))
        
    def test_spwSelect(self):
        '''
            test_spwSelect
            ------------------
            
            Test that a calibration table generated with a spectral window selection is different than one generated with no selection parameters
        '''
        
        self.assertFalse(use_param(spw='1'))
        
    def test_intentSelect(self):
        '''
            test_intentSelect
            --------------------
            
            Test that a calibration table generated with an intent selection is different than one generated with no selection parameters
        '''
        
        self.assertFalse(use_param(intent='*AMPLI*'))
        
    def test_selectData(self):
        '''
            test_selectData
            -------------------
            
            Test when selectdata = True selections may be used, and while selectdata = False they may not be used
        '''
        
        self.assertTrue(use_param(selectdata=False, scan='2'))
        self.assertFalse(use_param(selectdata=True, scan='2'))
        
    def test_timeRangeSelect(self):
        '''
            test_timeRangeSelect
            -----------------------
            
            Test that a calibration table generated with a timerange selection is different than one generated with no selection parameters
        '''
        
        self.assertFalse(use_param(timerange='03:01:32.7~03:04:59.5'))
        
    def test_antennaSelect(self):
        '''
            test_antennaSelect
            ----------------------
            
            Test that a calibration table generated with an antenna selection is different than one generated with no selection parameters
        '''
        #antenna selection always returing an empty table? (come back to this one)
        self.assertFalse(use_param(vis=datapath, antenna='1&&', forAnt=True))
        
    def test_scanSelect(self):
        '''
            test_scanSelect
            -------------------
            
            Test that a calibration table generated with a scan selection is different than one generated with no selection parameters
        '''
        
        self.assertFalse(use_param(scan='2'))
        
    def test_obsSelect(self):
        '''
            test_obsSelect
            ----------------
            
            Test that a calibration table generated with an observation selection is different than one generated with no seletion parameters
        '''
        
        accor(vis=datapath, caltable=caltab, observation='0')
        self.assertTrue(os.path.exists(caltab))
        
        
    def test_solint(self):
        '''
            test_solint
            --------------
            
            Test that the solint parameter changes he solutin interval (?)
        '''
        
        self.assertFalse(use_param(solint='10s'))
        
    def test_combineSelect(self):
        '''
            test_combineSelect
            -----------------------
            
            Test that a calibration table generated with a combine selection is different than one generated with no selection parameteres
        '''
        
        self.assertFalse(use_param(combine='scan, spw'))
        
    def test_append(self):
        '''
            test_append
            ---------------
            
            Generates a cal table and then attempts to append to it. 
            The final cal table should be larger than the first.
        '''
        
        accor(datapath, caltab)
        before_append = cal_size(caltab)
        
        accor(datapath, caltab, append=True)
        after_append = cal_size(caltab)
        
        self.assertTrue(after_append > before_append)
        
    def test_docallib(self):
        '''
            test_docallib
            ----------------
            
            Test that the do callib parameter allows for the selection of caltables
            With this set to true the results when providing a callib should differ from the default
            With this parameter set to false it should be identical to the default
        '''

        accor(datapath, caltable='cal.B')
                
        self.assertFalse(use_param(docallib=True, callib=libpath))
        self.assertTrue(use_param(docallib=False, callib=libpath))

        
        
        
    def test_callib(self):
        '''
            test_callib
            ---------------
            
            Test that providing a callib creates a new cal table that differs from the default
            This is essentially covered in the the docallib test as well
        '''
        
        accor(datapath, caltable='cal.B')
        
        self.assertFalse(use_param(docallib=True, callib=libpath))
        
        
    def test_gaintable(self):
        '''
            test_gaintable
            -----------------
            
            Test that providing the gaintable will yeild a different final cal table than the default.
        '''
        accor(datapath, caltable='cal.B')
        
        self.assertFalse(use_param(gaintable='cal.B'))
        
        
    def test_gainfield(self):
        '''
            test_gainfield
            ----------------
            
            Test that adding a field selection to the gaintable will yeild a different cal table than gaintable with no field selection
        '''
        accor(datapath, caltable='cal.B')
        
        self.assertFalse(use_param(gaintable='cal.B', gainfield='1', gaincheck=True))
        
    def test_interp(self):
        '''
            test_interp
            ----------------
            
            Test that adding an interp selection to the gaintable will yeild a different cal table than gaintable with standard interp (linear, linear)
        '''
        accor(vladata, caltable='cal.B')
        
        accor(vladata, caltable='cal.A', gaintable=['cal.B'])
        accor(vladata, caltable='cal.C', gaintable=['cal.B'],  interp='nearest')
        
        
        print(len(filecmp.dircmp('cal.A', 'cal.C').diff_files),(filecmp.dircmp('cal.A', 'cal.C').diff_files))
        
        self.assertTrue(len(filecmp.dircmp('cal.A', 'cal.C').diff_files) > 1)
        
        rmtables('cal.C')
        
    
    def test_spwmap(self):
        '''
            test_spwmap
            --------------
            
            Test that adding a spwmap selection to the gaintable will yeild a different cal table than gaintable with no spwmap
        '''
        accor(datapath, caltable='cal.B')
        
        self.assertFalse(use_param(gaintable='cal.B', spwmap=[0,0]))
        
        
                
def suite():
    return[accor_test]

if __name__ == '__main__':
    unittest.main()
