#########################################################################
# test_tool_sdm.py
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
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.calibrater.html
#
# Testing of methods frommms and summarystr
#
##########################################################################
import os
import sys
import shutil
import unittest

from casatools import ctsys, sdm, ms
ctsys_resolve = ctsys.resolve

from casatestutils import testhelper as th

datapath = ctsys_resolve('unittest/sdmtool/')

_ms = ms( )

# Class of tests for sdm.summarystr
class exportasdm_test(unittest.TestCase):
    
    vis_b = 'ngc4826_bima_7fields_7spw.ms'
    vis_c = 'M100-X220-shortened.ms'
    vis_d = 'ngc4826.tutorial.ngc4826.ll.5.ms'
    vis_e = 'g19_d2usb_targets_line-shortened.ms'
    vis_f = 'Itziar.ms'
    vis_g = 'M51.ms'
    vis_h = 'xosro2ref.ms'
    vis_i = 'asdm.ms'
    out = 'exportasdm-output.asdm'
    rval = False
    
    def setUp(self):    
        self.rval = False
        if(not os.path.exists(self.vis_b)):
            # renamed test.ms to ngc4826_bima_7fields_7spw.ms
            os.system('cp -RH '+datapath+'ngc4826_bima_7fields_7spw.ms'+' .')
        if(not os.path.exists(self.vis_c)):
            os.system('cp -RH '+datapath+'M100-X220-shortened.ms'+' .')
        if(not os.path.exists(self.vis_d)):
            _ms.fromfits( self.vis_d, datapath+'/ngc4826.ll.fits5' )
            _ms.close( )
        if(not os.path.exists(self.vis_e)):
            os.system('cp -RH '+datapath+'g19_d2usb_targets_line-shortened.ms'+' .')
        if(not os.path.exists(self.vis_f)):
            os.system('cp -RH '+datapath+'Itziar.ms'+' .')
        if(not os.path.exists(self.vis_g)):
            os.system('cp -RH '+datapath+'M51.ms'+' .')
        if(not os.path.exists(self.vis_h)):
            os.system('ln -sf '+datapath+'X_osro_013.55979.93803716435'+' .')
            mysdm = sdm('X_osro_013.55979.93803716435')
            mysdm.toms('xosro2ref.ms',process_flags=False,scans='0:2',ocorr_mode='co',with_pointing_correction=True)
        if(not os.path.exists(self.vis_i)):
            os.system('ln -sf '+datapath+'uid___A002_X72bc38_X000'+' .')
            mysdm = sdm('uid___A002_X72bc38_X000')
            mysdm.toms('asdm.ms', scans='0:2')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.vis_c,ignore_errors=True)
        shutil.rmtree(cls.vis_b,ignore_errors=True)
        shutil.rmtree(cls.vis_d,ignore_errors=True)
        shutil.rmtree(cls.vis_e,ignore_errors=True)
        shutil.rmtree(cls.vis_f,ignore_errors=True)
        shutil.rmtree(cls.vis_g,ignore_errors=True)
        shutil.rmtree(cls.vis_h,ignore_errors=True)
        shutil.rmtree(cls.vis_i,ignore_errors=True)
        os.system('rm -rf test*exportasdm*.asdm')

    def tearDown(self):
        os.system('rm -rf myinput.ms')
        os.system('rm -rf ' + self.out)

    def verify_asdm(self,asdmname, withPointing):
        print("Verifying asdm %s" % asdmname)
        if(not os.path.exists(asdmname)):
            print("asdm %s doesn't exist." % asdmname)
            raise Exception
        # test for the existence of all obligatory tables
        allTables = [ "Antenna.xml",
                      "ASDM.xml",
                     # "CalData.xml",
                     # "CalDelay.xml",
                     # "CalReduction.xml",
                      "ConfigDescription.xml",
                      "CorrelatorMode.xml",
                      "DataDescription.xml",
                      "ExecBlock.xml",
                      "Feed.xml",
                      "Field.xml",
                     #"FocusModel.xml",
                     #"Focus.xml",
                      "Main.xml",
                      "PointingModel.xml",
                      "Polarization.xml",
                      "Processor.xml",
                      "Receiver.xml",
                      "SBSummary.xml",
                      "Scan.xml",
                      "Source.xml",
                      "SpectralWindow.xml",
                      "State.xml",
                      "Station.xml",
                      "Subscan.xml",
                      "SwitchCycle.xml"
                      ]
        isOK = True
        # test if xmllint is available and can be used in the following
        xmllint_ok = (os.system('xmllint --version') == 0)
        for fileName in allTables:
            filePath = asdmname+'/'+fileName
            if(not os.path.exists(filePath)):
                print("ASDM table file %s doesn't exist." % filePath)
                isOK = False
            elif(xmllint_ok):
                # test if well formed
                rval = os.system('xmllint --noout '+filePath)
                if(rval !=0):
                    print("Table %s is not a well formed XML document." % filePath)
                    isOK = False
        if(isOK and not xmllint_ok):
            print("Note: Test of XML well-formedness not possible since xmllint not available.")
        else:
            print("Note: xml validation not possible since ASDM DTDs (schemas) not yet online.")
            
        if(not os.path.exists(asdmname+"/ASDMBinary")):
            print("ASDM binary directory %s/ASDMBinary doesn't exist." % asdmname)
            isOK = False
    
        if(withPointing and not os.path.exists(asdmname+"/Pointing.bin")):
            print("ASDM binary file %s/Pointing.bin doesn't exist." % asdmname)
            isOK = False
    
        if (not isOK):
            raise Exception

# Test cases    
    def test1(self):
        '''Test 1: Testing default'''
        myvis = self.vis_b
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = False
        try:
            # this should throw an exception
            mysdm = sdm(self.out)
            mysdm.fromms( )
            exportasdm()
        except:
            self.rval = True
        self.assertTrue(self.rval)

    def test2(self):
        '''Test 2: small input MS, default output'''
        myvis = self.vis_b
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        mysdm = sdm(self.out)
        self.rval = mysdm.fromms('myinput.ms', archiveid="S1", verbose=True, apcorrected=False)

        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(2)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, False)

    def test3(self):
        '''Test 3: simulated input MS, default output'''
        myvis = self.vis_f
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        mysdm = sdm(self.out)
        self.rval = mysdm.fromms('myinput.ms',archiveid="S1")
        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(3)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, True)

    def test4(self):
        '''Test 4: real input MS, default output'''
        myvis = self.vis_d
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        mysdm = sdm(self.out)
        self.rval = mysdm.fromms('myinput.ms', archiveid="S1", apcorrected=False)
        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(4)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, False)

    def test5(self):
        '''Test 5: real input MS, MS has several SPWs observed in parallel - not supported, expected error'''
        myvis = self.vis_e
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        mysdm = sdm(self.out)
        self.rval = mysdm.fromms('myinput.ms', archiveid="S1", apcorrected=False)
        self.assertFalse(self.rval)

    def test6(self):
        '''Test 6: simulated input MS with pointing table, default output'''
        myvis = self.vis_g
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        mysdm = sdm(self.out)
        self.rval = mysdm.fromms('myinput.ms', archiveid="S002", apcorrected=False)
        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(6)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, True)

    def test7(self):
        '''Test 7: simulated input MS, default output'''
        myvis = self.vis_f
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        mysdm = sdm(self.out)
        self.rval = mysdm.fromms('myinput.ms', archiveid="S1")
        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(7)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, True)

    def test8(self):
        '''Test 8: real input MS, default output'''
        myvis = self.vis_d
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        mysdm = sdm(self.out)
        self.rval = mysdm.fromms('myinput.ms', archiveid="S1", apcorrected=False)
        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(8)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, False)

    def test9(self):
        '''Test 9: simulated input MS with pointing table, default output'''
        myvis = self.vis_g
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        mysdm = sdm(self.out)
        self.rval = mysdm.fromms('myinput.ms', archiveid="S002", apcorrected=False)
        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(9)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, True)

    @unittest.skip("test disabled")
    def test10(self):
        '''Test 10: ALMA input MS with pointing table and various shortcomings, default output'''
        myvis = self.vis_c
##        os.system('rm -rf myinput.ms')
##        os.system('cp -R ' + myvis + ' myinput.ms')
##         self.rval = exportasdm(
##             vis = 'myinput.ms',
##             asdm = self.out,
##             archiveid="S002",
##             apcorrected=False,
##             )
        self.assertNotEqual(self.rval,False)
##        omsname = "test"+str(10)+self.out
##        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
##        self.verify_asdm(omsname, True)

    def test11(self):
        '''Test 11: EVLA MS from X_osro_013.55979.93803716435 scan 2, full pol!'''
        myvis = self.vis_h
        os.system('rm -rf xosro2ref-reimp.ms xosro2asdm')
        mysdm = sdm('xosro2asdm')
        self.rval = mysdm.fromms(myvis,apcorrected=False,verbose=True)
        # mirroring the sdm.toms used to create the MS that was jsut exported
        self.rval = self.rval and mysdm.toms('xosro2ref-reimp.ms',process_flags=False,verbose=True,ocorr_mode='co',with_pointing_correction=True)
        self.rval = self.rval and th.compmsmainnumcol(myvis, 'xosro2ref-reimp.ms', 1E-5)
        self.rval = self.rval and th.compmsmainboolcol(myvis, 'xosro2ref-reimp.ms')

        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(11)+self.out
        os.system('rm -rf '+omsname+'; mv  xosro2asdm '+omsname)

    # test12 is skipped because exportasdm does additional steps to the MS before
    # sdm.fromms can be run on the result. That test doesn't add any additional value
    # in this tool test and so it is skipped here.
    @unittest.skip("This is a task test, not appropriate for the sdm tool alone.")
    def test12(self):
        '''Test 12: ALMA MS from uid___A002_X72bc38_X000 scan 2, only XX and YY'''
        myvis = self.vis_i
        os.system('rm -rf asdmasdm asdm-reimp.ms')
        mysdm = sdm('asdmasdm')
        self.rval = mysdm.fromms(myvis, apcorrected=False, verbose=True)
        self.rval = self.rval and mysdm.toms('asdm-reimp.ms', verbose=True)

        self.rval = self.rval  and th.compmsmainnumcol(myvis, 'asdm-reimp.ms', 1E-5)
        self.rval = self.rval and th.compmsmainboolcol(myvis, 'asdm-reimp.ms')

        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(12)+self.out
        os.system('rm -rf '+omsname+'; mv  asdm '+omsname)

# Class of tests for sdm.summarystr
class sdm_summarystr_test(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_alma_asdm(self):
        ''' ALMA M51 data'''
        # used in test_importasdm, test_importasdm_mms.
        mysdm = sdm(datapath + 'uid___X5f_X18951_X1')
        # self.doASDMSummary(asdmpath,174)
        summary = mysdm.summarystr().splitlines()
        self.assertTrue(len(summary) == 166, summary[0] if len(summary) == 1 else None)

    def test_vla_asdm(self):
        '''VLA data'''
        # used in test_importevla, test_importasdm_mms, test_importasdm
        mysdm = sdm(datapath + 'X_osro_013.55979.93803716435')
        summary = mysdm.summarystr().splitlines()
        self.assertTrue(len(summary) == 246, summary[0] if len(summary) == 1 else None)

    def test_aca_asdm(self):
        '''ACA with mixed pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        mysdm = sdm(datapath + 'uid___A002_X72bc38_X000')
        summary = mysdm.summarystr().splitlines()
        self.assertTrue(len(summary) == 2513, summary[0] if len(summary) == 1 else None)

    def test_12m_asdm(self):
        ''' 12m with mixedl pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        mysdm = sdm(datapath + 'uid___A002_X71e4ae_X317_short')
        summary = mysdm.summarystr().splitlines()
        self.assertTrue(len(summary) == 1017, summary[0] if len(summary) == 1 else None)

    def test_bogus_file(self):
        ''' 12m with mixedl pol/channelisation'''
        # used in test_importasdm_mms, test_importasdm
        passes = False
        try:
            mysdm = sdm('/tmp/some_nonexistent_directory_for_our_test')
            summary = mysdm.summarystr()
        except:
            passes = True
        self.assertTrue(passes, "non-existent file fails to throw an exception")

if __name__ == '__main__':
    unittest.main()
