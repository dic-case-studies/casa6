# unit test for the exportasdm task

import os
import sys
import shutil

import unittest
from CASAtools import ctsys, ms
from CASAtasks import exportasdm
### for testhelper import
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import testhelper as th

_ms = ms( )
class exportasdm_test(unittest.TestCase):
    
    #vis_a = 'ngc4826.ms'
    vis_b = 'test.ms'
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
            os.system('cp -R '+ctsys.resolve('regression/fits-import-export/input/test.ms')+' .')
        if(not os.path.exists(self.vis_c)):
            os.system('cp -R '+ctsys.resolve('regression/exportasdm/input/M100-X220-shortened.ms')+' .')
        if(not os.path.exists(self.vis_d)):
            _ms.fromfits( self.vis_d, 'regression/ngc4826/fitsfiles/ngc4826.ll.fits5' )
            _ms.close( )
        if(not os.path.exists(self.vis_e)):
            os.system('cp -R '+ctsys.resolve('regression/cvel/input/g19_d2usb_targets_line-shortened.ms')+' .')
        if(not os.path.exists(self.vis_f)):
            os.system('cp -R '+ctsys.resolve('regression/exportasdm/input/Itziar.ms')+' .')
        if(not os.path.exists(self.vis_g)):
            os.system('cp -R '+ctsys.resolve('regression/exportasdm/input/M51.ms')+' .')

        #> to be factored back in when sdm.toms( ) is implemented
        #> ---- --- --- ---- --- --- ---- --- --- ---- --- --- ---- --- --- ---- --- --- ---- --- ---
        #> if(not os.path.exists(self.vis_i)):
        #>     os.system('ln -sf '+os.environ['CASAPATH'].split()[0]+'/data/regression/asdm-import/input/uid___A002_X72bc38_X000')
        #>     importasdm('uid___A002_X72bc38_X000', vis = 'asdm.ms', scans='0:2')

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
            self.rval = exportasdm()
        except:
            self.rval = True
        self.assertTrue(self.rval)

    def test2(self):
        '''Test 2: small input MS, default output, v3'''
        myvis = self.vis_b
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = exportasdm(
                vis = 'myinput.ms',
                asdm = self.out,
                archiveid="S1",
                verbose=True,
                apcorrected=False,
                useversion='v3')

        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(2)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, False)

    def test3(self):
        '''Test 3: simulated input MS, default output, v3'''
        myvis = self.vis_f
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = exportasdm(vis = 'myinput.ms',asdm = self.out,archiveid="S1", useversion='v3')
        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(3)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, True)

    def test4(self):
        '''Test 4: real input MS, default output, v3'''
        myvis = self.vis_d
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = exportasdm(
            vis = 'myinput.ms',
            asdm = self.out,
            archiveid="S1",
            apcorrected=False,
            useversion='v3'
            )

        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(4)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, False)

    def test5(self):
        '''Test 5: real input MS, MS has several SPWs observed in parallel, v3 - not supported, expected error'''
        myvis = self.vis_e
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = exportasdm(
            vis = 'myinput.ms',
            asdm = self.out,
            archiveid="S1",
            apcorrected=False,
            useversion = 'v3'
            )

        self.assertFalse(self.rval)

    def test6(self):
        '''Test 6: simulated input MS with pointing table, default output, v3'''
        myvis = self.vis_g
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = exportasdm(
            vis = 'myinput.ms',
            asdm = self.out,
            archiveid="S002",
            apcorrected=False,
            useversion = 'v3'
            )

        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(6)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, True)

    def test7(self):
        '''Test 7: v3, simulated input MS, default output'''
        myvis = self.vis_f
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = exportasdm(vis = 'myinput.ms',asdm = self.out,archiveid="S1", useversion='v3')
        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(7)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, True)

    def test8(self):
        '''Test 8: v3, real input MS, default output'''
        myvis = self.vis_d
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = exportasdm(
            vis = 'myinput.ms',
            asdm = self.out,
            archiveid="S1",
            apcorrected=False,
            useversion='v3'
            )

        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(8)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, False)

    def test9(self):
        '''Test 9: v3, simulated input MS with pointing table, default output'''
        myvis = self.vis_g
        os.system('rm -rf myinput.ms')
        os.system('cp -R ' + myvis + ' myinput.ms')
        self.rval = exportasdm(
            vis = 'myinput.ms',
            asdm = self.out,
            archiveid="S002",
            apcorrected=False,
            useversion='v3'
            )

        self.assertNotEqual(self.rval,False)
        omsname = "test"+str(9)+self.out
        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
        self.verify_asdm(omsname, True)

    @unittest.skip("test disabled")
    def test10(self):
        '''Test 10: v3, ALMA input MS with pointing table and various shortcomings, default output'''
        myvis = self.vis_c
##        os.system('rm -rf myinput.ms')
##        os.system('cp -R ' + myvis + ' myinput.ms')
##         self.rval = exportasdm(
##             vis = 'myinput.ms',
##             asdm = self.out,
##             archiveid="S002",
##             apcorrected=False,
##             useversion='v3'
##             )
        self.assertNotEqual(self.rval,False)
##        omsname = "test"+str(10)+self.out
##        os.system('rm -rf '+omsname+'; mv exportasdm-output.asdm '+omsname)
##        self.verify_asdm(omsname, True)

    @unittest.skip("waiting for sdm.toms(...)")
    def test12(self):
        '''Test 12: v3, ALMA MS from uid___A002_X72bc38_X000 scan 2, only XX and YY'''
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



class exportasdm_test2(unittest.TestCase):

    vis_b = 'test.ms'
    vis_c = 'M100-X220-shortened.ms'
    vis_d = 'ngc4826.tutorial.ngc4826.ll.5.ms'
    vis_e = 'g19_d2usb_targets_line-shortened.ms'
    vis_f = 'Itziar.ms'
    vis_g = 'M51.ms'
    vis_h = 'xosro2ref.ms'
    vis_i = 'asdm.ms'
    
    def setUp(self):  
        pass

    def tearDown(self):
        shutil.rmtree(self.vis_c,ignore_errors=True)
        shutil.rmtree(self.vis_b,ignore_errors=True)
        shutil.rmtree(self.vis_d,ignore_errors=True)
        shutil.rmtree(self.vis_e,ignore_errors=True)
        shutil.rmtree(self.vis_f,ignore_errors=True)
        shutil.rmtree(self.vis_g,ignore_errors=True)
        shutil.rmtree(self.vis_h,ignore_errors=True)
        shutil.rmtree(self.vis_i,ignore_errors=True)
        os.system('rm -rf test*exportasdm*.asdm')
    
    def test1a(self):
        '''Exportasdm: Cleanup'''
        pass

def suite():
    return [exportasdm_test,exportasdm_test2]

if __name__ == '__main__':
    unittest.main()
