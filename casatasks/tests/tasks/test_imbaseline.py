import random
import os
import numpy
import shutil
import unittest
import math

from casatools import ctsys, image, regionmanager, componentlist, table, quanta
from casatasks import imbaseline, casalog, imsubimage
from casatasks.private.task_imbaseline import ImBaselineVals as blv, do_imsmooth

_ia = image()
_rg = regionmanager()
_tb = table()
_qa = quanta()
ctsys_resolve = ctsys.resolve


class imsmooth_test(unittest.TestCase):

    image_names=['g192_a2.image', 'g192_a2.image-2.rgn']
    datapath = ctsys_resolve('unittest/imbaseline/')
    targetres_im = "imsmooth_targetres.fits"
    tiny = "tiny.im"

    def setUp(self):
        if os.path.exists(self.image_names[0]):
            for file in self.image_names:
                if os.path.exists(file): 
                    if os.path.isdir(file):
                        shutil.rmtree(file)
                    else:
                        os.unlink(file)
        
        for file in self.image_names:
            os.system('cp -RH '+ os.path.join(self.datapath,file)+' ' + file)
        self.ia = image()
        for f in [self.targetres_im, self.tiny]:
            if os.path.exists(f): 
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.unlink(f)
            os.system('cp -RH '+ os.path.join(self.datapath,f)+' ' + f)

    def tearDown(self):
        self.assertTrue(len(_tb.showcache()) == 0)
        # make sure directory is clean as per verification test requirement
        cwd = os.getcwd()
        for filename in os.listdir(cwd):
            file_path = os.path.join(cwd, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    # CASA 5 tests need this directory
                    if filename != 'xml':
                        shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e)) 
    
    def test_input(self):
        '''Imsmooth: Testing INPUT/OUTPUT tests'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        casalog.post( "Starting imsmooth INPUT/OUTPUT tests.", 'NORMAL2' )

        #######################################################################
        # Testing the imagename parameter.
        #    1. Bad file name should throw and exception
        #    2. Good file name, a file should be
        #######################################################################
        casalog.post( "The IMAGENAME parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        major = "2.5arcsec"
        minor = "2arcsec"
        pa = "0deg"
        result = None   
        beam = {"major": major, "minor": minor, "pa": pa}

        try:
            blv(imagename='g192')
        except:
            no_op='noop'
        else:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Badfile, 'g192', was not reported as missing."

        vals = blv(imagename = self.tiny, dirkernel='gaussian', major="2.5arcsec", minor="2arcsec", pa="0deg")
        vals.imsmooth_output = 'input_test1'
        do_imsmooth(vals)
        self.assertTrue(os.path.exists(vals.imsmooth_output))

        #######################################################################
        # Testing the outfile parameter.
        #    1. Bad file, file already exists, exception should be thrown
        #    2. Good file name, a file should be
        #######################################################################
        casalog.post( "The OUTFILE parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
        vals = blv(imagename = self.tiny, dirkernel='gaussian', major="2.5arcsec", minor="2arcsec", pa="0deg")
        vals.imsmooth_output = 'input_test2'
        do_imsmooth(vals)
        self.assertTrue(os.path.exists(vals.imsmooth_output))

        vals = blv(imagename = self.tiny, dirkernel='gaussian', major="2.5arcsec", minor="2arcsec", pa="0deg")
        vals.imsmooth_output = 'input_test2'
        do_imsmooth(vals)
        self.assertTrue(os.path.exists(vals.imsmooth_output))

        #######################################################################
        # Testing KERNEL parameter, valid values 0 and greater
        #    1. Below invalid range: junk, ''
        #    3. valid: gaussian, boxcar, tophat, 
        #######################################################################
        casalog.post( "The KERNEL parameter tests will cause errors to occur, do not be alarmed", 'WARN' )

        try:
            vals = blv(imagename = self.tiny, dirkernel='junk', major="2.5arcsec", minor="2arcsec", pa="0deg")
            vals.imsmooth_output = 'input_test3'
            do_imsmooth(vals)
        except:
            pass
        else:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                +"\nError: No exception thrown for bad kernel value, 'junk'"

        try:
            vals = blv(imagename = self.tiny, dirkernel='gaussian', major="2.5arcsec", minor="2arcsec", pa="0deg")
            vals.imsmooth_output = 'input_test4'
            do_imsmooth(vals)
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Gaussian smoothing failed."
        if ( not os.path.exists( vals.imsmooth_output ) ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +f"\nError: output file {vals.imsmooth_output} was NOT created."

        try:
            vals = blv(imagename = self.tiny, dirkernel='boxcar', major="2arcsec", minor="2arcsec", pa="0deg")
            vals.imsmooth_output = 'input_test5'
            do_imsmooth(vals)
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Boxcar smoothing failed. " \
                     +str(err)
        if ( not os.path.exists( vals.imsmooth_output )): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +f"\nError: output file {vals.imsmooth_output} was NOT created."

        #######################################################################
        # Testing MAJOR parameter
        # Expects a numerical value: 1, '2pix', '0.5arcsec'
        # Tests include: invalid values, valid values, major < minor 
        #######################################################################
        casalog.post( "The MAJOR parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
        try:
            vals = blv(imagename = self.tiny, dirkernel='boxcar', major="bad", minor="2arcsec", pa="0deg")
            vals.imsmooth_output = 'input_test8'
            do_imsmooth(vals)
        except:
            no_op='noop'
        else:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Bad major value, 'bad', was not reported."

        try:
            vals = blv(imagename = self.tiny, dirkernel='boxcar', major="-5", minor="2arcsec", pa="0deg")
            vals.imsmooth_output = 'input_test9'
            do_imsmooth( vals )
        except:
            no_op='noop'
        else:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Bad major value, '-5', was not reported."
    