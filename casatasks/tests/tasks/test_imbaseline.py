import random
import os
import numpy
import shutil
import unittest
import math
from scipy import signal

from casatools import ctsys, image, regionmanager, componentlist, table, quanta
from casatasks import imbaseline, casalog, imsubimage
from casatasks.private.task_imbaseline import ImBaselineVals as blv, do_imsmooth

_ia = image()
_rg = regionmanager()
_tb = table()
_qa = quanta()
ctsys_resolve = ctsys.resolve

### imsmooth ###

class imsmooth_test(unittest.TestCase):

    image_names=['g192_a2.image', 'g192_a2.image-2.rgn']
    datapath = ctsys_resolve('unittest/imsmooth/')
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
            vals = blv(imagename = self.tiny, dirkernel='boxcar', major="-5arcsec", minor="2arcsec", pa="0deg")
            vals.imsmooth_output = 'input_test9'
            do_imsmooth( vals )
        except:
            no_op='noop'
        else:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Bad major value, '-5arcsec', was not reported."

        vals = blv(imagename = self.tiny, dirkernel='gaussian', major="2arcsec", minor="1arcsec", pa="0deg")
        vals.imsmooth_output = 'input_test11'
        do_imsmooth(vals)
        self.assertTrue(os.path.exists(vals.imsmooth_output))
            
        vals = blv(imagename = self.tiny, dirkernel='gaussian', major="2pix", minor="1pix", pa="0deg")
        vals.imsmooth_output = 'input_test12'
        do_imsmooth(vals)
        self.assertTrue(os.path.exists(vals.imsmooth_output))

        vals = blv(imagename = self.tiny, dirkernel='gaussian', major="1.5arcsec", minor="1arcsec", pa="0deg")
        vals.imsmooth_output = 'input_test13'
        do_imsmooth(vals)
        self.assertTrue(os.path.exists(vals.imsmooth_output))

        try:
            vals = blv(imagename = self.tiny, dirkernel='gaussian', major="0.5arcsec", minor="2arcsec", pa="0deg")
            vals.imsmooth_output = 'input_test14'
            do_imsmooth(vals)
        except:
            pass
        else:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Bad major value less than minor value was not reported."


    def test_smooth(self):
        '''Imsmooth: Testing correctness'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        casalog.post( "Starting imsmooth CORRECTNESS tests.", 'NORMAL2' )
    
        # First step get rid of all the old test files!
        for file in os.listdir( '.' ):
            if file.startswith( 'smooth_test' ):
                shutil.rmtree( file )
    
        if os.path.exists( 'smooth.pointsrc.image' ):
            shutil.rmtree( 'smooth.pointsrc.image' )
    
        # Second step is to create a file with a single point
        # source so that we can check the correctness.  The
        # resulting convolution should be the same shape as
        # the kernel that is used if done correctly.  Also the
        # area under the kernel should be equivalent to the value
        # our point source.
        #
        # We use the the coordinate system from the g192 image
        # and make our image the same size.  In theory it could
        # be anything, it's nice having a coordinate system for
        # the image.
        try:
            # Get the coordinate system and size of the image
            _ia.open( 'g192_a2.image' )
            csys = _ia.coordsys()
            bb = _ia.boundingbox()
            shape = bb['bbShape']
            _ia.done()
    
            # Create an array of zero's, then set position 212,220,0,20
            # to 100 (our point source).
            #
            # Note that 
            inputArray = numpy.zeros( (shape[0], shape[1], shape[2], shape[3]), 'float' )
            inputArray[212,220,0,20] = 100
    
            # Now make the image!
            _ia.fromarray( pixels=inputArray, csys=csys.torecord(), \
                          outfile='smooth.pointsrc.image' )
            _ia.done()
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Unable to create point source image."\
                     +"\n\t REULTS: "+str(err)        
    
        # Do a Gaussian smoothing with major axis of 50, and minor of 25
        # pixels.  We expect the resulting image to have a Gussian shape,
        # with    max at:    212,220,0,20
        #         1st sd:    from
        #         2nd sd:    from
    
        try:
            vals = blv(imagename='smooth.pointsrc.image', dirkernel='gaussian', major="50arcsec", minor="25arcsec", pa="0deg")
            vals.imsmooth_output = 'smooth_test1'
            do_imsmooth( vals )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: boxcar smooth failed on point source file."
    
            
        if ( not os.path.exists( 'smooth_test1' ) ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Gaussian smoothfailed on point source file."
        else:
            # Now that we know something has been done lets check the results!
            #      1. Check that the sum of the values under the curve is 100
            #      2. Check that the max is at 212, 220, 0 , 20
            allowedError = 0.009
            
            _ia.open('smooth_test1')
            stats = _ia.statistics()
            sum = stats['sum'][0]
            if ( ( sum < 100 and sum < ( 100-allowedError ) )
                 or ( sum > 100 and sum > ( 100+allowedError) ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Sum under Gaussian is "+str(stats['sum'][0])\
                    +" expected 100."
    
            maxpos=stats['maxpos'].tolist()
            if ( maxpos[0]!=212 or maxpos[1]!=220 or \
                 maxpos[2]!=0 or maxpos[3]!=20 ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Max position found at "+str(maxpos)\
                    +" expected it to be at 212,220,0,20."            
            
            _ia.done()
                
    
        # Do a box car smooth and verify expected results as follows:
        #
        try:
            vals = blv(imagename='smooth.pointsrc.image', dirkernel='boxcar', major="20arcsec", minor="10arcsec", pa="0deg")
            vals.imsmooth_output = 'smooth_test2'
            do_imsmooth( vals )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: boxcar smooth failed on point source file."
            
        if ( not os.path.exists( 'smooth_test2' ) ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: output file 'smooth_test2' was NOT created."
        else:
            # Now that we know something has been done lets check the results!
            #        1. Check that the sum of the points is 100
            #        2. That the points in the box are 0.125=(100/((10+10)*(20+20))
            _ia.open('smooth_test2')
            stats = _ia.statistics()
            if ( ( sum < 100 and sum < ( 100-allowedError ) )
                 or ( sum > 100 and sum > ( 100+allowedError) ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Sum under Gaussian is "+str(stats['sum'][0])\
                    +" expected 100."

            pixels = list(map( lambda pv: _ia.pixelvalue(pv)['value']['value'],
                          [ [204,200,0,20], [222,236,0,20], [204,239,0,20],
                            [222,201,0,20], [212,220,0,20] ] ))

            for value in pixels:
                if ( not (value>(0.125-allowedError) and value<(0.125+allowedError)) and
                     not (value>-3.65746967e-09 and value<3.65746967e-09) ):
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']\
                        +"\nError: Values in the smoothed box are not all 0.125"\
                        +" found value of "+str(value)
            _ia.done()
    
        self.assertTrue(retValue['success'],retValue['error_msgs'])

    def test_image_kernel(self):
        """Test image as kernel, CAS-5844"""
        imagename = os.path.join(self.datapath,"point.im")
        kimage = os.path.join(self.datapath,"bessel.im")

        vals = blv(imagename=imagename, dirkernel="image", kimage=kimage)
        vals.imsmooth_output = "point_c_bessel.im"
        do_imsmooth(vals)
        self.assertTrue(os.path.exists(vals.imsmooth_output))
        myia = self.ia
        myia.open(kimage)
        bessel = myia.getchunk()
        myia.open(vals.imsmooth_output)
        conv = myia.getchunk()
        myia.done()
        ratio = conv/bessel
        print("max",abs(ratio - ratio[50,50]).max())
        self.assertTrue((abs(ratio - ratio[50,50]) < 2e-4).all())

        vals.imsmooth_scale = 1
        vals.imsmooth_output = "point_c_bessel_scale1.im"
        do_imsmooth(vals)
        self.assertTrue(os.path.exists(vals.imsmooth_output))
        myia.open(vals.imsmooth_output)
        conv = myia.getchunk()
        myia.done()
        diff = conv - bessel
        self.assertTrue(abs(diff).max() < 2e-7)


def suite():
    return [imsmooth_test]    

if __name__ == '__main__':
    unittest.main()