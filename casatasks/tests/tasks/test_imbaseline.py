import random
import os
import numpy
import re
import shutil
import sys
import unittest
import math
from scipy import signal

from casatools import ctsys, image, regionmanager, componentlist, table, quanta ,ms
from casatasks import imbaseline, casalog, imsubimage
from casatasks.private.task_imbaseline import ImBaselineVals as blv, do_imsmooth, do_sdbaseline, do_sdsmooth
from casatasks.private.sdutil import table_manager, calibrater_manager

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


### sdsmooth ###

def gaussian_kernel(nchan, kwidth):
    sigma = kwidth / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    g = signal.gaussian(nchan, sigma, False)
    g /= g.sum()
    g0 = g[0]
    g[:-1] = g[1:]
    g[-1] = g0
    return g

class sdsmooth_test_base(unittest.TestCase):
    """
    Base class for sdsmooth unit test.
    The following attributes/functions are defined here.

        datapath
        decorators (invalid_argument_case, exception_case)
    """
    datapath=ctsys.resolve('unittest/sdsmooth/')
    infile_data = 'tsdsmooth_test.ms'
    infile_float = 'tsdsmooth_test_float.ms'

    # task execution result
    result = None

    @property
    def outfile(self):
        return self.infile.rstrip('/') + '_out'

    # decorators
    @staticmethod
    def invalid_argument_case(func):
        """
        Decorator for the test case that is intended to fail
        due to invalid argument.
        """
        import functools
        @functools.wraps(func)
        def wrapper(self):
            func(self)
            self.assertFalse(self.result, msg='The task must return False')
        return wrapper

    @staticmethod
    def exception_case(exception_type, exception_pattern):
        """
        Decorator for the test case that is intended to throw
        exception.

            exception_type: type of exception
            exception_pattern: regex for inspecting exception message
                               using re.search
        """
        def wrapper(func):
            import functools
            @functools.wraps(func)
            def _wrapper(self):
                self.assertTrue(len(exception_pattern) > 0, msg='Internal Error')
                with self.assertRaises(exception_type) as ctx:
                    func(self)
                    self.fail(msg='The task must throw exception')
                the_exception = ctx.exception
                message = str(the_exception)
                self.assertIsNotNone(re.search(exception_pattern, message), msg='error message \'%s\' is not expected.'%(message))
            return _wrapper
        return wrapper

    @staticmethod
    def weight_case(func):
        import functools
        @functools.wraps(func)
        def wrapper(self):
            with table_manager(self.infile) as tb:
                for irow in range(tb.nrows()):
                    self.assertTrue(tb.iscelldefined('WEIGHT_SPECTRUM', irow))

            # weight mode flag
            self.weight_propagation = True

            func(self)

        return wrapper

    def exec_sdsmooth(self, infile, outfile, vals):
        vals.temporary_vis = infile
        vals.sdsmooth_output = outfile
        do_sdsmooth(vals)
    
    def set_blv_params(self, vals, args):
        for key, val in args.items():
            if key in vals.__dict__:
                vals.__dict__[key] = val

    def run_test(self, *args, **kwargs):
        datacol_name = self.datacolumn.upper()
        weight_mode = hasattr(self, 'weight_propagation') and getattr(self, 'weight_propagation') is True

        if 'kwidth' in kwargs:
            kwidth = kwargs['kwidth']
        else:
            kwidth = 5

        vals = blv(imagename=self.infile, spkernel='gaussian', kwidth=kwidth)
        vals.datacolumn = datacol_name
        self.exec_sdsmooth(self.infile, self.outfile, vals)

        # sanity check
        self.assertTrue(os.path.exists(self.outfile), msg='Output file is not properly created.')

        if 'spw' in kwargs:
            spw = kwargs['spw']
        else:
            spw = ''
        dd_selection = None
        if len(spw) == 0:
            expected_nrow = 2
            with table_manager(self.infile) as tb:
                data_in = tb.getvarcol(datacol_name)
                flag_in = tb.getvarcol('FLAG')
                if weight_mode:
                    weight_in = tb.getvarcol('WEIGHT_SPECTRUM')
        else:
            myms = ms()
            a = myms.msseltoindex(self.infile, spw=spw)
            spw_selection = a['spw']
            dd_selection = a['dd']
            expected_nrow = len(spw_selection)
            with table_manager(self.infile) as tb:
                try:
                    tsel = tb.query('DATA_DESC_ID IN %s'%(dd_selection.tolist()))
                    data_in = tsel.getvarcol(datacol_name)
                    flag_in = tsel.getvarcol('FLAG')
                    if weight_mode:
                        weight_in = tsel.getvarcol('WEIGHT_SPECTRUM')
                finally:
                    tsel.close()

        with table_manager(self.outfile) as tb:
            nrow = tb.nrows()
            data_out = tb.getvarcol(datacol_name)
            flag_out = tb.getvarcol('FLAG')
            if weight_mode:
                weight_out = tb.getvarcol('WEIGHT_SPECTRUM')

        # verify nrow
        self.assertEqual(nrow, expected_nrow, msg='Number of rows mismatch (expected %s actual %s)'%(expected_nrow, nrow))

        # verify data
        eps = 1.0e-6
        for key in data_out.keys():
            row_in = data_in[key]
            flg_in = flag_in[key]
            row_in[numpy.where(flg_in == True)] = 0.0
            row_out = data_out[key]
            self.assertEqual(row_in.shape, row_out.shape, msg='Shape mismatch in row %s'%(key))

            npol, nchan, _ = row_out.shape
            kernel_array = gaussian_kernel(nchan, kwidth)
            expected = numpy.convolve(row_in[0,:,0], kernel_array, mode='same')
            output = row_out[0,:,0]
            zero_index = numpy.where(numpy.abs(expected) <= eps)
            self.assertTrue(all(numpy.abs(output[zero_index]) < eps), msg='Failed to verify zero values: row %s'%(key))
            nonzero_index= numpy.where(numpy.abs(expected) > eps)
            diff = numpy.abs((output[nonzero_index] - expected[nonzero_index]) / expected[nonzero_index].max())
            #print diff
            #print output[nonzero_index]
            #print expected[nonzero_index]
            self.assertTrue(all(diff < eps), msg='Failed to verify nonzero values: row %s'%(key))
            #print 'row_in', row_in[0,:,0].tolist()
            #print 'gaussian', kernel_array.tolist()
            #print 'expected', expected.tolist()
            #print 'result', row_out[0,:,0].tolist()

            # weight check if this is weight test
            if weight_mode:
                #print 'Weight propagation test'
                wgt_in = weight_in[key]
                wgt_out = weight_out[key]
                wkwidth = int(kwidth + 0.5)
                wkwidth += (1 if wkwidth % 2 == 0 else 0)
                half_width = wkwidth // 2
                peak_chan = kernel_array.argmax()
                start_chan = peak_chan - half_width
                wkernel = kernel_array[start_chan:start_chan+wkwidth].copy()
                wkernel /= sum(wkernel)
                weight_expected = wgt_in.copy()
                for ichan in range(half_width, nchan-half_width):
                    s = numpy.zeros(npol, dtype=float)
                    for jchan in range(wkwidth):
                        s += wkernel[jchan] * wkernel[jchan] / wgt_in[:,ichan-half_width+jchan,0]
                    weight_expected[:,ichan,0] = 1.0 / s
                #print weight_expected[:,:10]
                diff = numpy.abs((wgt_out - weight_expected) / weight_expected)
                self.assertTrue(all(diff.flatten() < eps), msg='Failed to verify spectral weight: row %s'%(key))

    def _setUp(self, files):
        for f in files:
            if os.path.exists(f):
                shutil.rmtree(f)
            shutil.copytree(os.path.join(self.datapath, f), f)

    def _tearDown(self, files):
        for f in files:
            if os.path.exists(f):
                shutil.rmtree(f)

    def setUp(self):
        self._setUp([self.infile])

    def tearDown(self):
        self._tearDown([self.infile, self.outfile])

class sdsmooth_test_fail(sdsmooth_test_base):
    """
    Unit test for task sdsmooth.

    The list of tests:
    test_sdsmooth_fail01 --- default parameters (raises an error)
    test_sdsmooth_fail02 --- invalid kernel type
    test_sdsmooth_fail03 --- invalid selection (empty selection result)
    """
    invalid_argument_case = sdsmooth_test_base.invalid_argument_case
    exception_case = sdsmooth_test_base.exception_case
    infile = sdsmooth_test_base.infile_data

    @invalid_argument_case
    def test_sdsmooth_fail01(self):
        """test_sdsmooth_fail01 --- default parameters (raises an error)"""
        self.assertRaises(Exception, do_sdsmooth)

    @invalid_argument_case
    def test_sdsmooth_fail02(self):
        """test_sdsmooth_fail02 --- invalid kernel type"""
        self.assertRaises(ValueError, blv, imagename=self.infile, spkernel='normal') # must move to blv test part

    @exception_case(RuntimeError, 'Spw Expression: No match found for 3')
    def test_sdsmooth_fail03(self):
        """test_sdsmooth_fail03 --- invalid selection (empty selection result)"""
        vals = blv(imagename=self.infile, spkernel='gaussian')
        vals.sdsmooth_spw = '3'
        self.exec_sdsmooth(self.infile, self.outfile, vals)


class sdsmooth_test_complex(sdsmooth_test_base):
    """
    Unit test for task sdsmooth. Process MS having DATA column.

    The list of tests:
    test_sdsmooth_complex_fail01 --- non-existing data column (FLOAT_DATA)
    test_sdsmooth_complex_gauss01 --- gaussian smoothing (kwidth 5)
    test_sdsmooth_complex_gauss02 --- gaussian smoothing (kwidth 3)
    """
    exception_case = sdsmooth_test_base.exception_case
    infile = sdsmooth_test_base.infile_data
    datacolumn = 'data'

    @exception_case(RuntimeError, 'Desired column \(FLOAT_DATA\) not found in the input MS')
    def test_sdsmooth_complex_fail01(self):
        """test_sdsmooth_complex_fail01 --- non-existing data column (FLOAT_DATA)"""
        vals = blv(imagename=self.infile, spkernel='gaussian')
        vals.datacolumn = 'float_data'
        self.exec_sdsmooth(self.infile, self.outfile, vals)

    def test_sdsmooth_complex_gauss01(self):
        """test_sdsmooth_complex_gauss01 --- gaussian smoothing (kwidth 5)"""
        self.run_test(kwidth=5)

    def test_sdsmooth_complex_gauss02(self):
        """test_sdsmooth_complex_gauss02 --- gaussian smoothing (kwidth 3)"""
        self.run_test(kwidth=3)


class sdsmooth_test_float(sdsmooth_test_base):
    """
    Unit test for task sdsmooth. Process MS having FLOAT_DATA column.

    The list of tests:
    test_sdsmooth_float_fail01 --- non-existing data column (DATA)
    test_sdsmooth_float_gauss01 --- gaussian smoothing (kwidth 5)
    test_sdsmooth_float_gauss02 --- gaussian smoothing (kwidth 3)
    """
    exception_case = sdsmooth_test_base.exception_case
    infile = sdsmooth_test_base.infile_float
    datacolumn = 'float_data'

    @exception_case(RuntimeError, 'Desired column \(DATA\) not found in the input MS')
    def test_sdsmooth_float_fail01(self):
        """test_sdsmooth_complex_fail01 --- non-existing data column (DATA)"""
        vals = blv(imagename=self.infile, spkernel='gaussian')
        vals.datacolumn = 'data'
        self.exec_sdsmooth(self.infile, self.outfile, vals)

    def test_sdsmooth_float_gauss01(self):
        """test_sdsmooth_float_gauss01 --- gaussian smoothing (kwidth 5)"""
        self.run_test(kwidth=5)

    def test_sdsmooth_float_gauss02(self):
        """test_sdsmooth_float_gauss02 --- gaussian smoothing (kwidth 3)"""
        self.run_test(kwidth=3)

class sdsmooth_test_weight(sdsmooth_test_base):
    """
    Unit test for task sdsmooth. Verify weight propagation.

    The list of tests:
    test_sdsmooth_weight_gauss01 --- gaussian smoothing (kwidth 5)
    test_sdsmooth_weight_gauss02 --- gaussian smoothing (kwidth 3)
    """
    weight_case = sdsmooth_test_base.weight_case
    infile = sdsmooth_test_base.infile_data
    datacolumn = 'data'

    def setUp(self):
        super(sdsmooth_test_weight, self).setUp()

        # initialize WEIGHT_SPECTRUM
        with calibrater_manager(self.infile) as cb:
            cb.initweights()

    @weight_case
    def test_sdsmooth_weight_gauss01(self):
        """test_sdsmooth_weight_gauss01 --- gaussian smoothing (kwidth 5)"""
        self.run_test(kwidth=5)

    @weight_case
    def test_sdsmooth_weight_gauss02(self):
        """test_sdsmooth_weight_gauss02 --- gaussian smoothing (kwidth 3)"""
        self.run_test(kwidth=3)


class sdsmooth_test_boxcar(sdsmooth_test_base):
    """
    Unit test for checking boxcar smoothing.

    The input data (sdsmooth_delta.ms) has data with the following features:
      in row0, pol0: 1 at index 100, 0 elsewhere,
      in row0, pol1: 1 at index 0 and 2047(i.e., at both ends), 0 elsewhere,
      in row1, pol0: 1 at index 10 and 11, 0 elsewhere,
      in row1, pol1: 0 throughout.
    If input spectrum has delta-function-like feature, the
    expected output spectrum will be smoothing kernel itself.
    As for the data at [row0, pol0], the output data will be:
      kwidth==1 -> spec[100] = 1
      kwidth==2 -> spec[100,101] = 1/2 (=0.5)
      kwidth==3 -> spec[99,100,101] = 1/3 (=0.333...)
      kwidth==4 -> spec[99,100,101,102] = 1/4 (=0.25)
      kwidth==5 -> spec[98,99,100,101,102] = 1/5 (=0.2)
      and so on.
    """

    infile = 'tsdsmooth_delta.ms'
    datacolumn = 'float_data'
    centers = {'00': [100], '01': [0,2047], '10': [10,11], '11':[]}

    def _getLeftWidth(self, kwidth):
        assert(0 < kwidth)
        return (2-kwidth)//2

    def _getRightWidth(self, kwidth):
        assert(0 < kwidth)
        return kwidth//2

    def _checkResult(self, spec, kwidth, centers, tol=5.0e-06):
        sys.stdout.write('testing kernel_width = '+str(kwidth)+'...')
        for i in range(len(spec)):
            count = 0
            for j in range(len(centers)):
                lidx = centers[j] + self._getLeftWidth(kwidth)
                ridx = centers[j] + self._getRightWidth(kwidth)
                if (lidx <= i) and (i <= ridx): count += 1
            value = count/float(kwidth)
            self.assertTrue(((spec[i] - value) < tol), msg='Failed.')
        sys.stdout.write('OK.\n')

    def setUp(self):
        super(sdsmooth_test_boxcar, self).setUp()

    def test000(self):
        # testing kwidth from 1 to 5.
        for kwidth in range(1,6):
            vals = blv(imagename=self.infile, spkernel='boxcar', kwidth = kwidth)
            vals.datacolumn = self.datacolumn
            self.exec_sdsmooth(self.infile, self.outfile, vals)
            with table_manager(self.outfile) as tb:
                for irow in range(tb.nrows()):
                    spec = tb.getcell(self.datacolumn.upper(), irow)
                    for ipol in range(len(spec)):
                        center = self.centers[str(irow)+str(ipol)]
                        self._checkResult(spec[ipol], kwidth, center)

    def test000_datacolumn_uppercase(self):
        # testing kwidth from 1 to 5.
        datacolumn = "FLOAT_DATA"
        for kwidth in range(1,6):
            vals = blv(imagename=self.infile, spkernel='boxcar', kwidth = kwidth)
            vals.datacolumn = self.datacolumn
            self.exec_sdsmooth(self.infile, self.outfile, vals)
            with table_manager(self.outfile) as tb:
                for irow in range(tb.nrows()):
                    spec = tb.getcell(datacolumn.upper(), irow)
                    for ipol in range(len(spec)):
                        center = self.centers[str(irow)+str(ipol)]
                        self._checkResult(spec[ipol], kwidth, center)


class sdsmooth_selection(sdsmooth_test_base, unittest.TestCase):
    infile = "analytic_type1.sm.ms"
    outfile = "smoothed.ms"
    # common_param = dict(infile=infile, outfile=outfile,
    #                     kernel='boxcar', kwidth=5)
    selections=dict(sdsmooth_intent=("CALIBRATE_ATMOSPHERE#OFF*", [1]),
                    sdsmooth_antenna=("DA99", [1]),
                    sdsmooth_field=("M1*", [0]),
                    sdsmooth_spw=(">6", [1]),
                    sdsmooth_timerange=("2013/4/28/4:13:21",[1]),
                    sdsmooth_scan=("0~8", [0]),
                    sdsmooth_pol=("YY", [1]))
    verbose = False

    def _get_selection_string(self, key):
        if key not in self.selections.keys():
            raise ValueError("Invalid selection parameter %s" % key)
        return {key: self.selections[key][0]}

    def _get_selected_row_and_pol(self, key):
        if key not in self.selections.keys():
            raise ValueError("Invalid selection parameter %s" % key)
        pols = [0,1]
        rows = [0,1]
        if key == 'sdsmooth_pol':  #self.selection stores pol ids
            pols = self.selections[key][1]
        else: #self.selection stores row ids
            rows = self.selections[key][1]
        return (rows, pols)

    def _get_reference(self, nchan, row_offset, pol_offset, datacol):
        if datacol.startswith("float"):
            col_offset = 10
        elif datacol.startswith("corr"):
            col_offset = 50
        else:
            raise ValueError("Got unexpected datacolumn.")
        spike_chan = col_offset + 20*row_offset + 10*pol_offset
        reference = numpy.zeros(nchan)
        reference[spike_chan-2:spike_chan+3] = 0.2
        if self.verbose: print("reference=%s" % str(reference))
        return reference

    def run_test(self, sel_param, datacolumn, reindex=True):
        inparams = self._get_selection_string(sel_param)
        vals = blv(imagename=self.infile, spkernel='boxcar', kwidth=5)
        vals.datacolumn = datacolumn
        vals.sdsmooth_reindex = reindex
        self.set_blv_params(vals, inparams)
        self.exec_sdsmooth(self.infile, self.outfile, vals)
        self._test_result(self.outfile, sel_param, datacolumn)

    def _test_result(self, msname, sel_param, dcol, atol=1.e-5, rtol=1.e-5):
        # Make sure output MS exists
        self.assertTrue(os.path.exists(msname), "Could not find output MS")
        # Compare output MS with reference (nrow, npol, and spectral values)
        (rowids, polids) = self._get_selected_row_and_pol(sel_param)
        if dcol.startswith("float"):
            testcolumn = "FLOAT_DATA"
        else: #output is in DATA column
            testcolumn = "DATA"
        _tb.open(msname)
        try:
            self.assertEqual(_tb.nrows(), len(rowids), "Row number is wrong %d (expected: %d)" % (_tb.nrows(), len(rowids)))
            for out_row in range(len(rowids)):
                in_row = rowids[out_row]
                sp = _tb.getcell(testcolumn, out_row)
                self.assertEqual(sp.shape[0], len(polids), "Number of pol is wrong in row=%d:  %d (expected: %d)" % (out_row,len(polids),sp.shape[0]))
                nchan = sp.shape[1]
                for out_pol in range(len(polids)):
                    in_pol = polids[out_pol]
                    reference = self._get_reference(nchan, in_row, in_pol, dcol)
                    if self.verbose: print("data=%s" % str(sp[out_pol]))
                    self.assertTrue(numpy.allclose(sp[out_pol], reference,
                                                   atol=atol, rtol=rtol),
                                    "Smoothed spectrum differs in row=%d, pol=%d" % (out_row, out_pol))
        finally:
            _tb.close()


    def testIntentF(self):
        """Test selection by intent (float_data)"""
        self.run_test("sdsmooth_intent", "float_data")

    def testIntentC(self):
        """Test selection by intent (corrected)"""
        self.run_test("sdsmooth_intent", "corrected")

    def testAntennaF(self):
        """Test selection by antenna (float_data)"""
        self.run_test("sdsmooth_antenna", "float_data")

    def testAntennaC(self):
        """Test selection by antenna (corrected)"""
        self.run_test("sdsmooth_antenna", "corrected")

    def testFieldF(self):
        """Test selection by field (float_data)"""
        self.run_test("sdsmooth_field", "float_data")

    def testFieldC(self):
        """Test selection by field (corrected)"""
        self.run_test("sdsmooth_field", "corrected")

    def testSpwF(self):
        """Test selection by spw (float_data)"""
        self.run_test("sdsmooth_spw", "float_data")

    def testSpwC(self):
        """Test selection by spw (corrected)"""
        self.run_test("sdsmooth_spw", "corrected")

    def testTimerangeF(self):
        """Test selection by timerange (float_data)"""
        self.run_test("sdsmooth_timerange", "float_data")

    def testTimerangeC(self):
        """Test selection by timerange (corrected)"""
        self.run_test("sdsmooth_timerange", "corrected")

    def testScanF(self):
        """Test selection by scan (float_data)"""
        self.run_test("sdsmooth_scan", "float_data")

    def testScanC(self):
        """Test selection by scan (corrected)"""
        self.run_test("sdsmooth_scan", "corrected")

    def testPolF(self):
        """Test selection by pol (float_data)"""
        self.run_test("sdsmooth_pol", "float_data")

    def testPolC(self):
        """Test selection by pol (corrected)"""
        self.run_test("sdsmooth_pol", "corrected")

    def testReindexSpw(self):
        """Test reindex =T/F in spw selection"""
        for datacol in ['float_data', 'corrected']:
            print("Test: %s" % datacol.upper())
            for (reindex, ddid, spid) in zip([True, False], [0, 1], [0,7]):
                print("- reindex=%s" % str(reindex))
                self.run_test("sdsmooth_spw", datacol, reindex=reindex)
                _tb.open(self.outfile)
                try:
                    self.assertEqual(ddid, _tb.getcell('DATA_DESC_ID', 0),
                                     "comparison of DATA_DESCRIPTION_ID failed.")
                finally: _tb.close()
                _tb.open(self.outfile+'/DATA_DESCRIPTION')
                try:
                    self.assertEqual(spid, _tb.getcell('SPECTRAL_WINDOW_ID', ddid),
                                     "comparison of SPW_ID failed.")
                finally: _tb.close()
                shutil.rmtree(self.outfile)

    def testReindexIntent(self):
        """Test reindex =T/F in intent selection"""
        for datacol in ['float_data', 'corrected']:
            print("Test: %s" % datacol.upper())
            for (reindex, idx) in zip([True, False], [0, 4]):
                print("- reindex=%s" % str(reindex))
                self.run_test("sdsmooth_intent", datacol, reindex=reindex)
                _tb.open(self.outfile)
                try:
                    self.assertEqual(idx, _tb.getcell('STATE_ID', 0),
                                     "comparison of state_id failed.")
                finally: _tb.close()
                shutil.rmtree(self.outfile)


### sdbaseline ###

### Utilities for reading blparam file
class FileReader(object):
    def __init__(self, filename):
        self.__filename = filename
        self.__data = None
        self.__nline = None

    def read(self):
        if self.__data is None:
            f = open(self.__filename, 'r')
            self.__data = f.readlines()
            f.close()
            self.__nline = len(self.__data)
        return

    def nline(self):
        self.read()
        return self.__nline

    def index(self, txt, start):
        return self.__data[start:].index(txt) + 1 + start

    def getline(self, idx):
        return self.__data[idx]


class BlparamFileParser(FileReader):
    def __init__(self, blfile):
        FileReader.__init__(self, blfile)
        self.__nrow = None
        self.__coeff = None
        self.__rms = None
        self.__ctxt = 'Baseline parameters\n'
        self.__rtxt = 'Results of baseline fit\n'

    def nrow(self):
        self.read()
        if self.__nrow is None:
            return self._nrow()
        else:
            return self.__nrow

    def coeff(self):
        self.read()
        if self.__coeff is None:
            self.parseCoeff()
        return self.__coeff

    def rms(self):
        self.read()
        if self.__rms is None:
            self.parseRms()
        return self.__rms

    def _nrow(self):
        self.__nrow = 0
        for i in range(self.nline()):
            if self.getline(i) == self.__ctxt:
                self.__nrow += 1
        return self.__nrow

    def parse(self):
        self.read()
        self.parseCoeff()
        self.parseRms()
        return
        
    def parseCoeff(self):
        self.__coeff = []
        nrow = self.nrow()
        idx = 0
        while (len(self.__coeff) < nrow):
            try:
                idx = self.index(self.__ctxt, idx)
                coeffs = []
                while(self.getline(idx) != self.__rtxt):
                    coeff = self.__parseCoeff(idx)
                    coeffs += coeff
                    idx += 1
                self.__coeff.append(coeffs)
            except:
                break
        return

    def parseRms(self):
        self.__rms = []
        nrow = self.nrow()
        idx = 0
        while (len(self.__rms) < nrow):
            try:
                idx = self.index(self.__rtxt, idx)
                self.__rms.append(self.__parseRms(idx))
            except:
                break   
        return

    def __parseCoeff(self, idx):
        return parseCoeff(self.getline(idx))

    def __parseRms(self, idx):
        return parseRms(self.getline(idx))

def parseCoeff(txt):
    clist = txt.rstrip('\n').split(',')
    ret = []
    for c in clist:
        ret.append(float(c.split('=')[1]))
    return ret


def parseRms(txt):
    t = txt.lstrip().rstrip('\n')[6:]
    return float(t)


def remove_single_file_dir(filename):
    """
    Remove a single file or a single directory.
    For filename, '.' and those end with '..' (namely, '..', '../..' etc.)
    are not allowed.
    """
    if filename == '.' or filename[-2:] == '..':
        raise Exception("Caution! Attempting to remove '" + filename + "'!!")
    
    if os.path.exists(filename):
        if os.path.isdir(filename):
            shutil.rmtree(filename)
        else: # file or symlink
            os.remove(filename)


def remove_files_dirs(filename):
    """
    Remove files/directories/symlinks 'filename*'.
    For filename, '', '.' and those end with '..' (namely, '..', '../..' etc.)
    are not allowed.
    """
    if filename == '.' or filename[-2:] == '..':
        raise Exception("Caution! Attempting to remove '" + filename + "*'!!")
    elif filename == '':
        raise Exception("The parameter 'filename' must not be a null string.")
    
    import glob
    filenames = glob.glob('{}*'.format(filename.rstrip('/')))

    for filename in filenames:
        remove_single_file_dir(filename)


class sdbaseline_unittest_base(unittest.TestCase):
    """
    Base class for sdbaseline unit test
    """
    # Data path of input/output
    datapath = ctsys_resolve('unittest/sdbaseline/')
    taskname = 'imbaseline_sdbaseline'
    verboselog = False

    #complist = ['max','min','rms','median','stddev']

    blparam_order = ['row', 'pol', 'mask', 'nclip', 'cthre',
                     'uself', 'lthre', 'ledge', 'redge', 'chavg',
                     'btype', 'order', 'npiec', 'nwave']
    blparam_dic = {}
    blparam_dic['row']   = [0, 0, 1, 1, 2, 2, 3, 3]
    blparam_dic['pol']   = [0, 1, 0, 1, 0, 1, 0, 1]
    #blparam_dic['mask']  = ['0~4000;6000~8000']*3 + ['']*5
    blparam_dic['mask']  = ['500~2500;5000~7500']*8
    blparam_dic['nclip'] = [0]*8
    blparam_dic['cthre'] = ['3.']*8
    blparam_dic['uself'] = ['false']*4 + ['true'] + ['false']*3
    blparam_dic['lthre'] = ['0.']*4 + ['3.', '', '', '0.']
    blparam_dic['ledge'] = [0]*4 + [10, 50, '', 0]
    blparam_dic['redge'] = [0]*4 + [10, 50, '', 0]
    blparam_dic['chavg'] = [0]*4 + [4, '', '', 0]
    blparam_dic['btype'] = ['poly'] + ['chebyshev']*2 + ['poly', 'chebyshev', 'poly'] + ['cspline']*2
    blparam_dic['order'] = [0, 0, 1, 1, 2, 2, '', '']
    blparam_dic['npiec'] = [0]*6 + [1]*2
    blparam_dic['nwave'] = [[]]*3 + ['']*2 + [[]]*3

    ### helper functions for tests ###
    def _createBlparamFile(self, file, param_order, val, option=''):
        nspec = 8
        f = open(file, 'w')
        assert(len(param_order) == len(val.keys()))
        for key in val.keys():
            assert(len(val[key]) == nspec)
        for i in range(nspec):
            do_write = True
            s = ''
            for key in param_order:
                v = val[key][i]
                if key == 'nwave':
                    if v != '':
                        s += ','
                        s += str(v)
                else:
                    s += str(v)
                    if key != 'npiec': s += ','
            s += '\n'
            if (option == 'r2p1less') and (val['row'][i] == 2) and (val['pol'][i] == 1):
                do_write = False
            if (option == 'r2p1cout') and (val['row'][i] == 2) and (val['pol'][i] == 1):
                s = '#' + s
            if do_write:
                f.write(s)
        f.close()

    def _checkfile(self, name, fail=True):
        """
        Check if the file exists.
        name : the path and file name to test
        fail : if True, Error if the file does not exists.
               if False, return if the file exists
        """
        isthere = os.path.exists(name)
        if fail:
            self.assertTrue(isthere, msg='Could not find, %s'%(name))
        else:
            return isthere

    def _remove(self, names):
        """
        Remove a list of files and directories from disk
        """
        for name in names:
            remove_single_file_dir(name)

    def _copy(self, names, from_dir=None, dest_dir=None):
        """
        Copy a list of files and directories from a directory (from_dir) to
        another (dest_dir) in the same name.
        
        names : a list of files and directories to copy
        from_dir : a path to directory from which search and copy files
                   and directories (the default is the current path)
        to_dir   : a path to directory to which copy files and directories
                   (the default is the current path)
        NOTE: it is not allowed to specify 
        """
        # Check for paths
        if not from_dir and not dest_dir:
            raise ValueError("The value of from_dir or dest_dir is empty.")
        from_path = os.path.abspath("." if not from_dir else from_dir.rstrip("/"))
        to_path = os.path.abspath("." if not dest_dir else dest_dir.rstrip("/"))
        if from_path == to_path:
            raise ValueError("Can not copy files to exactly the same path.")
        # Copy a list of files and directories
        for name in names:
            from_name = from_path + "/" + name
            to_name = to_path + "/" + name
            if os.path.exists(from_name):
                if os.path.isdir(from_name):
                    shutil.copytree(from_name, to_name)
                else:
                    shutil.copyfile(from_name, to_name)
                if self.verboselog:
                    casalog.post("Copying '%s' FROM %s TO %s" % (name, from_path, to_path))
            else:
                casalog.post("Could not find '%s'...skipping copy" % from_name, 'WARN')
    
    def _getUniqList(self, val):
        """Accepts a python list and returns a list of unique values"""
        if not isinstance(val, list):
            raise Exception('_getUniqList: input value must be a list.')
        return list(set(val))

    def _getListSelection(self, val):
        """
        Converts input to a list of unique integers
        Input: Either comma separated string of IDs, an integer, or a list of values.
        Output: a list of unique integers in input arguments for string and integer input.
                In case the input is a list of values, output will be a list of unique values.
        """
        if isinstance(val, str):
            val_split = val.split(',')
            val_sel = []
            for j in range(len(val_split)):
                val_sel.append(int(val_split[j]))
        elif isinstance(val, int):
            val_sel = [val]
        elif isinstance(val, list) or isinstance(val, tuple):
            val_sel = val.copy()
        else:
            raise Exception('_getListSelection: wrong value ' + str(val) + ' for selection.')
        return self._getUniqList(val_sel)
    
    def _getListSelectedRowID(self, data_list, sel_list):
        """
        Returns IDs of data_list that contains values equal to one in
        sel_list.
        The function is used to get row IDs that corresponds to a
        selected IDs. In that use case, data_list is typically a list
        of values in a column of an MS (e.g., SCAN_NUMBER) and sel_list is
        a list of selected (scan) IDs.

        data_list : a list to test and get IDs from
        sel_list  : a list of values to look for existance in data_list
        """
        res = []
        for i in range(len(data_list)):
            if data_list[i] in sel_list:
                #idx = sel_list.index(data_list[i])
                res.append(i)
        return self._getUniqList(res)
    
    def _getEffective(self, spec, mask):
        """
        Returns an array made by selected elements in spec array.
        Only the elements in the ID range in mask are returned.

        spec : a data array
        mask : a mask list of the channel ranges to use. The format is
               [[start_idx0, end_idx0], [start_idx1, end_idx1], ...]
        """
        res = []
        for i in range(len(mask)):
            for j in range(mask[i][0], mask[i][1]):
                res.append(spec[j])
        return numpy.array(res)

    def _getStats(self, filename=None, spw=None, pol=None, colname=None, mask=None):
        """
        Returns a list of statistics dictionary of selected rows in an MS.

        filename : the name of MS
        spw      : spw ID selection (default: all spws in MS)
        pol      : pol ID selection (default: all pols in MS)
        colname  : the name of data column (default: 'FLOAT_DATA')
        mask     : a mask list of the channel ranges to use. The format is
                   [[start_idx0, end_idx0], [start_idx1, end_idx1], ...]
        
        The order of output list is in the ascending order of selected row IDs.
        The dictionary in output list has keys:
        'row' (row ID in MS), 'pol' (pol ID), 'rms', 'min', 'max', 'median',
        and 'stddev'
        """
        # Get selected row and pol IDs in MS. Also get spectrumn in the MS
        if not spw: spw = ''
        select_spw = (spw not in ['', '*'])
        if select_spw: spw_sel = self._getListSelection(spw)
        if not pol: pol = ''
        select_pol = (pol not in ['', '*'])
        if select_pol: pol_sel = self._getListSelection(pol)
        if not colname: colname='FLOAT_DATA'
        self._checkfile(filename)
        with table_manager(filename) as tb:
            data = tb.getcol(colname)
            ddid = tb.getcol('DATA_DESC_ID')
        with table_manager(filename + '/DATA_DESCRIPTION') as tb:
            spwid = tb.getcol('SPECTRAL_WINDOW_ID').tolist()
        if not select_spw: spw_sel = spwid
        # get the selected DD IDs from selected SPW IDs.
        dd_sel = self._getListSelectedRowID(spwid, spw_sel)
        # get the selected row IDs from selected DD IDs
        row_sel = self._getListSelectedRowID(ddid, dd_sel)
        if not select_spw: row_sel = range(len(ddid))
        if not select_pol: pol_sel = range(len(data))

        res = []
        for irow in row_sel:
            for ipol in pol_sel:
                spec = data[ipol,:,irow]
                res_elem = self._calc_stats_of_array(spec, mask=mask)
                res_elem['row'] = irow
                res_elem['pol'] = ipol
                
                res.append(res_elem)

        return res

    def _calc_stats_of_array(self, data, mask=None):
        """
        """
        if mask is not None:
            spec = self._getEffective(data, mask)
        else:
            spec = numpy.array(data)
        res_elem = {}
        res_elem['rms'] = numpy.sqrt(numpy.var(spec))
        res_elem['min'] = numpy.min(spec)
        res_elem['max'] = numpy.max(spec)
        spec_mea = numpy.mean(spec)
        res_elem['median'] = numpy.median(spec)
        res_elem['stddev'] = numpy.std(spec)
        return res_elem
        

    def _convert_statslist_to_dict(self, stat_list):
        """
        Returns a disctionary of statistics of selected rows in an MS.

        stat_list: a list of stats dictionary (e.g., return value of _getStats)

        The output dictionary is in form:
        {'max': [max0, max1, max2, ...], 'min': [min0, min1,...], ...}
        The order of elements are in ascending order of row and pol IDs pair, i.e.,
        (row0, pol0), (row0, pol1), (row1, pol0), ....
        """
        #if len(stat_list)==0: raise Exception, "No row selected in MS"
        keys=stat_list[0].keys()
        stat_dict={}
        for key in keys:
            stat_dict[key] = []
        for stat in stat_list:
            for key in keys:
                stat_dict[key].append(stat[key])
        return stat_dict

    def _compareStats(self, currstat, refstat, rtol=1.0e-2, atol=1.0e-5, complist=None):
        """
        Compare statistics results (dictionaries) and test if the values are within
        an allowed tolerance.

        currstat : the statistic values to test (either an MS name or
                   a dictionary)
        refstat  : the reference statistics values (a dictionary)
        rtol   : tolerance of relative difference
        atol   : tolerance of absolute difference
        complist : statistics to compare (default: keys in refstat)
        """
        # test if the statistics of baselined spectra are equal to
        # the reference values
        printstat = False #True
        # In case currstat is filename
        if isinstance(currstat, str) and os.path.exists(currstat):
            #print "calculating statistics from '%s'" % currstat
            currstat = self._getStats(currstat)

        self.assertTrue(isinstance(currstat,dict) and \
                        isinstance(refstat, dict),\
                        "Need to specify two dictionaries to compare")
        if complist:
            keylist = complist
        else:
            keylist = refstat.keys()
            #keylist = self.complist
        
        for key in keylist:
            self.assertTrue(key in currstat,\
                            msg="%s is not defined in the current results." % key)
            self.assertTrue(key in refstat,\
                            msg="%s is not defined in the reference data." % key)
            refval = refstat[key]
            currval = currstat[key]
            # Quantum values
            if isinstance(refval, dict):
                if 'unit' in refval and 'unit' in currval:
                    if printstat:
                        print("Comparing unit of '%s': %s (current run), %s (reference)" %\
                              (key,currval['unit'],refval['unit']))
                    self.assertEqual(refval['unit'],currval['unit'],\
                                     "The units of '%s' differs: %s (expected: %s)" % \
                                     (key, currval['unit'], refval['unit']))
                    refval = refval['value']
                    currval = currval['value']
                else:
                    raise Exception("Invalid quantum values. %s (current run) %s (reference)" %\
                                    (str(currval),str(refval)))
            currval = self._to_list(currval)
            refval = self._to_list(refval)
            if printstat:
                print("Comparing '%s': %s (current run), %s (reference)" %\
                      (key,str(currval),str(refval)))
            self.assertTrue(len(currval)==len(refval),"Number of elemnets in '%s' differs." % key)
            if isinstance(refval[0],str):
                for i in range(len(currval)):
                    if isinstance(refval[i],str):
                        self.assertTrue(currval[i]==refval[i],\
                                        msg="%s[%d] differs: %s (expected: %s) " % \
                                        (key, i, str(currval[i]), str(refval[i])))
            else:
                # numpy.allclose handles almost zero case more properly.
                self.assertTrue(numpy.allclose(currval, refval, rtol=rtol, atol=atol),
                                msg="%s differs: %s" % (key, str(currval)))
            del currval, refval

    def _to_list(self, input):
        """
        Convert input to a list
        If input is None, this method simply returns None.
        """
        import numpy
        listtypes = (list, tuple, numpy.ndarray)
        if input == None:
            return None
        elif type(input) in listtypes:
            return list(input)
        else:
            return [input]

    def _compareBLparam(self, out, reference):
        # test if baseline parameters are equal to the reference values
        # currently comparing every lines in the files
        # TO DO: compare only "Fitter range" and "Baseline parameters"
        self._checkfile(out)
        self._checkfile(reference)
        
        blparse_out = BlparamFileParser(out)
        blparse_out.parse()
        coeffs_out = blparse_out.coeff()
        rms_out = blparse_out.rms()
        blparse_ref = BlparamFileParser(reference)
        blparse_ref.parse()
        coeffs_ref = blparse_ref.coeff()
        rms_ref = blparse_ref.rms()
        allowdiff = 0.01
        print('Check baseline parameters:')
        for irow in range(len(rms_out)):
            print('Row %s:'%(irow))
            print('   Reference rms  = %s'%(rms_ref[irow]))
            print('   Calculated rms = %s'%(rms_out[irow]))
            print('   Reference coeffs  = %s'%(coeffs_ref[irow]))
            print('   Calculated coeffs = %s'%(coeffs_out[irow]))
            r0 = rms_ref[irow]
            r1 = rms_out[irow]
            rdiff = (r1 - r0) / r0
            self.assertTrue((abs(rdiff)<allowdiff),
                            msg='row %s: rms is different'%(irow))
            c0 = coeffs_ref[irow]
            c1 = coeffs_out[irow]
            for ic in range(len(c1)):
                rdiff = (c1[ic] - c0[ic]) / c0[ic]
                self.assertTrue((abs(rdiff)<allowdiff),
                                msg='row %s: coefficient for order %s is different'%(irow,ic))
        print('')


class sdbaseline_basic_test(sdbaseline_unittest_base):
    """
    Basic unit tests for task sdbaseline. No interactive testing.

    List of tests:
    test000 --- default values for all parameters
    test001 --- polynominal baselining with no mask (maskmode = 'list'). spw and pol specified.
    test002 --- Chebyshev polynominal baselining with no mask (maskmode = 'list'). spw and pol specified.
    test003 --- cubic spline baselining with no mask (maskmode = 'list'). spw and pol specified.
    test004 --- sinusoidal baselining with no mask (maskmode = 'list'). spw and pol specified.
    test050 --- existing file as outfile with overwrite=False (raises an exception)
    test051 --- no data after selection (raises an exception)
    test060 --- blparam file (infile+'_blparam.txt') should be removed if it exists

    Note: The input data 'OrionS_rawACSmod_calave.ms' is generated
          from a single dish regression data 'OrionS_rawACSmod' as follows:
          
          default(sdcal)
          sdcal(infile='OrionS_rawACSmod',scanlist=[20,21,22,23],
                calmode='ps',tau=0.09,outfile='temp.asap')
          default(sdcal)
          sdcal(infile='temp.asap',timeaverage=True,
                tweight='tintsys',outfile='temp2.asap')
          sdsave(infile='temp2.asap',outformat='MS2',
                 outfile='OrionS_rawACSmod_calave.ms')
    """
    # Input and output names
    infile = 'OrionS_rawACSmod_calave.ms'
    outroot = sdbaseline_unittest_base.taskname+'_basictest'
    blrefroot = os.path.join(sdbaseline_unittest_base.datapath,'refblparam')
    tid = None

    def setUp(self):
        if os.path.exists(self.infile):
            shutil.rmtree(self.infile)
        shutil.copytree(os.path.join(self.datapath,self.infile), self.infile)

        if os.path.exists(self.infile+'_blparam.txt'):
            os.remove(self.infile+ '_blparam.txt')
        if os.path.exists(self.infile+'_blparam.csv'):
            os.remove(self.infile+ '_blparam.csv')
        if os.path.exists(self.infile+'_blparam.btable'):
            shutil.rmtree(self.infile+ '_blparam.btable')

    def tearDown(self):
        remove_files_dirs(self.infile)
        remove_files_dirs(self.outroot)

    def test000(self):
        """Basic Test 000: default values for all parameters"""
        tid = '000'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'float_data'
        
        vals = blv(imagename=infile)
        vals.datacolumn = datacolumn
        vals.sdsmooth_output = infile
        vals.sdbaseline_output = outfile
        do_sdbaseline(vals)
        # uncomment the next line once blparam file can be output
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid)
        row = 3
        pol = 1
        results = self._getStats(outfile, '')
        theresult = None
        for i in range(len(results)):
            if ((results[i]['row'] == int(row)) and (results[i]['pol'] == int(pol))):
                theresult = results[i]
        reference = {'rms': 0.16677055621054496,
                     'min': -2.5817961692810059,
                     'max': 1.3842859268188477,
                     'median': -0.00086212158203125,
                     'stddev': 0.16677055621054496,
                     }
        self._compareStats(theresult, reference)

    def test001(self):
        """Basic Test 001: simple successful case: blfunc = 'poly', maskmode = 'list' and masklist=[] (no mask)"""
        tid = '001'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'float_data'
        maskmode = 'list'
        blfunc = 'poly'
        spw = '3'
        pol = 'LL'
        vals = blv(imagename=infile, maskmode=maskmode, blfunc=blfunc)
        vals.datacolumn = datacolumn
        vals.sdbaseline_spw = spw
        vals.sdbaseline_pol = pol
        vals.sdsmooth_output = infile
        vals.sdbaseline_output = outfile
        do_sdbaseline(vals)
        # uncomment the next line once blparam file can be output
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid)
        results = self._getStats(outfile)
        print(self._getStats(outfile))
        theresult = None
        for i in range(len(results)):
            theresult = results[i]

        reference = {'rms': 0.16677055621054496,
                     'min': -2.5817961692810059,
                     'max': 1.3842859268188477,
                     'median': -0.00086212158203125,
                     'stddev': 0.16677055621054496,
                     }

        self._compareStats(theresult, reference)

    def test001_uppercase_params(self):
        """Basic Test 001: simple successful case: blfunc = 'poly', maskmode = 'list' and masklist=[] (no mask)"""
        tid = '001'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'FLOAT_DATA'
        maskmode = 'LIST'
        blfunc = 'POLY'
        spw = '3'
        pol = 'LL'
        vals = blv(imagename=infile, maskmode=maskmode, blfunc=blfunc)
        vals.datacolumn = datacolumn
        vals.sdbaseline_spw = spw
        vals.sdbaseline_pol = pol
        vals.sdsmooth_output = infile
        vals.sdbaseline_output = outfile
        do_sdbaseline(vals)
        # uncomment the next line once blparam file can be output
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid)
        results = self._getStats(outfile)
        print(self._getStats(outfile))
        theresult = None
        for i in range(len(results)):
            theresult = results[i]

        reference = {'rms': 0.16677055621054496,
                     'min': -2.5817961692810059,
                     'max': 1.3842859268188477,
                     'median': -0.00086212158203125,
                     'stddev': 0.16677055621054496,
                     }

        self._compareStats(theresult, reference)

    def test002(self):
        """Basic Test 002: simple successful case: blfunc = 'chebyshev', maskmode = 'list' and masklist=[] (no mask)"""
        tid = '002'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'float_data'
        maskmode = 'list'
        blfunc = 'chebyshev'
        spw = '3'
        pol = 'LL'
        vals = blv(imagename=infile, maskmode=maskmode, blfunc=blfunc)
        vals.datacolumn = datacolumn
        vals.sdbaseline_spw = spw
        vals.sdbaseline_pol = pol
        vals.sdsmooth_output = infile
        vals.sdbaseline_output = outfile
        do_sdbaseline(vals)

        # uncomment the next line once blparam file can be output
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid)
        results = self._getStats(outfile)
        print(self._getStats(outfile))
        theresult = None
        for i in range(len(results)):
            theresult = results[i]

        reference = {'rms': 0.16677055621054496,
                     'min': -2.5817961692810059,
                     'max': 1.3842859268188477,
                     'median': -0.00086212158203125,
                     'stddev': 0.16677055621054496,
                     }

        self._compareStats(theresult, reference)
    
    
    def test003(self):
        """Basic Test 003: simple successful case: blfunc = 'cspline', maskmode = 'list' and masklist=[] (no mask)"""
        print("")

        tid = '003'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        datacolumn = 'float_data'  
        maskmode = 'list'
        blfunc = 'cspline'
        npiece = 3
        spw='3'
        pol='LL'
        vals = blv(imagename=infile, maskmode=maskmode, blfunc=blfunc)
        vals.datacolumn = datacolumn
        vals.sdbaseline_spw = spw
        vals.sdbaseline_pol = pol
        vals.sdbaseline_npiece = npiece
        vals.sdsmooth_output = infile
        vals.sdbaseline_output = outfile
        do_sdbaseline(vals)
        
        # uncomment the next line once blparam file can be output
        #self._compareBLparam(outfile+"_blparam.txt",self.blrefroot+tid) 
        results = self._getStats(outfile)
        print(self._getStats(outfile))
        
        theresult = None
        for i in range(len(results)):
            theresult = results[i]

        reference = {'rms': 0.16685959517745799,
                     'min': -2.5928177833557129,
                     'max': 1.3953156471252441,
                     'median': -0.00089824199676513672,
                     'stddev': 0.16685959517745766,
                    }

        self._compareStats(theresult, reference)

        #***
        #*** check if baseline is subtracted ***
        #***
        # Output MS only has the selected pol, LL
        in_pol=1
        out_pol=0

        # open the original MS
        _tb.open(infile)
        orig_pol1_value = numpy.array(_tb.getcell('FLOAT_DATA', int(spw))[in_pol,:])
        _tb.close()
        variance_orig_pol1 = numpy.var(orig_pol1_value)
        
        # open the MS after sdbaseline
        _tb.open(outfile)
        pol1_value = numpy.array(_tb.getcell('FLOAT_DATA', 0)[out_pol,:])
        _tb.close()
        variance_pol1 = numpy.var(pol1_value)

        # assert pol1_value < orig_pol1_value
        self.assertTrue((pol1_value<orig_pol1_value).all())
        
        # assert variance of pol1_value < variance of orig_pol1_value
        self.assertLess(variance_pol1**0.5, variance_orig_pol1**0.5)

        #print '1sigma before cspline (pol1)', variance_orig_pol1**0.5 
        #print '1sigma after cspline (pol1)',  variance_pol1**0.5 

    def test051(self):
        """Basic Test 051: failure case: no data after selection"""
        tid = '051'
        infile = self.infile
        outfile = self.outroot+tid+'.ms'
        spw = '10' # non-existent IF value
        mode = 'list'
        try:
            vals = blv(imagename=infile, maskmode=mode)
            vals.sdsmooth_output = infile
            vals.sdbaseline_output = outfile
            vals.sdbaseline_spw = spw
            do_sdbaseline(vals)
        except Exception as e:
            self.assertIn('Spw Expression: No match found for 10,', str(e))



def suite():
    return [imsmooth_test, 
            sdsmooth_test_fail, sdsmooth_test_complex, sdsmooth_test_float,
            sdsmooth_test_weight, sdsmooth_test_boxcar, sdsmooth_selection,
            sdbaseline_basic_test]    

if __name__ == '__main__':
    unittest.main()