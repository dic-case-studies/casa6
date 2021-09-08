import random
import os
import numpy
import re
import shutil
import unittest
import math
from scipy import signal

from casatools import ctsys, image, regionmanager, componentlist, table, quanta ,ms
from casatasks import imbaseline, casalog, imsubimage
from casatasks.private.task_imbaseline import ImBaselineVals as blv, do_imsmooth, do_sdsmooth
from casatasks.private.sdutil import table_manager

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

    def run_test(self, *args, **kwargs):
        datacol_name = self.datacolumn.upper()
        weight_mode = hasattr(self, 'weight_propagation') and getattr(self, 'weight_propagation') is True

        if 'kwidth' in kwargs:
            kwidth = kwargs['kwidth']
        else:
            kwidth = 5

        vals = blv(imagename=self.infile, spkernel='gaussian', kwidth=kwidth)
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
    test_sdsmooth_fail04 --- outfile exists (overwrite=False)
    test_sdsmooth_fail05 --- empty outfile
    test_sdsmooth_fail06 --- invalid data column name
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
    test_sdsmooth_complex_select --- data selection (spw)
    test_sdsmooth_complex_overwrite --- overwrite existing outfile (overwrite=True)
    """
    exception_case = sdsmooth_test_base.exception_case
    infile = sdsmooth_test_base.infile_data
    datacolumn = 'data'

    @exception_case(RuntimeError, 'Desired column \(FLOAT_DATA\) not found in the input MS')
    def test_sdsmooth_complex_fail01(self):
        """test_sdsmooth_complex_fail01 --- non-existing data column (FLOAT_DATA)"""
        vals = blv(imagename=self.infile, spkernel='gaussian')
        vals.datacolumn = datacolumn='float_data'
        self.exec_sdsmooth(self.infile, self.outfile, vals)

    def test_sdsmooth_complex_gauss01(self):
        """test_sdsmooth_complex_gauss01 --- gaussian smoothing (kwidth 5)"""
        self.run_test(kwidth=5)

    def test_sdsmooth_complex_gauss02(self):
        """test_sdsmooth_complex_gauss02 --- gaussian smoothing (kwidth 3)"""
        self.run_test(kwidth=3)

def suite():
    return [imsmooth_test, sdsmooth_test_fail, sdsmooth_test_complex]    

if __name__ == '__main__':
    os.chdir("/work/dev/shimada/casa6.13520/tmp")
    unittest.main()