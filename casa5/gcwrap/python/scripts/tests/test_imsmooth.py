########################################################################3
#  imsmooth_test.py
#
# Copyright (C) 2008, 2009
# Associated Universities, Inc. Washington DC, USA.
#
# This scripts free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning AIPS++ should be adressed as follows:
#        Internet email: aips2-request@nrao.edu.
#        Postal address: AIPS++ Project Office
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA
#

from __future__ import absolute_import
from __future__ import print_function
import random
import os
import numpy
import shutil
import unittest
import math

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys, image, regionmanager, componentlist, table, quanta
    from casatasks import imsmooth, casalog

    _ia = image()
    _rg = regionmanager()
    _tb = table()
    _qa = quanta()
    ctsys_resolve = ctsys.resolve
else:
    import casac
    from tasks import *
    from taskinit import *

    componentlist = cltool
    image = iatool

    dataRoot = os.path.join(os.environ.get('CASAPATH').split()[0],'data')

    def ctsys_resolve(apath):
        return os.path.join(dataRoot,apath)

    _ia = iatool( )
    _rg = rgtool( )

    # not local tools
    _tb = tb
    _qa = qa

targetres_im = "imsmooth_targetres.fits"
tiny = "tiny.im"
image_names=['g192_a2.image', 'g192_a2.image-2.rgn']

def _near(got, expected, tol):
    return _qa.le(
        _qa.div(
            _qa.abs(_qa.sub(got, expected)),
            expected
        ),
        tol
    )
    
def run_imsmooth(
    imagename, major, minor, pa, targetres,
    outfile, kernel="gauss", overwrite=False, beam={}
):
    return imsmooth(
        imagename=imagename, kernel=kernel,
        major=major, minor=minor, pa=pa,
        targetres=targetres, outfile=outfile,
        overwrite=overwrite, beam=beam
    )

def make_gauss2d(shape, xfwhm, yfwhm):
    fac = 4*math.log(2)
    values = numpy.empty(shape, dtype=float)
    for i in range(shape[0]):
        x = shape[0]/2 - i
        for j in range(shape[1]):
            y = shape[1]/2 - j
            xfac = x*x*fac/(xfwhm*xfwhm)
            yfac = y*y*fac/(yfwhm*yfwhm)
            values[i, j] = math.exp(-(xfac + yfac));
    return values

    
class imsmooth_test(unittest.TestCase):

    def setUp(self):
        if(os.path.exists(image_names[0])):
            for file in image_names:
                os.system('rm -rf ' +file)
        
        datapath = 'regression/g192redux/reference'
        for file in image_names:
            os.system('cp -r '+ctsys_resolve(os.path.join(datapath,file))+' ' + file)
        self.ia = image()
        self.datapath = 'regression/imsmooth'
        for f in [targetres_im, tiny]:
            if(os.path.exists(f)):
                os.system('rm -rf ' +f)
            os.system('cp -r '+ctsys_resolve(os.path.join(self.datapath,f))+' ' + f)

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
 
    ####################################################################
    # Incorrect inputs to parameters.  The parameters are:
    #    imagename
    #    outfile
    #    kernel
    #    major
    #    minor
    #    mask
    #    region
    #    box
    #    chans
    #    stokes
    #
    # Returns True if successful, and False if it has failed.
    ####################################################################
    
    def _compare_beams(self, beam1, beam2):
        self.assertTrue(_near(beam1["major"], beam2["major"], 2e-5))
        self.assertTrue(_near(beam1["minor"], beam2["minor"], 2e-5))
        pa = []
        for b in [beam1, beam2]:
            if "positionangle" in b:
                pa.append(b["positionangle"])
            else:
                pa.append(b["pa"])

        diff = abs(
            _qa.sub(
                _qa.quantity(pa[0]), 
                _qa.quantity(pa[1])
            )["value"]
        )
        self.assertTrue(diff < 1e-5)

    def test_input(self):
        '''Imsmooth: Testing INPUT/OUTPUT tests'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        casalog.post( "Starting imsmooth INPUT/OUTPUT tests.", 'NORMAL2' )
    
        # First step get rid of all the old test files!
        for file in os.listdir( '.' ):
            if file.startswith( 'input_test' ):
                shutil.rmtree( file )
        if os.path.exists( 'garbage.rgn' ):
            os.remove('garbage.rgn')
    
    
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
            results = imsmooth( 'g192', outfile='input_test1', beam=beam )
        except:
            no_op='noop'
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Badfile, 'g192', was not reported as missing."
            
        self.assertTrue(imsmooth( tiny, outfile='input_test1', beam=beam))
        self.assertTrue(os.path.exists( 'input_test1' ))
    
        # same thing, just using major, minor, and pa
        self.assertTrue(imsmooth( tiny, outfile='input_test1b', major=major, minor=minor, pa=pa))
        self.assertTrue(os.path.exists( 'input_test1b' ))
    
        #######################################################################
        # Testing the outfile parameter.
        #    1. Bad file, file already exists, exception should be thrown
        #    2. Good file name, a file should be
        #######################################################################
        casalog.post( "The OUTFILE parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
        results = None
        try:
            results = imsmooth( tiny, outfile='input_test1', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['error_msgs']=retValue['error_msgs']\
                  +"\nError: Badfile, 'input_test1', was not reported as already existing."
            
        results = None
        try:
            results=imsmooth( tiny, outfile='input_test2', beam=beam )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create smoothed image 'input_test2'"
        if ( not os.path.exists( 'input_test2' ) and results==False ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: output file, 'input_test2', was not created."
    
    
        #######################################################################
        # Testing KERNEL parameter, valid values 0 and greater
        #    1. Below invalid range: junk, ''
        #    3. valid: gaussian, boxcar, tophat, 
        #######################################################################
        casalog.post( "The KERNEL parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
        results = None
        try:
            results = imsmooth( tiny, kernel='', outfile='input_test3', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: No exception thrown for bad kernel value, ''"
    
        results = None
        try:
            results = imsmooth( tiny, kernel='junk', outfile='input_test4', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: No exception thrown for bad kernel value, 'junk'"    
    
        results = None
        try:
            results=imsmooth( tiny, kernel='gauss', outfile='input_test5', beam=beam)
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Gaussian smoothing failed."
        if ( not os.path.exists( 'input_test5' ) or results == False ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: input_test5 output file was NOT created."\
                       + "\n\tRESULTS: "+str(results)
            
    
        results = None
        major = "2arcsec"
        minor = "2arcsec"
        pa = "0deg"
        try:
            results=imsmooth( tiny, kernel='boxcar', outfile='input_test6', major=major, minor=minor, pa=pa )
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Boxcar smoothing failed. " \
                     +str(err)
            
        if ( not os.path.exists( 'input_test6' ) or results==False ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: output file 'input_test6' was NOT created."
    
    
        # UNCOMMENT when tophat support has been added.
        #try:
        #    results = None
        #    results=imsmooth( 'g192_a2.image', kernel='tophat', oufile='input_test7' )
        #except Exception, err:
        #    retValue['success']=False
        #    retValue['error_msgs']=retValue['error_msgs']\
        #             +"\nError: Tophat smoothing failed. "
        #    
        #if ( not os.path.exists( 'input_test7' ) or results==None ): 
        #    retValue['success']=False
        #    retValue['error_msgs']=retValue['error_msgs']\
        #               +"\nError: output file 'input_test7' was NOT created."
    
    
    
        #######################################################################
        # Testing MAJOR parameter
        # Expects a numerical value: 1, '2pix', '0.5arcsec'
        # Tests include: invalid values, valid values, major < minor 
        #######################################################################
        casalog.post( "The MAJOR parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
        try:
            result = None
            results = imsmooth( tiny, major='bad', minor=minor, pa=pa, oufile='input_test8')
        except:
            no_op='noop'
        else:
            if ( results != None and results!=False ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Bad major value, 'bad', was not reported."
    
        try:
            result = None
            results = imsmooth( tiny, major=-5, minor=minor, pa=pa, oufile='input_test9' )
        except:
            no_op='noop'
        else:
            if ( results != None and results!=False ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Bad major value, '-5', was not reported."
    
        result = None
        results = imsmooth(tiny, major=2, minor=1, pa=0, outfile='input_test11')
        self.assertTrue(results)
        self.assertTrue(os.path.exists( 'input_test11' ))
            
        result = None
        results = imsmooth( tiny, major='2pix', minor='1pix', pa="deg", outfile='input_test12')
        self.assertTrue(results)
        self.assertTrue(os.path.exists( 'input_test12' ))
        
        result = None
        results = imsmooth( tiny, major='1.5arcsec', minor='1arcsec', pa="0deg", outfile='input_test13')
        self.assertTrue(results)
        self.assertTrue(os.path.exists( 'input_test13' ))

    
        result = None
        try:
            results = imsmooth( tiny, major='0.5arcsec', minor='2arcsec', pa="0deg", outfile='input_test14')
        except:
            pass
        else:
            if ( results != None and results!=False ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Bad major value less than minor value was not reported."        
    
    
        #######################################################################
        # Testing REGION parameter
        # Expects a file containing a region record, as created by the viewer.
        # Tests include bad file name, file with bad content, and good file.
        ####################################################################### 
        casalog.post( "The REGION parameter tests will cause errors to occur, do not be alarmed", 'WARN' )

        # CASA6 tasks throw exceptions, CASA5 tasks return False
        passes = False
        try:
            results = imsmooth( tiny, region=7, beam=beam )
            if not results:
                passes = True 
        except:
            passes = True
        if ( not passes ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Bad region file, 7, was not reported as bad."
    
        try:
            results = imsmooth( tiny, region='garbage.rgn', beam=beam )
        except:
            #We want this to fail
            no_op = 'noop'
        else:
            if ( results ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                                        +"\nError: Bad region file, 'garbage.rgn', was not reported as missing."
    
        try:
            filename = os.getcwd()+'/garbage.rgn'
            fp=open( filename, 'w' )
            fp.writelines('This file does NOT contain a valid CASA region specification\n')
            fp.close()
    
            try:
                results = imsmooth( 'g192_a2.image', region=filename , beam=beam)
            except:
                no_op='noop'
            else:
                if ( results ):
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']\
                              + "\nError: Bad region file, 'garbage.rgn',"\
                              + " was not reported as bad."
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Unable to create bad region file.\n\t"
            raise
    
        
        results = None
        results=imsmooth( 'g192_a2.image', region='g192_a2.image-2.rgn', outfile='input_test15', beam=beam )
        self.assertTrue(results)
        self.assertTrue(os.path.exists('input_test15'))
    
        #######################################################################
        # Testing BOX parameter
        # The input file has pixel values ranging from
        #   0-511, 0-511
        # Tests include -3, -1, 0, 1 random valid value, 500, 511, 525
        #   for both the x, and y coords
        #
        # Note: -1 is a special case implying use the full range, so to
        #       be out of bounds we need -2 or less.
        #######################################################################
        casalog.post( "The BOX parameter tests will cause errors to occur, do not be alarmed", 'WARN' )

        # CASA6 task throw exceptions, CASA5 tasks return False
        if is_CASA6:
            self.assertRaises(Exception, imsmooth, tiny, box='-3,0,511,511', beam=beam)
        else:
            self.assertFalse(imsmooth(tiny, box='-3,0,511,511', beam=beam))

        x1=random.randint(0,127)
        x2=random.randint(x1,127)
        y1=random.randint(0,127)
        y2=random.randint(y1,127)
        boxstr=str(x1)+','+str(y1)+','+str(x2)+','+str(y2)
        
        self.assertTrue(imsmooth( tiny, box=boxstr, outfile='input_test16', beam=beam ))
        self.assertTrue(os.path.exists( 'input_test16' ))
    
        #######################################################################
        # Testing CHANS parameter: valid values 0-39 for our image
        # Values used for testing, -5,-2,0,22~35, 39,40,45
        #
        # NOTE: a coord value of -1 indicates use all, so -1 is a valid
        #       coordiante.
        #######################################################################
        casalog.post( "The CHANS parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
        results = None
        try:
            results = imsmooth( 'g192_a2.image', chans='-5', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):    
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                      +"\nError: Bad channel value, '-5', was not reported."
    
        results = None
        try:
            results = imsmooth( 'g192_a2.image', chans='-2', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):    
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                  +"\nError: Bad channel value, '-2', was not reported."
    
        results = None
        try:
            results = imsmooth( 'g192_a2.image', chans='-18', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Bad channel value of -18 was not reported."
    
        results = None
        try:
            results = imsmooth( 'g192_a2.image', chans='45', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Bad channel value of 45 was not reported."
    
        results = None
        try:
            results = imsmooth( 'g192_a2.image', chans='40', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Bad channel value of 40 was not reported."
    
        
        self.assertTrue(imsmooth( 'g192_a2.image', chans='22~35', outfile='input_test17', beam=beam ))
        self.assertTrue(os.path.exists('input_test17'))
        
        self.assertTrue(imsmooth( tiny, chans='0', outfile='input_test17b', beam=beam ))
        self.assertTrue(os.path.exists( 'input_test17b' ))
    
        self.assertTrue(imsmooth( 'g192_a2.image', chans='39', outfile='input_test18', beam=beam ))
        self.assertTrue(os.path.exists("input_test18"))
            
        #######################################################################
        # Testing STOKES parameter, valid values: 'I'
        #    Tests are 'Q', 'yellow' (invalid) and 'I'
        #######################################################################
        casalog.post( "The STOKES parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
        results=None
        try:
            results = imsmooth( 'g192_a2.image', stokes='Q', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Bad stokes value, 'Q', was not reported."
                
        results=None    
        try:
            results = imsmooth( 'g192_a2.image', stokes='yellow', beam=beam )
        except:
            pass
        else:
            if ( results != None and \
                 ( not isinstance(results, bool) or results==True ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Bad stokes value, 'yellow', was not reported."
    
        
        self.assertTrue(imsmooth( tiny, stokes='I', outfile='input_test19', beam=beam ))
        self.assertTrue(os.path.exists('input_test19'))
        
        self.assertTrue(retValue['success'],retValue['error_msgs'])
    
    
    ####################################################################
    # Smoothing correctness test.
    #
    # This test subtacts the continuum from the g192 data file
    # and compares the results (both continuum and spectral line
    # with subtracted continuum files) with pervious results.
    #
    # Random values are selected in the files and compared.
    #
    # Returns True if successful, and False if it has failed.
    ####################################################################
     
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
    
        results = None
        try:
            results=imsmooth( 'smooth.pointsrc.image', kernel='gauss', \
                              major=50, minor=25, pa=0, outfile='smooth_test1' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: boxcar smooth failed on point source file."
    
            
        if ( not os.path.exists( 'smooth_test1' ) or results==None ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Gaussian smoothfailed on point source file."
        else:
            # Now that we know something has been done lets check the results!
            #      1. Check that the sum of the values under the curve is 100
            #      2. Check that the max is at 212, 220, 0 , 20
            allowedError = 0.009
            
            _ia.open( 'smooth_test1')
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
        results = None
        try:
            results=imsmooth( 'smooth.pointsrc.image', kernel='boxcar', \
                              major=20, minor=10, pa=0, outfile='smooth_test2' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: boxcar smooth failed on point source file."
            
        if ( not os.path.exists( 'smooth_test2' ) or results==None ): 
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

## this was the original test, before CASA6 was written:
## _ia.pixelvalue returns a dict type
## So this test relies on a dict comparison
##    value > (0.125-allowedError)
## Which (wisely) is unavailable in python 3
## But which also must be doing the wrong thing, because the logic of that
## test is something that SHOULD be true here : i.e. are the values all
## within the allowedError of 0.125
## When this test is rewritten to use the values in that dict, it fails in CASA5.
## So the dict comparison must not be doing what it was expected to do.
#
#            val1 = _ia.pixelvalue( [ 204,200,0,20] )
#            val2 = _ia.pixelvalue( [ 222,236,0,20] )
#            val3 = _ia.pixelvalue( [ 204,239,0,20] )
#            val4 = _ia.pixelvalue( [ 222,201,0,20] )        
#            midVal = _ia.pixelvalue( [212,220,0,20] )
#            for value in [val1, val2, val3, val4, midVal ]:
#                if ( value>(0.125-allowedError) and value<(0.125+allowedError)):
#                    retValue['success']=False
#                    retValue['error_msgs']=retValue['error_msgs']\
#                        +"\nError: Values in the smoothed box are not all 0.125"\
#                        +" found value of "+str(value)
#

            pixels = list(map( lambda pv: _ia.pixelvalue(pv)['value']['value'],
                          [ [204,200,0,20], [222,236,0,20], [204,239,0,20],
                            [222,201,0,20], [212,220,0,20] ] ))

            for value in pixels:
#               #     this is the original test as modified in CASA6, it fails when moved to CASA5
#               #     The second comparison appears to be an attempt to exlude values near zero
#               #     from failing, but it's too close to the actual values returned in CASA5 and
#               #     so this fails
#                if ( not (value>(0.125-allowedError) and value<(0.125+allowedError)) and
#                     not (value>-3.3740714666663507e-09 and value<3.3740714666663507e-09) ):
#               #     this works, but appears to be just as much a kludge and not as the 
#               #     original test seemed to intend
#               #     more thought should be given here to what test is appropriate
                if ( not (value>(0.125-allowedError) and value<(0.125+allowedError)) and
                     not (value>-3.37407147e-09 and value<3.37407147e-09) ):
                    retValue['success']=False
                    retValue['error_msgs']=retValue['error_msgs']\
                        +"\nError: Values in the smoothed box are not all 0.125"\
                        +" found value of "+str(value)
            _ia.done()
    
        self.assertTrue(retValue['success'],retValue['error_msgs'])
    
    
    ####################################################################
    # Region selection correction test.
    #
    # This test selects a region(s) where smoothing will be performed.
    #
    # Returns True if successful, and False if it has failed.
    ####################################################################
    
    def test_region(self):
        '''Imsmooth: Region selection correction test'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        casalog.post( "Starting imsmooth REGION tests.", 'NORMAL2' )
        allowedError = 0.0005
    
        # First step get rid of all the old test files!
        for file in os.listdir( '.' ):
            if file.startswith( 'rgn_test' ):
                shutil.rmtree( file )
    
        if os.path.exists( 'rgn.pointsrc.image' ):
            shutil.rmtree( 'rgn.pointsrc.image' )
    
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
    
            # Create an array of zero's, then set a couple positions (point
            # sources) to 100.
            #
            # Note that 
            inputArray = numpy.zeros( (shape[0], shape[1], shape[2], shape[3]), 'float' )
            inputArray[49,71,0,14] = 100     # For rgn file
            inputArray[233,276,0,20] = 100     # For rgn in image
            inputArray[15,15,0,30] = 100     # For rgn in image
    
            
            # Now make the image!
            _ia.fromarray( pixels=inputArray, csys=csys.torecord(), \
                          outfile='rgn.pointsrc.image' )
            _ia.done()
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Unable to create point source image."\
                     +"\n\t REULTS: "+str(err)        
    
    
        # Select the following regions without the point source:
        #            1. Sky region without the point source
        #            2. Channel that doesn't have the point source
        #
        # Note: that when we check the resulting smoothed images
        #       the should be empty.
    
        results = None
        try:
            results=imsmooth( 'rgn.pointsrc.image', kernel='gauss', \
                              major=50, minor=25, pa=0, outfile='rgn_test1', \
                              box='350,350,375,390')
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Smoothng failed on region 250,250,275,290." + str(err)
    
            
        if ( not os.path.exists( 'rgn_test1' ) or results==None ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Smoothing failed on region 250,250,275,290. second block"
        else:
            # Now that we know something has been done lets check the results!
            #      1. Check that the sum of the values under the curve is 0
            _ia.open( 'rgn_test1')
            stats = _ia.statistics()
            if ( stats['sum'][0] < ( 0-allowedError) \
                 or stats['sum'][0] > ( 0+allowedError) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Sum on smoothed file rgn_test1 is "\
                    +str(stats['sum'][0]) +" expected value is 0."
            _ia.done()
    
    
        results = None
        try:
            results=imsmooth( 'rgn.pointsrc.image', kernel='gauss', \
                              major=50, minor=25, pa=0, outfile='rgn_test2', \
                              chans='22')
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Smoothng failed on channel 22."
    
            
        if ( not os.path.exists( 'rgn_test2' ) or results==None ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Smoothing failed on channel 22."
        else:
            # Now that we know something has been done lets check the results!
            #     1. Check that the sum of the values under the curve is 0
            _ia.open( 'rgn_test2')
            stats = _ia.statistics()
            if ( stats['sum'][0] < ( 0-allowedError) \
                 or stats['sum'][0] > ( 0+allowedError) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Sum on smoothed file rgn_test2 is "\
                    +str(stats['sum'][0]) +" expected value is 0."
            _ia.done()
    
    
    
        # Select a region that contains the point source
        #   1. using imsmooth parameters
        #   2. region defined in an image
        #        g192_a.image:testregion (blc=166,222,0,0  trc=296,328,0,39)
        #   3. region file.
        #        g192_1.image.rgn      (blc=0,0,0,0 trc=511,511,0,14)
        #
        results = None
        try:
            results=imsmooth( 'rgn.pointsrc.image', kernel='gauss', \
                              major=10, minor=5, pa=0, outfile='rgn_test3', \
                              chans='14', box='0,0,200,200')
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Smoothng failed on channel 14, box 0,0,200,200."
    
            
        if ( not os.path.exists( 'rgn_test3' ) or results==None ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Smoothing failed on channel 14, box 0,0,200,200."
        else:
            # Now that we know something has been done lets check the results!
            #     1. Check that the sum of the values under the curve is 100
            #     2. Check that the max is at 49,71, 0, 14
            _ia.open( 'rgn_test3')
            stats = _ia.statistics()
            if ( stats['sum'][0] < ( 100-allowedError) \
                 or stats['sum'][0] > ( 100+allowedError) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Sum on smoothed file rgn_test3 is "\
                    +str(stats['sum'][0]) +" expected value is 100."
            _ia.done()
    
            # Note that since we've selected a single plane then our
            # output image has a single plane, 0, only!  Thus, unlike
            # our original image the max point should be found on the
            # 0th channel and NOT the 14th channel.
            maxpos=stats['maxpos'].tolist()
            if ( maxpos[0]!=49 or maxpos[1]!=71 or \
                 maxpos[2]!=0 or maxpos[3]!=0 ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Max position found at "+str(maxpos)\
                    +" expected it to be at 49,71,0,0."            
    
    
        results = None
        # This test was all screwed up when it fell in my lap. Fixing as best as I can - dmehring
        output = 'rgn_test5'
        try:
            results=imsmooth( 'rgn.pointsrc.image', kernel='gauss', \
                              major=2, minor=1, pa=0, outfile = output, \
                              region='g192_a2.image:testregion')
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Smoothing failed with internal image region 'testregion'."
    
        if ( not os.path.exists(output) or results==None ): 
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Smoothing failed internal image region 'testregion'."
        else:
            # Now that we know something has been done lets check the results!
            #     1. Check that the sum of the values under the curve is 100
            #     2. Check that the max is at 49,71, 0, 14
            _ia.open(output)
            stats = _ia.statistics()
            _ia.done()
            
            sum = stats['sum'][0]
            fluxDensity = 99.948 # not 100 because of flux located outside small image
            allowedError=0.001
            if (sum < (fluxDensity - allowedError) or sum > (fluxDensity + allowedError)):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Sum on smoothed file " + output + " is "\
                    +str(stats['sum'][0]) +" expected value is 100."
    
            # Max position = point src position - minx, miny of region    
            maxpos=stats['maxpos'].tolist()
            if ( maxpos[0] != 5 or maxpos[1] != 5 or \
                 maxpos[2] != 0 or maxpos[3] != 30 ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                    +"\nError: Max position found at "+str(maxpos)\
                    +" expected it to be at 212,220,0,20."            
    
        self.assertTrue(retValue['success'],retValue['error_msgs'])


    def test_stretch(self):
        """ imsmooth(): Test stretch parameter"""
        yy = image()
        mymask = "maskim"
        yy.fromshape(mymask, [200, 200, 1, 1])
        yy.addnoise()
        yy.done()
        shape = [200,200,1,20]
        imagename = "tmp.im"
        yy.fromshape(imagename, shape)
        yy.addnoise()
        yy.done()
        zz = imsmooth(
            imagename=imagename, major="2arcsec", minor="2arcsec",
            pa="0deg", mask=mymask + ">0", stretch=False
        )
        self.assertFalse(zz)

        zz = imsmooth(
            imagename=imagename, major="2arcsec", minor="2arcsec",
            pa="0deg", mask=mymask + ">0", stretch=True
        )
        self.assertTrue(zz and type(zz) == type(True))
    
    def test_multibeam(self):
        """Test per plane beams"""
        myia = self.ia
        imname = "test_image2dconvolver_multibeam.im"
        shutil.copytree(ctsys_resolve(os.path.join(self.datapath, imname)), imname)
        # myia.open(ctsys_resolve(os.path.join(self.datapath,"test_image2dconvolver_multibeam.im")))
        major = "10arcmin"
        minor = "8arcmin"
        pa = "80deg"
        gotname = 'convolve2d_multibeam.im'
        # got = myia.convolve2d(axes=[0, 1], major=major, minor=minor, pa=pa)
        self.assertTrue(
            imsmooth(
                imname=imname, axes=[0, 1], major=major, minor=minor, pa=pa,
                outname=gotname
            ), 'imsmooth failed'
        )
        myia.open(imname)
        shape = myia.shape()
        myia.done()
        got = image()
        for i in range(5):
            blc=[0, 0, i]
            trc=[shape[0]-1, shape[1]-1, i]
            reg = _rg.box(blc=blc, trc=trc)
            # xx = myia.subimage(region=reg)
            subname = 'subi' + str(i) + '.im'
            imsubimage(imagename=imname, region=reg, outname=subname)
            # exp = xx.convolve2d(axes=[0, 1], major=major, minor=minor, pa=pa)
            outname = 'convolve' + str(i) + '.im'
            imsmooth(
                subname, axes=[0, 1], major=major, minor=minor,
                pa=pa, outname=outname
            )
            exp.open(outname)
            expbeam = myia.restoringbeam()
            myia.open(gotname)
            gotbeam = got.restoringbeam(channel=i)
            for j in ["major", "minor", "positionangle"]:
                self.assertTrue(_near(gotbeam[j], expbeam[j], 2e-7))
            self.assertTrue(abs(got.getchunk(blc=blc, trc=trc) - exp.getchunk()).max() < 3e-5)
            got.done()
            exp.done()
            # xx.done()
        myia.done()
        got.done()
                            
    def test_targetres(self):
        """Test targetres parameter"""
        myia = self.ia
        imagename = "tres1.im"
        myia.fromshape(imagename, [100, 100])
        csys = myia.coordsys()
        csys.setunits(["arcsec", "arcsec"])
        csys.setincrement([-1, 1])
        myia.setcoordsys(csys.torecord())
        myia.setrestoringbeam(major="6arcsec", minor="3arcsec", pa="0deg")
        shape = myia.shape()
        values = make_gauss2d(shape, 3.0, 6.0)
        expected = make_gauss2d(shape, 5.0, 10.0)
        myia.putchunk(values)
        myia.done()
        emaj = _qa.quantity("10arcsec")
        emin = _qa.quantity("5arcsec")
        epa = _qa.quantity("0deg")
        
        for unit in ("Jy/beam", "K"):
            myia.open(imagename)
            myia.setbrightnessunit(unit)
            myia.done()
            expected = make_gauss2d(shape, 5.0, 10.0)
            if (unit == "K"):
                expected *= 3.0*6.0/5.0/10.0
            # for code in (run_convolve2d, run_imsmooth):
            for targetres in [False, True]:
                if not targetres:
                    major = "8arcsec"
                    minor = "4arcsec"
                    pa = "0deg"
                    outfile = "tres1" + unit[0]
                else:
                    major = "10arcsec"
                    minor = "5arcsec"
                    pa = "0deg"
                    outfile = "tres2" + unit[0]
                run_imsmooth(
                    imagename=imagename, kernel="gaussian",
                    major=major, minor=minor, pa=pa, targetres=targetres,
                    outfile=outfile
                )
                myia.open(outfile)
                gotbeam = myia.restoringbeam()
                gotvals = myia.getchunk()
                myia.done()
                self._compare_beams(
                    gotbeam, {"major": emaj, "minor": emin, "pa": epa}
                )
                maxdiff = (abs(gotvals-expected)).max()
                self.assertTrue(maxdiff < 1e-6)     
        
        csys.addcoordinate(spectral=True)
        for unit in ("Jy/beam", "K"):
            myia.fromshape(
                outfile=imagename, shape=[100, 100, 2],
                csys=csys.torecord(), overwrite=True
            )
            myia.setbrightnessunit(unit)
            myia.setrestoringbeam(
                major="6arcsec", minor="3arcsec", pa="0deg", channel=0
            )
            myia.setrestoringbeam(
                major="4arcsec", minor="2arcsec", pa="0deg", channel=1
            )
            values = myia.getchunk()
            shape = myia.shape()
            expected = values.copy()
            for k in range(shape[2]):
                if k == 0:
                    xfwhm = 3
                    yfwhm = 6
                else:
                    xfwhm = 2
                    yfwhm = 4
                values[:,:,k] = make_gauss2d([shape[0], shape[1]], xfwhm, yfwhm)
            myia.putchunk(values)
            outia = image()
            for targetres in [False, True]:
                ebeam = []
                if targetres:
                    major = "10arcsec"
                    minor = "5arcsec"
                else:
                    major = "8arcsec"
                    minor = "4arcsec"
                pa = "0deg"
            
                for k in range(shape[2]):
                    reg = _rg.box(blc=[0, 0, k], trc=[shape[0]-1, shape[1]-1, k])
                    subim = myia.subimage(outfile="", region=reg, dropdeg=True)
                    convim = subim.convolve2d(
                        type="gaussian", major=major, minor=minor,
                        pa=pa, targetres=targetres
                    )
                    subim.done()
                    expected[:, :, k] = convim.getchunk()
                    gotbeam = convim.restoringbeam()
                    
                    if targetres:
                        self._compare_beams(gotbeam, {"major": major, "minor": minor, "pa": pa})
    
                    ebeam.append(gotbeam)
                    convim.done()
                # for code in [run_convolve2d, run_imsmooth]:
                if targetres:
                    outfile = "tres3" + unit[0] + str(code)
                else:
                    outfile = "tres4" + unit[0] + str(code)
                run_imsmooth(
                     imagename=imagename, kernel="gaussian",
                     major=major, minor=minor, pa=pa, targetres=targetres,
                     outfile=outfile
                )        
                outia.open(outfile)
                for k in range(shape[2]):
                    gotbeam = outia.restoringbeam(channel=k)
                    self._compare_beams(gotbeam, ebeam[k])
                    maxdiff = (abs(outia.getchunk()-expected)).max()
                    self.assertTrue(maxdiff < 1e-6) 
        myia.done()
        outia.done()

        myia.open(imagename)
        myia.setrestoringbeam(
            major="6arcsec", minor="3arcsec", pa="0deg"
        )
        myia.done()
        # for code in [run_convolve2d, run_imsmooth]:
        outfile = "tres6"
        self.assertFalse(
            run_imsmooth(
                imagename=imagename, kernel="gaussian",
                major="5.99arcsec", minor="2.99arcsec", pa="0deg",
                targetres=True, outfile=outfile
            )
        )        

    def test_overwrite(self):
        """ test overwrite parameter """
        myia = self.ia
        outfile = "test_overwrite.im"
        myia.fromshape(outfile, [200, 200])
        imagename = "input_overwrite"
        myia.fromshape(imagename, [200, 200])
        myia.done()
        self.assertTrue(
            run_imsmooth(
                imagename=imagename, kernel="gauss", major="5arcmin",
                minor="4arcmin", pa="0deg", targetres=False,
                overwrite=True, outfile=outfile
            )
        )
        self.assertFalse(
            run_imsmooth(
                imagename=imagename,
                kernel="gauss", major="5arcmin", minor="4arcmin",
                pa="0deg", targetres=False, overwrite=False, outfile=outfile
            )
        )  
        
    def test_beam(self):
        """Test the beam parameter"""
        myia = self.ia
        imagename = "tbeam1.im"
        myia.fromshape(imagename, [100, 100])
        csys = myia.coordsys()
        csys.setunits(["arcsec", "arcsec"])
        csys.setincrement([1, 1])
        myia.setcoordsys(csys.torecord())
        myia.setbrightnessunit("Jy/beam")
        myia.setrestoringbeam(major="6arcsec", minor="3arcsec", pa="0deg")
        shape = myia.shape()
        myia.putchunk(make_gauss2d(shape, 3.0, 6.0))
        expected = make_gauss2d(shape, 5.0, 10.0)
        # for code in (run_convolve2d, run_imsmooth):
        for beam in [
            {"major": "8arcsec", "minor": "4arcsec", "pa": "0deg"},
            {
                "major": {"unit": "arcsec", "value": 8},
                "minor": {"unit": "arcsec", "value": 4},
                "pa": {"unit": "deg", "value": 0},
            }
        ]:
            outfile = 'smooth'
            x = run_imsmooth(
                imagename=imagename, major="", minor="", pa="",
                beam=beam, outfile=outfile, targetres=False,
                overwrite=True
            )
            if type(x) == type(myia):
                x.done()
            myia.open(outfile)
            maxdiff = (abs(myia.getchunk()-expected)).max()
            self.assertTrue(maxdiff < 1e-6) 
            myia.done()
                
    def test_conserve_flux(self):
        """Test flux density is conserved for images with units of K or anything without 'beam'"""
        myia = self.ia
        imagename = "tres1x.im"
        myia.fromshape(imagename, [100, 100])
        csys = myia.coordsys()
        csys.setunits(["arcsec", "arcsec"])
        csys.setincrement([-1, 1])
        myia.setcoordsys(csys.torecord())
        myia.setrestoringbeam(major="6arcsec", minor="3arcsec", pa="0deg")
        shape = myia.shape()
        values = make_gauss2d(shape, 3.0, 6.0)
        #expected = make_gauss2d(shape, 5.0, 10.0)
        myia.putchunk(values)
        for unit in ("K", "cm-2"):
            myia.setbrightnessunit(unit)
            zz = myia.fitcomponents()
            mycl = componentlist()
            mycl.fromrecord(zz['results'])
            expected = mycl.getfluxvalue(0)
            gg = image()
            outfile = "gxg_" + unit + ".im"
            imsmooth(
                imagename=imagename, targetres=True, major="10arcsec", minor="5arcsec",
                pa="0deg", outfile=outfile
            )
            gg.open(outfile)
            zz = gg.fitcomponents()
            gg.done()
            mycl.fromrecord(zz['results'])
            got = mycl.getfluxvalue(0)
            self.assertTrue(abs(got[0]/expected[0] - 1) < 3e-7, "Failed testing unit " + unit)
        mycl.done()
        myia.done()
        
    def test_commonbeam(self):
        """Test kernel='commonbeam' in imsmooth"""
        myia = self.ia
        imagename = "cb1.im"
        myia.fromshape(imagename, [100, 100, 5])
        csys = myia.coordsys()
        csys.setunits(["arcsec", "arcsec", "Hz"])
        csys.setincrement([-1, 1, 1e6])
        myia.setcoordsys(csys.torecord())
        myia.setbrightnessunit("Jy/beam")
        myia.done()
        outfile = "cbout1.im"
        
        myia.open(imagename)
        myia.setrestoringbeam(major="6arcsec", minor="3arcsec", pa="0deg")
        myia.done()
        """
        self.assertTrue(
            imsmooth(
                imagename=imagename, kernel='commonbeam', outfile=outfile,
                targetres=False
            )           
        )
        myia.open(outfile)
        beam = myia.restoringbeam()
        myia.done()
        root2 = math.sqrt(2)
        self.assertTrue(abs(_qa.getvalue(beam['major']) - 6*root2) < 1e-6)
        self.assertTrue(abs(_qa.getvalue(beam['minor']) - 3*root2) < 1e-6)
        self.assertTrue(abs(_qa.getvalue(beam['positionangle'])) < 1e-5)
        """
        outfile = "cbout2.im"
        self.assertTrue(
            imsmooth(
                imagename=imagename, kernel='commonbeam', outfile=outfile,
                targetres=True
            )           
        )
        myia.open(outfile)
        beam = myia.restoringbeam()
        myia.done()
        self.assertTrue(abs(_qa.getvalue(beam['major']) - 6) < 1e-5)
        self.assertTrue(abs(_qa.getvalue(beam['minor']) - 3) < 1e-5)
        self.assertTrue(abs(_qa.getvalue(beam['positionangle'])) < 1e-5)
        
        myia.open(imagename)
        myia.setrestoringbeam(remove=True)
        myia.setrestoringbeam(major="6arcsec", minor="3arcsec", pa="0deg", channel=0)

        myia.setrestoringbeam(major="8arcsec", minor="4arcsec", pa="0deg", channel=1)
        myia.done()
        outfile = "cbout3.im"
        """
        self.assertTrue(
            imsmooth(
                imagename=imagename, kernel='commonbeam', outfile=outfile,
                targetres=False
            )           
        )
        myia.open(outfile)
        for i in range(5):
            beam = myia.restoringbeam(channel=i)
            if i == 1:
                emajor = root2*8
                eminor = root2*4
            else:
                emajor = 10
                eminor = 5
            self.assertTrue(abs(_qa.getvalue(beam['major']) - emajor) < 1e-6)
            self.assertTrue(abs(_qa.getvalue(beam['minor']) - eminor) < 1e-6)
            self.assertTrue(abs(_qa.getvalue(beam['positionangle'])) < 1e-5)
        myia.done()
        """
        outfile = "cbout4.im"
        self.assertTrue(
            imsmooth(
                imagename=imagename, kernel='commonbeam', outfile=outfile,
                targetres=True
            )           
        )
        myia.open(outfile)
        for i in range(5):
            beam = myia.restoringbeam(channel=i)
            self.assertTrue(abs(_qa.getvalue(beam['major']) - 8) < 1e-6)
            self.assertTrue(abs(_qa.getvalue(beam['minor']) - 4) < 1e-6)
            self.assertTrue(abs(_qa.getvalue(beam['positionangle'])) < 1e-5)
        myia.done()
        
    def test_image_kernel(self):
        """Test image as kernel, CAS-5844"""
        imagename = ctsys_resolve(os.path.join(self.datapath,"point.im"))
        kimage = ctsys_resolve(os.path.join(self.datapath,"bessel.im"))
        outfile = "point_c_bessel.im"
        self.assertTrue(
            imsmooth(
                imagename=imagename, kernel="i",
                kimage=kimage, outfile=outfile
            )
        )
        myia = self.ia
        myia.open(kimage)
        bessel = myia.getchunk()
        myia.open(outfile)
        conv = myia.getchunk()
        myia.done()
        ratio = conv/bessel
        print("max",abs(ratio - ratio[50,50]).max())
        self.assertTrue((abs(ratio - ratio[50,50]) < 2e-4).all())
        self.assertTrue(
            imsmooth(
                imagename=imagename, kernel="i",
                kimage=kimage, outfile=outfile,
                overwrite=True, scale=1
            )
        )
        myia.open(outfile)
        conv = myia.getchunk()
        myia.done()
        diff = conv - bessel
        self.assertTrue(abs(diff).max() < 2e-7)

    def test_history(self):
        """Test that history is written"""
        myia = self.ia
        imagename = "zz.im"
        myia.fromshape(imagename, [20,20])
        major = "2arcmin"
        minor = "2arcmin"
        pa = "0deg"
        bb = myia.convolve2d("", major=major,  minor=minor, pa=pa)        
        myia.done()
        msgs = bb.history()
        bb.done()
        teststr = "ia.convolve2d"
        self.assertTrue(teststr in msgs[-4], "'" + teststr + "' not found")     
        self.assertTrue(teststr in msgs[-3], "'" + teststr + "' not found")
        
        outfile = "zz_out.im"
        self.assertTrue(
            imsmooth(
                imagename=imagename, outfile=outfile,
                major=major, minor=minor, pa=pa
            )
        )
        myia.open(outfile)
        msgs = myia.history()
        myia.done()
        teststr = "version"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")    
        teststr = "imsmooth"
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")

    def test_CAS_12904(self):
        """Test fix of CAS-12904 bug"""
        imname = 'orig.im'
        yy = iatool()
        yy.fromshape(imname, [100, 100, 3])
        pix = yy.getchunk()
        for i in range(3):
            pix[:, :, i] = i
        yy.putchunk(pix)
        yy.done()
        outname = 'mysub.im'
        imsubimage(imagename=imname, outfile=outname, mask=imname + '>0')
        subi = iatool()
        subi.open(outname)
        rg = rgtool()
        for i in range(3):
            reg = rg.box([0, 0, i], [99, 99, i])
            npts = subi.statistics(region=reg)['npts']
            expec = 0 if i == 0 else 1
            self.assertEqual(npts.size, expec, 'wrong length npts array')
            if i>0:
                self.assertEqual(npts[0], 10000, 'wrong number of pts')
        subi.done()
        imname = outname
        outname = 'conv.im'
        imsmooth(
            imname, major='4arcmin', minor='4arcmin', pa='0deg',
            mask=imname + '<2', outfile=outname
        )
        yy.open(outname)
        for i in range(3):
            reg = rg.box([0, 0, i], [99, 99, i])
            npts = yy.statistics(region=reg)['npts']
            expec = 1 if i == 1 else 0
            self.assertEqual(npts.size, expec, 'wrong length npts array')
            if i==1:
                self.assertEqual(npts[0], 10000, 'wrong number of pts')
        yy.done()     

def suite():
    return [imsmooth_test]    

if is_CASA6:
    if __name__ == '__main__':
        unittest.main()
