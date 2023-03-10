#########################################################################
# test_task_immoments.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.immoments.html
#
##########################################################################
import random
import sys
import os
import shutil
import numpy
import math
import unittest
from casatools import ctsys, regionmanager, image, table
from casatasks import casalog, imstat, imhistory, imval, immath, immoments

_rg = regionmanager( )
_ia = image( )

    #
    # To make things a little more interesting, I've used the
    # lyrics to song for debug msgs, guess which song and you win
    # a prize! If only I could get it to play at the same time!.
    #
    # To turn them off just set debugMsgs=False
    
debug_msgs={}
debug_msgs[0]= "This time I wonder what it feels like"
debug_msgs[1]= "To find the one in this life"
debug_msgs[2]= "The one we all dream of"
debug_msgs[3]= "but dreams just aren't enough"
debug_msgs[4]= "So I'll be waiting for the real thing"
debug_msgs[5]= "I'll know it by the feeling"
debug_msgs[6]= "The moment when we're meeting"
debug_msgs[7]= "Will play out like a scene"
debug_msgs[8]= "right off the silver screen."
debug_msgs[9]= "So I'll be holdin' my breath"
debug_msgs[10]= 'Until that moment when'
debug_msgs[11]= "I find the one that I'll spend forever with!"

debug_msgs[12]= 'Cause nobody wants to be the last one there'
debug_msgs[13]= "'Cause everyone wants to feel like somone cares"
debug_msgs[14]= 'Somone to love with my life in their hands'
debug_msgs[15]= "There's gotta be somebody for me like that"

debug_msgs[16]= " 'Cause nobody wants to go it on their own"
debug_msgs[17]= "And everyone wants to know they're not alone"
debug_msgs[18]= "There's somebody else that feels the same somewhere"
debug_msgs[19]= "There's gotta be somebody for me out there"

debug_msgs[20]= "Tonight, out on the street out in the moonlight"
debug_msgs[21]= "And dammit this feels too right"
debug_msgs[22]= "It's just like Deja Vu"
debug_msgs[23]= "Me standin$B!G(B here with you"
debug_msgs[24]= "So I'll be holdin`my own breath"
debug_msgs[25]= "Could this be the end?"
debug_msgs[26]= "Is it that moment when"
debug_msgs[27]= "I find the one that I'll spend forever with?"

debug_msgs[28]= " 'Cause nobody wants to be the last one there"
debug_msgs[29]= " 'Cause everyone wants to feel like someone cares."
debug_msgs[30]= "Someone to love with my life in their hands."
debug_msgs[31]= "There's gotta be somebody for me like that."

debug_msgs[32]= "'Cause nobody wants to do it on their own"
debug_msgs[33]= "And everyone wants to know theyN4re not alone."
debug_msgs[34]= "There's somebody else that feels the same somewhere"
debug_msgs[35]= "There`s gotta be somebody for me out there."

debug_msgs[36]= "You can't give up!"
debug_msgs[37]= "Lookin' for that diamond in the rough"
debug_msgs[38]= "You never know but when it shows up"
debug_msgs[39]= "Make sure you're holdin` on"
debug_msgs[40]= "'Cause it could be the one, the one you're waiting on"

debug_msgs[41]= "'Cause nobody wants to be the last one there."
debug_msgs[42]= "And everyone wants to feel like someone cares."
debug_msgs[43]= "Someone to love with my life in their hands."
debug_msgs[44]= "There has gotta be somebody for me"
debug_msgs[45]= "Ohhhhhh."

debug_msgs[46]= "Nobody wants to do it on their own"
debug_msgs[47]= "And everyone wants to know they're not alone."
debug_msgs[48]= "Is there somebody else that feels the same somewhere?"
debug_msgs[49]= "There`s gotta be somebody for me out there."

debug_msgs[50]= "Nobody wants to be the last one there"
debug_msgs[51]= "'Cause everyone wants to feel like someone cares."
debug_msgs[52]= "Is there somebody else that feels the same somewhere?"
debug_msgs[53]= "There has gotta be somebody for me out there."

debugMsgs = False
def _momentTest_debug_msg( msgNum=0 ):
    if ( not debugMsgs ):
        return
    idx = msgNum % 54
    print(str(msgNum)+":  "+debug_msgs[idx])
    return

datapath = ctsys.resolve('unittest/immoments/')

# input files
list1=['n1333_both.image','n1333_both.image.rgn']
list2=['n1333_both.image','n1333_both.src.tmom0.all','n1333_both.image.rgn', 
       'immoment_image', 'first_moment.im']
#list=['n1333_both.image', 'n1333_both.src.tmom0.blu', 'n1333_both.src.tmom0.red', 
#       'n1333_both.src.tmom0.all', 'n1333_both.src.tmom1.all', 'n1333_both.image.rgn', 
#       'immoment_image', 'first_moment.im']

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

####################################################################
# Incorrect inputs to parameters.  The parameters are:
#    imagename
#    moments
#    axis
#    region
#    box
#    chans
#    stokes
#    mask
#    includepix
#    excludepix
#    outfile
#
# Returns True if successful, and False if it has failed.
####################################################################
class immoment_test1(unittest.TestCase):

    garbagergn = 'garbage.rgn'
    
    def setUp(self):
        self.tb = table( )
        if(os.path.exists(list1[0])):
            for file in list1:
                os.system('rm -rf ' +file)
        
        for file in list1:
            os.system('cp -RL ' +os.path.join(datapath,file)+' ' + file)


    def tearDown(self):
        for file in list1:
            os.system('rm -rf ' +file)
            os.system('rm -rf input_test*')
            os.system('rm -rf moment_test*')
        os.system('rm -rf ' + self.garbagergn)
        self.assertTrue(len(self.tb.showcache()) == 0)
        self.tb.done( )
        
    def test_input(self):
        '''Immoment: Test input/output parameters'''
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
        casalog.post( "Starting immoments INPUT/OUTPUT tests.", 'NORMAL2' )
    
    
        #######################################################################
        # Testing the imagename parameter.
        #    1. Bad file name should throw and exception
        #    2. Good file name, a file should be
        #######################################################################
#        _momentTest_debug_msg( 5 )
        results = None
        success = False
        try:
            results = immoments( 'n1333_both', moments=[0], outfile='input_test_1' )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Badfile, 'n1333_both', was not reported as bad."
            
#        _momentTest_debug_msg( 6 )
        try:
            immoments( 'n1333_both.image', moments=[0], outfile='input_test_1' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment 0 on n133_both.image"
        if ( not os.path.exists( 'input_test_1' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Moment file, 'input_test_1', was not created."
            
        #######################################################################
        # Testing MOMENTS parameter, valid values are -1 to 11
        #    1. Below valid range: -2, and -10
        #    2. Above valid range: 12 and 21
        #    3. Within range: -1,5, 11
        #######################################################################
        casalog.post( "The moment parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
#        _momentTest_debug_msg( 7 )
        results=None
        try:
            results= immoments( 'n1333_both.image', moments=[-2], outfile='moment_test' )
        except:
            no_op='noop'
        else:
            if ( results != None ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: No exception thrown for bad moment value (-2)"
                 
#        _momentTest_debug_msg( 8 )
        results=None    
        try:
            results=immoments( 'n1333_both.image', moments=[-10], outfile='moment_test' )
        except:
            no_op='noop'
        else:
            if ( results != None and results!=True ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: No exception thrown for bad moment value (-10)"
                
#        _momentTest_debug_msg( 9 )
        results=None    
        try:
            results=immoments( 'n1333_both.image', moments=[12], outfile='moment_test' )
        except:
            no_op='noop'
        else:
            if ( results != None and results!=True ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: No exception thrown for bad moment value (12)"\
                     +"\n\t REULTS: "+str(results)
                
                
    
#        _momentTest_debug_msg( 10 )
        results=None
        try:
            results=immoments( 'n1333_both.image', moments=[21], outfile='moment_test' )
        except:
            no_op='noop'
        else:
            if ( results != None ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: No exception thrown for bad moment value (21)  "+str(type(results))
    
    
        # The remaining tests should succeed.
#        _momentTest_debug_msg( 11 )
        try:    
            immoments( 'n1333_both.image', moments=[-1], axis='spec', outfile='moment_test_2_1' )
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment -1 on n133_both.image\n"\
                       +str(err)
        if ( not os.path.exists( 'moment_test_2_1' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Moment file, 'moment_test_2_1', was not created."
            
#        _momentTest_debug_msg( 12 )
        self.assertTrue(len(self.tb.showcache()) == 0)

        try:    
            immoments( 'n1333_both.image', moments=[5], axis='spec', outfile='moment_test_2_5' )
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment 5 on n133_both.image\n"\
                       +str(err)
        if ( not os.path.exists( 'moment_test_2_5' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Moment file, 'moment_test_2_5', was not created."

#        _momentTest_debug_msg( 13 )
        try:    
            immoments( 'n1333_both.image', moments=[11], axis='spec', outfile='moment_test_2_11' )
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment 11 on n133_both.image\n"\
                       +str(err)
        if ( not os.path.exists( 'moment_test_2_11' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Moment file, 'moment_test_2_11', was not created."
            
        #######################################################################
        # Testing AXIS parameter, valid values are spec, stokes, ra, dec,
        #    as well as 0,1,2,3 ... (depending on the number of axes)
        #######################################################################
#        _momentTest_debug_msg( 14 )
        try:    
            immoments( 'n1333_both.image', moments=[0], axis='ra', outfile='input_test_axis_ra' )
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment 0 on axis ra on n133_both.image\n"\
                       +str(err)
        if ( not os.path.exists( 'input_test_axis_ra' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Moment file, 'input_test_axis_ra', was not created."

#        _momentTest_debug_msg( 15 )
        try:
            immoments( 'n1333_both.image', moments=[0], axis=1, outfile='input_test_axis_dec' )
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       + "\nError: Unable to create moment 0 on axis spec on n133_both.image\n"\
                       +str(err)
    
        if ( not os.path.exists( 'input_test_axis_dec' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Moment file, 'input_test_axis_dec', was not created."

#        _momentTest_debug_msg( 16 )
        try:    
            immoments( 'n1333_both.image', moments=[0], axis='spec', outfile='input_test_axis_spec' )
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment 0 on axis spec on n133_both.image\n"\
                       +str(err)
        if ( not os.path.exists( 'input_test_axis_spec' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Moment file, 'input_test_axis_spec', was not created."

#        _momentTest_debug_msg( 17 )
        try:    
            immoments( 'n1333_both.image', moments=[0], axis='stokes', outfile='input_test_axis_stokes' )
        except Exception as err:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment 0 on axis stokes on n133_both.image\n"\
                       +str(err)
    #    if ( not os.path.exists( 'input_test_axis_stokes' ) or not isinstance( results, object) ):
    #        retValue['success']=False
    #        retValue['error_msgs']=retValue['error_msgs']\
    #                   +"\nError: Moment file, 'input_test_axis_stokes', was not created."

        casalog.post( "The axis parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
#        _momentTest_debug_msg( 18 )
        try:
            immoments( 'n1333_both.image', moments=[0], axis=-1, outfile='input_test_bad_axis' )
        except Exception:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Bad axis value, '-1', was not reported as bad."

#        _momentTest_debug_msg( 19 )
        success = False
        try:
            results = immoments( 'n1333_both.image', moments=[0], axis=4, outfile='input_test_bad_axis' )
        except:
            success = True
        self.assertTrue(success)
        success = False
        try:
            results = immoments( 'n1333_both.image', moments=[0], axis='whatever', outfile='input_test_bad_axis_unique' )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Bad axis value, 'whatever', was not reported as bad."
    
        #######################################################################
        # Testing REGION parameter
        # Expects a file containing a region record, as created by the viewer.
        # Tests include bad file name, file with bad content, and good file.
        ####################################################################### 
        casalog.post( "The axis parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
#        _momentTest_debug_msg( 21 )
        self.assertTrue(len(self.tb.showcache()) == 0)

        success = False
        try:
            immoments( 'n1333_both.image', region=3, outfile='input_test_bad_rgn' )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Bad region file, 3, was not reported as bad."
    
#        _momentTest_debug_msg( 22 )

        success = False
        try:
            immoments(
                'n1333_both.image', region=self.garbagergn,
                outfile='input_test_bad_rgn'
            )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +f"\nError: Bad region file, '{self.garbagergn}', was not reported as missing."
    
#        _momentTest_debug_msg( 23 )
        with open(self.garbagergn,'w') as fp:
            fp.write('This file does NOT contain a valid CASA region specification')

        success = False
        try:
            immoments(
                'n1333_both.image', region=self.garbagergn,
                outfile='input_test_bad_rgn'
            )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +f"\nError: Bad region file, '{self.garbagergn}', was not reported as bad."
#        _momentTest_debug_msg( 24 )

        try:
            immoments( 'n1333_both.image', moments=[0], region='n1333_both.image.rgn', outfile='input_test_rgn_1' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment 0 on n133_both.image with region file n133_both.image.rgn"

    
        #######################################################################
        # Testing BOX parameter
        # The input file has pixel values ranging from
        #   0-799, 0-799
        # Tests include -3, -1, 0, 1 random valid value, 799, 800, 820
        #   for both the x, and y coords
        #######################################################################
        casalog.post( "The BOX parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
#        _momentTest_debug_msg( 34 )
        results = None
        success = False
        try:
            results = immoments( 'n1333_both.image', box='-3,0,799,799', outfile='input_test_bad_box' )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Bad box value, 'x=-3', was not reported as missing."
    
#        _momentTest_debug_msg( 26 )
        results = None
        success = False
        try:
            results = immoments( 'n1333_both.image', box='0,-3,799,799', outfile='input_test_bad_box' )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Bad box value, 'y=-3', was not reported as missing."
    
    
#        _momentTest_debug_msg( 27 )
        self.assertTrue(len(self.tb.showcache()) == 0)

        results = None
        success = False
        try:
            results = immoments( 'n1333_both.image', box='-2,0,798,798', outfile='input_test_bad_box' )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Bad box value, 'x=-2', was not reported."
    
#        _momentTest_debug_msg( 28 )
        results = None
        success = False
        try:
            results = immoments( 'n1333_both.image', box='0,-2,799,799', outfile='input_test_bad_box' )
        except:
            success = True
        if ( not success ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                 +"\nError: Bad box value, 'y=-2', was not reported as missing."
    
#        _momentTest_debug_msg( 29 )
        self.assertRaises(
            Exception, immoments, 'n1333_both.image', box='0,0,800,799',
            outfile='input_test_bad_box'
        )
    
        self.assertTrue(len(self.tb.showcache()) == 0)

   
        x1=random.randint(0,799)
        x2=random.randint(x1,799)
        y1=random.randint(0,799)
        y2=random.randint(y1,799)
        boxstr=str(x1)+','+str(y1)+','+str(x2)+','+str(y2)
#        _momentTest_debug_msg( 33 )
        self.assertTrue(len(self.tb.showcache()) == 0)

        try:
            immoments( 'n1333_both.image', box=boxstr, outfile='input_test_box_1' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment with box="+boxstr

        self.assertTrue(len(self.tb.showcache()) == 0)

        #######################################################################
        # Testing CHANS parameter: valid values 0-17 for our image
        # Values used for testing, -5,-1,0,2~5, 17,18,32
        #######################################################################
        casalog.post( "The CHANS parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        
        self.assertRaises(
            Exception, immoments, 'n1333_both.image', chans='-5',
            outfile='input_test_bad_chans'
        )
        self.assertTrue(len(self.tb.showcache()) == 0)
        self.assertRaises(
            Exception, immoments, 'n1333_both.image', chans='18',
            outfile='input_test_bad_chans'
        )
        self.assertTrue(len(self.tb.showcache()) == 0)

        try:
            immoments( 'n1333_both.image', chans='0~17', outfile='input_test_chans_1' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment wit chans=0~17"
        self.assertTrue(len(self.tb.showcache()) == 0)

        try:
            results = immoments( 'n1333_both.image', chans='2~5', outfile='input_test_chans_2' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment wit chans=2~5"
        if ( not os.path.exists( 'input_test_chans_2' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Moment file, 'input_test_chan_2', was not created."
    
            
        #######################################################################
        # Testing STOKES parameter, valid values: 'I'
        #    Tests are 'Q', 'yellow' (invalid) and 'I'
        #######################################################################
        casalog.post( "The STOKES parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
        self.assertRaises(
            Exception, immoments, 'n1333_both.image', stokes='Q',
            outfile='input_test_bad_stokes'
        )
        self.assertTrue(len(self.tb.showcache()) == 0)
        results = None
        try:
            results = immoments( 'n1333_both.image', stokes='I', outfile='input_test_stokes_1' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment with stokes=Q"
        if ( not os.path.exists( 'input_test_stokes_1' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Moment file, 'input_test_stokes_1', was not created."\
                     +"\n RESULTS: "+str(results)

        self.assertTrue(len(self.tb.showcache()) == 0)       
    
        #######################################################################
        # Testing MASK parameter:
        # There are many, many valid mask values.
        # See AIPS++ Note 223
        #    http://www.astron.nl/aips++/docs/notes/223/223.html
        # The test file n1333 has a mask alread, 'mask0' we will use this
        # for testing.
        # Invalid tests: 'blarg', 'n1333_both.image:mask1', 'bad_file:mask0'
        # Valid tests: n1333_both.image:nomask, n133_both.image:mask0, '<0.5'
        #######################################################################
        casalog.post( "The MASK parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
    
        self.assertRaises(
            Exception, immoments, 'n1333_both.image', mask='blarg',
            outfile='input_test_bad_mask'
        )
        self.assertTrue(len(self.tb.showcache()) == 0)
        try:
            immoments( 'n1333_both.image', mask='mask(n1333_both.image:mask0)', outfile='input_test_mask_2' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment with mask='mask0'"

    
        #######################################################################
        # Testing includepix parameter: Valid values are a vector of two or
        #    1 float values, indicating a single value or a range of values.
        #    values in our test image range from -0.04201228 to 0.04867625
        # Invalid tests: ['bad']
        # Valid test:    [-0.1,0.1], [-6,-5] (this is valid, but creates an
        #                                     array of NaNs)
        #######################################################################
        casalog.post( "The INCLUDEPIX parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
    
        self.assertRaises(
            Exception, immoments, 'n1333_both.image', includepix='bad',
            outfile='input_test_bad_incpix'
        )
        immoments( 'n1333_both.image', includepix=[-0.1,0.1], outfile='input_test_incpix_2' )
        self.assertTrue(os.path.exists('input_test_incpix_2'))
        #######################################################################
        # Testing excludepix parameter: Valid values are a vector of two or
        #    1 float values, indicating a single value or a range of values.
        #    values in our test image range from -0.04201228 to 0.04867625
        # Invalid tests: ['badpix']
        # Valid test:    [0.1,-0.1], [6,5] (this is valid, but creates an
        #                                     array of NaNs)
        #######################################################################
        casalog.post( "The EXCLUDEPIX parameter tests will cause errors to occur, do not be alarmed", 'WARN' )
    
        self.assertRaises(
            Exception, immoments, 'n1333_both.image', excludepix='badpix',
            outfile='input_test_bad_expix'
        )
        immoments( 'n1333_both.image', excludepix=[0.1,-0.1], outfile='input_test_expix_2')
        self.assertTrue(os.path.exists('input_test_expix_2'))
        #######################################################################
        # Testing outfile parameter
        # Expects a file containing a region record, as created by the viewer.
        # Tests include bad file name (one writing in a read only aread and
        # one rewriting on an existing file., and good file.
        #######################################################################
        self.assertTrue(len(self.tb.showcache()) == 0)
        immoments( 'n1333_both.image', outfile='input_test_outfile_1')
        self.assertTrue(os.path.exists('input_test_outfile_1'))
        self.assertRaises(
            Exception, immoments, 'n1333_both.image', outfile='input_test_outfile_1'
        )
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        self.assertTrue(len(self.tb.showcache()) == 0)

class immoment_test2(unittest.TestCase):    

    myim = 'myim.im'
    mymask = 'mask.im'
    mymask1 = 'mask1.im'
    mymask2 = 'mask2.im'
    mymask3 = 'mask3.im'
    cas5287 = 'cas5287.im'
    myexp = 'expected.im'
    expected_integrated = myexp + '.integrated'
    expected_median = myexp + '.median'
    mygot = 'got.im'
    temporary = 'Temporary_Image.'
    temporary_integrated = temporary + 'integrated'
    temporary_median = temporary + 'median'
    cas8570 = 'cas8570.im'
    zz = 'zz.im'
    zz_out = 'zz_out.im'
    zz_out_wc = zz_out + '.weighted_coord'
    zz_out_avg = zz_out + '.average'
    zz_out_wdc = zz_out + '.weighted_dispersion_coord'
    median_coord = 'median_coord.im'
    median_coord_out = 'median_coord_out.im'
    myim2 = 'myim2.im'
    coord = 'coord.im'
    coord_max = coord + '.max'
    coord_min = coord + '.min'
    coord_max_out = coord_max + '.out'
    coord_min_out = coord_min + '.out'
    bb = 'bb.im'
    jj = 'jj.fits'
    myout1 = 'myout1.im'
    myout2 = 'myout2.im'
    testim = 'test.image'

    def setUp(self):
        self.tb = table( )
        if(os.path.exists(list2[1])):
            for file in list2:
                os.system('rm -rf ' +file)
        for file in list2:
            os.system('cp -RL ' +os.path.join(datapath,file)+' ' + file)


    def tearDown(self):
        for file in list2:
            os.system('rm -rf ' +file)
            os.system('rm -rf mask_test*')
            os.system('rm -rf moment_test*')
        list4 = [
            self.myim, self.mymask, self.mymask1, self.mymask2, self.mymask3,
            self.cas5287, self.myexp, self.mygot, self.expected_integrated,
            self.expected_median, self.temporary_integrated,
            self.temporary_median, self.cas8570, self.zz, self.zz_out,
            self.zz_out_wc, self.zz_out_avg, self.zz_out_wdc, self.median_coord,
            self.median_coord_out, self.myim2, self.coord, self.coord_max,
            self.coord_min, self.coord_max_out, self.coord_min_out, self.bb,
            self.jj, self.myout1, self.myout2, self.testim
        ]
        for f in list4:
            if os.path.exists(f):
                if os.path.isdir(f):
                    shutil.rmtree(f)
                else:
                    os.remove(f)
        self.assertTrue(len(self.tb.showcache()) == 0)
        self.tb.done( )

    
    ####################################################################
    # Testing calculation of each type of moment
    #
    # Notes: Size of input image is 800 x 800 (sky pixels), 1 stokes
    #        value, 18 channels.
    #
    # TODO: moments -1, 2-11
    ####################################################################
    def test_moments(self):
        '''Immoment: Testing calculation of each type of moment'''
        casalog.post( "Starting MOMENT CALCULATION tests", 'NORMAL2' )
        retValue = {'success': True, \
                    'msgs': '', \
                    'error_msgs': "" }
        
        # Calculate the "0" moment and compare with previous results
        # using imval to check a few data points in the image.
        #
        # Note: that the images we are comparing with come from the
        #       ngc1333_regression.py so when this changes, these
        #       test will also need to change.
#        _momentTest_debug_msg( 56 )
        try:
            immoments( 'n1333_both.image', moments=[0], axis='spec', chans='2~15', includepix=[0.003,100.0],excludepix=[-1],outfile = 'moment_test.mom0' )
        except Exception as e:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to very moment 0.\n"\
                       + str( e )
        else:
            if ( not os.path.exists( 'moment_test.mom0' ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                          +"\nError: Moment file, 'moment_test.mom0', was not created."

        # Output some stats on the moment file created.  Just for
        # debugging purposes.  This part could be commented out if
        # its felt to be not needed.
#        _momentTest_debug_msg( 57 )
        stats = imstat( 'moment_test.mom0' )
        casalog.post( str(stats), 'NORMAL3' )
    
        # Check 1% of the data points in the image to see if
        # they match.
        err_margin=0.17  # Needs to take into acct. the noise of the image
                         # no bigger then the value of 1 sigma on the source.
        for i in range( int( (800+800+1+1)*.01 ) ):
            x = random.randint( 0, 799 )
            y = random.randint( 0, 799 )
            sky=str(x)+','+str(y)
            stokes = 'I'
            chan = 0
            current = imval( 'moment_test.mom0', box=sky, stokes=stokes, chans=str(chan) )
            orig    = imval( 'n1333_both.src.tmom0.all', box=sky, stokes=stokes, chans=str(chan) )
            if ( abs(current['data'][0]-orig['data'][0]) > err_margin ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                          +"\nError: Moment test 0 values differ at "+str(x)\
                          +','+str(y)+',I,0.\nThe values are '+str(current)\
                          +" and "+str(orig)
    
        # Find the residual of the original 0th moment and the one
        # we just calculated.  Checking the min/max values to see
        # if they are within our error margin.
        immath( outfile='moment_test.resid0', expr='"moment_test.mom0"-"n1333_both.src.tmom0.all"' )
        resid0_stats=imstat( 'moment_test.resid0' )
        resid0_max   = float(resid0_stats['max'][0])
        resid0_min   = float(resid0_stats['min'][0])
        if ( (resid0_max > 0 and resid0_max > err_margin ) or \
             (resid0_min < 0 and resid0_min < (-1*err_margin) ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                          +"\nError: Moment test 0 residual file varies too"\
                          +" much. Values range from "\
                          +str(float(resid0_stats['min'][0]))+" to "\
                          +str(float(resid0_stats['max'][0]))                                                                                     
                
        # Calculate the 1st moment and compare with previous results
        # using imval to check a few data points in the image.
        #
#        _momentTest_debug_msg( 58 )
        err_margin=0.95  # Needs to take into acct. the noise of the image
                         # no bigger then the value of 1 sigma on the source.
        infile = "immoment_image"
        expected = "first_moment.im"
        got = "moment_test.mom1"
        difference = "first_moment_diff"
        try:
            immoments(
                infile, moments=[1], axis='spec', chans='2~15',
                includepix=[0.02,100.0],excludepix=[-1], outfile=got
            )
        except Exception as e:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to very moment 1.\n"\
                       + str( e )
        else:
            if ( not os.path.exists( got ) ):
                retValue['success']=False
                retValue['error_msgs']=retValue['error_msgs']\
                          +"\nError: Moment file, " + got + ", was not created."

        immath( outfile=difference, expr='"' + got + '"-"' + expected + '"' )
        myia = image()
        myia.open(difference)
        stats = myia.statistics(list=True, verbose=True)
        myia.close()
        casalog.post("moment 1 difference image stats " + str(stats))
        self.assertTrue(stats['sumsq'][0] == 0)

        # Cleanup
        os.system('rm -rf first_moment_diff')
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        self.assertTrue(len(self.tb.showcache()) == 0)

    
    ####################################################################
    # Testing the correctness of the include/exclude pix parameters
    ####################################################################
    def test_pixel_set(self):
        '''Immoment: Testing the correctness of the include/exclude pix parameters'''
        print("starting PIXEL selection tests")
        retValue = {'success': True, \
                    'msgs': '', \
                    'error_msgs': "Pixel selection test NOT implemented yet." }
        
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        self.assertTrue(len(self.tb.showcache()) == 0)

    
    ####################################################################
    # Testing the correctness of the region selection with moment calcs.
    ####################################################################
    def test_region(self):
        '''Immoment: Testing the correctness of the region selection with moment calcs'''
        print("starting REGION selection tests")
        retValue = {'success': True, \
                    'msgs': '', \
                    'error_msgs': "region selection test NOT implemented yet." }
    
    
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        self.assertTrue(len(self.tb.showcache()) == 0)

    
    ####################################################################
    # Testing the correctness of mask usage with moment creation
    #
    # We've already tested valid input and output to the maks
    # parameter so we do not duplicate this testing here. In
    # this function we are more interested in the correctness of
    # the tests, ie. checking the data produced and also testing
    # simple to complex mask specifications.
    #
    # Note: http://www.astron.nl/aips++/docs/notes/223/223.html
    #       contains the legal syntax for masks
    #
    # Test mask is set to and all on n1333_both.image:
    #    n1333_both.image:nomask   - use no mask
    #    n1333_both.image:mask0    - only mask defined in the image.
    #    n1333_both.image>0.1      - all points >0.1
    #    test_image>0.1            - all points >0.1 i test image.
    #
    # TODO create orig images to compare the results too.
    #
    ####################################################################
    def test_mask(self):
        '''Immoment: Testing the correctness of mask usage with moment creation'''
        casalog.post( "Starting MASK tests", 'NORMAL2' )
        retValue = {'success': True, 'msgs': "", 'error_msgs': '' }
    
        try:
            immoments( 'n1333_both.image', mask='mask(n1333_both.image:nomask)', outfile='mask_test_1' )
        except:
            retValue['success']=False
    
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment 0 with mask set as no mssk'"
        if ( not os.path.exists( 'mask_test_1' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Mask test failed, 'mask_test_1`', was not created."
    
        self.assertTrue(len(self.tb.showcache()) == 0)

        try:
            immoments( 'n1333_both.image', mask='mask(n1333_both.image:mask0)',  outfile='mask_test_2' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moments with mask0'"
        if ( not os.path.exists( 'mask_test_2' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Mask test failed, 'mask_test_2', was not created."
                     

        try:
            immoments( 'n1333_both.image', mask='n1333_both.image>0.1', outfile='mask_test_3' )
        except:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moments with mask-'"
        if ( not os.path.exists( 'mask_test_3' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Mask test failed, 'mask_test_3', was not created."
             
        self.assertTrue(len(self.tb.showcache()) == 0)

        try:
            _ia.open( 'n1333_both.image' )
            csys=_ia.coordsys()
            _ia.done()
            _ia.fromshape(
                self.testim, shape=[800,800,1,18], csys=csys.torecord(),
                overwrite=True, log=True
            )
            _ia.addnoise()
            stats=_ia.statistics(list=True, verbose=True)
            _ia.done()
            # pick a place that is slightly bigger then the mid value
            # range for doing the mask, just for fun.
            maskPt=float((stats['max'][0]+stats['min'][0])/2.0)-1.5
            maskStr=self.testim + '>' + str(maskPt)
                                      
            immoments( 'n1333_both.image', mask=maskStr, outfile='mask_test_4' )
        except Exception as e:
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                       +"\nError: Unable to create moment with mask from a 2nd file"\
                       +"\n"+str(e)
        if ( not os.path.exists( 'mask_test_4' ) ):
            retValue['success']=False
            retValue['error_msgs']=retValue['error_msgs']\
                     +"\nError: Moment file, 'mask_test_4', was not created."

        casalog.post( "Done MASK tests", 'NORMAL2' )
        self.assertTrue(retValue['success'],retValue['error_msgs'])
        self.assertTrue(len(self.tb.showcache()) == 0)

    
    def test_CAS2943(self):
        """Test the stretch parameter"""
        myia = image()
        myim = self.myim
        myia.fromshape(myim, [10, 20, 4, 40])
        myia.done()
        myia.fromshape(self.mymask1, [10, 20, 4, 40])
        myia.done()
        myia.fromshape(self.mymask2, [10, 20, 1, 1])
        myia.done()
        myia.fromshape(self.mymask3, [10, 20, 4, 2])
        myia.done()
        self.assertTrue(len(self.tb.showcache()) == 0)
        for i in range(4):
            outfile = ""
            mask = ""
            stretch = False
            exp = []
            atype = True
            if (i == 0):
                outfile = "out1.im"
                mask = "mask1.im > 5"
                stretch = False
                exp = [10, 20, 4, 1]
                atype = True
            if (i == 1):
                outfile = "out2.im"
                mask="mask2.im > 5"
                stretch=False
                atype = False
            if (i == 2):
                outfile = "out3.im"
                mask = "mask2.im > 5"
                stretch = True
                exp = [10, 20, 4, 1]
                atype = True
            if (i == 3):
                outfile = "out4.im"
                mask = "mask3.im > 5"
                stretch = False
                atype = False
            ret = True
            func = ""
            if (atype):
                immoments(
                    imagename=myim, outfile=outfile,
                    mask=mask, stretch=stretch
                )
                self.assertTrue(myia.open(outfile))
                self.assertTrue(myia.shape().tolist() == exp)
                myia.done()
            else:
                self.assertRaises(
                    Exception, immoments, imagename=myim, outfile=outfile,
                    mask=mask, stretch=stretch
                )

            if os.path.exists(outfile):
                shutil.rmtree(outfile)
        self.assertTrue(len(self.tb.showcache()) == 0)

    def test_CAS_4400(self):
        """Verify feature addtion to preconvolve with largest beam when multiple varying beams."""
        myia = image()
        shape = [100, 100, 1, 3]
        myia.fromshape("", shape)
        myia = myia.subimage(self.myexp, dropdeg=True)
        aa = make_gauss2d([shape[0], shape[1]], 10, 10)
        for i in range(shape[3]):
            myia.putchunk(pixels=aa, blc=[0, 0, i])
        myia.setbrightnessunit("Jy/pixel")
        moments = [0, 3]
        ret = myia.moments(moments=moments, axis=2)
        ret.done()
        myia.done()
        for im in ["integrated", "median"]:
            myia.open(self.myexp + '.' + im)
            self.assertTrue(myia.getchunk(getmask=True).all(), "bad mask for " + myia.name())
            myia.done()		
        myia.fromshape("", shape)
        myia = myia.subimage(self.mygot, dropdeg=True)
        myia.setbrightnessunit("Jy/beam")
        for i in range(shape[3]):
            fwhm = 6 + 2*i
            aa = make_gauss2d([shape[0], shape[1]], fwhm, fwhm)
            myia.putchunk(pixels=aa, blc=[0, 0, i])
            f = str(fwhm) + "arcmin"
            myia.setrestoringbeam(major=f, minor=f, pa="0deg", channel=i)
        ret = myia.moments(moments=moments, axis=2)
        ret.done()
        myia.done()
        got = image()
        exp = image()
        mask = image()
        mask.open(self.myexp)
        mask = mask.subimage(
            self.mymask, region=_rg.box(blc=[0, 0, 0], trc=[99, 99,0]),
            dropdeg=True
        )
        cc = mask.getchunk()
        mask.done()
        self.assertTrue(len(self.tb.showcache()) == 0)
        for i in range(shape[0]):
            for j in range(shape[1]):
                if cc[i,j] < 0.01:
                    cc[i,j] = 0
                else:
                    cc[i,j] = 1
        for im in ["integrated", "median"]:
            self.assertTrue(got.open(self.temporary + im))
            self.assertTrue(exp.open(self.myexp + '.' + im))
            shape = got.shape()
            gotpix = got.getchunk() * cc
            exppix = exp.getchunk() * cc
            got.done()
            exp.done()
            epsilon = 1e-5
            
            mymax = abs(gotpix - exppix).max()
            self.assertTrue(mymax < epsilon)
        self.assertTrue(len(self.tb.showcache()) == 0)


    def test_CAS4526(self):
        """Verify CAS-4526 fix"""
        myia = image()
        myia.fromshape("", [100, 100, 1, 10])
        myia.addnoise()
        myia.setbrightnessunit("Jy/beam")
        for i in range(10):
            bmaj = str(4+i) + "arcmin"
            bmin = str(3+i) + "arcmin"
            bpa = "60deg"
            myia.setrestoringbeam(major=bmaj, minor=bmin, pa=bpa, channel=i, polarization=0)
        self.assertTrue(myia.moments(moments=[0], axis=3))
        self.assertTrue(myia.moments(moments=[0]))
        self.assertTrue(len(self.tb.showcache()) == 0)

    def test_CAS5287(self):
        """ Verify fix of CAS-5287"""
        outfile = self.cas5287
        immoments(
            imagename=os.path.join(datapath,"38chans.im"),moments=[8],axis="spectral",
            outfile=outfile
        )
        myia = image()
        myia.open(outfile)
        self.assertTrue((myia.shape() == [1,1,1,1]).all())
        myia.done()
        
    def test_region_selected_in_fits(self):
        """Test that region selection happens in non-paged images, CAS-5278"""
        myia = image()
        test1 = self.bb
        myia.fromshape(test1 ,[2, 2, 5])
        bb = myia.getchunk()
        for i in range(5):
            bb[:,:,i] = 5
        myia.putchunk(bb)
        fits = self.jj
        myia.tofits(fits)
        out1 = self.myout1
        immoments(imagename=test1, outfile=out1, chans="0~2")
        myia.open(out1)
        aa = myia.getchunk()
        out2 = self.myout2
        immoments(imagename=fits, outfile=out2, chans="0~2")
        myia.open(out2)
        bb = myia.getchunk()
        myia.done()
        self.assertTrue((aa == bb).all())
        
    def test_minmax_coord(self):
        """Verify CAS-5376, test min/max coords"""
        myia = image()
        for minmax, im in zip(("max", "min"),(self.coord_max, self.coord_min)):
            imagename = im
            myia.fromshape(imagename, [10, 10, 10])
            vels = []
            bb = myia.getchunk()
            val = 1
            if (minmax == "min"):
                val = -1
            for i in range(10):
                world = myia.toworld([0, 0, i])['numeric']
                vels.append(myia.coordsys().frequencytovelocity(world[2])[0])
                for j in range(10):
                    bb[i, j, i] = val
            myia.putchunk(bb)
            moments = 9
            if minmax == "min":
                moments = 11
            myia.done()
            outfile = imagename + ".out"
            immoments(imagename, outfile=outfile, moments=moments)
            myia.open(outfile)
            self.assertTrue(myia.brightnessunit() == "km/s")
            bb = myia.getchunk()
            for i in range(10):
                for j in range(10):
                    self.assertTrue(abs(bb[i,j] - vels[i]) < 1e-4)
            myia.done()
            
    def test_median_coord(self):
        """Verify CAS-5570 fix for median"""
        myia = image()
        imagename = self.median_coord
        myia.fromshape(imagename, [11, 1, 11])
        cc = myia.getchunk()
        vels = []
        for i in range(11):
            for j in range(11):
                cc[i, 0, j] = (j - i + 5) % 11 + 1
            world = myia.toworld([0, 0, i])['numeric']
            vels.append(myia.coordsys().frequencytovelocity(world[2])[0])
            
        myia.putchunk(cc)
        myia.done()
        # FIXME this should be a task, not tool method, call for casa6
        outfile = self.median_coord_out
        immoments(imagename, outfile=outfile, moments=4, includepix=[0,12])
        myia.open(outfile)
        bb = myia.getchunk()
        myia.done()
        for i in range(11):
            self.assertTrue(abs(bb[i,0] - vels[i]) < 10)
        
    def test_CAS7850(self):
        """Verify support for rest frequency=0"""
        myia = image()
        myia.fromshape("", [20, 20, 20])
        csys = myia.coordsys()
        csys.setrestfrequency("0Hz")
        myia.setcoordsys(csys.torecord())
        csys.done()
        zz = myia.moments(moments=0)
        self.assertTrue(zz)
        zz = myia.moments(moments=1)
        self.assertTrue(zz)
        zz.done()
        myia.done()
                
    def test_history(self):
        """Verify that history is written"""
        myia = image()
        imagename = self.zz
        myia.fromshape(imagename, [20, 20, 20])
        bb = myia.moments()
        myia.done()
        msgs = bb.history()
        bb.done()
        teststr = "ia.moments"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
        outfile = self.zz_out
        moments = [-1, 1, 2]
        axis = 2
        immoments(imagename=imagename, outfile=outfile, axis=axis, moments=moments)
        for s in [".weighted_coord", ".average", ".weighted_dispersion_coord"]:
            im = outfile + s
            myia.open(im)
            msgs = myia.history()
            myia.done()
            teststr = "version"
            self.assertTrue(teststr in msgs[-2], f'"{teststr}" not found')
            teststr = "immoments"
            self.assertTrue(teststr in msgs[-1], f'"{teststr}" not found')

    def test_flush(self):
        """CAS-8570: Ensure moments images are flushed to disk"""
        myia = image()
        myia.fromshape("", [20, 20, 20])
        outfile = self.cas8570
        bb = myia.moments(outfile=outfile, moments=8)
        # this running correctly indicates issue resolved
        imhistory(outfile)
        myia.done()
        bb.done()
    
if __name__ == '__main__':
    unittest.main()
