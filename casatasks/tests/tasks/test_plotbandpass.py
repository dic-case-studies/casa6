##########################################################################
# test_plotbandpass.py
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
# [https://open-jira.nrao.edu/browse/CAS-9912]
#
#
##########################################################################

CASA6 = False
try:
    import casatools
    from casatasks import plotbandpass
    CASA6 = True
except ImportError:
    from __main__ import *
    from tasks import *
    from taskinit import *

import struct

def pngWidthHeight(filename):
    """
    Reads the width and height of a png image (in pixels).
    -Todd Hunter
    """
    if (os.path.exists(filename) == False):
        print("Cannot find file = ", filename)
        return
    f = open(filename, 'rb')
    data = f.read()
    f.close()
    if CASA6:
        if (data[12:16].decode("utf-8") == 'IHDR'):
            w, h = struct.unpack('>LL', data[16:24])
            width = int(w)
            height = int(h)
    else:
        if (data[:8] == '\211PNG\r\n\032\n'and (data[12:16] == 'IHDR')):
            w, h = struct.unpack('>LL', data[16:24])
            width = int(w)
            height = int(h)
    #else:
    #    raise Exception('not a png image')
    
    return width, height
    
import sys
import os
import unittest
import shutil
import numpy as np
from filecmp import dircmp
datapath = "/export/home/jarvis_2/awells/JIRA/CAS-9912/data/"
datapath ="/lustre/cv/users/awells/CAS-9912/data/test_plotbandpass/"
figdir= os.getcwd() + '/'

class plotbandpass_1_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        if not CASA6:
            default(plotbandpass)
        os.symlink(datapath+'Band7multi_april22.ms', os.getcwd() + '/Band7multi_april22.ms')
            
    def tearDown(self):
        os.unlink(os.getcwd() + '/Band7multi_april22.ms')

    @classmethod
    def tearDownClass(cls):
        pass

    # 0
    def test_createImage_regression00(self):

        #regression00.pdf
        #regression00.spw00.t00.png
        #regression00.spw01.t01.png
        #regression00.spw02.t02.png

        plotbandpass(datapath + 'bandpass.bcal',showtsky=False,xaxis='freq',yaxis='amp',overlay='antenna',spw='',field='0', interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(0),debug=False)

        width,height = pngWidthHeight(os.getcwd()+'/regression00.spw02.t02.png')
        self.assertTrue(width == 864,"Observed: {}, Expected: {}".format(width, 864) )

    #1 
    def test_createImage_regression01(self):

        #regression01.DV04.spw00.t00.png
        #regression01.DV07.spw00.t01.png
        #regression01.DV08.spw00.t02.png
        #regression01.DV10.spw00.t00.png
        #regression01.pdf
        #regression01.PM01.spw00.t01.png
        #regression01.PM02.spw00.t02.png

        plotbandpass(datapath+'bandpass.bcal',showtsky=False,xaxis='freq',yaxis='amp',overlay='baseband',spw='',field='3c279', interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(1),debug=False)

    # 2
    def test_createImage_regression02(self):
        #regression02.pdf
        #regression02.spw00.t00.png
        #regression02.spw01.t01.png
        #regression02.spw02.t02.png

        plotbandpass(datapath+'bandpass.bcal',showtsky=False,xaxis='freq',yaxis='ampdb',overlay='antenna',spw='',
             field='!Titan,!TW Hya,!J1147-382=QSO',
             interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(2),debug=False)

    # 3
    def test_createImage_regression03(self):
        '''Test 3'''
        """
        regression03.DV04.spw00.t00.png
        regression03.DV06.spw00.t02.png
        regression03.DV06.spw01.t01.png
        regression03.DV06.spw02.t00.png
        regression03.DV06.spw02.t02.png
        regression03.DV06.spw03.t01.png
        regression03.DV07.spw00.t00.png
        regression03.DV07.spw00.t02.png
        regression03.DV07.spw01.t01.png
        regression03.DV07.spw02.t00.png
        regression03.DV07.spw02.t02.png
        regression03.DV07.spw03.t01.png
        regression03.DV08.spw00.t00.png
        regression03.DV08.spw00.t02.png
        regression03.DV08.spw01.t01.png
        regression03.DV08.spw02.t00.png
        regression03.DV08.spw02.t02.png
        regression03.DV08.spw03.t01.png
        regression03.DV09.spw00.t00.png
        regression03.DV09.spw00.t02.png
        regression03.DV09.spw01.t01.png
        regression03.DV09.spw02.t00.png
        regression03.DV09.spw02.t02.png
        regression03.DV09.spw03.t01.png
        regression03.DV10.spw00.t00.png
        regression03.DV10.spw00.t02.png
        regression03.DV10.spw01.t01.png
        regression03.DV10.spw02.t00.png
        regression03.DV10.spw02.t02.png
        regression03.DV10.spw03.t01.png
        regression03.pdf
        regression03.PM01.spw00.t00.png
        regression03.PM01.spw00.t02.png
        regression03.PM01.spw01.t01.png
        regression03.PM01.spw02.t00.png
        regression03.PM01.spw02.t02.png
        regression03.PM01.spw03.t01.png
        regression03.PM02.spw00.t00.png
        regression03.PM02.spw00.t02.png
        regression03.PM02.spw01.t01.png
        regression03.PM02.spw02.t00.png
        regression03.PM02.spw02.t02.png
        regression03.PM02.spw03.t01.png
        regression03.PM03.spw00.t00.png
        regression03.PM03.spw00.t02.png
        regression03.PM03.spw01.t01.png
        regression03.PM03.spw02.t00.png
        regression03.PM03.spw03.t02.png
        """
        plotbandpass(datapath+'bandpass.bcal',showtsky=False,xaxis='freq',yaxis='both',phase=[-180,180],plotrange=[0,0,0,2],
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(3),debug=False)

    # 4
    def test_createImage_regression04(self):
        #regression04.pdf
        #regression04.spw00.t00.png
        #regression04.spw01.t01.png
        #regression04.spw02.t02.png

        plotbandpass(datapath+'bandpass.bcal',showtsky=False,xaxis='freq',yaxis='amp',overlay='antenna',spw='',field='',
             chanrange='1200~2000',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(4),
             debug=False, chanrangeSetXrange=True)

    # 5
    def test_createImage_regression05(self):
        plotbandpass(datapath+'bandpass.bcal',showatm=True,xaxis='chan',yaxis='amp',overlay='antenna',spw='',field='',
                  plotrange=[1200,2000,0,0],interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(5),debug=False)

    # 6
    def test_createImage_regression06(self):
        plotbandpass(datapath+'bandpass.bcal',showatm=True,xaxis='chan',yaxis='amp',spw='',field='',plotrange=[1200,3840,0,0],
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(6),debug=False)

    # 7
    def test_createImage_regression07(self):
        plotbandpass(datapath+'bandpass.bcal',showtsky=False,xaxis='chan',yaxis='amp',overlay='antenna',spw='',field='',showatm=True,
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(7),debug=False)

    # 8
    def test_createImage_regression08(self):
        plotbandpass(datapath+'bandpass.bcal',showtsky=False,xaxis='freq',yaxis='phase',overlay='antenna',spw='',field='',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(8),debug=False)

    # 9
    def test_createImage_regression09(self):
        plotbandpass(datapath+'bandpass.bcal',showtsky=False,xaxis='freq',yaxis='phase',overlay='antenna',spw='',field='',
             chanrange='1200~1800',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(9),
             debug=False, chanrangeSetXrange=True)

    # 10
    def test_createImage_regression10(self):
        plotbandpass(datapath+'bandpass.bcal',overlay='antenna',yaxis='amp',field='0~1,4',xaxis='freq',showtsky=True,
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(10),debug=False)

class plotbandpass_2_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        if not CASA6:
            default(plotbandpass)
        os.symlink(datapath+'X3c1.ms', os.getcwd() + '/X3c1.ms')
            
    def tearDown(self):
        os.unlink(os.getcwd() + '/X3c1.ms')


    @classmethod
    def tearDownClass(cls):
        pass

    #11
    def test_b_skipspw19high_regression11(self):
        plotbandpass(datapath + 'bandpass_b_skipspw19high.bcal',yaxis='amp',field='0',xaxis='freq',
                  caltable2=datapath + 'bandpass_bpoly_skipspw19high.bcal',showpoints=True,spw='0,1',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(11),debug=False)

    # 12
    def test_b_skipspw19high_regression12(self):
        plotbandpass(datapath +'bandpass_b_skipspw19high.bcal',yaxis='phase',field='0',xaxis='freq',
                  caltable2= datapath + 'bandpass_bpoly_skipspw19high.bcal',showpoints=True,spw=0,
                  caltable3= datapath + 'bandpass_bpoly_skipspw19high.bcal',interactive=False,
                  buildpdf=True,figfile=figdir+'regression%02d'%(12),debug=False)


    # 13
    def test_b_skipspw19high_regression13(self):
        plotbandpass(datapath + 'bandpass_b_skipspw19high.bcal',yaxis='both',field='0',xaxis='freq',
                  caltable2= datapath + 'bandpass_bpoly_skipspw19high.bcal',showpoints=True,spw=0,
                  caltable3= datapath + 'bandpass_bpoly_skipspw19high.bcal',interactive=False,
                  buildpdf=True,figfile=figdir+'regression%02d'%(13),debug=False)

class plotbandpass_X3c1_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        if not CASA6:
            default(plotbandpass)
        os.symlink(datapath+'X3c1.ms', os.getcwd() + '/X3c1.ms')
            
    def tearDown(self):
        os.unlink(os.getcwd() + '/X3c1.ms')


    @classmethod
    def tearDownClass(cls):
        pass

    #14
    def test_X3c1_tsys_fdm_regression14(self):
        plotbandpass(datapath + 'X3c1.tsys.fdm',overlay='antenna',yaxis='amp',field='1',xaxis='chan',
                  showtsky=True,interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(14)) 

    # 15
    def test_X3c1_tsys_fdm_regression15(self):
        plotbandpass(datapath + 'X3c1.tsys.fdm',overlay='antenna',yaxis='amp',field='1',xaxis='chan',
                  poln='y',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(15),debug=False)


    # 16
    def test_X3c1_tsys_regression16(self):
        plotbandpass(datapath + 'X3c1.tsys',overlay='antenna',yaxis='amp',field='0~1,4',xaxis='chan',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(16),debug=False)

    # 17
    def test_X3c1_tsys_regression17(self):
        plotbandpass(datapath + 'X3c1.tsys',overlay='time',yaxis='amp',field='2',xaxis='chan',interactive=False,
                  buildpdf=True,figfile=figdir+'regression%02d'%(17),debug=False)

    # 18
    def test_X3c1_tsys_regression18(self):
        plotbandpass(datapath + 'X3c1.tsys',overlay='',yaxis='amp',field='',xaxis='freq',showfdm=True,interactive=False,
                  buildpdf=True,figfile=figdir+'regression%02d'%(18),debug=False)

    # 19
    def test_X3c1_tsys_regression19(self):
        plotbandpass(datapath + 'X3c1.tsys',overlay='time',yaxis='amp',field='',xaxis='freq',showfdm=True,interactive=False,
                  buildpdf=True,figfile=figdir+'regression%02d'%(19),debug=False)

    # 20
    def test_X3c1_tsys_regression20(self):
        plotbandpass(datapath + 'X3c1.tsys',overlay='',yaxis='amp',field='2',xaxis='freq',chanrange='45~65',
             interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(20),debug=False,
             chanrangeSetXrange=True)

class plotbandpass_smooth_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        if not CASA6:
            default(plotbandpass)
        os.symlink(datapath+'Band7multi_april22.ms', os.getcwd() + '/Band7multi_april22.ms')
            
    def tearDown(self):
        os.unlink(os.getcwd() + '/Band7multi_april22.ms')

    @classmethod
    def tearDownClass(cls):
        pass

    # 21
    def test_b_skipspw19high_regression21(self):
        plotbandpass(datapath + 'bandpass_bpoly_skipspw19high.bcal',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(21),debug=False,xaxis='freq')


    # 22
    def test_b_smooth_regression22(self):
        plotbandpass(datapath +'bandpass.bcal',caltable2=datapath +'bandpass.bcal_smooth',xaxis='freq',interactive=False,
                  buildpdf=True,figfile=figdir+'regression%02d'%(22),debug=False) 
  
    # 23
    def test_b_smooth_regression23(self):
        plotbandpass(datapath +'bandpass.bcal',caltable2=datapath +'bandpass.bcal_smooth',xaxis='freq',yaxis='ampdb',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(23),debug=False) 

    # 24
    def test_b_smooth_regression24(self):
        plotbandpass(datapath +'bandpass.bcal',caltable2=datapath +'bandpass.bcal_smooth',yaxis='phase', xaxis='freq',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(24),debug=False)
  
    # 25
    def test_b_smooth_regression25(self):
        plotbandpass(datapath +'bandpass.bcal',caltable2=datapath +'bandpass.bcal_smooth',xaxis='freq',chanrange='1000~3000',
             interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(25),debug=False,
             chanrangeSetXrange=True)
  
    # 26
    def test_b_smooth_regression26(self):
        plotbandpass(datapath +'bandpass.bcal',caltable2=datapath +'bandpass.bcal_smooth',xaxis='freq',showtsky=True,
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(26),debug=False)
  
    # 27
    def test_b_smooth_regression27(self):
        plotbandpass(datapath +'bandpass.bcal',caltable2=datapath +'bandpass.bcal_smooth',xaxis='freq',poln='x',
                  showflagged=True,interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(27),debug=False) 
  
    # 28
    def test_b_smooth_regression28(self):
        plotbandpass(datapath +'bandpass.bcal',caltable2=datapath +'bandpass.bcal_smooth',xaxis='freq',poln='x',
                  showflagged=True, showtsky=True,interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(28),debug=False)

class plotbandpass_X3c1_2_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        if not CASA6:
            default(plotbandpass)
        os.symlink(datapath+'X3c1.ms', os.getcwd() + '/X3c1.ms')
            
    def tearDown(self):
        os.unlink(os.getcwd() + '/X3c1.ms')

    #29
    def test_X3c1_tsys_regression29(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',caltable2=datapath +'X3c1.tsys',    yaxis='amp',xaxis='freq',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(29),debug=False) 

    # 30
    def test_X3c1_tsys_regression30(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',caltable2=datapath +'X3c1.tsys',    yaxis='amp',xaxis='freq',
                  showtsky=True, interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(30),debug=False) 
  
    # 31
    def test_X3c1_tsys_regression31(self):
        plotbandpass(caltable=datapath +'X3c1.tsys',    caltable2=datapath +'X3c1.tsys.fdm',yaxis='amp',xaxis='freq',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(31),debug=False)

    # 32
    def test_X3c1_tsys_regression32(self):
        plotbandpass(caltable=datapath +'X3c1.tsys',    caltable2=datapath +'X3c1.tsys.fdm',yaxis='amp',xaxis='freq',
                  showtsky=True,interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(32),debug=False) 

    # 33  
    def test_X3c1_tsys_regression33(self):
        plotbandpass(caltable=datapath +'X3c1.tsys',caltable2=datapath +'X3c1.tsys.fdm',    yaxis='both',xaxis='freq',
             chanrange='10~118',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(33),
             debug=False, chanrangeSetXrange=True)

    # 34
    def test_X3c1_tsys_regression34(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',caltable2=datapath +'X3c1.tsys',    yaxis='amp',xaxis='freq',
             chanrange='1000~3000',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(34),
             debug=False, chanrangeSetXrange=True)
  
    # 35
    def test_X3c1_tsys_regression35(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',caltable2=datapath +'X3c1.tsys',    yaxis='amp',xaxis='freq',
             chanrange='1000~3000',showtsky=True,interactive=False,buildpdf=True,
             figfile=figdir+'regression%02d'%(35),debug=False, chanrangeSetXrange=True) 

    # 36
    def test_X3c1_tsys_regression36(self):
        plotbandpass(caltable=datapath +'X3c1.tsys',    caltable2=datapath +'X3c1.tsys.fdm',yaxis='both',xaxis='freq',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(36),debug=False) 

  # 37
    def test_X3c1_tsys_regression37(self):
        plotbandpass(caltable=datapath +'X3c1.tsys',    caltable2=datapath +'X3c1.tsys.fdm',yaxis='both',xaxis='freq',
                  showtsky=True,interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(37),debug=False)

    # 38
    def test_X3c1_tsys_regression38(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',caltable2=datapath +'X3c1.tsys',    yaxis='both',xaxis='freq',
                  zoom='intersect',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(38),debug=False)

    # 39
    def test_X3c1_tsys_regression39(self):
        plotbandpass(caltable=datapath +'X3c1.tsys',    caltable2=datapath +'X3c1.tsys.fdm',yaxis='both',xaxis='freq',
                  zoom='intersect',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(39),debug=False)
  
    # 40
    def test_X3c1_tsys_regression40(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',caltable2=datapath +'X3c1.tsys',yaxis='amp',xaxis='freq',poln='XX',
                  interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(40),debug=False) 

    # 41
    def test_X3c1_tsys_regression41(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',caltable2=datapath +'X3c1.tsys',yaxis='amp',field='1',xaxis='freq',
                  poln='YY',zoom='intersect',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(41),debug=False)

    # 42
    def test_X3c1_tsys_regression42(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',caltable2=datapath +'X3c1.tsys',yaxis='amp',field='1',xaxis='freq',
                  poln='YY',zoom='intersect',showatm=True,interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(42),debug=False)

    # 43 
    def test_X3c1_tsys_regression43(self):
        plotbandpass(caltable=datapath +'X3c1.tsys.fdm',yaxis='amp',field='1',xaxis='chan',poln='YY',zoom='intersect',
                  showatm=True,showimage=True,interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(43),debug=False)

class plotbandpass_multi_field_bandpass_test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        if not CASA6:
            default(plotbandpass)
        os.symlink(datapath+'Band7multi_april22.ms', os.getcwd() + '/Band7multi_april22.ms')
            
    def tearDown(self):
        os.unlink(os.getcwd() + '/Band7multi_april22.ms')


    #tests for multi-field bandpass solution overlay: in SciVer/TWHya
    # 44
    def test_multi_field_bandpass_regression44(self):
        plotbandpass(datapath + 'band7multi_a6p7_titan.bcal',caltable2=datapath +'band7multi_b.bcal',xaxis='freq',
             yaxis='amp',chanrange='',interactive=False,buildpdf=True,figfile=figdir+'regression%02d'%(44),
             debug=False)
  
    # 45
    def test_multi_field_bandpass_regression45(self):
        plotbandpass(datapath +'band7multi_a6p7_titan.bcal',caltable2=datapath +'band7multi_b.bcal',xaxis='freq',
             yaxis='amp',chanrange='',showflagged=True,interactive=False,buildpdf=True,
             figfile=figdir+'regression%02d'%(45),debug=False)
  
    # 46
    def test_multi_field_bandpass_regression46(self):
        plotbandpass(caltable=datapath +'band7multi_b.bcal',caltable3=datapath +'band7multi_bpoly_a6p7_titan.bcal',
                  caltable2=datapath +'band7multi_bpoly.bcal',xaxis='freq',yaxis='both',interactive=False,
                  buildpdf=True,figfile=figdir+'regression%02d'%(46),debug=False)
def suite():
    return[plotbandpass_1_test,plotbandpass_2_test, plotbandpass_X3c1_test, plotbandpass_smooth_test]

if __name__ == '__main__':
    unittest.main()

