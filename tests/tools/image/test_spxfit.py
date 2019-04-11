##########################################################################
# imfit_test.py
#
# Copyright (C) 2008, 2009
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
# <author>
# Dave Mehringer
# </author>
#
# <summary>
# Test suite for the CASA task spxfit and tool method ia.fitprofile
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <prerequisite>
# <ul>
#   <li> <linkto class="task_specfit.py:description">spxfit</linkto> 
# </ul>
# </prerequisite>
#
# <etymology>
# Test for the spxfit task and ia.fitprofile() tool method.
# </etymology>
#
# <synopsis>
# Test the spxfit task and the ia.fitprofile() method upon which it is built.
# </synopsis> 
#
# <example>
#
# This test runs as part of the CASA python unit test suite and can be run from
# the command line via eg
# 
# `echo $CASAPATH/bin/casa | sed -e 's$ $/$'` --nologger --log2term -c `echo $CASAPATH | awk '{print $1}'`/code/xmlcasa/scripts/regressions/admin/runUnitTest.py test_spxfit[test1,test2,...]
#
# </example>
#
# <motivation>
# To provide a test standard for the spxfit task to ensure
# coding changes do not break the associated bits 
# </motivation>
#

###########################################################################
import shutil
import unittest
import numpy
import os
from math import log

from casatools import image as iatool
from casatools import functional as fntool
from casatools import table

nanvalue = 4.53345345

datapath='regression/unittest/spxfit/'
myia = iatool()
fn = fntool()

def run_fitprofile (
    imagename, box, region, chans, stokes,
    axis, mask, ngauss, poly, multifit, model="",
    residual="", amp="", amperr="", center="", centererr="",
    fwhm="", fwhmerr="", integral="", integralerr="",
    estimates="", logresults=True, pampest="", pcenterest="", pfwhmest="",
    pfix="", gmncomps=0, gmampcon="", gmcentercon="",
    gmfwhmcon="", gmampest=[0], gmcenterest=[0],
    gmfwhmest=[0], gmfix="", logfile="", pfunc="",
    goodamprange=[0.0], goodcenterrange=[0.0], goodfwhmrange=[0.0],
    sigma="", outsigma=""
):
    myia = iatool()
    myia.open(imagename)
    if (not myia.isopen()):
        myia.done()
        raise Exception
    res = myia.fitprofile(
        box=box, region=region, chans=chans,
        stokes=stokes, axis=axis, mask=mask,
        ngauss=ngauss, poly=poly, estimates=estimates,
        multifit=multifit,
        model=model, residual=residual, amp=amp,
        amperr=amperr, center=center, centererr=centererr,
        fwhm=fwhm, fwhmerr=fwhmerr, integral=integral,
        integralerr=integralerr, logresults=logresults, pampest=pampest,
        pcenterest=pcenterest, pfwhmest=pfwhmest, pfix=pfix,
        gmncomps=gmncomps, gmampcon=gmampcon,
        gmcentercon=gmcentercon, gmfwhmcon=gmfwhmcon,
        gmampest=gmampest, gmcenterest=gmcenterest,
        gmfwhmest=gmfwhmest, gmfix=gmfix, logfile=logfile,
        pfunc=pfunc, goodamprange=goodamprange,
        goodcenterrange=goodcenterrange,
        goodfwhmrange=goodfwhmrange, sigma=sigma, outsigma=outsigma    )
    myia.close()
    myia.done()
    return res

class spxfit_test(unittest.TestCase):
    
    def setUp(self):
        pass

    def tearDown(self):
        myia.done()
        tb = table( )
        self.assertTrue(len(tb.showcache()) == 0)
        tb.done( )

    def checkArray(self, gotArray, expectedArray):
        mytype = type(gotArray.ravel()[0])
        self.assertTrue(gotArray.shape == expectedArray.shape)
        if mytype == numpy.float64:
            newgot = gotArray.copy()
            newexp = expectedArray.copy()
            for arr in [newgot, newexp]:
                if (numpy.isnan(arr).any()):
                    arr2 = arr.ravel()
                    for i in range(len(arr2)):
                        if isnan(arr2[i]):
                            arr2[i] = nanvalue
                    arr = arr2
            diffData = newgot - newexp
            self.assertTrue(abs(diffData).max() < 2e-7)
        else:
            self.assertTrue((gotArray == expectedArray).all())

    def checkImage(self, gotImage, expectedName):
        expected = iatool()                                
        self.assertTrue(expected.open(expectedName))
        got = iatool()
        if type(gotImage) == str:
            self.assertTrue(got.open(gotImage))
        else:
            got = gotImage
        self.assertTrue((got.shape() == expected.shape()).all())
        gotchunk = got.getchunk()
        expchunk = expected.getchunk()
        self.checkArray(gotchunk, expchunk)
        self.checkArray(got.getchunk(getmask=True), expected.getchunk(getmask=True))
        gotCsys = got.coordsys()
        expectedCsys = expected.coordsys()
        diffPixels = gotCsys.referencepixel()['numeric'] - expectedCsys.referencepixel()['numeric']
        self.assertTrue(abs(diffPixels).max() == 0)
        
        diffRef = gotCsys.referencevalue()['numeric'] - expectedCsys.referencevalue()['numeric']
        # fracDiffRef = (diffRef)/expectedCsys.referencevalue()['numeric'];
        self.assertTrue(abs(diffRef).max() == 0)
        got.close()
        got.done()
        expected.close()
        expected.done()

    def test_exceptions(self):
        """spxfit: Test various exception cases"""
        myia.fromshape("", [1,1,10])
        self.assertRaises(Exception, myia.fitprofile, poly=2, plpest=[1,2])
        
    def test_plpfit(self):
        """ Test fitting a power logarithmic polynomial"""
        imagename = "spxfit.im"
        myia.fromshape(imagename,[2, 2, 100])
        csys = myia.coordsys()
        inc = csys.increment()['numeric']
        inc[2] = 1e7
        csys.setincrement(inc)
        myia.setcoordsys(csys.torecord())
        zz = myia.getchunk()
        plpest = [0.5, 2]
        myfn = fn.powerlogpoly(plpest)
        for i in range(zz.shape[2]):
            world = myia.toworld([0,0,i])['numeric'][2]
            zz[:,:,i] = myfn.f(world/1e9)
        myia.putchunk(zz)
        
        rec = myia.fitprofile(ngauss=0, spxtype="plp", spxest=plpest)
        sols = rec['plp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/plpest) < 0.1e-7).all())
        if i == 1:
            rec = myia.fitprofile(ngauss=0, spxtype="plp", spxest=[0.4, 3])
        sols = rec['plp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/plpest) < 0.1e-7).all())
            
        myia.addnoise(pars=[0, 0.001])
        plpestoff = [0.4, 3]
        rec = myia.fitprofile(ngauss=0, spxtype="plp", spxest=plpestoff)
        sols = rec['plp']['solution'].ravel()
        print("*** i %d" % i)
        print("** max %s" % (sols/plpest))
        self.assertTrue((abs(1 - sols/plpest) < 4e-2).all())
        plpsol = "plpsol.im"
        plperr = "plperr.im"
        plpestoff = [0.4, 2.2]
        plpfix = [False, True]
        rec = myia.fitprofile(
            ngauss=0, spxtype="plp", spxest=plpestoff, spxfix=plpfix,
            multifit=True, spxsol=plpsol, spxerr=plperr
        )
        sols = rec['plp']['solution']
        self.assertTrue((sols[:,:,:,1] == 2.2).all())
        for j in [0, 1]:
            myia.open(plpsol + "_" + str(j))
            self.assertTrue(
                (abs(myia.getchunk()/sols[:,:,:,j] - 1) < 1e-7).all()
            )
            myia.done(remove=True)

            myia.open(plperr + "_" + str(j))
            self.assertTrue(
                (abs(myia.getchunk() - rec['plp']['error'][:, :, :, j]) < 1e-8).all()
            )
            myia.done(remove=True)
            
    def test_ltpfit(self):
        """ Test fitting a logarithmic transformed polynomial"""
        imagename = "ltpfit.im"
        myia.fromshape(imagename,[2, 2, 100])
        csys = myia.coordsys()
        inc = csys.increment()['numeric']
        inc[2] = 1e7
        csys.setincrement(inc)
        myia.setcoordsys(csys.torecord())
        zz = myia.getchunk()
        plpest = [0.5, 2]
        myfn = fn.powerlogpoly(plpest)
        for i in range(zz.shape[2]):
            world = myia.toworld([0,0,i])['numeric'][2]
            zz[:,:,i] = myfn.f(world/1e9)
        myia.putchunk(zz)
        ltpest = plpest
        ltpest[:] = plpest
        ltpest[0] = log(plpest[0])
        rec = myia.fitprofile(ngauss=0, spxtype="ltp", spxest=ltpest)
        print(str(rec))
        sols = rec['ltp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/ltpest) < 0.1e-7).all())
        rec = myia.fitprofile(ngauss=0, spxtype="ltp", spxest=[0.4, 3])
        sols = rec['ltp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/ltpest) < 0.1e-7).all())
        print('*** xUnit %s' % rec['xUnit'])
        self.assertTrue(rec['xUnit'] == "Hz")
        
        myia.addnoise(pars=[0, 0.001])
        ltpestoff = [0.4, 3]
        rec = myia.fitprofile(ngauss=0, spxtype="ltp", spxest=ltpestoff)
        sols = rec['ltp']['solution'].ravel()
        self.assertTrue(
            (abs(1 - sols/ltpest) < 3e-2).all(),
            "sols " + str(sols) + " ltpest " + str(ltpest)
        )
        spxsol = "ltpsol.im"
        spxerr = "ltperr.im"
        ltpestoff = [0.4, 2.2]
        ltpfix = [False, True]
        rec = myia.fitprofile(
            ngauss=0,  spxtype="ltp", spxest=ltpestoff, spxfix=ltpfix,
            multifit=True, spxsol=spxsol, spxerr=spxerr
        )
        sols = rec['ltp']['solution']
        self.assertTrue((sols[:,:,:,1] == 2.2).all())
        for j in [0, 1]:
            self.assertTrue(myia.open(spxsol + "_" + str(j)))
            self.assertTrue(
                (abs(myia.getchunk()/sols[:,:,:,j] - 1) < 1e-7).all()
            )
            myia.done(remove=True)

            self.assertTrue(myia.open(spxerr + "_" + str(j)))
            self.assertTrue(
                (abs(myia.getchunk() - rec['ltp']['error'][:,:,:,j]) < 1e-8).all()
            )
            myia.done(remove=True)
            
    def test_ltpfit_with_negative_values(self):
        """Test fitting of ltp when y values are negative returns something reasonable because of auto masking"""
        global myia
        myia.fromshape("", [1, 1, 10])
        f = []
        for i in range(10):
            f.append(myia.toworld([0,0,i])['numeric'][2])
        f = numpy.array(f)
        vals = (f/f[0])**(3.0)
        # add a negative value to the array to test that it's not used in the fit
        vals[0] = -1
        bb = myia.getchunk()
        bb[0, 0, :] = vals
        myia.putchunk(bb)
        res = myia.fitprofile(spxtype="ltp",spxest=[0, 3],div=f[0])
        myia.done()
        self.assertTrue(abs(res['ltp']['solution'][0][0][0][1] - 3) < 0.01 )
        
def suite():
    return [spxfit_test]

if __name__ == '__main__':
    unittest.main()
