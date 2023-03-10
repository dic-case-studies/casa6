#########################################################################
# test_task_spxfit.py
#
# Copyright (C) 2018
# Associated Universities, Inc. Washington DC, USA
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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.spxfit.html
#
#
##########################################################################
import shutil
import unittest
import numpy
import math
import os
import glob
import shutil

from casatools import ctsys, image, functional, table
from casatasks import spxfit

nanvalue = 4.53345345

datapath=ctsys.resolve('unittest/spxfit/')

myia = image()
_fn = functional()

class spxfit_test(unittest.TestCase):

    def setUp(self):
        self._tb = table( )

    def tearDown(self):
        myia.done()
        files = glob.glob('cubeApF.*')
        files.extend(glob.glob('test_sol.im_*'))
        files.extend(glob.glob('concat*.im'))
        files.extend(['spxfit.im', 'spxfit.log', 'model.im'])
        for f in files:
            if os.path.islink(f) or os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)

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
                        if numpy.isnan(arr2[i]):
                            arr2[i] = nanvalue
                    arr = arr2
            diffData = newgot - newexp
            self.assertTrue(abs(diffData).max() < 2e-7)
        else:
            self.assertTrue((gotArray == expectedArray).all())

    def checkImage(self, gotImage, expectedName):
        expected = image()
        self.assertTrue(expected.open(expectedName))
        got = image()
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
        myfn = _fn.powerlogpoly(plpest)
        for i in range(zz.shape[2]):
            world = myia.toworld([0,0,i])['numeric'][2]
            zz[:,:,i] = myfn.f(world/1e9)
        myia.putchunk(zz)

        rec = spxfit(imagename=imagename, spxtype="plp", spxest=plpest)
        sols = rec['plp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/plpest) < 0.1e-7).all())
        rec = spxfit(imagename=imagename, spxtype="plp", spxest=[0.4, 3])
        sols = rec['plp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/plpest) < 0.1e-7).all())

        myia.addnoise(pars=[0, 0.001])
        plpestoff = [0.4, 3]
        rec = spxfit(imagename=imagename, spxtype="plp", spxest=plpestoff)
        sols = rec['plp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/plpest) < 4e-2).all())
        plpsol = "plpsol.im"
        plperr = "plperr.im"
        plpestoff = [0.4, 2.2]
        plpfix = [False, True]
        rec = spxfit(
            imagename=imagename, spxtype="plp", spxest=plpestoff, spxfix=plpfix,
            multifit=True, spxsol=plpsol, spxerr=plperr
        )
        myia.done(remove=True)
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
        myfn = _fn.powerlogpoly(plpest)
        for i in range(zz.shape[2]):
            world = myia.toworld([0,0,i])['numeric'][2]
            zz[:,:,i] = myfn.f(world/1e9)
        myia.putchunk(zz)
        ltpest = plpest
        ltpest[:] = plpest
        ltpest[0] = math.log(plpest[0])
        rec = spxfit(imagename=imagename, spxtype="ltp", spxest=ltpest)
        sols = rec['ltp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/ltpest) < 0.1e-7).all())
        rec = spxfit(imagename=imagename, spxtype="ltp", spxest=[0.4, 3])
        sols = rec['ltp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/ltpest) < 0.1e-7).all())
        self.assertTrue(rec['xUnit'] == "Hz")

        myia.addnoise(pars=[0, 0.001])
        ltpestoff = [0.4, 3]
        rec = spxfit(imagename=imagename, spxtype="ltp", spxest=ltpestoff)
        sols = rec['ltp']['solution'].ravel()
        self.assertTrue(
            (abs(1 - sols/ltpest) < 3e-2).all(),
            "sols " + str(sols) + " ltpest " + str(ltpest)
        )
        spxsol = "ltpsol.im"
        spxerr = "ltperr.im"
        ltpestoff = [0.4, 2.2]
        ltpfix = [False, True]
        rec = spxfit(
            imagename=imagename, spxtype="ltp", spxest=ltpestoff, spxfix=ltpfix,
            multifit=True, spxsol=spxsol, spxerr=spxerr
        )
        myia.done(remove=True)
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

    def test_multi_image(self):
        """Test multi image support"""
        imagename1 = "concat1.im"
        myia.fromshape(imagename1,[2, 2, 100])
        csys = myia.coordsys()
        inc = csys.increment()['numeric']
        inc[2] = 1e7
        csys.setincrement(inc)
        myia.setcoordsys(csys.torecord())

        imagename2 = "concat2.im"
        myia.fromshape(imagename2,[2, 2, 100])
        csys = myia.coordsys()
        inc = csys.increment()['numeric']
        inc[2] = 1e7
        csys.setincrement(inc)
        inc[2] = 1e7
        refval = csys.referencevalue()['numeric']
        refval[2] = 3e9
        csys.setreferencevalue(refval)
        myia.setcoordsys(csys.torecord())
        plpest = [0.5, 2]
        for imagename in [imagename1, imagename2]:
            myia.open(imagename)
            zz = myia.getchunk()
            myfn = _fn.powerlogpoly(plpest)
            for i in range(zz.shape[2]):
                world = myia.toworld([0,0,i])['numeric'][2]
                zz[:,:,i] = myfn.f(world/1e9)
                myia.putchunk(zz)
        rec = spxfit(imagename=[imagename1, imagename2], spxtype="plp", spxest=plpest, model="model.im")
        sols = rec['plp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/plpest) < 1e-9).all())
        rec = spxfit(imagename=[imagename1, imagename2], spxtype="plp", spxest=[0.4, 3])
        sols = rec['plp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/plpest) < 1e-9).all())

        myia.open(imagename1)
        myia.addnoise(pars=[0, 0.001])
        myia.open(imagename2)
        myia.addnoise(pars=[0, 0.001])

        plpestoff = [0.4, 3]
        rec = spxfit(imagename=[imagename1, imagename2], spxtype="plp", spxest=plpestoff)
        sols = rec['plp']['solution'].ravel()
        self.assertTrue((abs(1 - sols/plpest) < 1e-2).all())

        plpsol = "plpsol.im"
        plperr = "plperr.im"
        plpestoff = [0.4, 2.2]
        plpfix = [False, True]
        rec = spxfit(
            imagename=[imagename1, imagename2], spxtype="plp", spxest=plpestoff,
            spxfix=plpfix, multifit=True, spxsol=plpsol, spxerr=plperr
        )
        myia.done(remove=True)
        sols = rec['plp']['solution']
        self.assertTrue((sols[:,:,:,1] == 2.2).all())
        for j in [0, 1]:
            myia.open(plpsol + "_" + str(j))
            self.assertTrue(
                (
                 abs(myia.getchunk()/sols[:,:,:,j] - 1) < 1e-7
                ).all()
            )
            myia.done(remove=True)

            myia.open(plperr + "_" + str(j))
            self.assertTrue(
                (
                 abs(myia.getchunk() - rec['plp']['error'][:,:,:,j]) < 1e-8
                ).all()
            )

    def test_output_mask(self):
        """ Test the the output solution image mask is correct"""
        im45 = os.path.join(datapath,"small_42GHz_map.image")
        im690 = os.path.join(datapath,"small_690GHz_map.image")
        outfile = "test_sol.im"
        res = spxfit(
            imagename=[im45, im690],multifit=True,
            spxest=[1e-3,-3],spxtype='ltp',div='100GHz',spxsol=outfile
        )
        self.assertTrue(res)
        for i in (0, 1):
            myia.open(outfile + "_" + str(i))
            mask = myia.getchunk(getmask=True)
            self.assertTrue((mask.shape == myia.shape()).all())
            self.assertTrue(mask.all())
            myia.done()

    def test_mask_and_pixels(self):
        """Test that the mask and pixels of the output are correct"""
        image = "cubeApF.im"
        os.symlink(os.path.join(datapath,image), image)
        mask = mask = 'mask("'+image +':spxmask")'
        sol = 'cubeApF.spx'
        err = 'cubeApF.spxerr'
        model = 'cubeApF.spxmodel'
        resid = 'cubeApF.spxresidual'
        logfile = "spxfit.log"
        res = spxfit(
            imagename = image, spxsol=sol,
            spxerr=err, model=model,
            residual=resid, spxest=[0.5,0.0],
            multifit=True,mask=mask,logresults=False,
            logfile=logfile
        )
        for im in (
            sol + "_0", sol + "_1", err + "_0",
            err + "_1", model, resid
        ):
            print("checking image product " + im)
            self.checkImage(im, os.path.join(datapath,im))
            # check history
            myia.open(im)
            msgs = myia.history()
            teststr = "version"
            self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")
            teststr = "spxfit"
            self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
            myia.done()
        for im in (sol, err):
            key = 'solution'
            if (im == err):
                key = 'error'
            for j in [0, 1]:
                myname = im + "_" + str(j)
                myia.open(myname)
                print("checking " + key + " " + str(j) + " array")
                self.checkArray(res['plp'][key][:,:,:,:,j], myia.getchunk())
                myia.done()
        got = open(logfile).readlines()
        expec = open(os.path.join(datapath,logfile)).readlines()
        for i in [9, 6]:
            del got[i]
            del expec[i]
        for i in range(len(got)):
            self.assertTrue(got[i] == expec[i]);

        for pixel in ([340, 265], [508,378]):
            box = str(pixel[0]) + "," + str(pixel[1])
            box = box + "," + box
            res2 = spxfit(
                imagename = image, box=box,
                spxest=[0.5,0.0],
                multifit=True,mask=mask,
                logfile=logfile
            )
            for t in ("solution", "error"):
                got = res2['plp'][t][0,0,0,0]
                expec = res['plp'][t][pixel[0],pixel[1],0,0]
                if (got.max() != 0 and expec.max() != 0):
                    diff = abs((got - expec)/expec)
                    self.assertTrue(diff.max() < 0.003, "got " + str(got) + " exp " + str(expec))
            got = res2['direction'][0,0,0,0]
            expec = res['direction'][pixel[0],pixel[1],0,0]
            self.assertTrue(got == expec)

if __name__ == '__main__':
    unittest.main()
