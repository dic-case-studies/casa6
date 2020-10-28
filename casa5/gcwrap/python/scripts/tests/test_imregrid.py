from __future__ import print_function
import sys
import traceback
import os
import shutil
import random
import re
import time
import numpy as np
import glob
import struct
import unittest

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import ctsys, image, regionmanager, table, measures, componentlist
    from casatasks import immath, casalog
    _ia = image()
    _rg = regionmanager()
    _cs = coordsys()
    _tb = table()
    _me = measures()
    _cl = componentlist()
    ctsys_resolve = ctsys.resolve
else:
    import casac
    from tasks import *
    from taskinit import *
    _ia = iatool()
    _rg = rgtool()
    _cs = cstool()
    _tb = tbtool()
    _me = metool()
    _cl = cltool()
    dataRoot = os.path.join(os.environ.get('CASAPATH').split()[0],'data')
    from casa_stack_manip import stack_frame_find
    casa_stack_rethrow = stack_frame_find().get('__rethrow_casa_exceptions', False)

    def ctsys_resolve(apath):
        return os.path.join(dataRoot,apath)

sep = os.sep
datapath = ctsys_resolve(os.path.join('regression','unittest','imregrid'))

IMAGE = 'image.im'
gim = "gaussian_source.im"

total = 0
fail  = 0
current_test =""
stars = "*************"

def alleqnum(x,num,tolerance=0):
    if len(x.shape)==1:
        for i in range(x.shape[0]):
            if not (abs(x[i]-num) < tolerance):
                print("x[",i,"]=", x[i])
                return False
    if len(x.shape)==2:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                if not (abs(x[i][j]-num) < tolerance):
                    print("x[",i,"][",j,"]=", x[i][j])
                    return False
    if len(x.shape)==3:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                for k in range(x.shape[2]):
                    if not (abs(x[i][j][k]-num) < tolerance):
                        print("x[",i,"][",j,"][",k,"]=", x[i][j][k])
                        return False
    if len(x.shape)==4:
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                for k in range(x.shape[2]):
                    for l in range(x.shape[3]):
                        if not (abs(x[i][j][k][l]-num) < tolerance):
                            print("x[",i,"][",j,"][",k,"][",l,"]=", x[i][j][k])
                            return False
    if len(x.shape)>4:
        stop('unhandled array shape in alleq')
    return True

def test_start(msg):
    global total, current_test
    total += 1
    print()
    print(stars + " Test " + msg + " start " + stars)
    current_test = msg
    
def test_end(condition, error_msg):
    global total, fail
    status = "OK"
    if not condition:
        print(error_msg, file=sys.stderr)
        fail += 1
        status = "FAIL"
    print(stars + " Test " + current_test + " " + status + " " + stars)
        
out1 = 'regridded'
out2 = 'bigger_image'
out3 = 'shifted_image'
out4 = 'back_to_image'
out5 = 'template'
out6 = 'gal_coords.im'

class imregrid_test(unittest.TestCase):

    def setUp(self):
        pass
    
    def tearDown(self):
        _ia.done()
        
        for i in (IMAGE, out1, out2, out3, out4, out5, out6):
            if (os.path.exists(i)):
                os.system('rm -rf ' + i)
        
        self.assertTrue(len(_tb.showcache()) == 0)
        
    def test1(self):    
        myia = _ia  
        myia.maketestimage(outfile = IMAGE)
        
        outim=imregrid(
            imagename = IMAGE,
            template = IMAGE,
            output = out1
        )
        im1 = myia.newimage(IMAGE)
        im2 = myia.newimage(out1)
        
        im1.statistics()
        im2.statistics()
        
        rec1 = im1.torecord()
        print('*************')
        print(rec1['shape'])
        print('*************')
        shape = im1.shape()
        print(shape)
        checked = 0
        for x in range(shape[0]):
            for y in range(shape[1]):
                p1 = im1.pixelvalue([x, y])
                p2 = im2.pixelvalue([x, y])
                if p1['mask'] != p2['mask']:
                    raise Exception(p1['mask'] + ' != ' + p2['mask'])
                if p1['value']['value'] != p2['value']['value']:
                    raise Exception(p1['value']['value'] + ' != ' + p2['value']['value'])
                if p1['value']['unit'] != p2['value']['unit']:
                    raise Exception(p1['value']['unit'] + ' != ' + p2['value']['unit'])
                checked += 3
        
        im2.done()
        
        print(str(checked) + ' values checked')
        
        # rescale by factors 3 x 2
        rec1 = im1.torecord()
        print('*************')
        print("shape before " + str(rec1['shape']))
        print('*************')
        rec1['shape'] = numpy.array([3*rec1['shape'][0], 2*rec1['shape'][1]], numpy.int32)
        print("shape after " + str(rec1['shape']))

        rec1['coordsys']['coordsys']['direction0']['cdelt'] = [
            rec1['coordsys']['coordsys']['direction0']['cdelt'][0]/3.0,
            rec1['coordsys']['coordsys']['direction0']['cdelt'][1]/2.0]
        rec1['coordsys']['coordsys']['direction0']['crpix'] = [
            rec1['coordsys']['coordsys']['direction0']['crpix'][0]*3.0,
            rec1['coordsys']['coordsys']['direction0']['crpix'][1]*2.0]
        print(rec1)
        
        myia.fromrecord(rec1, out2)
        
        # First we need to remove the output file.
        if (  os.path.exists(out1) ):
              shutil.rmtree( out1 )
        outim=imregrid(
            imagename=IMAGE, template=out2,
            output=out1, shape=rec1["shape"]
        )
        s1 = imstat(IMAGE)
        s2 = imstat(out1)
        _ia.open(out1)
        print("out shape " + str(_ia.shape()))
        print("S1: ", s1)
        print(" ")
        print(" ")
        print("S2: ", s2)
        
        if s1['maxpos'][0]*3 != s2['maxpos'][0]:
            raise Exception(str(s1['maxpos'][0]*3) + ' != ' + str(s2['maxpos'][0]))
        if s1['maxpos'][1]*2 != s2['maxpos'][1]:
            raise Exception(str(s1['maxpos'][1]*2) + ' != ' + str(s2['maxpos'][1]))
        
        
        
        # shift by -13, 1 pixels

        rec1 = im1.torecord()
        rec1['coordsys']['coordsys']['direction0']['crpix'] = [
            rec1['coordsys']['coordsys']['direction0']['crpix'][0]-13,
            rec1['coordsys']['coordsys']['direction0']['crpix'][1]+1]
        
        myia.fromrecord(rec1, out3)
        myia.close()
        # First we need to remove the output file.
        if (  os.path.exists(out1 ) ):
              shutil.rmtree( out1)
        outim = imregrid(imagename = IMAGE,
                 template = out3,
                 output = out1)
        s1 = imstat(IMAGE)
        s2 = imstat(out1)
        if s1['maxpos'][0]-13 != s2['maxpos'][0]:
            raise Exception(str(s1['maxpos'][0]-13) + ' != ' + str(s2['maxpos'][0]))
        if s1['maxpos'][1]+1 != s2['maxpos'][1]:
            raise Exception(str(s1['maxpos'][1]+1) + ' != ' + str(s2['maxpos'][1]))
        
        
        # Shift back to original
        rec1['coordsys']['coordsys']['direction0']['crpix'] = [
            rec1['coordsys']['coordsys']['direction0']['crpix'][0]+13,
            rec1['coordsys']['coordsys']['direction0']['crpix'][1]-1]
        if (  os.path.exists(out3 ) ):
            shutil.rmtree( out3)
        myia.fromrecord(rec1, out3)
        myia.close()
        outim = imregrid(imagename = IMAGE,
                 template = out3,
                 output = out4)
        s1 = imstat(IMAGE)
        s2 = imstat(out4)
        print(s1)
        print(s2)
        for stat in ['rms', 'medabsdevmed', 'minpos',
                     'min', 'max', 'sum', 'minposf',
                     'median', 'flux', 'sumsq', 'maxposf',
                     'trcf', 'quartile', 'npts', 'maxpos',
                     'mean', 'sigma', 'trc', 'blc', 'blcf']:
            if type(s1[stat]) == type('a string'):
                print("Checking string", stat, s1[stat])
                if s1[stat] != s2[stat]:
                    raise Exception
            else:
                for i in range(len(s1[stat])):
                    print("Checking", stat, "[", i, "]", s1[stat][i])
                    if s1[stat][i] != s2[stat][i]:
                        # Note:  == comparison of floating point values,
                        # it works right now on this computer but might need to get fixed...
                        raise Exception
        
        
        # Exercise various reference codes (no check on output)
        codes = _cs.newcoordsys(direction=True).referencecode('dir', True)
        rec1 = im1.torecord()
        im1.done()
        for ref in codes:
            print("Regrid to", ref)
            if ref not in ['JMEAN', 'JTRUE', 'APP',
                           'BMEAN', 'BTRUE', 'HADEC',
                           'AZEL', 'AZELSW', 'AZELNE',
                           'AZELGEO',
                           'AZELSWGEO', 'AZELNEGEO',
                           'JNAT',
                           'MECLIPTIC', 'TECLIPTIC',
                           'ITRF', 'TOPO']:
                rec1['coordsys']['coordsys']['direction0']['conversionSystem'] = ref
                if ref == "COMET":
                    ### To regrid to COMET frame you have to specify a COMET table or ephemeris
                    ### to measures. As this is not possible in this interface,
                    ### avoiding this regridding (CAS-11403)
                    self.assertRaises(Exception, myia.fromrecord, rec1, out5)
                    continue
                myia.fromrecord(rec1, out5)
                myia.close()
                if (  os.path.exists(out1 ) ):
                    shutil.rmtree( out1 )
                    try:
                        imregrid(
                            imagename=IMAGE, template=out5,
                            output=out1
                        )
                    except RuntimeError as exc:
                        self.assertTrue(
                            'All output pixels are masked' in str(exc),
                            'Wrong error'
                        )
        self.assertTrue(len(_tb.showcache()) == 0)

    def test_axes(self):
        imagename = "test_axes.im"
        templatename = "test_axes.tmp"
        output = "test_axes.out"
        myia = _ia
        myia.fromshape(imagename, [10, 10, 10])
        exp = myia.coordsys().increment()["numeric"]
        myia.fromshape(templatename, [10, 10, 10])
        mycsys = myia.coordsys()
        mycsys.setincrement(mycsys.increment()["numeric"]/2)
        myia.setcoordsys(mycsys.torecord())
        exp[2] = mycsys.increment()["numeric"][2]
        zz = imregrid(imagename, template=templatename, output=output, axes=2, decimate=3)
        myia.open(output)
        got = myia.coordsys().increment()["numeric"]
        self.assertTrue((got == exp).all())
        myia.done()
        self.assertTrue(len(_tb.showcache()) == 0)

    def test_ref_code_preserves_position(self):
        """Test that regridding to new refcode preserves source positions"""
        shutil.copytree(ctsys.resolve(os.path.join(datapath,gim)), gim)
        orig = _ia
        myme = _me
        for rah in (0, 4, 8, 12, 16, 20):
            # image has axis units of arcmin
            ra = rah*60*15
            for decdeg in (-80, -60, -40, -20, 0, 20, 40, 68, 80):
                # image has axis units of arcmin
                dec = decdeg*60
                orig.open(gim)
                csys = orig.coordsys()
                csys.setreferencevalue([ra, dec])
                orig.setcoordsys(csys.torecord())
                csys.done()
                orig.done()
                ctype_range = ["J2000"]
                if decdeg == 0 and rah == 0:
                    ctype_range = ["J2000", "GALACTIC"]
                for ctype in ctype_range:
                    if ctype == "GALACITC":
                        orig.open(gim)
                        csys = orig.coordsys()
                        csys.setconversiontype(ctype)
                        orig.setcoordsys(csys.torecord())
                        csys.done()
                        orig.done()
                    image_id = (
                        "_ra" + str(rah) + "_dec" + str(decdeg) 
                        + "_origconv" + ctype
                    )
                    gal = "mygalactic" + image_id + ".im"
                    imregrid(gim,template="GALACTIC", output=gal)
                    orig.open(gim)
                    ofit = orig.fitcomponents(box="850,150,950,250")
                    orig.done()
                    ocl = _cl
                    ocl.add(ofit['results']['component0'])
                    orefdir = ocl.getrefdir(0)
                    galtool = _ia
                    galtool.open(gal)
                    #gfit = galtool.fitcomponents(box="1120,520,1170,570")
                    gfit = galtool.fitcomponents(mask= "'" + gal + "'> 0.001")
                    galtool.done()
                    gcl = _cl
                    gcl.add(gfit['results']['component0'])
                    grefdir = gcl.getrefdir(0)
                    print(
                        "diff", _qa.getvalue(
                            _qa.convert(
                                myme.separation(orefdir, grefdir), "arcsec"
                            )
                        )
                    )
                    self.assertTrue(
                        _qa.getvalue(
                            _qa.convert(
                                myme.separation(orefdir, grefdir), "arcsec"
                            )
                        ) < 0.003
                    )
                    rev = "back_to_J2000" + image_id + ".im"
                    imregrid(gal,template="J2000", output=rev)
                    revtool = _ia
                    revtool.open(rev)
                    rfit = revtool.fitcomponents(box="850,150,950,250")
                    revtool.done()
                    rcl = _cl
                    rcl.add(rfit['results']['component0'])
                    rrefdir = rcl.getrefdir(0)
                    rcl.done()
                    self.assertTrue(
                        _qa.getvalue(
                            _qa.convert(
                                myme.separation(orefdir, rrefdir), "arcsec"
                            )
                        ) < 1e-2
                    )
            
    def test_get(self):
        """Test using template='get' works"""
        tempfile = "xyz.im"
        myia = _ia 
        myia.fromshape(tempfile,[20,20,20])
        dicttemp = imregrid(tempfile, template="get")
        dicttemp['csys']['direction0']['crpix'] = [2.5, 2.5]
        output = "out.im"
        self.assertRaises(
            RuntimeError, imregrid, tempfile, template=dicttemp, output=output
        )
        
    def test_interpolate(self):
        """Test interpolation parameter is recognized"""
        imagename = "zzx.im"
        myia = _ia
        myia.fromshape(imagename, [30, 30])
        csys = myia.coordsys()
        incr = csys.increment()['numeric']
        incr[0] = incr[0]*0.9
        incr[1] = incr[1]*0.9
        csys.setincrement(incr)
        template = {}
        template['csys'] = csys.torecord()
        template['shap'] = myia.shape()
        myia.done()
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output="blah", interpolation="x"
        )
        output = 'blah3'
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, interpolation="cubic"
        )
        self.assertTrue(os.path.exists(output))
        
    def test_default_shape(self):
        """ Verify default shape is what users have requested, CAS-4959"""
        myia = _ia
        imagename = "myim.im"
        myia.fromshape(imagename,[20,20,20])
        template = "mytemp.im"
        myia.fromshape(template,[10,10,10])
        csys = myia.coordsys()
        csys.setreferencepixel([5,5,5])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "cas_4959_0"
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=3
        )
        myia.open(output)
        self.assertTrue((myia.shape() == [10, 10, 10]).all())
        output = "CAS_4959_1"
        imregrid(
            imagename=imagename, template=template,
            output=output, axes=[0,1], decimate=3
        )
        myia.open(output)
        self.assertTrue((myia.shape() == [10, 10, 20]).all())
        output = "CAS_4959_2"
        imregrid(
            imagename=imagename, template=template,
            output=output, axes=[2], decimate=3
        )
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 10]).all())
        
    def test_axis_recognition(self):
        """Test that imregrid recognizes axis by type, not position"""
        myia = _ia
        target = "target.im"
        myia.fromshape(target, [4,4,2,30])
        template = "template.im"
        myia.fromshape(template, [6, 6, 36, 2])
        outfile = "myout.im"
        # returns None upon success
        imregrid(imagename=target, template=template, output=outfile)
        self.assertTrue(os.path.exists(outfile))
        myia.open(outfile)
        self.assertTrue((myia.shape() == [6, 6, 2, 36]).all())
        myia.done()
        outfile = "myout1.im"
        # returns None upon success
        imregrid(
            imagename=target, template=template,
            output=outfile, axes=[0, 1], decimate=2
        )
        self.assertTrue(os.path.exists(outfile))
        myia.open(outfile)
        self.assertTrue((myia.shape() == [6, 6, 2, 30]).all())
        myia.done()

    def test_no_output_stokes(self):
        """Test rule that if input image has no stokes and template image has stokes, output image has no stokes"""
        myia = _ia
        imagename = "aa.im"
        myia.fromshape(imagename, [20, 20, 20], overwrite=True)
        template = "aa_temp.im"
        myia.fromshape(template, [20, 20, 2, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "aa.out.im"
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertFalse(myia.coordsys().findcoordinate("stokes")['return'])
        myia.done()
        
    def test_no_template_stokes(self):
        """Test rule that if input image has stokes and template image does not have stokes, output image has stokes"""
        myia = _ia
        imagename = "ab.im"
        myia.fromshape(imagename, [20, 20, 2, 20])
        template = "ab_temp.im"
        myia.fromshape(template, [20, 20, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ab.out.im"
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        stokes_info = myia.coordsys().findcoordinate("stokes")
        self.assertTrue(stokes_info['return'])
        exp_axis = 2
        self.assertTrue(stokes_info['pixel'][0] == exp_axis)
        self.assertTrue(myia.shape()[exp_axis] == 2)
        self.assertTrue(myia.coordsys().stokes() == ['I','Q'])
        myia.done()
        # specifying an output stokes length other than the input stokes length
        # is not allowed
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True,
            shape=[20, 20, 1, 20]
        )
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True,
            shape=[20, 20, 3, 20]
        )
        # specifying an output stokes length equal to the input stokes length
        # is allowed
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True,
            shape=[20, 20, 2, 20]
        )
        self.assertTrue(os.path.exists(output))

    def test_degenerate_template_stokes_axis_and_input_stokes_length_gt_0(self):
        """Verify correct behavior for the template image having a degenerate stokes axis"""
        myia = _ia
        imagename = "ac.im"
        myia.fromshape(imagename, [20, 20, 2, 20])
        template = "ac_temp.im"
        myia.fromshape(template, [20, 20, 1, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ac.out.im"
        # all input stokes in output if shape and axes not specified
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue(myia.shape()[2] == 2)
        myia.done()
        # not allowed if output stokes length different from input stokes length
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, shape=[20, 20, 1, 20], overwrite=True
        )
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, shape=[20, 20, 3, 20], overwrite=True
        )
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, shape=[20, 20, 2, 20],
            overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        
    def test_template_stokes_length_gt_1_and_input_stokes_length_gt_0(self):
        """Verify correct behavior for the template image having a stokes axis of length > 1"""
        myia = _ia
        imagename = "ad.im"
        myia.fromshape(imagename, [20, 20, 4, 20])
        template = "ad_temp.im"
        myia.fromshape(template, [20, 20, 4, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ad.out.im"
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(template)
        csys = myia.coordsys()
        csys.setstokes('XX RL LR YY')
        myia.setcoordsys(csys.torecord())
        myia.done()
        # no match between input and template stokes => not allowed
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        csys.setstokes("XX I LL RR")
        myia.open(template)
        myia.setcoordsys(csys.torecord())
        myia.done()
        # specified output stokes axis length != number of common stokes => not allowed
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True, shape=[20, 20, 3, 20]
        )
        # no output shape and number of common stokes > 0 => allowed
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        expec = ["I"]
        self.assertTrue(myia.coordsys().stokes() == expec)
        myia.done()
        
        csys.setstokes("XX I U RR")
        myia.open(template)
        myia.setcoordsys(csys.torecord())

        myia.done()
        # no output shape and number of common stokes > 0 => allowed
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        expec = ["I", "U"]
        self.assertTrue(myia.coordsys().stokes() == expec)
        self.assertTrue((myia.shape() == [20, 20, 2, 20]).all())
        myia.done()
    
    def test_no_input_spectral(self):
        """Verify if input image has no spectral axis, output will not have spectral axis"""
        myia = _ia
        imagename = "ae.im"
        myia.fromshape(imagename, [20, 20, 4])
        template = "ae_temp.im"
        myia.fromshape(template, [20, 20, 4, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ae.out.im"
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4]).all())
        myia.done()
        
    def test_no_template_spectral_axis(self):
        """Verify behavior for when template has no spectral axis, but input does"""
        myia = _ia
        imagename = "af.im"
        myia.fromshape(imagename, [20, 20, 4, 20])
        template = "af_temp.im"
        myia.fromshape(template, [20, 20, 4])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "af.out.im"
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4, 20]).all())
        myia.done()
        # Cannot explicitly specify to regrid spectral axis if template has no such axis
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True, axes=[0, 1, 3]
        )
        
    def test_degenerate_template_spectral_axis(self):
        """Verify correct behavior for when template has a degenerate spectral axis"""
        myia = _ia
        imagename = "ag.im"
        myia.fromshape(imagename, [20, 20, 4, 20])
        template = "ag_temp.im"
        myia.fromshape(template, [20, 20, 4, 1])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 900])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ag.out.im"
        # input spectral axis copied to output
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4, 20]).all())
        # should fail the image to be overwritten is open
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True, axes=[0, 1, 3]
        )
        # so close the image and then try to overwrite it
        myia.done() 
        # the spectral axis is removed from the list of axes, a warning is emitted
        # that it cannot be regridded, and the input spectral axis is copied to
        # the ouptut image
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True,
            axes=[0, 1, 3]
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4, 20]).all())
        myia.done()
    
    def test_degenerate_input_spectral_axis(self):
        """Verify correct behavior for when input has a degenerate spectral axis"""
        myia = _ia
        imagename = "ah.im"
        myia.fromshape(imagename, [20, 20, 4, 1])
        template = "ah_temp.im"
        myia.fromshape(template, [20, 20, 4, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 900])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ah.out.im"
        # when spectral axis not specified, input spectral axis is copied to
        # output spectral axis
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4, 1]).all())
        # should fail because image to be overwritten is open in the table cache
        self.assertRaises(
            RuntimeError, imregrid, imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True, axes=[0, 1, 3]
        )
        # so close it and try again
        myia.done()
        # if explicitly specified in the axis parameter, the template spectral
        # axis is copied to the output and the output's spectral axis length as
        # the same as the template's spectral axis length. The output pixel values
        # are replicated from the input image, all spectral hyperplanes in the output
        # will have identical pixel value arrays.
        # returns None upon success
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True,
            axes=[0, 1, 3]
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4, 20]).all())
        got = myia.coordsys().increment()['numeric']
        expec = csys.increment()['numeric']
        self.assertTrue((got == expec).all())
        myia.done()
    
    def test_bad_shape(self):
        """ Verify that bad shape specification results in exception"""
        myia = _ia
        imagename = "aj.im"
        myia.fromshape(imagename, [20, 20, 1, 1])
        template = "aj_temp.im"
        myia.fromshape(template, [20, 20, 1, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 900])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "aj.out.im"
        self.assertFalse(
            imregrid(
                imagename=imagename,
                template=template, output=output, decimate=5,
                overwrite=True, shape=[20, 20, 20, 1]
            )
        )
    
    def test_nested_image(self):
        """ Verify that one image which lies completely inside the other will not cause failure"""
        myia = _ia
        imagename = "ak.im"
        myia.fromshape(imagename, [20, 20])
        csys = myia.coordsys()
        csys.setreferencevalue([1800, 1800])
        myia.setcoordsys(csys.torecord())
        template = "ak_temp.im"
        myia.fromshape(template, [4, 4])
        csys = myia.coordsys()
        csys.setreferencevalue([1800, 1800])
        csys.setincrement([-0.9, 0.9])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ak.out.im"
        self.assertTrue(
            imregrid(
                imagename=imagename,
                template=template, output=output, decimate=5,
                overwrite=True
            )
        )
        self.assertTrue(
            imregrid(
                imagename=template,
                template=imagename, output=output, decimate=5,
                overwrite=True
            )
        )
    
    def test_regrid_galactic(self):
        """Verify fix for CAS-5534"""
        myia = _ia
        myia.open(datapath + "ngc5921.clean.image")
        csys = myia.coordsys()
        csys.setreferencecode('GALACTIC', type='direction', adjust=True)
        zz = myia.regrid(outfile='gal_regrid.image', shape=[300, 300, 1, 46], csys=csys.torecord(), overwrite=True)  
        myia.open(datapath + "gal_regrid.image")
        self.assertTrue(numpy.max(numpy.abs(zz.getchunk() - myia.getchunk())) < 1e-8)
        self.assertTrue((zz.getchunk(getmask=True) == myia.getchunk(getmask=True)).all())
        myia.done()
        zz.done()

    def test_linear_overlap(self):
        """Test that overlapping linear coordinates works, CAS-5767"""
        myia = _ia
        myia.open(datapath + "lin_template.im")
        csys = myia.coordsys().torecord()
        shape = myia.shape()
        myia.done()
        myia.open(datapath + "lin_source.im")
        xx = myia.regrid(axes=[0,1], csys=csys, shape=shape)
        myia.done()
        self.assertTrue((xx.shape() == shape).all())
        xx.done()
        
    def test_template_stokes_length_and_input_stokes_length_gt_1(self):
        myia = _ia
        my_image = numpy.zeros([128,128,32,4])
        os.system("rm -rf fake.image")
        myia.fromarray(
            outfile="fake.image",
            pixels=my_image,
            overwrite=True
        )
        myia.close()

        my_image = numpy.zeros([128,128,32,3])
        os.system("rm -rf fake_2.image")
        myia.fromarray(
            outfile="fake_2.image",
            pixels=my_image,
            overwrite=True
        )
        myia.close()

        output = "dummy.image"
        os.system("rm -rf " + output)
        imregrid(
            imagename = "fake.image",
            template="fake_2.image",
            output=output,
            axes=[0,1,2]
        )
        myia.open(output)
        self.assertTrue(myia.shape()[3] == 4)
        myia.done()
        
    def test_decimate(self):
        """imregrid, test too high a value for decimate throws exception - CAS-5313"""
        myia = _ia
        myia.fromshape("", [10, 10, 10])
        self.assertRaises(
            Exception, myia.regrid, axes=[0, 1], csys=myia.coordsys().torecord()
        )
        regridded = myia.regrid(
            axes=[0, 1], csys=myia.coordsys().torecord(), decimate=3
        )
        self.assertTrue(regridded)
        regridded.done()
        # decimate doen't matter for non-direction axis regridding
        regridded = myia.regrid(
            axes=[2], csys=myia.coordsys().torecord()
        )
        self.assertTrue(regridded)
        regridded.done()
        
    def test_complex(self):
        """Test regridding a complex image, CAS-1390"""
        myia = _ia
        for i in (0, 1):
            myia.open(datapath + "real1.im")
            if i == 1:
                gg = myia.getchunk()
                myia.fromarray("", gg, myia.coordsys().torecord(), type='d')
            realpart = myia.getchunk()
            csys = myia.coordsys()
            csys.setincrement([-0.9, 0.9])
            rrg = myia.regrid(csys=csys.torecord())
            rrgpart = rrg.getchunk()
            rrg.done()
            myia.open(datapath + "imag1.im")
            if i == 1:
                gg = myia.getchunk()
                myia.fromarray("", gg, myia.coordsys().torecord(), type='d')
            imagpart = myia.getchunk()
            irg = myia.regrid(csys=csys.torecord())
            irgpart = irg.getchunk()
            irg.done()
            compltype = 'c'
            if i == 1:
                compltype = 'cd'
            myia.fromshape("", myia.shape(), type=compltype)
            comp = myia.getchunk()
            comp = realpart + imagpart*1j
            myia.putchunk(comp)
            crg = myia.regrid(csys=csys.torecord())
            crgpart = crg.getchunk()
            crg.done()
            myia.done()
            self.assertTrue((crgpart == rrgpart + irgpart*1j).all())        
        
    def test_multibeam(self):
        """test multibeams cannot be regridded"""
        myia = _ia
        myia.fromshape("",[100,100,20])
        myia.setrestoringbeam("20arcsec", "20arcsec", "0deg", channel=0)
        myia.setrestoringbeam("30arcsec", "30arcsec", "0deg", channel=1)
        self.assertRaises(
            Exception, myia.regrid, "", shape=[10, 10, 10],
            csys=myia.coordsys().torecord()
        )
        self.assertRaises(
            Exception, myia.regrid, "", shape=[10, 10, 10],
            csys=myia.coordsys().torecord(), axes=2
        )
        bb = myia.regrid(
            "", shape=[50, 50, 20],
            csys=myia.coordsys().torecord(), axes=[0, 1]
        )
        bb.done()
        myia.done()
        
    def test_CAS_8345(self):
        """verify fix to CAS-8345, channels not replicating properly"""
        myia = _ia
        iname = "CAS_8345.im"
        myia.fromshape(iname, [20, 20, 10])
        myia.addnoise()
        sname = "8345.sub"
        sub = myia.subimage(sname, region=_rg.box([0,0,5], [19,19,5]), mask=iname + ">0")
        expec = sub.getchunk()
        expm = sub.getchunk(getmask=True)
        myia.done()
        sub.done()
        out = "8345.out"
        imregrid(imagename=sname, template=iname, output=out, axes=[2], decimate=0)
        myia.open(out)
        self.assertTrue((myia.shape() == [20, 20, 10]).all(), "Shape error")
        got = myia.getchunk()
        gotm = myia.getchunk(getmask=True)
        for i in xrange(20):
            for j in xrange(20):
                self.assertTrue((got[i,j,:] == expec[i,j]).all(), "incorrect values" )
                self.assertTrue((gotm[i,j,:] == expm[i,j]).all(), "incorrect mask" )

        myia.done()
                
    def test_history(self):
        """Test history writing"""
        myia = _ia
        imagename = "zz.im"
        myia.fromshape(imagename, [20, 20, 10])
        myia.done()
        outfile = "zz_out.im"
        self.assertTrue(
            imregrid(
                imagename=imagename, output=outfile,
                template="GALACTIC", decimate=5
            )
        )
        myia.open(outfile)
        msgs = myia.history()
        myia.done()
        teststr = "version"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")
        teststr = "imregrid"
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
        
def suite():
    return [imregrid_test]
    
