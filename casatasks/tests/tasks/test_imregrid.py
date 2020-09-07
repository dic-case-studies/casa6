import os
import shutil
import numpy
import unittest

from casatools import ctsys, image, regionmanager, coordsys, measures, componentlist, table, quanta
from casatasks import imregrid, imstat

_tb = table( )
_ia = image( )
_rg = regionmanager( )
_cs = coordsys( )
_qa = quanta( )

IMAGE = 'image.im'
gim = "gaussian_source.im"

total = 0
fail  = 0
current_test =""
stars = "*************"

datapath = ctsys.resolve('unittest/imregrid/')

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
    print( )
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
        self._myia = image( )
    
    def tearDown(self):
        self._myia.done()
        
        for i in (IMAGE, out1, out2, out3, out4, out5, out6):
            if (os.path.exists(i)):
                os.system('rm -rf ' + i)
        
        self.assertTrue(len(_tb.showcache()) == 0)
        
    def test1(self):    
        myia = self._myia  
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
                if p1['value']['value'] != p2['value']['value']: raise Exception(p1['value']['value'] + ' != ' + p2['value']['value'])
                if p1['value']['unit'] != p2['value']['unit']: raise Exception(p1['value']['unit'] + ' != ' + p2['value']['unit'])
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
                        outim=imregrid(imagename = IMAGE,
                                       template = out5,
                                       output = out1)
                    except RuntimeError as exc:
                        if not 'All output pixels are masked' in str(exc):
                            raise

        self.assertTrue(len(_tb.showcache()) == 0)

    def test_axes(self):
        imagename = "test_axes.im"
        templatename = "test_axes.tmp"
        output = "test_axes.out"
        myia = self._myia
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
        shutil.copytree(os.path.join(datapath,gim), gim)
        orig = image( )
        myme = measures( )
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
                    ocl = componentlist( )
                    ocl.add(ofit['results']['component0'])
                    orefdir = ocl.getrefdir(0)
                    galtool = image( )
                    galtool.open(gal)
                    #gfit = galtool.fitcomponents(box="1120,520,1170,570")
                    gfit = galtool.fitcomponents(mask= "'" + gal + "'> 0.001")
                    galtool.done()
                    gcl = componentlist( )
                    gcl.add(gfit['results']['component0'])
                    grefdir = gcl.getrefdir(0)
                    print("diff", _qa.getvalue(
                            _qa.convert(
                                myme.separation(orefdir, grefdir), "arcsec"
                            )
                        ))
                    self.assertTrue(
                        _qa.getvalue(
                            _qa.convert(
                                myme.separation(orefdir, grefdir), "arcsec"
                            )
                        ) < 0.003
                    )
                    rev = "back_to_J2000" + image_id + ".im"
                    imregrid(gal,template="J2000", output=rev)
                    revtool = image( )
                    revtool.open(rev)
                    rfit = revtool.fitcomponents(box="850,150,950,250")
                    revtool.done()
                    rcl = componentlist( )
                    rcl.add(rfit['results']['component0'])
                    rrefdir = rcl.getrefdir(0)
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
        myia = self._myia 
        myia.fromshape(tempfile,[20,20,20])
        dicttemp = imregrid(tempfile, template="get")
        dicttemp['csys']['direction0']['crpix'] = [2.5, 2.5]
        output = "out.im"
        with self.assertRaises(RuntimeError):
            imregrid (tempfile, template=dicttemp, output=output)
        
    def test_interpolate(self):
        """Test interpolation parameter is recognized"""
        imagename = "zzx.im"
        myia = self._myia
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
        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output="blah", interpolation="x"
            )
        output = "blah3"
        imregrid(
            imagename=imagename, template=template,
            output=output, interpolation="cubic"
        )
        self.assertTrue(os.path.exists(output))
        
    def test_default_shape(self):
        """ Verify default shape is what users have requested, CAS-4959"""
        myia = self._myia
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
        myia = self._myia
        target = "target.im"
        myia.fromshape(target, [4,4,2,30])
        template = "template.im"
        myia.fromshape(template, [6, 6, 36, 2])
        outfile = "myout.im"
        imregrid(imagename=target, template=template, output=outfile)
        self.assertTrue(os.path.exists(outfile))
        myia.open(outfile)
        self.assertTrue((myia.shape() == [6, 6, 2, 36]).all())
        myia.done()
        outfile = "myout1.im"
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
        myia = self._myia
        imagename = "aa.im"
        myia.fromshape(imagename, [20, 20, 20], overwrite=True)
        template = "aa_temp.im"
        myia.fromshape(template, [20, 20, 2, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "aa.out.im"
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
        myia = self._myia
        imagename = "ab.im"
        myia.fromshape(imagename, [20, 20, 2, 20])
        template = "ab_temp.im"
        myia.fromshape(template, [20, 20, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ab.out.im"
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
        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output=output, decimate=5, overwrite=True,
                shape=[20, 20, 1, 20]
            )

        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output=output, decimate=5, overwrite=True,
                shape=[20, 20, 3, 20]
            )

        # specifying an output stokes length other than the input stokes length
        # is allowed
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True,
            shape=[20, 20, 2, 20]
        )
        self.assertTrue(os.path.exists(output))
        
    def test_degenerate_template_stokes_axis_and_input_stokes_length_gt_0(self):
        """Verify correct behavior for the template image having a degenerate stokes axis"""
        myia = self._myia
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
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue(myia.shape()[2] == 2)
        myia.done()
        # not allowed if output stokes length different from input stokes length
        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output=output, decimate=5, shape=[20, 20, 1, 20],
                overwrite=True
            )

        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output=output, decimate=5, shape=[20, 20, 3, 20],
                overwrite=True
            )

        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, shape=[20, 20, 2, 20],
            overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        
    def test_template_stokes_length_gt_1_and_input_stokes_length_gt_0(self):
        """Verify correct behavior for the template image having a stokes axis of length > 1"""
        myia = self._myia
        imagename = "ad.im"
        myia.fromshape(imagename, [20, 20, 4, 20])
        template = "ad_temp.im"
        myia.fromshape(template, [20, 20, 4, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ad.out.im"
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
        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output=output, decimate=5, overwrite=True
            )

        csys.setstokes("XX I LL RR")
        myia.open(template)
        myia.setcoordsys(csys.torecord())
        myia.done()
        # specified output stokes axis length != number of common stokes => not allowed
        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output=output, decimate=5, overwrite=True,
                shape=[20, 20, 3, 20]
            )
        # no output shape and number of common stokes > 0 => allowed
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
        myia = self._myia
        imagename = "ae.im"
        myia.fromshape(imagename, [20, 20, 4])
        template = "ae_temp.im"
        myia.fromshape(template, [20, 20, 4, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 1500])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "ae.out.im"
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
        myia = self._myia
        imagename = "af.im"
        myia.fromshape(imagename, [20, 20, 4, 20])
        template = "af_temp.im"
        myia.fromshape(template, [20, 20, 4])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "af.out.im"
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4, 20]).all())
        myia.done()
        # Cannot explicitly specify to regrid spectral axis if template has no such axis
        self.assertRaises(Exception, imregrid,
                imagename=imagename, template=template,
                output=output, decimate=5, overwrite=True,
                axes=[0, 1, 3]
        )
        
    def test_degenerate_template_spectral_axis(self):
        """Verify correct behavior for when template has a degenerate spectral axis"""
        myia = self._myia
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
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4, 20]).all())
        # should fail the image to be overwritten is open
        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output=output, decimate=5, overwrite=True,
                axes=[0, 1, 3]
            )

        # so close the image and then try to overwrite it
        myia.done() 
        # the spectral axis is removed from the list of axes, a warning is emitted
        # that it cannot be regridded, and the input spectral axis is copied to
        # the ouptut image
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
        myia = self._myia
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
        imregrid(
            imagename=imagename, template=template,
            output=output, decimate=5, overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        myia.open(output)
        self.assertTrue((myia.shape() == [20, 20, 4, 1]).all())
        # should fail because image to be overwritten is open in the table cache
        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename, template=template,
                output=output, decimate=5, overwrite=True,
                axes=[0, 1, 3]
            )

        # so close it and try again
        myia.done()
        # if explicitly specified in the axis parameter, the template spectral
        # axis is copied to the output and the output's spectral axis length as
        # the same as the template's spectral axis length. The output pixel values
        # are replicated from the input image, all spectral hyperplanes in the output
        # will have identical pixel value arrays.
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
        myia = self._myia
        imagename = "aj.im"
        myia.fromshape(imagename, [20, 20, 1, 1])
        template = "aj_temp.im"
        myia.fromshape(template, [20, 20, 1, 20])
        csys = myia.coordsys()
        csys.setincrement([-0.9, 0.9, 1, 900])
        myia.setcoordsys(csys.torecord())
        myia.done()
        output = "aj.out.im"
        with self.assertRaises(RuntimeError):
            imregrid(
                imagename=imagename,
                template=template, output=output, decimate=5,
                overwrite=True, shape=[20, 20, 20, 1]
            )
    
    def test_nested_image(self):
        """ Verify that one image which lies completely inside the other will not cause failure"""
        myia = self._myia
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
        imregrid(
            imagename=imagename,
            template=template, output=output, decimate=5,
            overwrite=True
        )
        self.assertTrue(os.path.exists(output))
        imregrid(
            imagename=template,
            template=imagename, output=output, decimate=5,
            overwrite=True
        )
        self.assertTrue(os.path.exists(output))

    def test_template_stokes_length_and_input_stokes_length_gt_1(self):
        myia = self._myia
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
        
    def test_CAS_8345(self):
        """verify fix to CAS-8345, channels not replicating properly"""
        myia = self._myia
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
        for i in range(20):
            for j in range(20):
                self.assertTrue((got[i,j,:] == expec[i,j]).all(), "incorrect values" )
                self.assertTrue((gotm[i,j,:] == expm[i,j]).all(), "incorrect mask" )

        myia.done()
                
    def test_history(self):
        """Test history writing"""
        myia = self._myia
        imagename = "zz.im"
        myia.fromshape(imagename, [20, 20, 10])
        myia.done()
        outfile = "zz_out.im"
        imregrid(
            imagename=imagename, output=outfile,
            template="GALACTIC", decimate=5
        )
        self.assertTrue(os.path.exists(outfile))
        myia.open(outfile)
        msgs = myia.history()
        myia.done()
        teststr = "version"
        self.assertTrue(teststr in msgs[-2], "'" + teststr + "' not found")
        teststr = "imregrid"
        self.assertTrue(teststr in msgs[-1], "'" + teststr + "' not found")
        
def suite():
    return [imregrid_test]

if __name__ == '__main__':
    unittest.main()
