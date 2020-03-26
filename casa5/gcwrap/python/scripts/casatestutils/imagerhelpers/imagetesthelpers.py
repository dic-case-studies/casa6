
##################################       Imports        ##################################
############################################################################################

import os
import sys
import time
from functools import wraps
import fnmatch
import logging
import filecmp
import unittest
import pickle
import numbers
import operator
import subprocess
import numpy
import six

casa5 = False
casa6 = False
__bypass_parallel_processing = 0

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:

    # CASA 6
    logging.debug("Importing CASAtools")
    import casatools
    logging.debug("Importing CASAtasks")
    import casatasks
    _tb = casatools.table()
    _tbt = casatools.table()
    _ia  = casatools.image()
    from casatasks import casalog
    casa6 = True

else:

    # CASA 5
    logging.debug("Import casa6 errors. Trying CASA5...")
    from __main__ import default
    from taskinit import tbtool, mstool, iatool
    from taskinit import *
    from casa_stack_manip import stack_find, find_casa
    try:
        from mpi4casa.MPIEnvironment import MPIEnvironment
        if not MPIEnvironment.is_mpi_enabled:
            __bypass_parallel_processing = 1
    except ImportError:
        print("MPIEnvironment not Enabled")
    _tb = tbtool()
    _tbt = tbtool()
    _ia = iatool()
    casa = find_casa()
    if casa.has_key('state') and casa['state'].has_key('init_version') and casa['state']['init_version'] > 0:
        casaglobals=True
        casac = stack_find("casac")
        casalog = stack_find("casalog")
    casa5 = True

############################################################################################
##################################       imagerhelpers       ###############################
############################################################################################
class TestHelpers:
    def check_model(self, msname=""):
        """Check hasmodcol, modsum, hasvirmod"""
        logging.debug("Executing: check_model(msname={})".format(msname))
        hasmodcol = False
        modsum=0.0
        hasvirmod = False
        _tb.open( msname )
        hasmodcol = (  (_tb.colnames()).count('MODEL_DATA')>0 )
        if hasmodcol:
            model_data = _tb.getcol('MODEL_DATA')
            modsum = model_data.sum()
        _tb.close()
        _tb.open( msname+'/SOURCE' )
        keys = _tb.getkeywords()
        if len(keys)>0:
            hasvirmod=True
        _tb.close()
        _tb.open( msname )
        keys = _tb.getkeywords()
        for key in keys:
            if key.count("model_")>0:
                hasvirmod=True
        _tb.close()
        logging.info("MS Name: {}, modelcol= {},  modsum = {}, virmod = {}".format( msname, hasmodcol, modsum, hasvirmod ))
        return hasmodcol, modsum, hasvirmod

    def get_max(self, imname):
        """Get Image max"""
        logging.debug("Executing: get_max(imname={})".format(imname))
        _ia.open(imname)
        stat = _ia.statistics()
        _ia.close()
        logging.debug("stat['max'] = {}".format(stat['max']))
        logging.debug("stat['maxpos'] = {}".format(stat['maxpos']))
        return stat['max'],stat['maxpos']

    def get_pix(self, imname,pos):
        """Get Image val"""
        _ia.open(imname)
        apos = _ia.pixelvalue(pos)
        _ia.close()
        if apos == {}:
            return None
        else:
            return apos['value']['value']

    def get_pixmask(self, imname,pos):
        """Get Image Mask val"""
        _ia.open(imname)
        apos = _ia.pixelvalue(pos)
        _ia.close()
        if apos == {}:
            return None
        else:
            return apos['mask']

    def check_beam_compare(self, image1, image2, op=operator.le):
        """Compare all plane of cube beam image1 operator op than image1"""
        _ia.open(image1)
        nchan = _ia.shape()[3]
        beam1 = numpy.zeros(nchan)
        for k in range(nchan):
            beam1[k]= _ia.beamarea(k,0)['arcsec2']
        _ia.close()
        _ia.open(image2)
        if(nchan != _ia.shape()[3]):
            return False
        beam2 = numpy.zeros(nchan)
        for k in range(nchan):
            beam2[k] = _ia.beamarea(k,0)['arcsec2']
        _ia.close()
        return numpy.alltrue(op(beam1, beam2))

    def image_exists(self, imname):
        """ Image exists """
        return os.path.exists(imname)

    def get_peak_res(self, summ):
        """Get Peak Res"""
        if summ.has_key('summaryminor'):
            reslist = summ['summaryminor'][1,:]
            peakres = reslist[ len(reslist)-1 ]
        else:
            peakres = None
        return peakres

    def check_peak_res(self, summ,correctres, epsilon=0.05):
        """Check Peak Res"""
        peakres = get_peak_res(summ)
        out = True
        if correctres == None and peakres != None:
            out = False
            return out,peakres
        if correctres != None and peakres == None:
            out = False
            return out,peakres
        if out==True and peakres != None:
            if abs(correctres - peakres)/abs(correctres) > epsilon:
                out=False
                return out,peakres
        return out,peakres

    def get_mod_flux(self, summ):
        """Get Mod Flux"""
        if summ.has_key('summaryminor'):
            modlist = summ['summaryminor'][2,:]
            modflux = modlist[ len(modlist)-1 ]
        else:
            modflux = None
        return modflux

    def check_mod_flux(self, summ,correctmod, epsilon=0.05):
        """Check Mod Flux"""
        modflux = get_mod_flux(summ)
        out = True
        if correctmod == None and modflux != None:
            out = False
            return out,modflux
        if correctmod != None and modflux == None:
            out = False
            return out,modflux
        if out==True and modflux != None:
            if abs(correctmod - modflux)/abs(correctmod) > epsilon:
                out=False
                return out,modflux
        return out,modflux

    def get_iter_done(self, summ):
        """Get Iterdone"""
        if summ.has_key('iterdone'):
            iters = summ['iterdone']
        else:
            iters = None
        return iters

    def verdict(self, boolval):
        """Return the string 'Pass' if boolean is True, Else return string 'Fail'"""
        return "Pass" if boolval else "Fail"

    ############################################################################################
    ##############################       imagerhelpers: Checks       ###########################
    ############################################################################################

    def check_ret(self, summ,correctres,correctmod,epsilon = 0.05, testname ="check_ret"):
        pstr = ''
        retres, peakres = TestHelpers().check_peak_res(summ, correctres, epsilon)
        retmod, modflux = TestHelpers().check_mod_flux(summ, correctmod, epsilon)
        pstr_peak =  "[ {} ] PeakRes is  {}  ( {} : should be  {}, Epsilon: {} + )\n".format(testname, str(peakres), TestHelpers().verdict(retres) , str(correctres), str(epsilon))
        pstr_mod  =  "[ {} ] Modflux is  {}  ( {} : should be  {}, Epsilon: {} + )\n".format(testname, str(modflux), TestHelpers().verdict(retmod) , str(correctmod), str(epsilon))
        pstr =  pstr_peak + pstr_mod
        logging.info(pstr)
        if retres == False or retmod == False:
            return False, pstr
        else:
            return True, pstr

    def check_val(self, val, correctval, valname='Value', exact=False, epsilon=0.05, testname = "check_val"):
        pstr = ''
        out = True
        if numpy.isnan(val) or numpy.isinf(val):
            out = False
        if correctval == None and val != None:
            out = False
        if correctval != None and val == None:
            out = False
        if out==True and val != None:
            if exact==True:
                if correctval != val:
                    out = False
            else:
                if abs(correctval - val)/abs(correctval) > epsilon:
                    out=False
        if exact == True:
            pstr = "[ {} ] {} is {} ( {} : should be {}, Exact: True )\n".format(testname, valname, str(val), TestHelpers().verdict(out), str(correctval) )
        else:
            pstr = "[ {} ] {} is {} ( {} : should be {}, Epsilon: {})\n".format(testname, valname, str(val), TestHelpers().verdict(out), str(correctval), str(epsilon) )
        logging.info(pstr)
        return out, pstr

    def check_val_less_than(self, val, bound, valname='Value',testname ="check_val_less_than"):
        pstr = ''
        out = True
        if numpy.isnan(val) or numpy.isinf(val):
            out = False
        if bound == None and val != None:
            out = False
        if bound != None and val == None:
            out = False
        if out == True and val != None:
            if val > bound:
                out = False
        pstr = "[ {} ] {} is {} ( {} : should be less than {} )\n".format(testname, valname, str(val), TestHelpers().verdict(out), str(bound))
        logging.info(pstr)
        return out, pstr

    def check_val_greater_than(self, val, bound, valname='Value',testname ="check_val_greater_than"):
        pstr = ''
        out = True
        if numpy.isnan(val) or numpy.isinf(val):
            out = False
        if bound == None and val != None:
            out = False
        if bound != None and val == None:
            out = False
        if out==True and val != None:
            if val < bound:
                out = False
        pstr = "[ {} ] {} is {} ( {} : should be greater than {} )\n".format(testname, valname, str(val), TestHelpers().verdict(out), str(bound))
        logging.info(pstr)
        return out, pstr

    def check_ims(self, imlist, truth, testname="check_ims"):
        pstr = ''
        imex = []
        out = True
        for imname in imlist:
            ondisk = TestHelpers().image_exists(imname)
            imex.append(ondisk)
            if ondisk != truth:
                out = False
        pstr = "[ {} ] Image made : {} =  {} ( {} : should all be {} )\n".format(testname, str(imlist), str(imex), TestHelpers().verdict(out), str(truth))
        logging.info(pstr)
        return pstr

    def check_keywords(self, imlist, testname="check_keywords"):
        """
        Keyword related checks (presence/absence of records and entries in these records,
        in the keywords of the image table).
        :param imlist: names of the images produced by a test execution.
        :returns: the usual (test_imager_helper) string with success/error messages.
        """
        # Keeping the general approach. This is fragile!
        # accumulator of error strings
        pstr = ''
        for imname in imlist:
            if os.path.exists(imname):
                issues = TestHelpers().check_im_keywords(imname, check_misc=True, check_extended=True)
                if issues:
                    pstr += '[{0}] {1}: {2}'.format(testname, imname, issues)
        if not pstr:
            pstr += 'All expected keywords in imageinfo, miscinfo, and coords found.\n'
        return pstr

    def check_im_keywords(self, imname, check_misc=True, check_extended=True):
        """
        Checks several lists of expected and forbidden keywords and entries of these
        keywords.
        Forbidden keywords lists introduced with CAS-9231 (prevent duplication of
        TELESCOP and OBJECT).
        Note that if imname is the top level of a refconcat image, there's no table to open
        to look for its keywords. In these cases nothing is checked. We would not have the
        'imageinfo' keywords, only the MiscInfo that goes in imageconcat.json and I'm not
        sure yet how that one is supposed to behave.
        Tests should check the 'getNParts() from imname' to make sure the components of
        the refconcat image exist, have the expected keywords, etc.
        :param imname: image name (output image from tclean)
        :param check_misc: whether to check miscinfo in addition to imageinfo'
        :param check_extended: can leave enabled for images other than .tt?, .alpha, etc.
        :returns: the usual (test_imager_helper) string with success/error messages.
        Errors marked with '(Fail' as per self.verdict().
        """
        try:
            _tbt.open(imname)
            keys = _tbt.getkeywords()
        except RuntimeError as exc:
            if os.path.isfile(os.path.join(os.path.abspath(imname), 'imageconcat.json')):
                # Looks like a refconcat image, nothing to check
                #return ''
                # make a bit more informative
                pstr = 'Looks like it is a refconcat image. Skipping the imageinfo keywords check.\n'
                return pstr
            else:
                pstr = 'Cannot open image table to check keywords: {0}\n'.format(imname)
                return pstr
        finally:
            _tbt.close()
        pstr = ''
        if len(keys) <= 0:
            pstr += ('No keywords found ({0})\n'.format(TestHelpers().verdict(False)))
            return pstr
        # Records that need to be present
        imageinfo = 'imageinfo'
        miscinfo = 'miscinfo'
        coords = 'coords'
        mandatory_recs = [imageinfo, coords]
        if check_misc:
            mandatory_recs.append(miscinfo)
        for rec in mandatory_recs:
            if rec not in keys:
                pstr += ('{0} record not found ({1})\n'.format(rec, TestHelpers().verdict(False)))
        if len(pstr) > 0:
            return pstr
        mandatory_imageinfo = ['objectname', 'imagetype']
        pstr += TestHelpers().check_expected_entries(mandatory_imageinfo, imageinfo, keys)
        if check_misc:
            if check_extended:
                mandatory_miscinfo = ['INSTRUME', 'distance']
                pstr += TestHelpers().check_expected_entries(mandatory_miscinfo, miscinfo, keys)
            forbidden_miscinfo = ['OBJECT', 'TELESCOP']
            pstr += TestHelpers().check_forbidden_entries(forbidden_miscinfo, miscinfo, keys)
        mandatory_coords = ['telescope']
        pstr += TestHelpers().check_expected_entries(mandatory_coords, coords, keys)
        return pstr
    def check_expected_entries(self, entries, record, keys):
        pstr = ''
        for entry in entries:
            if entry not in keys[record]:
                pstr += ('entry {0} not found in record {1} ({2})\n'.format(entry, record, TestHelpers().verdict(False)))
            else:
                # TODO: many tests leave 'distance' empty. Assume that's acceptable...
                if entry != 'distance' and not keys[record][entry]:
                    pstr += ('entry {0} is found in record {1} but it is empty ({2})\n'.format(entry, record, TestHelpers().verdict(False)))
        return pstr

    def check_forbidden_entries(self, entries, record, keys):
        pstr = ''
        for entry in entries:
            if entry in keys[record]:
                pstr += ('entry {0} should not be in record {1} ({2})\n'.format(entry, record, TestHelpers().verdict(False)))
        return pstr

    def check_pix_val(self, imname, theval=0, thepos=[0, 0, 0, 0], exact=False, epsilon=0.05, testname="check_pix_val"):
        pstr = ''
        readval = TestHelpers().get_pix(imname, thepos)
        res = True
        if readval == None:
            res = False
        elif numpy.isnan(readval) or numpy.isinf(readval):
            res = False
        else:
            if abs(theval) > epsilon:
                if exact == False:
                    if abs(readval - theval)/abs(theval) > epsilon:
                        res = False
                    else:
                        res = True
                else:
                    if abs(readval - theval) > 0.0:
                        res = False
                    else:
                        res = True
            else:  ## this is to guard against exact zero... sort of.
                if abs(readval - theval) > epsilon:
                    res = False
                else:
                    res = True
        pstr = "[ {} ] {} : Value is {} at {} ( {} : should be {} , Epsilon: {})\n".format(testname, imname, str(readval), str(thepos), TestHelpers().verdict(res), str(theval), str(epsilon))
        logging.info(pstr)
        return pstr

    def check_pixmask(self, imname, theval=True, thepos=[0, 0, 0, 0], testname="check_pixmask"):
        pstr = ''
        readval = get_pixmask(imname, thepos)
        res = True
        if readval == None:
            res = False
        elif numpy.isnan(readval) or numpy.isinf(readval) or type(readval) != bool:
            res = False
        else:
            if readval == theval:
                res = True
            else:
                res = False
        pstr = "[ {} ] {} : Mask is {} at {} ( {} : should be {} )\n".format(testname, imname, str(readval), str(thepos), TestHelpers().verdict(res), str(theval))
        logging.info(pstr)
        return pstr

    def check_ref_freq(self, imname, theval=0, epsilon=0.05, testname="check_ref_freq"):
        pstr = ''
        retres = True
        _ia.open(imname)
        csys = _ia.coordsys()
        _ia.close()
        reffreq = csys.referencevalue()['numeric'][3]
        if  abs(reffreq - theval)/theval > epsilon:
            retres = False
        else:
            retres = True
        pstr = "[ {} ] Ref-Freq is {} ( {} : should be {} )\n".format(testname, str(reffreq), TestHelpers().verdict(retres), str(theval))
        logging.info(pstr)
        return pstr

    def check_imexist(self, imgexist):
        pstr = ''
        if imgexist != None:
            if type(imgexist) == list:
                pstr += TestHelpers().check_ims(imgexist, True)
                print("pstr after checkims = {}".format(pstr))
                pstr += TestHelpers().check_keywords(imgexist)
                print("pstr after check_keywords = {}".format(pstr))
        return pstr

    def check_imexistnot(self, imgexistnot):
        pstr = ''
        if imgexistnot != None:
            if type(imgexistnot) == list:
                pstr += TestHelpers().check_ims(imgexistnot, False)
        return pstr

    def check_imval(self, imgval, epsilon=0.05):
        pstr = ''
        if imgval != None:
            if type(imgval) == list:
                for ii in imgval:
                    if type(ii) == tuple and len(ii) == 3:
                        pstr += TestHelpers().check_pix_val(ii[0], ii[1], ii[2], epsilon=epsilon)
        return pstr

    def check_imvalexact(self, imgvalexact, epsilon=0.05):
        pstr = ''
        if imgvalexact != None:
            if type(imgvalexact) == list:
                for ii in imgvalexact:
                    if type(ii) == tuple and len(ii) == 3:
                        pstr += TestHelpers().check_pix_val(ii[0], ii[1], ii[2], exact=True, epsilon=epsilon)
        return pstr

    def check_immask(self, imgmask):
        pstr = ''
        if imgmask != None:
            if type(imgmask) == list:
                for ii in imgmask:
                    if type(ii) == tuple and len(ii) == 3:
                        pstr += TestHelpers().check_pixmask(ii[0], ii[1], ii[2])
        return pstr

    def check_tabcache(self, tabcache, testname="check_tabcache"):
        pstr = ''
        if tabcache == True:
            opentabs = _tb.showcache()
            if len(opentabs) > 0:
                pstr += "["+ testname +"] " + TestHelpers().verdict(False) + ": Found open tables after run \n"
        return pstr

    def check_stopcode(self, stopcode, ret, testname="check_stopcode"):
        pstr = ''
        if stopcode != None:
            if type(stopcode) == int:
                stopstr = "["+ testname +"] Stopcode is " + str(ret['stopcode']) + " (" + TestHelpers().verdict(ret['stopcode'] == stopcode)  +  " : should be " + str(stopcode) + ")\n"
                print(stopstr)
                pstr += stopstr
        return pstr

    def check_reffreq(self, reffreq):
        pstr = ''
        if reffreq != None:
            if type(reffreq) == list:
                for ii in reffreq:
                    if type(ii) == tuple and len(ii) == 2:
                        pstr += TestHelpers().check_ref_freq(ii[0], ii[1])
        return pstr

    def checkall(self, ret=None, peakres=None, modflux=None, iterdone=None, nmajordone=None, imgexist=None, imgexistnot=None, imgval=None, imgvalexact=None, imgmask=None, tabcache=True, stopcode=None, reffreq=None, epsilon=0.05):
        """
            ret=None,
            peakres=None, # a float
            modflux=None, # a float
            iterdone=None, # an int
            nmajordone=None, # an int
            imgexist=None,  # list of image names
            imgexistnot=None, # list of image names
            imgval=None,  # list of tuples of (imagename,val,pos)
            imgvalexact=None, # list of tuples of (imagename,val,pos)
            imgmask=None,  #list of tuples to check mask value
            tabcache=True,
            stopcode=None,
            reffreq=None # list of tuples of (imagename, reffreq)
        """
        pstr = "[ checkall ] \n"
        if ret != None and type(ret) == dict:
            try:
                if peakres != None:
                    pstr += TestHelpers().check_val(val=get_peak_res(ret), correctval=peakres, valname="peak res")
                if modflux != None:
                    pstr += TestHelpers().check_val(val=get_mod_flux(ret), correctval=modflux, valname="mod flux")
                if iterdone != None:
                    pstr += TestHelpers().check_val(val=ret['iterdone'], correctval=iterdone, valname="iterdone", exact=True)
                if nmajordone != None:
                    pstr += TestHelpers().check_val(val=ret['nmajordone'], correctval=nmajordone, valname="nmajordone", exact=True)
            except Exception as e:
                logging.info(ret)
                raise
        logging.info("Epsilon: {}".format(epsilon))
        pstr += TestHelpers().check_imexist(imgexist)
        pstr += TestHelpers().check_imexistnot(imgexistnot)
        pstr += TestHelpers().check_imval(imgval, epsilon=epsilon)
        pstr += TestHelpers().check_imvalexact(imgvalexact, epsilon=epsilon)
        pstr += TestHelpers().check_immask(imgmask)
        pstr += TestHelpers().check_tabcache(tabcache)
        pstr += TestHelpers().check_stopcode(stopcode, ret)
        pstr += TestHelpers().check_reffreq(reffreq)
        return pstr

    def check_final(self, pstr=""):

        if not isinstance(pstr, six.string_types):
            return False
        casalog.post(pstr, 'INFO')
        if pstr.count("Fail") > 0:
            return False
        return True

