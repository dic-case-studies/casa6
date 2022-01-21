
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

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:

    # CASA 6
    logging.debug("Importing CASAtools")
    import casatools
    logging.debug("Importing CASAtasks")
    import casatasks
    _cb = casatools.calibrater()
    _tb = casatools.table()
    _tbt = casatools.table()
    _ia  = casatools.image()
    _cb = casatools.calibrater()
    from casatasks import casalog

    casampi_imported = False
    import importlib
    _casampi_spec = importlib.util.find_spec('casampi')
    if _casampi_spec:
        # don't catch import error from casampi if it is found in the system modules
        from casampi.MPIEnvironment import MPIEnvironment
        casampi_imported = True
    else:
        casalog.post('casampi not available - not testing MPIEnvironment stuff', 'WARN')

    def tclean_param_names():
        from casatasks.tclean import _tclean_t
        return _tclean_t.__code__.co_varnames[:_tclean_t.__code__.co_argcount]
    def decon_param_names():
        from casatasks.deconvolve import _deconvolve_t
        return _deconvolve_t.__code__.co_varnames[:_deconvolve_t.__code__.co_argcount]
    def sdint_param_names():
        from casatasks.sdintimaging import _sdintimaging_t
        return _sdintimaging_t.__code__.co_varnames[:_sdintimaging_t.__code__.co_argcount]

    casa6 = True

else:

    # CASA 5
    logging.debug("Import casa6 errors. Trying CASA5...")
    from __main__ import default
    from taskinit import tbtool, mstool, iatool, cbtool
    from taskinit import *
    from casa_stack_manip import stack_find, find_casa
    from mpi4casa.MPIEnvironment import MPIEnvironment
    casampi_imported = True

    _tb = tbtool()
    _tbt = tbtool()
    _ia = iatool()
    _cb = cbtool()
    casa = find_casa()
    if casa.has_key('state') and casa['state'].has_key('init_version') and casa['state']['init_version'] > 0:
        casaglobals=True
        casac = stack_find("casac")
        casalog = stack_find("casalog")

    def tclean_param_names():
        # alternatively could use from tasks import tclean; tclean.parameters
        from task_tclean import tclean
        return tclean.__code__.co_varnames[:tclean.__code__.co_argcount]
    def decon_param_names():
        from tasks import deconvolve
        return deconvolve.parameters.keys()
    def sdint_param_names():
        from tasks import sdintimaging
        return sdintimaging.parameters.keys()

    casa5 = True

############################################################################################
##################################       imagerhelpers       ###############################
############################################################################################
class TestHelpers:

    # For comparison with keywords added by tclean in its output images
    if casampi_imported:
        num_mpi_procs = 1 + len(MPIEnvironment.mpi_server_rank_list())
    else:
        num_mpi_procs = 1

    def delmodels(self,msname="",modcol='nochange'):
       TestHelpers().delmodkeywords(msname) ## Get rid of extra OTF model keywords that sometimes persist...
       if modcol=='delete':
           TestHelpers().delmodelcol(msname) ## Delete model column
       if modcol=='reset0':
           TestHelpers().resetmodelcol(msname,0.0)  ## Set model column to zero
       if modcol=='reset1':
           TestHelpers().resetmodelcol(msname,1.0)  ## Set model column to one

    def delmodkeywords(self,msname=""):
        #delmod(msname)
        _tb.open( msname+'/SOURCE', nomodify=False )
        keys = _tb.getkeywords()
        for key in keys:
            _tb.removekeyword( key )
        _tb.close()

    def resetmodelcol(self,msname="",val=0.0):
        _tb.open( msname, nomodify=False )
        hasmodcol = (  (_tb.colnames()).count('MODEL_DATA')>0 )
        if not hasmodcol:
            _cb.open(msname)
            _cb.close()
        hasmodcol = (  (_tb.colnames()).count('MODEL_DATA')>0 )
        if hasmodcol:
            dat = _tb.getcol('MODEL_DATA')
            dat.fill( complex(val,0.0) )
            _tb.putcol('MODEL_DATA', dat)
        _tb.close();



    def delmodelcol(self,msname=""):
        _tb.open( msname, nomodify=False )
        hasmodcol = (  (_tb.colnames()).count('MODEL_DATA')>0 )
        if hasmodcol:
            _tb.removecols('MODEL_DATA')
        _tb.close()

    def get_coordsys(self,imname):
        try:
            _ia.open(imname)
            csys = _ia.coordsys()
            csys_rec = csys.torecord()
            csys.done()
        finally:
            _ia.close()

        return csys_rec
        
    def check_spec_frame(self,imname,frame, crval=0.0, cdelt=0.0):
        testname = "check_spec_frame"
        pstr = ""
        if os.path.exists(imname):
           res = True
           expcrval=""
           expcdelt=""
           thecval=""
           thecdelt=""
           coordsys = TestHelpers().get_coordsys(imname)
           baseframe = coordsys['spectral2']['system']
           basecrval = coordsys['spectral2']['wcs']['crval']
           basecdelt = coordsys['spectral2']['wcs']['cdelt']
           if baseframe != frame:
                res = False
           else:
                res = True
                if crval!=0.0:
                     if abs(basecrval - crval)/abs(crval) > 1.0e-6:
                          res = False
                     thecrval = " with crval " + str(basecrval)
                     expcrval = " with expected crval " + str(crval)
                else:
                     # skip the crval test
                     thecrval = ""
                     expcrval = ""
                if cdelt!=0.0:
                     if abs(basecdelt - cdelt)/abs(cdelt) > 1.0e-6:
                          res = False
                     thecdelt = " with cdelt " + str(basecdelt)
                     expcdelt = " with expected cdelt " + str(cdelt)
                else:
                     # skip the crval test
                     thecdelt = ""
           thecorrectans = frame +  expcrval + expcdelt
           pstr =  "[" + testname + "] " + imname + ": Spec frame is " +\
           str(baseframe) + thecrval + thecdelt + " (" +\
           TestHelpers().verdict(res) +" : should be " + thecorrectans +" )"
           print(pstr)
           pstr=pstr+"\n"
           
        #self.checkfinal(pstr)
        return pstr

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

    def exists(self, imname):
        # AW: This is a duplicate function, but it helps maintain uniformity to test_tclean.py
        """ Exists """
        return os.path.exists(imname)

    def get_peak_res(self, summ):
        """Get Peak Res"""
        # AW:  This can be reduced down for readability but putting in a fix for CAS-13182
        peakres = None
        if is_CASA6:
            if 'summaryminor' in summ:
                reslist = summ['summaryminor'][1,:]
                peakres = reslist[ len(reslist)-1 ]

        else:
            if summ.has_key('summaryminor'):
                reslist = summ['summaryminor'][1,:]
                peakres = reslist[ len(reslist)-1 ]
                
        return peakres

    def check_peak_res(self, summ,correctres, epsilon=0.05):
        """Check Peak Res"""
        peakres = TestHelpers().get_peak_res(summ)
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
        # AW: This can be reduced down for readability but putting in a fix for CAS-13182
        modflux = None
        if is_CASA6:
            if 'summaryminor' in summ:
                modlist = summ['summaryminor'][2,:]
                modflux = modlist[ len(modlist)-1 ]
        else:
            if summ.has_key('summaryminor'):
                modlist = summ['summaryminor'][2,:]
                modflux = modlist[ len(modlist)-1 ]
        return modflux

    def check_mod_flux(self, summ,correctmod, epsilon=0.05):
        """Check Mod Flux"""
        modflux = TestHelpers().get_mod_flux(summ)
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

    def check_chanvals(self,msname,vallist, epsilon = 0.05): # list of tuples of (channumber, relation, value) e.g. (10,">",1.0)
        testname = "check_chanvals"
        pstr = ""
        for val in vallist:
            if len(val)==3:
                thisval = TestHelpers().check_modelchan(msname,val[0])
                if val[1]==">":
                    ok = thisval > val[2]
                elif val[1]=="==":     
                    ok = abs( (thisval - val[2])/val[2] ) < epsilon
                elif val[1]=="<":     
                    ok = thisval < val[2]
                else:
                    ok=False
        pstr = pstr + "[" + testname + "] Chan "+ str(val[0]) + "  is " + str(thisval) + " ("+self.verdict(ok)+" : should be " + str(val[1]) + str(val[2]) + ")\n"

        print(pstr)
        return pstr

    def check_modelchan(self,msname="",chan=0):
        _tb.open( msname )
        hasmodcol = (  (_tb.colnames()).count('MODEL_DATA')>0 )
        modsum=0.0
        if hasmodcol:
            dat = _tb.getcol('MODEL_DATA')[:,chan,:]
            modsum=dat.mean()
        _tb.close()
        ##print(modsum)
        return modsum

    def get_iter_done(self, summ):
        """Get Iterdone"""
        # AW:  This can be reduced down for readability but putting in a fix for CAS-13182
        iters = None
        if is_CASA6:
            if 'iterdone' in summ:
                iters = summ['iterdone']
        else:
            if summ.has_key('iterdone'):
                iters = summ['iterdone']
        return iters

    def delmodkeywords(self,msname=""):
        """get rid of extra model keywords that sometimes persist"""

        _tb.open( msname+'/SOURCE', nomodify=False )
        keys = _tb.getkeywords()
        for key in keys:
            _tb.removekeyword( key )
        _tb.close()

    def resetmodelcol(self,msname="",val=0.0):
        """ set model column to val"""
        _tb.open( msname, nomodify=False )
        hasmodcol = (  (_tb.colnames()).count('MODEL_DATA')>0 )
        if not hasmodcol:
            _cb.open(msname)
            _cb.close()
        hasmodcol = (  (_tb.colnames()).count('MODEL_DATA')>0 )
        if hasmodcol:
            dat = _tb.getcol('MODEL_DATA')
            dat.fill( complex(val,0.0) )
            _tb.putcol('MODEL_DATA', dat)
        _tb.close();

    def delmodels(self,msname="",modcol='nochange'):
        """Get rid of OTF model and model column"""
         
        self.delmodkeywords(msname) ## Get rid of extra OTF model keywords that sometimes persist...

        if modcol=='delete':
            self.delmodelcol(msname) ## Delete model column
        if modcol=='reset0':
            self.resetmodelcol(msname,0.0)  ## Set model column to zero
        if modcol=='reset1':
            self.resetmodelcol(msname,1.0)  ## Set model column to one

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

    def check_list_vals(self, list1, list2, test, epsilon=None):
        """ compares 2 lists and returns if they are equivalent (within error) 
        """
        report = ''
        if len(list1) == len(list2):
            exact = (epsilon is None)
            i = 0
            while i < len(list1):
                result, pstr = self.check_val(list1[i], list2[i], \
                    valname=test+' index '+str(i), exact=exact, epsilon=epsilon)
                if result == False:
                    report += pstr
                    #report = pstr
                    #break
                i += 1
        else:
            result = False

        return result, report

    def check_dict_vals(self, exp_dict, act_dict, suffix, epsilon=0.01):
        """ Compares expected dictionary with actual dictionary.

            Parameters
            ----------
            exp_dict: dictionary
                Expected values, as key:value pairs.
                Values must be lists, where
                val[0] is "True" for an exact match, or a float for a non-exact
                epsilon for check_val(). Default is input epsilon.
                val[1] is either a value to check_val() against or a list to
                check_list_vals() against.
            act_dict: dictionary
                Actual values to compare to exp_dict (and just the values).
                Keys must match between exp_dict and act_dict.
            suffix: string
                For use with summary print statements.

            Notes
            -----
            Each exp_dict key:value pair does not need to match in
            exactness or listedness. One value could be exact, and
            the next relative to epsilon. One value could be a list,
            and the next an integer.

            Notes
            -----
            Example exp_dict
            {
                'end': [True, 2203765e5],
                'start': [False, 220257e5], # uses default epsilon value
                'start_delta': [0.01, 2202527e5],
            }
        """
        report = ''
        passed = True
        chans = 0
        for key in exp_dict:

            # convert the expected value in exp_dict[key] to its own dictionary
            exp_val = exp_dict[key]
            v = {
                'val': exp_val[1],
                'exact': False,
                'epsilon': epsilon,
            }
            if type(exp_val[0]) == bool:
                v['exact'] = exp_val[0]
            else:
                v['epsilon'] = exp_val[0]
            if v['val'] == 0.0: # special case for "0"
                v['exact'] = True
            if v['exact'] == True:
                v['epsilon'] = None

            # evaluate the expected value against the actual value
            if type(v['val']) == list:
                result, pstr = self.check_list_vals(act_dict[key], 
                    v['val'], test=suffix+' '+key, epsilon=v['epsilon'])
                report += self.check_val(result, True, \
                    valname=suffix+' '+key, exact=True)[1]
                report += pstr
            else:
                report += self.check_val(act_dict[key], \
                    v['val'], exact=v['exact'], epsilon=v['epsilon'],
                    valname=suffix+' '+key)[1]

        return report

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

    def check_keywords(self, imlist, testname="check_keywords", check_misc=True):
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
                issues = TestHelpers().check_im_keywords(imname, check_misc=check_misc)
                if issues:
                    pstr += '[{0}] {1}: {2}'.format(testname, imname, issues)
        if not pstr:
            pstr += ('All expected keywords in imageinfo, miscinfo, and coords found. '
                     '({})\n'.format(testname, TestHelpers().verdict(False)))
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
                # basic miscinfo and 'TcleanProcessingInfo' as per CAS-12204
                mandatory_miscinfo = ['INSTRUME', 'distance',
                                      'mpiprocs', 'chnchnks', 'memreq', 'memavail']
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

                # ensure mpiprocs is correct. Other keywords added in CAS-12204 have more
                # variable values (memavail, memreq, etc.) and cannot be compared here.
                if entry == 'mpiprocs':
                    if keys[record][entry] != self.num_mpi_procs:
                        pstr += ('mpiprocs is not as expected. It is {} but it should be {}, ({})'.
                                 format(keys[record][entry], self.num_mpi_procs, TestHelpers().verdict(False)))
        return pstr

    def check_forbidden_entries(self, entries, record, keys):
        pstr = ''
        for entry in entries:
            if entry in keys[record]:
                pstr += ('entry {0} should not be in record {1} ({2})\n'.format(entry, record, TestHelpers().verdict(False)))
        return pstr

    def check_history(self, imlist, testname="check_history"):
        """
        Checks presence of the logtable and rows with history information (task name,
        CASA version, all task parameters, etc.).

        :param imlist: names of the images produced by a test execution.
        :param testname: name to use in the checks report string
        :returns: the usual (test_imager_helper) string with success/error messages.
        """
        pstr = ''
        for imname in imlist:
            if os.path.exists(imname):
                issues = TestHelpers().check_im_history(imname)
                if issues:
                    pstr += '[{0}] {1}: {2}'.format(testname, imname, issues)
        if not pstr:
            pstr += ('[{}] All expected history entries found. ({})\n'.
                     format(testname, TestHelpers().verdict(True)))
        return pstr

    def check_im_history(self, imname):
        """
        Check the history records in an image, ensuring the taskname, CASA version, and
        full list of parameters is found (all the same number of times).

        :param imname: image name (output image from tclean)
        :returns: the usual (test_imager_helper) string with success/error messages.
        Errors are marked with the tag '(Fail' as per self.verdict().
        """
        import itertools

        ia_open = False
        try:
            _ia.open(imname)
            ia_open = True
            history = _ia.history(list=False)
        except RuntimeError as exc:
            pstr = ('Cannot retrieve history subtable from image: {}. Error: {}'.
                    format(imname, exc))
            return pstr
        finally:
            if ia_open:
                _ia.close()

        # build up a list of parameter names (to be evaluated as necessary)
        task_param_names = {
            "tclean": tclean_param_names,
            "deconvolve": decon_param_names,
            "sdintimaging": sdint_param_names
        }

        pstr = ''
        ncallsdict = {}
        ncallsdict['tclean']       = sum(line.startswith('taskname=tclean') for line in history)
        ncallsdict['deconvolve']   = sum(line == 'taskname=deconvolve' for line in history)
        ncallsdict['sdintimaging'] = sum(line == 'taskname=sdintimaging' for line in history)
        ncalls                     = ncallsdict['tclean'] + ncallsdict['deconvolve'] + ncallsdict['sdintimaging']
        nversions = sum(line.startswith('version:') for line in history)
        if ncalls < 1:
            pstr += ('No calls to cleaning tasks were found in history. ({})\n'.
                     format(TestHelpers().verdict(False)))
        if nversions < 1:
            pstr += ('No CASA version was found in history. ({})\n'.
                     format(TestHelpers().verdict(False)))
        # allow for impbcor history which puts one version line in some tests
        if ncalls != nversions and not nversions == ncalls+1:
            pstr += ('The number of taskname entries ({}) and CASA version entries ({}) do '
                     'not match. ({})\n'.format(ncalls, nversions,
                                                TestHelpers().verdict(False)))

        # check parameter names for cleaning tasks
        histories = TestHelpers().split_histories_by_task(history)
        for tname in ['tclean', 'deconvolve', 'sdintimaging']:
            ntcalls = ncallsdict[tname]
            if ntcalls == 0:
                continue

            # get task parameters and history to be checked
            if tname == 'tclean':
                tclean_keys = list(filter(lambda line: line.startswith("taskname=tclean"), histories.keys()))
                hist = list(itertools.chain.from_iterable([histories[k] for k in tclean_keys])) # concat all tclean* histories into one long list of lines
                pnames = tclean_param_names()
            else:
                hist = histories['taskname='+tname]
                if callable(task_param_names[tname]): # lazy evaluation of parameter names, since deconvolve doesn't exist in my branch yet
                    task_param_names[tname] = task_param_names[tname]()
                pnames = task_param_names[tname]

            # check parameters
            for param in pnames:
                nparval = sum('=' in line and line.split('=')[0].strip() == param for
                              line in hist)
                if nparval < 1:
                    pstr += ('No entries for {} parameter {} found in history. ({})'
                             '.'.format(tname, param, TestHelpers().verdict(False)))
                if nparval != ntcalls:
                    pstr += ("The number of history entries for parameter '{}' ({}) and task "
                             "calls ({}) do not match ({}).".
                             format(param, nparval, ntcalls, TestHelpers().verdict(False)))
        return pstr

    def split_histories_by_task(self, history):
        """
        Given the list of strings in history, split it wherever there is a "taskname=*" line.
        Return a dictionary where the keys are the "taskname=*" lines, and the values are all lines
        for that taskname.
        """
        ret = {}
        task_line = ""
        for line in history:
            if line.startswith("taskname="):
                task_line = line
            if task_line not in ret:
                ret[task_line] = []
            ret[task_line].append(line)
        return ret

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
        readval = TestHelpers().get_pixmask(imname, thepos)
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
        _ia.done()
        reffreq = csys.referencevalue()['numeric'][3]
        csys.done()
        if  abs(reffreq - theval)/theval > epsilon:
            retres = False
        else:
            retres = True
        pstr = "[ {} ] Ref-Freq is {} ( {} : should be {} )\n".format(testname, str(reffreq), TestHelpers().verdict(retres), str(theval))
        logging.info(pstr)
        return pstr

    def check_imexist(self, imgexist, check_keywords_misc=True):
        pstr = ''
        if imgexist != None:
            if type(imgexist) == list:
                pstr += TestHelpers().check_ims(imgexist, True)
                #print("pstr after checkims = {}".format(pstr))
                pstr += TestHelpers().check_keywords(imgexist, check_misc=check_keywords_misc)
                #print("pstr after check_keywords = {}".format(pstr))
                pstr += TestHelpers().check_history(imgexist)
                #print("pstr after check_history = {}".format(pstr))
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

    def check_reffreq(self, reffreq, epsilon=0.05):
        pstr = ''
        if reffreq != None:
            if type(reffreq) == list:
                for ii in reffreq:
                    if type(ii) == tuple and len(ii) == 2:
                        pstr += TestHelpers().check_ref_freq(ii[0], ii[1], epsilon=epsilon)
        return pstr

    def check_tfmask(self, tfmask, testname="check_tfmask"):
        pstr = ''
        if tfmask != None:
            if type(tfmask) == list:
                for ii in tfmask:
                    if type(ii) == tuple and len(ii) == 2:
                        _ia.open(ii[0])
                        mname = _ia.maskhandler('get')
                        mreport = "[" + testname + "]  T/F mask name for " + ii[0] +  " is : " + str(mname)
                        if mname==ii[1]:
                            mreport = mreport + " ("+TestHelpers().verdict(True) +" : should be " + str(ii[1]) + ") \n"
                        else:
                            mreport = mreport + " ("+TestHelpers().verdict(False) +" : should be " + str(ii[1]) + ") \n"
                        _ia.close()
                        pstr += mreport
            print(pstr)
        return pstr


    def get_log_length(self):
        return os.path.getsize(casalog.logfile())

    def check_logs(self, start, expected, testname="check_logs"):
        """ Test that there are log lines that match the expected regex strings (one line per expected string). """
        # Example usage:
        # logstart = test_helper.get_log_length()
        # tclean(...)
        # report = test_helper.check_logs(logstart, expected=[ r"-+ Run Minor Cycle Iterations  -+" ])
        import re

        # read all lines from the logfile that are relevant to this test
        lines = []
        with open(casalog.logfile(), 'r') as f:
            f.seek(start)
            lines = f.readlines()
        # casalog.post("Searching through " + str(len(lines)) + " log lines", "SEVERE") # debugging

        # find expected matches
        unmet = []
        for i in range(len(expected)):
            expectation = re.compile(expected[i])
            # casalog.post("expected["+str(i)+"]: "+expected[i], "SEVERE") # debugging
            found = -1
            for j in range(len(lines)):
                if (expectation.search(lines[j]) != None):
                    found = j
                    break
                else:
                    pass
                    # casalog.post("  X: " + lines[j].rstrip(), "SEVERE") # debugging
            if (found == -1):
                unmet.append(expected[i])
                # casalog.post("  expectation not met", "SEVERE") # debugging
            else:
                del lines[found]
                # casalog.post("  expectation met by line: " + lines[found].rstrip(), "SEVERE") # debugging

        # check that all expectations were met
        if (len(unmet) == 0):
            return "[ {} ]: found {} matching log lines (Pass)\n".format(testname, len(expected))
        else:
            return "[ {} ]: found {} out of {} matching log lines (Fail, unmet expectations: {})\n".format(testname, len(expected)-len(unmet), len(expected), ", ".join(unmet))

    def checkall(self, ret=None, peakres=None, modflux=None, iterdone=None, nmajordone=None, imgexist=None, imgexistnot=None, imgval=None, imgvalexact=None, imgmask=None, tabcache=True, stopcode=None, reffreq=None, firstcyclethresh=None, epsilon=0.05, tfmask=None, check_keywords_misc=True):
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
            firstcyclethresh = None, #a float 
            tfmask=None # list of tuples of (imagename, maskname). 
        """
        pstr = "[ checkall ] \n"
        if ret != None and type(ret) == dict:
            try:
                if peakres != None:
                    out, message = TestHelpers().check_val(val=TestHelpers().get_peak_res(ret), correctval=peakres, valname="peak res", epsilon=epsilon)
                    pstr = pstr + message
                if modflux != None:
                    out, message = TestHelpers().check_val(val=TestHelpers().get_mod_flux(ret), correctval=modflux, valname="mod flux", epsilon=epsilon)
                    pstr = pstr + message
                if iterdone != None:
                    out, message = TestHelpers().check_val(val=ret['iterdone'], correctval=iterdone, valname="iterdone", exact=True)
                    pstr = pstr + message
                if nmajordone != None:
                    out, message = TestHelpers().check_val(val=ret['nmajordone'], correctval=nmajordone, valname="nmajordone", exact=True)
                    pstr = pstr + message
                if firstcyclethresh != None:
                    out, message = TestHelpers().check_val(val=ret['summaryminor'][3][0], correctval=firstcyclethresh, valname='initial cyclethreshold', epsilon=epsilon)
                    pstr = pstr + message
            except Exception as e:
                logging.info(ret)
                raise
        logging.info("Epsilon: {}".format(epsilon))
        pstr += TestHelpers().check_imexist(imgexist, check_keywords_misc=check_keywords_misc)
        pstr += TestHelpers().check_imexistnot(imgexistnot)
        pstr += TestHelpers().check_imval(imgval, epsilon=epsilon)
        pstr += TestHelpers().check_imvalexact(imgvalexact, epsilon=epsilon)
        pstr += TestHelpers().check_immask(imgmask)
        pstr += TestHelpers().check_tabcache(tabcache)
        pstr += TestHelpers().check_stopcode(stopcode, ret)
        pstr += TestHelpers().check_reffreq(reffreq, epsilon=epsilon)
        pstr += TestHelpers().check_tfmask(tfmask)
        return pstr

    def check_final(self, pstr=""):

        import re
        casalog.post(pstr, 'INFO')
        if len(re.findall(r"\(.?Fail",pstr)) > 0:
            return False
        return True
        
    def write_file(self,filename,str_text):
        """Save the string in a text file"""
        inp = filename
        cmd = str_text
        # remove file first
        if os.path.exists(inp):
            os.remove(inp)
        # save to a file
        with open(inp, 'w') as f:
            f.write(cmd)
        f.close()
        return

    def mergeParaCubeResults(self,
                     ret=None,
                     parlist=[]
                     #peakres=None, # a float
                     #modflux=None, # a float
                     #iterdone=None, # an int
                     #nmajordone=None, # an int
                     #imexist=None,  # list of image names
                     #imexistnot=None, # list of image names
                     #imval=None,  # list of tuples of (imagename,val,pos)
                     #imvalexact=None, # list of tuples of (imagename,val,pos)
                     #immask=None,  #list of tuples to check mask value
                     #tabcache=True,
                     #stopcode=None,
                     #reffreq=None # list of tuples of (imagename, reffreq)
                     ):
        if ret!=None and isinstance(ret,dict):
            if list(ret.keys())[0].count('node'):
                mergedret={}
                nodenames = list(ret.keys())
                # must be parallel cube results
                if parlist.count('iterdone'):
                    retIterdone = 0
                    for inode in nodenames:
                        #print("ret[",inode,"]=",ret[inode])
                        #print("inode.strip = ", int(inode.strip('node')))
                        retIterdone+=ret[inode][int(inode.strip('node'))]['iterdone']
                    mergedret['iterdone']=retIterdone
                if parlist.count('nmajordone'):
                    retNmajordone = 0
                    for inode in nodenames:
                        retNmajordone = max(ret[inode][int(inode.strip('node'))]['nmajordone'],retNmajordone)
                    mergedret['nmajordone']=retNmajordone
                if parlist.count('peakres'):
                    #retPeakres = 0
                    #for inode in nodenames:
                        #tempreslist = ret[inode][int(inode.strip('node'))]['summaryminor'][1,:]
                        #if len(tempreslist)>0:
                        #    tempresval = tempreslist[len(tempreslist)-1]
                        #else:
                        #    tempresval=0.0
                        #retPeakres = max(tempresval,retPeakres)
                    #mergedret['summaryminor']=ret['node1'][1]['summaryminor']
                    if 'summaryminor' not in mergedret:
                        for inode in nodenames:
                            nodeid = int(inode.strip('node'))
                            if ret[inode][nodeid]['summaryminor'].size!=0:
                                lastnode = inode
                                lastid = nodeid
                           
                        mergedret['summaryminor']=ret[lastnode][lastid]['summaryminor']
                if parlist.count('modflux'):
                    #retModflux = 0
                    #for inode in nodenames:
                    #    tempmodlist = ret[inode][int(inode.strip('node'))]['summaryminor'][2,:]
                    #    print "tempmodlist for ",inode,"=",tempmodlist
                    #    if len(tempmodlist)>0:
                    #         tempmodval=tempmodlist[len(tempmodlist)-1]
                    #    else:
                    #         tempmodval=0.0
                    #    retModflux += tempmodval
                    #mergedret['modflux']=retModflux
                    if 'summaryminor' not in mergedret:
                        for inode in nodenames:
                            nodeid = int(inode.strip('node'))
                            if ret[inode][nodeid]['summaryminor'].size!=0:
                                lastnode = inode
                                lastid = nodeid
                           
                        mergedret['summaryminor']=ret[lastnode][lastid]['summaryminor']
                if parlist.count('stopcode'):
                    mergedret['stopcode']=ret['node1'][1]['stopcode']
            else:
                mergedret=ret

        return mergedret

