########################################################################3
#  task_immath.py
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
# <summary>
# CASA task for smoothing an image, by doing Forier-based convolution
# on a CASA image file.
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed
#
# <author>
# Shannon Jaeger, University of Calgary  (image math) 
# Takeshi Nakazato, National Radio Astronomy Obaservatory (polarization)
# </author>
#
# <prerequisite>
# </prerequisite>
#
# <etymology>
# immath stands for image mathematics
# </etymology>
#
# <synopsis>
#    This task evaluates mathematical expressions involving existing
#    image files. The results of the calculations are stored in the 
#    designated output file.  Options are available to specify mathematical 
#    expression directly or pre-defined expression for calculation of 
#    spectral index image, and polarization intensity and position angle 
#    images are available. The image file names imbedded in the expression or
#    specified in the imagename parameter for the pre-defined calculations may
#    be CASA images or FITS images.
#
#
#    NOTE: Index values for axes start at 0 for the box and chans
#          parameters, but at 1 when used with the indexin function
#          in expression. Use the imhead task to see the range of
#          values for each axes.
#    
#
#    Keyword arguments:
#    outfile -- The file where the results of the image calculations 
#                are stored.  Overwriting an existing outfile is not permitted.
#            Default: none;  Example: outfile='results.im'
#    mode -- mode for mathematical operation
#            Default: evalexpr
#            Options: 'evalexpr' : evalulate a mathematical expression defined in 'expr' 
#                     'spix' : spectalindex image 
#                     'pola' : polarization position angle image 
#                     'poli' : polarization intesity image 
#           mode expandable parameters
#            expr -- (for mode='evalexpr') A mathematical expression, with image file names.
#              Image file names MUST be enclosed in double quotes (&quot;)
#              Default: none 
#              Examples:
#                 Make an image that is image1.im - image2.im
#                   expr=' (&quot;image1.im&quot; - &quot;image2.im&quot; )'
#                 Clip an image below a value (0.5 in this case)
#                   expr = ' iif(&quot;image1.im&quot;>=0.5, &quot;image1.im&quot;, 0.0) '
#                         Note: iif (a, b, c)   a is the boolian expression
#                                               b is the value if true
#                                               c is the value if false
#                 Take the rms value of two images
#                   expr = ' sqrt(&quot;image1.im&quot; * &quot;image1.im&quot; + &quot;image2.im&quot; * &quot;image2.im&quot;) '
#                         Note: No exponentiaion available?
#                 Build an image pixel by pixel from the minimum of (image2.im, 2*image1.im)
#                   expr='min(&quot;image2.im&quot;,2*max(&quot;image1.im&quot;))'
#            imagename -- (for mode='spix','pola','poli') input image names        
#              Default: none;
#              Examples: mode='spix'; imagename=['image1.im','image2.im'] will calculate 
#                       an image of log(S1/S2)/log(f1/f2), where S1 and S2 are fluxes and 
#                       f1 and f2 are frequencies
#                       mode='pola'; imagename=['imageQ.im','imageU.im'] will calculate 
#                       an image of polarization angle distribution, where imageQ.im and 
#                       imageU.im are Stokes Q and U images, respectively. Calculate 0.5*arctan(U/Q).
#                       mode='poli'; imagename=['imageQ.im','imageU.im','imageV.im'] will calculate
#                       total polarization intensity image, where imageQ.im, imageU.im, imageV.im
#                       are Stokes Q, U, and V images, respectively.
#            sigma - (for mode='poli') standard deviation of noise of Stokes images with unit such as
#                    Jy/beam to correct for bias 
#              Default: '0.0Jy/beam' (= no debiasing)
#    mask -- Name of mask applied to each image in the calculation
#            Default '' means no mask;  Example: mask='orion.mask'.  
#    region -- File path to an ImageRegion file.
#            An ImageRegion file can be created with the CASA
#            viewer's region manager.  Typically ImageRegion files
#            will have the suffix '.rgn'.  If a region file is given
#            then the box, chans, and stokes selections whill be 
#            ignored.
#            Default: none
#            Example: region='myimage.im.rgn'
#    box --  A box region on the directional plane
#            Only pixel values acceptable at this time.
#            Default: none (whole 2-D plane);  Example: box='10,10,50,50'
#    chans -- channel numbers, velocity, and/or frequency
#            Only channel numbers acceptable at this time.
#            Default: none (all);  Example: chans='3~20'   
#    stokes -- Stokes parameters to image, may or may not be separated
#            by commas but best if you use commas.
#            Default: none (all); Example: stokes='IQUV';
#            Options: 'I','Q','U','V','RR','RL','LR','LL','XX','YX','XY','YY', ... 
#
#    Available functions in the <i>expr</i> and <i>mask</i> paramters:
#    pi(), e(), sin(), sinh(), asinh(), cos(), cosh(), tan(), tanh(),
#    atan(), exp(), log(), log10(), pow(), sqrt(), complex(), conj()
#    real(), imag(), abs(), arg(), phase(), aplitude(), min(), max()
#    round(), isgn(), floor(), ceil(), rebin(), spectralindex(), pa(), 
#    iif(), indexin(), replace(), ...
#
#    For a full description of the allowed syntax see the 
#    Lattice Expression Language (LEL) documentation on the at:
#    http://aips2.nrao.edu/docs/notes/223/223.html
#
#    NOTE: where indexing and axis numbering are used in the above
#    functions they are 1-based, ie. numbering starts at 1.
#
# </synopsis> 
#
# <example>
# <srcblock>
# </srcblock
#
# </example>
#
# <motivation>
# To provide a user-friendly task interface to imagecalc and ???
# as well as an more user-friendling math syntax then what is
# provided by the CASA Lattice Exprssion Language.
# </motivation>
#
# <todo>
#  Crystal wanted different masks for different inputs
#  but unlikely that its really needed.
#
#  Add an "overwrite" output file parameter
#
#  Add polygon and circle region selection 
# </todo>
########################################################################3

from __future__ import absolute_import
import os
import re
import shutil
import sys

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from casatools import image, imagepol, regionmanager, coordsys, quanta
    from casatasks import casalog
    from .ialib import write_image_history
else:
    from taskinit import *
    from ialib import write_image_history
    image = iatool
    imagepol = potool
    quanta = qatool
    regionmanager = rgtool
  
def immath(
    imagename, mode, outfile, expr, varnames, sigma,
    polithresh, mask, region, box, chans, stokes, stretch,
    imagemd, prec
):
    # Tell CASA who will be reporting
    casalog.origin('immath')
    tmpFilePrefix='_immath_tmp' + str(os.getpid()) + '_'
    try:
        myia = image()
        myia.dohistory(False)
        outia = None
        _immath_initial_cleanup(tmpFilePrefix)
        outfile = _immath_check_outfile(outfile)
        # Find the list of filenames in the expression
        # also do a quick check to see if all of the files
        # exist
        tmpfilenames = ''
        filenames = imagename
        if mode=='evalexpr':
            tmpfilenames = _immath_parse(expr)
        if isinstance(filenames, str):
            filenames = [filenames]
        varnames = _immath_varnames(varnames, filenames, tmpfilenames)
        filenames = _immath_filenames(filenames, tmpfilenames, varnames, mode)
        expr = expr.replace(' ', '')
        if mode == 'spix':
            expr = _immath_dospix(len(filenames), varnames)
        if mode == 'pola':
            _immath_new_pola(
                filenames, outfile, tmpFilePrefix, mask, region,
                box, chans, stokes, stretch, polithresh, myia
            )
            return True
        elif mode == 'poli' or mode == 'lpoli' or mode == 'tpoli':
            _immath_new_poli(
                filenames, outfile, tmpFilePrefix, mask, region,
                box, chans, stokes, stretch, sigma, myia, mode
            )
            return True
        elif mode == 'evalexpr' or mode == 'spix':
            if box or chans or stokes or region or mask:
                (subImages, file_map) = _immath_createsubimages(
                    box, chans, stokes, region, mask,
                    stretch, filenames, myia, tmpFilePrefix
                )
                if imagemd:
                    casalog.post(
                        "Specifying region, box, chan, or stokes will "
                        + "create smaller sub-images. The image "
                        + "metadata specified in imagemd will have to "
                        + "conform to the output, not the input image "
                        + "dimensions. Please check your output image "
                        + "for accurate header definition.", 'WARN'
                    )
                (expr, varnames, subImages) = _immath_updateexpr(
                    expr, varnames, subImages, filenames, file_map
                )
                outia = _immath_compute(
                    imagename, expr, outfile, imagemd, myia, prec
                )
            else:
                # If the user didn't give any region or mask information
                # then just evaluated the expression with the filenames in it.
                outia = _immath_dofull(
                    imagename, imagemd, outfile, mode, expr,
                    varnames, filenames, myia, prec
                )
        else:
            raise(Exception, "Unsupported mode " + str(mode))
        try:
            vars = locals()
            param_names = immath.__code__.co_varnames[:immath.__code__.co_argcount]
            param_vals = [vars[p] for p in param_names]
            write_image_history(
                outia, sys._getframe().f_code.co_name,
                param_names, param_vals, casalog
            )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')
        return True
    except Exception as error:
        if mode == 'evalexpr':
            casalog.post("Unable to process expression " + expr, 'SEVERE')
        else:
            casalog.post("Error running immath", 'SEVERE')
        casalog.post("Exception caught was: " + str(error), 'SEVERE')
        raise
    finally:
        if myia:
            myia.done()
        if outia:
            outia.done()
        _immath_cleanup(tmpFilePrefix)

def _immath_concat_stokes(filenames, target, myia):
    myia.open(filenames[0])
    stokes_axis = myia.coordsys().findaxisbyname("stokes")
    myia.done()
    casalog.post("Concatenating images along stokes axis")
    myia = myia.imageconcat(
        outfile=target, infiles=filenames, axis=stokes_axis
    )
    myia.done()
    
def _immath_getregion(region, box, chans, stokes, mode, myia, target):
    myreg = region
    if (type(region) != type({})):
        if stokes:
            casalog.post(
                "Ignoring stokes parameters selection for mode='"
                + mode + "'."
                ,'WARN' 
            )
            stokes=''
        myia.open(target)
        myrg = regionmanager()
        myreg = myrg.frombcs(
            csys=myia.coordsys().torecord(), shape=myia.shape(), box=box,
            chans=chans, stokes=stokes, stokescontrol="a", region=region
        )
        myia.done()
        myrg.done()
    return myreg

def _immath_new_pola(
    filenames, outfile, tmpFilePrefix, mask, region,
    box, chans, stokes, stretch, polithresh, myia
):
    target = filenames[0]
    if len(filenames) > 1:
        target = tmpFilePrefix + "_concat_for_pola"
        _immath_concat_stokes(filenames, target, myia)
    myreg = _immath_getregion(region, box, chans, stokes, "pola", myia, target)
    mypo = imagepol()
    myqa = quanta()
    if (polithresh):
        if (mask != ""):
            mask = ""
            casalog.post("Ignoring mask parameter in favor of polithresh parameter", 'WARN')
        if (myqa.getunit(polithresh) != ""):
            initUnit = myqa.getunit(polithresh)
            myia.dohistory(False)
            myia.open(filenames[0])
            bunit = myia.brightnessunit()
            polithresh = myqa.convert(polithresh, bunit)
            myia.done()
            if (myqa.getunit(polithresh) != bunit):
                raise Exception("Units of polithresh " + initUnit \
                + " do not conform to input image units of " + bunit \
                + " so cannot perform thresholding. Please correct units and try again.")
            polithresh = myqa.getvalue(polithresh)[0]
            lpol = tmpFilePrefix + "_lpol"
            mypo.open(target)
            myia = mypo.linpolint(debias=False, outfile=lpol, region=myreg)
            myia.done()
            mypo.done()
    mypo.open(target)
    myia = mypo.linpolposang(
        outfile=outfile, region=myreg, mask=mask, stretch=stretch
    )
    mypo.done()
    if (polithresh):
        myexpr = "'" + lpol + "' >= " + str(polithresh)
        myia.dohistory(False)
        myia.calcmask(name='mask0', mask=myexpr)
        casalog.post(
            'Calculated mask based on linear polarization threshold '
            + str(polithresh),
            'INFO'
        )
    myia.done()

def _immath_new_poli(
    filenames, outfile, tmpFilePrefix, mask, region,
    box, chans, stokes, stretch, sigma, myia, mode
):
    target = filenames[0]
    if len(filenames) > 1:
        target = tmpFilePrefix + "_concat_for_poli"
        _immath_concat_stokes(filenames, target, myia)
    if mode == 'tpoli':
        myia.open(target)
        csys = myia.coordsys()
        myia.done()
        stokes = csys.stokes()
        csys.done()
        for p in ["Q", "U", "V"]:
            if stokes.count(p) == 0:
                raise Exception(
                    "Stokes " + p + " is required for mode tpoli but was not found"
                )       
    debias = False
    newsigma = 0
    myqa = quanta()
    if sigma:
        qsigma = myqa.quantity(sigma)
        if myqa.getvalue(qsigma)[0] > 0:
            debias = True
            sigmaunit = myqa.getunit(qsigma)
            try:
                myia.open(filenames[0])
                iunit = myia.brightnessunit()
                myia.done()
            except:
                raise Exception('Unable to get brightness unit from image file ' + filenames[0])
            if sigmaunit != iunit:
                newsigma = myqa.convert(qsigma,iunit)
            else:
                newsigma = sigma
    myreg = _immath_getregion(region, box, chans, stokes, "poli", myia, target)
    mypo = imagepol()
    mypo.open(target)
    if is_CASA6:
        # for some annoying reason, qa.getvalue() returns an array in this context
        numeric_sigma = myqa.getvalue(myqa.quantity(newsigma))[0]
    else:
        numeric_sigma = myqa.getvalue(myqa.quantity(newsigma))
    if mode == 'tpoli' or mode == 'poli':
        myia = mypo.totpolint(
            debias=debias, sigma=numeric_sigma, outfile=outfile,
            region=myreg, mask=mask, stretch=stretch
        )
    elif mode == 'lpoli':
        myia = mypo.linpolint(
            debias=debias, sigma=numeric_sigma, outfile=outfile,
            region=myreg, mask=mask, stretch=stretch
        )
    else:
        raise Exception("Logic Error: Unhandled mode " + mode)
    myia.done()
    mypo.done()
    
def _immath_compute(
    imagename, expr, outfile, imagemd, myia, prec
):
    # Do the calculation
    res = myia.imagecalc(
        pixels=expr, outfile=outfile,
        imagemd=_immath_translate_imagemd(imagename, imagemd), prec=prec
    )
    res.dohistory(False)
    # modify stokes type for polarization intensity image
    return res

def _immath_updateexpr(expr, varnames, subImages, filenames, file_map):
    # Make sure no problems happened
    if len(filenames) != len(subImages):
        raise Exception(
            'Unable to create subimages for all image names given'
        )
    # because real file names also have to be mapped to a corresponding subimage, CAS-1830
    for k in file_map.keys():
        # we require actual image names to be in quotes when used in the expression
        varnames.extend(["'" + k + "'", '"' + k + '"'])
        subImages.extend(2 * [file_map[k]])
    # Put the subimage names into the expression
    try:
        expr = _immath_expr_from_varnames(expr, varnames, subImages)
    except Exception as e:
        casalog.post(
            "Unable to construct pixel expression aborting immath: " + str(e),
            'SEVERE'
        )
        raise
    return (expr, varnames, subImages)
    
def _immath_createsubimages(
    box, chans, stokes, region, mask,
    stretch, filenames, myia, tmpFilePrefix
):
    subImages = []
    file_map = {}
    i = 0
    for image in filenames:
        try:
            myia.open(image)
            myrg = regionmanager()
            reg = myrg.frombcs(csys=myia.coordsys().torecord(),
                shape=myia.shape(), box=box, chans=chans, stokes=stokes,
                stokescontrol="a", region=region
            )
            tmpFile = tmpFilePrefix + str(i)
            subim = myia.subimage(
                region=reg, mask=mask, outfile=tmpFile, stretch=stretch
            )
            subim.done()
            file_map[image] = tmpFile
            subImages.append(tmpFile)
            myia.done()
            i = i + 1
        except Exception as e:
            raise Exception(
                'Unable to apply region to file: ' + image
            )
        finally:
            myia.done()
    return (subImages, file_map)

def _immath_dofull(
    imagename, imagemd, outfile, mode, expr,
    varnames, filenames, myia, prec
):
    expr = _immath_expr_from_varnames(expr, varnames, filenames)    
    return _immath_compute(
        imagename, expr, outfile, imagemd, myia, prec
    )

def _immath_dospix(nfiles, varnames):
    if nfiles != 2:
        raise Exception('Requires two images at different frequencies')
    return 'spectralindex(' + varnames[0] + ', ' + varnames[1] + ')'

def _immath_filenames(filenames, tmpfilenames, varnames, mode):
    ignoreimagename = False
    if mode=='evalexpr':
        varnamesSet = set(varnames)
        count = 0
        for imname in tmpfilenames:
            # check if it is one of varnames, if not check the files in expr exist 
            if(not varnamesSet.issuperset(imname)):
               if( not os.path.exists(imname)):
                   raise Exception('Image data set not found - please verify ' + imname)
               else:
                   count = count + 1            
        if len(tmpfilenames) == count:
            ignoreimagename = True
            filenames = tmpfilenames
    if not ignoreimagename:
        for i in range(len(filenames)):
            if not os.path.exists(filenames[i]):
                casalog.post("Image data set not found - please verify " +filenames[i], "SEVERE")
                raise Exception('Image data set not found - please verify '+filenames[i])
    return filenames

def _immath_varnames(varnames, filenames, tmpfilenames):
    # Construct the variable name list.  We append to the list the
    # default variable names if the user hasn't supplied a full suite.
    if not isinstance(varnames, list):
        name0 = varnames
        varnames = []
        if name0:
            varnames.append(name0)
    nfile = max(len(filenames),len(tmpfilenames))
    for i in range(len(varnames), nfile):
        varnames.append('IM'+str(i))
    casalog.post( 'Variable name list is: '+str(varnames), 'DEBUG1' )
    return varnames

def _immath_initial_cleanup(tmpFilePrefix):
    try:
        _immath_cleanup(tmpFilePrefix)
    except Exception as e:
        casalog.post( 'Unable to cleanup working directory '+os.getcwd()+'\n'+str(e), 'SEVERE' )
        raise

def _immath_check_outfile(outfile):
    if not outfile:
        outfile = 'immath_results.im'
        casalog.post(
            "The outfile parameter is empty, consequently the "
            + "resultant image will be saved on disk in an image named "
            + outfile, 'WARN'
        )
    if (os.path.exists(outfile)):
        raise Exception(
            'Output file '+ outfile
            + ' exists. immath can not proceed, please '
            + 'remove it or change the output file name.'
        )
    return outfile

def _immath_cleanup(filePrefix):
    # Remove any old tmp files that may be left from
    # a previous run of immath
    fileList=os.listdir( os.getcwd() )
    for file in fileList:
        if ( file.startswith( filePrefix ) ):
            shutil.rmtree( file )

def _immath_parse( expr='' ):
        retValue=[]
        
        # Find out if the names are surrounded by single or double quotes
        quote=''
        if ( expr.find('"') > -1 ):
            quote='"'
        if ( expr.find("'") > -1 ):
            quote="'"

        current=expr;
        while( current.count(quote) > 1 ):
            start = current.index(quote)+1
            end   = current[start:].index(quote)+start
            if ( retValue.count( current[start:end] ) > 0 ) :
                # We already have this file name so we won't add it
                # to the list again.  This saves us work and disk space.
                current=current[end+1:]
                continue;
            
            retValue.append( current[start:end] )
            current=current[end+1:]

        return retValue

# it is important to sort the varnames in reverse order before doing
# the substitutions to assure the substitution set is performed correctly
# CAS-1678
def _immath_expr_from_varnames(expr, varnames, filenames):
        tmpfiles = {}
        for i in range(len(filenames)):
                tmpfiles[varnames[i]] = filenames[i]
        # python 3 requires explicit list conversion
        tmpvars = list(tmpfiles.keys())

        tmpvars.sort()
        tmpvars.reverse()
        for varname in tmpvars:
                expr = expr.replace(varname, '"' + tmpfiles[varname] + '"')
        return(expr)

def _immath_translate_imagemd(imagename, imagemd):
    # may IM0 etc is a real file name
    if os.path.exists(imagemd):
        return imagemd
    # match IM0, IM1, ... etc
    m = re.match("^IM(\d+)$", imagemd)
    if not m:
        return imagemd
    idx = int(m.groups()[0])
    if idx == 0 and type(imagename) == str:
        return imagename
    # if here, then imagename is an array of strings
    if idx >= len(imagename):
        # out of range
        return imagemd
    return imagename[idx]

