from __future__ import absolute_import
import os
import sys
import shutil

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import image, coordsys, quanta
    from casatasks import casalog
    from .ialib import write_image_history
else:
    from taskinit import *
    from ialib import write_image_history

    image = iatool
    coordsys = cstool
    quanta = qatool
    
def imregrid(
    imagename, template, output, asvelocity, axes, shape,
    interpolation, decimate, replicate, overwrite
):
    _myia = None
    outia = None
    csys = None
    try:
        casalog.origin('imregrid')
        if hasattr(template, 'lower') and not template.lower() == "get":
            # First check to see if the output file exists.  If it
            # does then we abort.  CASA doesn't allow files to be
            # over-written, just a policy.
            if len(output) < 1:
                output = imagename + '.regridded'
                casalog.post("output was not specified - defaulting to\n\t"
                     + output, 'INFO')
        _myia = image()
        _myia.dohistory(False)
        # Figure out what the user wants.
        if not isinstance(template, dict):
            if template.lower() == 'get':
                _myia.open(imagename)
                csys = _myia.coordsys()
                shape = _myia.shape()
                _myia.done()
                return {'csys': csys.torecord(), 'shap': shape}
            elif template.upper() in (
                'J2000', 'B1950', 'B1950_VLA', 'GALACTIC',
                'HADEC', 'AZEL', 'AZELSW', 'AZELNE', 'ECLIPTIC',
                'MECLIPTIC', 'TECLIPTIC', 'SUPERGAL'
            ):
                outia = _imregrid_to_new_ref_frame(
                    _myia, imagename, template, output, axes,
                    shape, overwrite, interpolation, decimate
                )
                try:
                    param_names = imregrid.__code__.co_varnames[:imregrid.__code__.co_argcount]
                    if is_python3:
                        vars = locals( )
                        param_vals = [vars[p] for p in param_names]
                    else:
                        param_vals = [eval(p) for p in param_names]   
                    write_image_history(
                        outia, sys._getframe().f_code.co_name,
                        param_names, param_vals, casalog
                    )
                except Exception as instance:
                    casalog.post(
                        "*** Error \'%s\' updating HISTORY" % (instance), 'WARN'
                    )
                outia.done()
                return
            else:
                if (
                    not os.path.isdir(template)
                    or not os.access(template, os.R_OK)
                ):
                    raise TypeError('Cannot read template image %s' % template)
                template_ia = image()
                template_ia.open(template)
                template_csys = template_ia.coordsys()
                image_ia = image()
                image_ia.open(imagename)
                image_csys = image_ia.coordsys()
                tempshape = template_ia.shape()
                imshape = image_ia.shape()
                axestoregrid = axes
                if (axes[0] < 0):
                    # default value of axes, need to determine actual axes to
                    # send to ia.regrid()
                    axestoregrid = []
                    image_ncoords = image_csys.ncoordinates()
                    for i in range(image_ncoords):
                        ctype = image_csys.coordinatetype(i)[0]
                        template_coord = template_csys.findcoordinate(ctype)
                        if ctype != 'Stokes' and template_coord["return"]:
                            # only regrid if not Stokes axis and coordinate
                            # exists in template
                            for template_pix_axis in template_coord['pixel']:
                                if tempshape[template_pix_axis] > 1:
                                    # only regrid if template axis is not
                                    # degenerate
                                    world_axes = image_csys.findcoordinate(
                                        ctype
                                    )['pixel']
                                    for world_pix_axis in world_axes:
                                        if imshape[world_pix_axis] > 1:
                                            # only regrid if the world axis is
                                            # not degenerate
                                            axestoregrid.append(world_pix_axis)
                    # eliminate dups
                    axestoregrid = list(set(axestoregrid))
                    if len(axestoregrid) == 0:
                        raise RuntimeError("Found no axes to regrid!")
                    axestoregrid.sort()
                if (len(shape) == 1 and shape[0] == -1):
                    shape = _imregrid_handle_default_shape(
                        imshape, image_csys, template_csys, 
                        axestoregrid, tempshape, axes
                    )
                template_ia.done()
                image_ia.done()
                csys = template_csys
        else:
            # csys and shape specified in dictionary generated by previous run
            # with template="get"
            csys = coordsys()
            csys.fromrecord(template['csys'])
            shape = template['shap']

        # The actual regridding.
        _myia.open(imagename)
        outia = _myia.regrid(
            outfile=output, shape=shape, csys=csys.torecord(),
            axes=axes, asvelocity=asvelocity,
            method=interpolation, decimate=decimate,
            replicate=replicate, overwrite=overwrite
        )
        try:
            param_names = imregrid.__code__.co_varnames[:imregrid.__code__.co_argcount]
            if is_python3:
                vars = locals( )
                param_vals = [vars[p] for p in param_names]
            else:
                param_vals = [eval(p) for p in param_names]   
            write_image_history(
                outia, sys._getframe().f_code.co_name,
                param_names, param_vals, casalog
            )
        except Exception as instance:
            casalog.post(
                "*** Error \'%s\' updating HISTORY" % (instance), 'WARN'
            )
        return
    finally:
        if _myia:
            _myia.done()
        if outia:
            outia.done()
        if csys:
            csys.done()

def _imregrid_to_new_ref_frame(
    _myia, imagename, template, output, axes,
    shape, overwrite, interpolation, decimate
):
    _myia.open(imagename)
    csys = _myia.coordsys()
    if len(shape) > 0 and shape != [-1]:
         casalog.post(
            "Specified shape parameter will be ignored when regridding to a "
            + "new reference frame", "WARN"
        )
    if len(axes) > 0 and axes != [-1]:
        casalog.post(
            "Specified axes parameter will be ignored when "
            + "regridding to a new reference frame",
            "WARN"
        )
    dirinfo = csys.findcoordinate("direction")
    if not dirinfo['return']:
        raise Exception("Image does not have a direction coordinate.")
    newrefcode = template.upper()
    oldrefcode = csys.referencecode("direction")[0]
    if oldrefcode == newrefcode:
        casalog.post(
            imagename + ' is already in ' + oldrefcode,
            'INFO'
        )
        casalog.post("...making a straight copy...", 'INFO')
        subi = _myia.subimage(output)
        _myia.done()
        csys.done()
        return subi
    if (csys.projection()['type'] == 'SFL'):
        raise Exception(
            "The direction coordinate of this image has a projection "
            "of SFL. Because of the constraints of this projection, "
            "this image cannot be easily rotated. You may wish to "
            "consider temporarily modifying the projection using "
            "cs.setprojection() to allow rotation of the image."
        )
    casalog.post(
        "Changing coordinate system from " + oldrefcode
        + " to " + newrefcode, 'INFO'
    )
    diraxes = dirinfo['pixel']
    if len(diraxes) != 2:
        raise Exception(
            "Unsupported number of direction axes. There must be exactly 2."
        )
    dirrefpix = csys.referencepixel("direction")["numeric"]
    shape = _myia.shape()
    centerpix = [int(shape[diraxes[0]]/2), int(shape[diraxes[1]]/2)]
    if centerpix[0] != dirrefpix[0] or centerpix[1] != dirrefpix[1]:
        casalog.post(
            "Center direction pixel and reference pixel are "
            + "different, making a temporary image and setting "
            + "the reference pixel equal to the center pixel. "
            + "The output image will have this modified coordinate system."
        )
        # so toworld() works in the correct frame
        csys.setconversiontype(oldrefcode)
        newrefpix = csys.referencepixel()["numeric"]
        newrefpix[diraxes[0]] = centerpix[0]
        newrefpix[diraxes[1]] = centerpix[1]
        newrefval = csys.toworld(newrefpix)["numeric"]
        csys.setreferencepixel(newrefpix)
        csys.setreferencevalue(newrefval)
        tsub = _myia.subimage()
        _myia.done()
        _myia = tsub
        _myia.dohistory(False)
        _myia.setcoordsys(csys.torecord())
    doref = (
        csys.referencecode("direction")[0] == csys.conversiontype("direction")
    )
    angle = csys.convertdirection(newrefcode)
    myqa = quanta()
    mysin = myqa.getvalue(myqa.sin(angle))
    mycos = myqa.getvalue(myqa.cos(angle))
    xnew = 0
    ynew = 0
    for xx in [-centerpix[0], centerpix[0]]:
        for yy in [-centerpix[1], centerpix[1]]:
            xnew = max(xnew, abs(xx*mycos - yy*mysin + 1))
            ynew = max(ynew, abs(xx*mysin + yy*mycos + 1))
    pad = int(max(xnew - shape[0]/2, ynew - shape[1]/2))
    if pad > 0:
        casalog.post(
            "Padding image by " + str(pad)
            + " pixels so no pixels are cut off in the regridding",
            "NORMAL"
        )
        _myia = _myia.pad("", pad, wantreturn=True)
        _myia.dohistory(False)
        shape = _myia.shape()
        newrefpix = csys.referencepixel()['numeric']
        newrefpix[diraxes[0]] = newrefpix[diraxes[0]] + pad
        newrefpix[diraxes[1]] = newrefpix[diraxes[1]] + pad
        csys.setreferencepixel(newrefpix)            
    regridded = _myia.regrid(
        outfile="",shape=shape, csys=csys.torecord(), axes=diraxes,
        method=interpolation, decimate=decimate ,doref=doref
    )
    regridded.dohistory(False)
    # beam is rotated counterclockwise
    angle_in_deg = format("%7.3f" % myqa.convert(angle, "deg")['value'])
    casalog.post(
        "Will rotate beams counterclockwise by " + angle_in_deg + " degrees, "
        + "if necessary, to account for angle between original and new frame "
        + "at the reference pixel", 'NORMAL'
    )
    regridded.rotatebeam(angle=myqa.mul(-1, angle))
    # now crop
    casalog.post("Cropping masked image boundaries", "NORMAL")
    cropped = regridded.crop(outfile=output, axes=diraxes, overwrite=overwrite) 
    regridded.done()
    _myia.done()
    return cropped

def _imregrid_handle_default_shape(
    imshape, image_csys, template_csys, 
    axestoregrid, tempshape, original_axes
):
    # CAS-4959, output shape should have template shape for axes being
    # regridded, input image shape for axes not being regridded, CAS-4960 in
    # cases where the input image and template both have multiple stokes, the
    # number of pixels on the output stokes axis is to be the number of stokes
    # the input and template have in common
    shape = imshape
    targetaxesnames = image_csys.names()
    template_spectral_info = template_csys.findcoordinate("Spectral")
    template_has_spectral = template_spectral_info['return']
    if template_has_spectral:
        template_spectral_axis = template_spectral_info['pixel'][0]
    atr = axestoregrid[:]
    for i in range(len(imshape)):
        atr_count = 0
        for j in atr:
            if i == j:
                # axis numbers may not correspond so have to look for the
                # template axis location by the axis name, CAS-4960
                template_axis = template_csys.findaxisbyname(targetaxesnames[i])
                template_axis_length = tempshape[template_axis]
                if (
                    template_has_spectral
                    and template_spectral_axis == template_axis
                    and template_axis_length == 1
                ):
                    casalog.post(
                        "You've specified that you want to regrid the spectral axis without specifying "
                        "the output shape. Normally the length chosen would be that of the corresponding "
                        "template axis, however, the template spectral axis is degenerate and one cannot regrid "
                        "an axis such that its output length is one. So, removing axis " + str(j)
                        + " from the axes list and just copying the input spectral information to the output image",
                        "WARN"
                    )
                    shape[i] = imshape[i]
                    del axestoregrid[atr_count]
                else:
                    shape[i] = tempshape[template_axis]
                 
                break
            atr_count += 1;
    if (
        template_csys.findcoordinate("stokes")['return']
        and image_csys.findcoordinate("stokes")['return']
        and len(template_csys.stokes()) > 1
        and len(image_csys.stokes()) > 1
    ):
        stokes_axis = image_csys.findaxisbyname("stokes")
        found = (
            len(original_axes) == 0
            or (len(original_axes) == 1 and original_axes[0] < 0)
        )
        if not found:
            for axis in axestoregrid:
                if axis == stokes_axis:
                    found = True
                    break
        if found:
            # adjust stokes length to be equal to number of common stokes
            common_stokes_count = 0
            for image_stokes in image_csys.stokes():
                for template_stokes in template_csys.stokes():
                    if image_stokes == template_stokes:
                        common_stokes_count += 1
                        break
            shape[stokes_axis] = common_stokes_count
        else:
            shape[stokes_axis] = imshape[stokes_axis]
    return shape
