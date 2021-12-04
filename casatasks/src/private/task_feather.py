from __future__ import absolute_import
import os

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
        from casatools import imager
        from casatasks import casalog
else:
        from taskinit import *
        imager = imtool

def feather(imagename=None,highres=None,lowres=None, sdfactor=None, effdishdiam=None, lowpassfiltersd=None):
        """ Feathering: Combine two images using the Fourier Addition:
        
        The algorithm converts each image to the gridded visibility
        plane, combines them, and reconverts them into an combined
        image.  The images must have a well-defined beam shape (clean beam),
        given in the image, in order for feathering to work well.
        The two images must be on the same grid and have the same
        flux density normalization scale.

        Keyword arguments:
        imagename -- Name of output feathered image
                default: none; example: imagename='orion_combined'
        highres -- Name of high resolution (interferometer) image
                default: none; example: imagename='orion_vla.im'
             This image is often a clean image obtained from synthesis
                observations.
        lowres -- Name of low resolution (single dish) image
                default: none; example: imagename='orion_gbt.im'
             This image is often a image from a single-dish observations
                or a clean image obtained from a lower resolution synthesis
                observations.

        """
        casalog.origin('feather')

        imFea=imager( )
        try:
                imFea.setvp(dovp=True)
                imFea.setsdoptions(scale=sdfactor)
                imFea.feather(image=imagename,highres=highres,lowres=lowres,
                              effdishdiam=effdishdiam,  lowpassfiltersd=lowpassfiltersd)
        finally:
                imFea.done( )
                del imFea
