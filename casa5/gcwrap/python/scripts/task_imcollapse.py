
##########################################################################
# task_imcollapse.py
#
# Copyright (C) 2008, 2009, 2010
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
# Task to collapse an image along a specified axis,
# computing a specified aggregate function of pixels along that axis
# </summary>
#
# <reviewed reviwer="" date="" tests="" demos="">
# </reviewed>
#
# <prerequisite>
# <ul>
#
# </ul>
# </prerequisite>
#
# <etymology>
# imtrans => im(age) collapse
# </etymology>
#
# <synopsis>
# imtrans collapses an image along a specified axis. It is built on top of ia.collapse()
# </synopsis> 
#
# <example>
# collapsed_image_tool = imcollapse(imagename="myim.im", outfile="collapsed.im", axis=2, function="variance", wantreturn=true)
#
# </example>
#
# <motivation>
# To make users happy (https://bugs.aoc.nrao.edu/browse/CAS-1222)
# and associated casacore class is prereq for specfit work
# </motivation>
#

###########################################################################

from __future__ import absolute_import
import sys

# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatools import image
    from casatasks import casalog
    from .ialib import write_image_history
else:
    from taskinit import *
    from ialib import write_image_history
    image = iatool

def imcollapse(
    imagename=None, function=None, axes=None, outfile=None, box=None,
    region=None, chans=None, stokes=None, mask=None,
    overwrite=None, stretch=None
):
    casalog.origin('imcollapse')
    myia = image()
    outia = None
    try :
        if len(outfile) == 0:
            raise RuntimeError("outfile must be specified")
        myia.dohistory(False)
        if (not myia.open(imagename)):
            raise RuntimeError("Cannot create image analysis tool using " + imagename)
        outia = myia.collapse(
            function, axes, outfile, region, box, chans,
            stokes, mask, overwrite, stretch
        )
        try:
            if is_CASA6:
                vars = locals()
                param_names = imcollapse.__code__.co_varnames[:imcollapse.__code__.co_argcount]
                param_vals = [vars[p] for p in param_names]
            else:
                param_names = imcollapse.func_code.co_varnames[:imcollapse.func_code.co_argcount]
                param_vals = [eval(p) for p in param_names]
            write_image_history(
                outia, sys._getframe().f_code.co_name,
                param_names, param_vals, casalog
            )
        except Exception as instance:
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')
        # returns None as per requirement
    finally:
        if myia:
            myia.done()
        if outia:
            outia.done()
