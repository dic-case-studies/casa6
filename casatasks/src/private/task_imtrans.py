##########################################################################
# task_imtrans.py
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
# Task to transpose an image
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
# imtrans => im(age) trans(pose)
# </etymology>
#
# <synopsis>
# imtrans transposes an image. It is built on top of ia.reorder()
# </synopsis> 
#
# <example>
# imtrans(imagename="myim.im", outfile="transposed.im", order="102")
#
# </example>
#
# <motivation>
# To make users happy, cf https://bugs.aoc.nrao.edu/browse/CAS-607
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

def imtrans(imagename, outfile, order):
    casalog.origin('imtrans')
    myia = image()
    myia.dohistory(False)
    outia = None
    try:
        if (not myia.open(imagename)):
            raise RuntimeError("Cannot create image analysis tool using " + imagename)
        if (len(outfile) == 0):
            raise ValueError("outfile parameter must be specified.")
        outia = myia.transpose(outfile=outfile, order=order)
        try:
            param_names = imtrans.__code__.co_varnames[:imtrans.__code__.co_argcount]
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
            casalog.post("*** Error \'%s\' updating HISTORY" % (instance), 'WARN')

    finally:
        if myia:
            myia.done()
        if outia:
            outia.done()
