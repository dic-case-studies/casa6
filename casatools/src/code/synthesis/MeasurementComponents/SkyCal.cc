//# SkyCal.cc:  this defines ClassName, which ...
//# Copyright (C) 2003
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

//# Includes
// synthesis module cannot depend on Sakura at this moment since it is required by
// wvrgcal which produces an executable. Once we move on to c++11 compiler, we
// can safely use Sakura functions inside synthesis.
//#include <libsakura/sakura.h>
#include <synthesis/MeasurementComponents/SkyCal.h>
#include <casacore/casa/BasicSL/Complex.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

template class SkyCal<Float, Float>;
template class SkyCal<Complex, Complex>;

} //# NAMESPACE CASA - END

