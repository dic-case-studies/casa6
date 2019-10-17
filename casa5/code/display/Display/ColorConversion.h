//# ColorConversion.h: utilities for converting between color spaces
//# Copyright (C) 1994,1995,1996,1997,1998,1999,2000
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

#ifndef TRIALDISPLAY_COLORCONVERSION_H
#define TRIALDISPLAY_COLORCONVERSION_H

// <summary>
// Utilities for converting between color spaces
// </summary>
//
// <synopsis>
// provides vectorized and single color conversion between RGB and
// HSV color spaces.  Operates using float type with the range of
// 0.0 to 1.0 for all R,G,B,H,S,V values.
// </synopsis>
//
// <motivation>
// want to be able to create something that acts like an "HSV" colormap
// </motivation>
//
// <group name="RGB/HSV casacore::Conversion Functions">


/*


#include <graphics/X11/X_enter.h>
#include <X11/Xlib.h>
#include <graphics/X11/X_exit.h>


*/
#include <casa/aips.h>
#include <casa/Arrays/Array.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// map a single hsv value to an rgb value
// assumes h, s, and v to be in range [0, 1]. Returns r, g, b in range [0, 1].
	void hsvToRgb(casacore::Float h, casacore::Float s, casacore::Float v,
	              casacore::Float &red, casacore::Float &green, casacore::Float &blue);

// map an array of hsv values to an array of rgb value
	void hsvToRgb(const casacore::Array<casacore::Float> & h, const casacore::Array<casacore::Float> & s,
	              const casacore::Array<casacore::Float> & v,
	              casacore::Array<casacore::Float> & r, casacore::Array<casacore::Float> & g, casacore::Array<casacore::Float> & b);

// map a single rgb value to its hsv value
//Assumes r, g, and b to be in range [0, 1]. Returns h, s, v in range [0, 1].
	void rgbToHsv(casacore::Float r, casacore::Float g, casacore::Float b,
	              casacore::Float &h, casacore::Float &s, casacore::Float &v);

// map an array of rgb values to an array of hsv values
	void rgbToHsv(const casacore::Array<casacore::Float> & r, const casacore::Array<casacore::Float> & g,
	              const casacore::Array<casacore::Float> & b,
	              casacore::Array<casacore::Float> & h, casacore::Array<casacore::Float> & s, casacore::Array<casacore::Float> & v);


// </group>


} //# NAMESPACE CASA - END

#endif
