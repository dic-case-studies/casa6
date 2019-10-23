//# DSPixelRectangle.h : Implementation of an absolute pixel DSRectangle
//# Copyright (C) 1998,1999,2000,2001,2002
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
//# $Id:

#ifndef TRIALDISPLAY_DSPIXELRECTANGLE_H
#define TRIALDISPLAY_DSPIXELRECTANGLE_H

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>

#include <display/DisplayShapes/DisplayShapeWithCoords.h>
#include <display/DisplayShapes/DSRectangle.h>

namespace casacore{

	class Record;
}

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// Implementation of an absolute pixel DSRectangle
// </summary>

	class DSWorldRectangle;
	class DSScreenRectangle;

	class DSPixelRectangle : public DSRectangle, public DisplayShapeWithCoords {

	public :
		DSPixelRectangle();
		DSPixelRectangle(const casacore::Record& settings);
		DSPixelRectangle(DSScreenRectangle& other);
		DSPixelRectangle(DSWorldRectangle& other);

		virtual ~DSPixelRectangle();

		virtual casacore::Bool setOptions(const casacore::Record& settings);
		virtual casacore::Record getOptions();

		virtual casacore::Record getRawOptions() {
			return DSRectangle::getOptions();
		}

	private:


	};

} //# NAMESPACE CASA - END

#endif
