//# DSPixelPolyLine.cc : Implementation of an absolute pixel DSPolyLine
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

#ifndef TRIALDISPLAY_DSPIXELPOLYLINE_H
#define TRIALDISPLAY_DSPIXELPOLYLINE_H

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>

#include <display/DisplayShapes/DSPolyLine.h>
#include <display/DisplayShapes/DisplayShapeWithCoords.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// Implementation of an absolute pixel DSPolyLine
// </summary>

	class PanelDisplay;
	class DSWorldPolyLine;
	class DSScreenPolyLine;

	class DSPixelPolyLine : public DSPolyLine , public DisplayShapeWithCoords {

	public:

		DSPixelPolyLine();
		DSPixelPolyLine(const casacore::Record& settings);
		DSPixelPolyLine(DSScreenPolyLine& other);
		DSPixelPolyLine(DSWorldPolyLine& other);

		virtual ~DSPixelPolyLine();

		virtual casacore::Bool setOptions(const casacore::Record& settings);
		virtual casacore::Record getOptions();

		virtual casacore::Record getRawOptions() {
			return DSPolyLine::getOptions();
		}

	private:

	};


} //# NAMESPACE CASA - END

#endif




