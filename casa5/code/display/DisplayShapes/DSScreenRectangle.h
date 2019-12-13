//# DSScreenRectangle : Implementation of a relative screen pos. DSRectangle
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
#ifndef TRIALDISPLAY_DSSCREENRECTANGLE_H
#define TRIALDISPLAY_DSSCREENRECTANGLE_H

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>

#include <display/DisplayShapes/DSRectangle.h>
#include <display/DisplayShapes/DisplayShapeWithCoords.h>

// <summary>
// Implementation of a relative screen pos. DSRectangle
// </summary>

namespace casa { //# NAMESPACE CASA - BEGIN

	class DSPixelRectangle;
	class DSWorldRectangle;
	class PixelCanvas;

	class DSScreenRectangle : public DSRectangle, public DisplayShapeWithCoords {

	public:

		DSScreenRectangle();
		DSScreenRectangle(const casacore::Record& setting, PixelCanvas* pc);
		DSScreenRectangle(DSPixelRectangle& other, PixelCanvas* pc);
		DSScreenRectangle(DSWorldRectangle& other);

		virtual ~DSScreenRectangle();

		virtual void recalculateScreenPosition();
		virtual casacore::Bool setOptions(const casacore::Record& settings);
		virtual casacore::Record getOptions();

		virtual void move(const casacore::Float& dX, const casacore::Float& dy);
		virtual void setCenter(const casacore::Float& xPos, const casacore::Float& yPos);
		virtual void rotate(const casacore::Float& angle);
		virtual void changePoint(const casacore::Vector<casacore::Float>&pos, const casacore::Int n);
		virtual void changePoint(const casacore::Vector<casacore::Float>& pos);
		// PLUS OTHERS! :)  change point etc

		virtual casacore::Record getRawOptions() {
			return DSRectangle::getOptions();
		}

		virtual void updateRelative();

	private:
		PixelCanvas* itsPC;
		casacore::Vector<casacore::Float> itsRelativeCenter;

		casacore::Float itsRelWidth;
		casacore::Float itsRelHeight;

	};

} //# NAMESPACE CASA - END

#endif
