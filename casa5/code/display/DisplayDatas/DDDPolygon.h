//# DDDPolygon.h: declaration of DDDPolygon
//# Copyright (C) 1999,2000,2001,2002
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


#ifndef TRIALDISPLAY_DDDPOLYGON_H
#define TRIALDISPLAY_DDDPOLYGON_H

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Containers/Record.h>
#include <casa/Containers/Block.h>
#include <casa/Logging/LogIO.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/Unit.h>
#include <display/Display/DisplayCoordinateSystem.h>
#include <display/DisplayDatas/DDDObject.h>
#include <display/DisplayDatas/DDDHandle.h>

namespace casacore{

	class LogIO;
}

namespace casa { //# NAMESPACE CASA - BEGIN


// <summary>
// Implementation of aan ellipse object for DrawingDisplayData class.
// </summary>

// <synopsis>
// This class implements the interface defined by DDDObject, to provide
// a polygon object for registration with DrawingDisplayData objects.
// </synopsis>
//
// <todo>
// Allow vertices to be edited.
// </todo>
//

	class DDDPolygon : public DDDObject {

	public:

		// Constructor taking a casacore::Record description.  Fields in the record,
		// on top of what is consumed by the DDDObject constructor, are:
		// Required : <src>x</src> (quantum of vector double) and
		// <src>y</src> (quantum of vector double).
		// Other, not required fields are <src>outline & editable</src>
		// which default to <src>true</src>.  When <src>editable=false</src>
		// the object can be moved but not reshaped, when <src>true</src>
		// it can be moved and reshaped.   When   <src>outline=false</src>
		// the ellipse is filled.
		DDDPolygon(const casacore::Record &description, DrawingDisplayData *owner);

		// Destructor.
		virtual ~DDDPolygon();

		// Draw this rectangle object for the given reason on the provided
		// WorldCanvas.
		virtual void draw(const Display::RefreshReason &reason,
		                  WorldCanvas *worldcanvas);

		// Return a record describing this object.  Presently returns
		// only construction description.
		virtual casacore::Record description();

		// Update this object based on the information in the provided
		// Record.
		virtual void setDescription(const casacore::Record &rec);

		// Event handlers.  The parent DrawingDisplayData will distribute
		// events as necessary to the various DDDObjects which comprise it.
		// <group>
		virtual void operator()(const WCRefreshEvent &/*ev*/) {
			;
		};
		virtual void operator()(const WCPositionEvent &ev);
		virtual void operator()(const WCMotionEvent &ev);
		// </group>




	protected:

		// (Required) default constructor.
		DDDPolygon();

		// (Required) copy constructor.
		DDDPolygon(const DDDPolygon &other);

		// (Required) copy assignment.
		void operator=(const DDDPolygon &other);

	private:

// Fill style
		casacore::Bool itsOutline;

// Define x and y
		casacore::Quantum<casacore::Vector<casacore::Double> > itsWorldX, itsWorldY;
		casacore::Vector<casacore::Double> itsPixelX, itsPixelY;
		casacore::uInt itsNPoints;

// Rotated Rectangle Corners to draw in screen pixel
		casacore::Matrix<casacore::Double> itsCorners;

// casacore::List of handles

		casacore::Block<DDDHandle> itsHandles;

// Mode.

		DDDObject::Mode itsMode;

// Store for movement bases.

		casacore::Double itsBaseMoveX, itsBaseMoveY;

// In function draw(), recompute all pixel coordinates when true
// else use what is currently set

		casacore::Bool itsRecompute;

// LOgger

		casacore::LogIO itsLogger;

// DisplayCoordinateSystem in screen pixel coordinates
		DisplayCoordinateSystem itsCoordinateSystem;

// pointer to world canvas (shallow copy)
		WorldCanvas* itsWorldCanvasPtr;

// Fractional indicators
		casacore::Bool itsFracX, itsFracY;

// COnvert parameters from world to pixel
		void convertToPixel ();

// Update private world parameters from current pixel values
		void updateWorldValues ();

// Decode record into private data
		void decode(const casacore::RecordInterface& description, casacore::Bool required);

// Encode private data into record
		void encode(casacore::RecordInterface& description) const;

// compute corners of rectangle of bounding box of polygon
		void createCorners();

// create handles, one per vertex
		void createHandles();


// Debugging routines

		void listWorld();
		void listPixel();

	};


} //# NAMESPACE CASA - END

#endif

