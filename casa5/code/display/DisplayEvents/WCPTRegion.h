//# WCPTRegion.h: WorldCanvas event-based polygon region drawer
//# Copyright (C) 1999,2000
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

#ifndef TRIALDISPLAY_WCPTREGION_H
#define TRIALDISPLAY_WCPTREGION_H

//# aips includes
#include <casa/aips.h>
#include <casa/Arrays/Vector.h>

//# trial includes

//# display library includes
#include <display/DisplayEvents/WCPolyTool.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// WorldCanvas event-based polygon region drawer
// </summary>

// <use visibility=export>

// <reviewed reviewer="" data="yyyy/mm/dd" tests="" demos="">
// </reviewed>

// <prerequisites>
//   <li> WCRectTool
// </prerequisites>

// <etymology>
// WCPTRegion stands for WorldCanvas Polygon Tool Region
// </etymology>

// <synopsis>
// This class finalises the functionality in WCRectTool to
// implement event-based drawing of polygonal regions on a
// WorldCanvas.  casacore::Double clicking inside the constructed polygon
// emits the region, whilst pressing Escape cancels the region.
// </synopsis>

// <example>
// </example>

// <motivation>
// Polygonal regions are the most common regions, and constructing
// them visually can be very efficient.
// </motivation>

// <todo asof="1999/02/10">
//   <li> Nothing known
// </todo>

	class WCPTRegion : public WCPolyTool {

	public:

		// Constructor
		WCPTRegion(WorldCanvas *wcanvas,
		           Display::KeySym keysym = Display::K_Pointer_Button1);

		// Destructor
		virtual ~WCPTRegion();

		// Functions specific to the region drawing - these are called by
		// lower-level event handlers in the base classes
		virtual void doubleInside();

		// This functio is called when a region is deliberately "finished"
		// by the user.  It can be implemented to do further actions, such
		// as report the region to the application
		virtual void regionReady() { };

	private:

		// get the linear coords of the current polygonal region
		//void getLinearCoords(casacore::Vector<casacore::Double> &x, casacore::Vector<casacore::Double> &y);

	};


} //# NAMESPACE CASA - END

#endif
