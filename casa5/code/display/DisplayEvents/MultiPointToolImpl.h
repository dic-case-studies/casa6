//# MultiPointToolImpl.h: WorldCanvas event-based point region drawer
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

#ifndef TRIALDISPLAY_MULTIPOINTTOOLIMPL_H
#define TRIALDISPLAY_MULTIPOINTTOOLIMPL_H

//# aips includes
#include <casa/aips.h>
#include <casa/Arrays/Vector.h>

//# trial includes

//# display library includes
#include <display/DisplayEvents/MultiPointTool.h>
#include <display/region/RegionSourceFactory.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// WorldCanvas event-based point region drawer
// </summary>

// <use visibility=export>

// <reviewed reviewer="" data="yyyy/mm/dd" tests="" demos="">
// </reviewed>

// <prerequisites>
//   <li> MWCPointTool
// </prerequisites>

// <etymology>
// MultiPointToolImpl stands for MultiWorldCanvas Point Tool Region
// </etymology>

// <synopsis>
// Why this "finalization" is required is beyond me... i.e. why not
// roll this functionality into the parent class??
// <drs> Thu Aug 18 20:13:53 UTC 2011
// ------------------------------------------------------------------------
// This class finalises the functionality in MWCPointTool to
// implement event-based drawing of elliptical regions on a
// WorldCanvas.  casacore::Double clicking inside the constructed point
// emits the region, whilst pressing Escape cancels the region.
// </synopsis>

// <example>
// </example>

// <motivation>
// Elliptical regions are the most common regions, and constructing
// them visually can be very efficient.
// </motivation>

// <todo asof="1999/02/10">
//   <li> Nothing known
// </todo>

	class MultiPointToolImpl : public MultiPointTool {

	public:

		// Constructor
		MultiPointToolImpl(viewer::RegionSourceFactory *rf, PanelDisplay* pd, Display::KeySym keysym = Display::K_Pointer_Button2);

		// Destructor
		virtual ~MultiPointToolImpl();

		// handle events, via new-style interface.  Currently just for reset event.
		virtual void handleEvent(DisplayEvent& ev);

	protected:

		// Handle double-click inside the point.  Invokes regionReady().
		virtual void doubleInside();

		// This function is called when a region is deliberately "finished"
		// by the user.  It can be implemented to do further actions, such
		// as report the region to the application
		virtual void regionReady() { };

		// get the world coords of the current elliptical region
		void getWorldCoords(casacore::Vector<casacore::Double> &blc, casacore::Vector<casacore::Double> &trc);

		// get the linear coords of the current elliptical region
		void getLinearCoords(casacore::Vector<casacore::Double> &blc, casacore::Vector<casacore::Double> &trc);

	};


} //# NAMESPACE CASA - END

#endif
