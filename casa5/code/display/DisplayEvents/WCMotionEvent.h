//# WCMotionEvent.h: class which stores WorldCanvas motion event information
//# Copyright (C) 1993,1994,1995,1996,1998,1999,2000
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

#ifndef TRIALDISPLAY_WCMOTIONEVENT_H
#define TRIALDISPLAY_WCMOTIONEVENT_H

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <display/DisplayEvents/WorldCanvasEvent.h>

namespace casa { //# NAMESPACE CASA - BEGIN

	class WorldCanvas;

// <summary>
// Class which stores WorldCanvas motion event information.
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>

// <etymology>
// "WCMotionEvent" is a contraction and concatentation of
// "WorldCanvas", "Motion" and "Event", and describes motion events
// occuring on WorldCanvases.
// </etymology>

// <synopsis>
// This class adds to the information stored in the <linkto
// class=WorldCanvasEvent>WorldCanvasEvent</linkto> class.  It adds
// informatino describing the current position of the mouse or
// pointing device, and the state of the keyboard modifiers (including
// the mouse buttons).
// </synopsis>

// <example>
// </example>

// <motivation>
// A compact way of passing motion event information around the
// WorldCanvas-oriented display classes was needed, with a functional
// but tight and efficient interface.
// </motivation>

// <thrown>
// None.
// </thrown>

// <todo asof="1999/10/15">
// None.
// </todo>

	class WCMotionEvent : public WorldCanvasEvent {

	public:

		// Constructor taking a pointer to the WorldCanvas for which the
		// event is valid, the state of the keyboard and pointer modifiers,
		// and the pixel, linear and world coordinates of the event.
		WCMotionEvent(WorldCanvas *wc, const casacore::uInt &modifiers,
		              const casacore::Int &pixX, const casacore::Int &pixY,
		              const casacore::Double &linX, const casacore::Double &linY,
		              const casacore::Vector<casacore::Double> &world);

		// copy constructor.
		WCMotionEvent(const WCMotionEvent &other);

		// Destructor.
		virtual ~WCMotionEvent();

		// The x and y pixel position of the pointer when the event
		// occurred.
		// <group>
		virtual casacore::Int pixX() const {
			return itsPixX;
		}
		virtual casacore::Int pixY() const {
			return itsPixY;
		}
		// </group>

		// The x and y linear coordinates of the event.
		// <group>
		virtual casacore::Double linX() const {
			return itsLinX;
		}
		virtual casacore::Double linY() const {
			return itsLinY;
		}
		// </group>

		// The world coordinates describing where the event occurred.
		virtual const casacore::Vector<casacore::Double> &world() const {
			return itsWorld;
		}

		// Return the state of the "modifiers": this is made up of mask bits
		// referring to various keys on the keyboard (eg. Control, Shift,
		// etc.) and the mouse buttons.
		virtual casacore::uInt modifiers() const {
			return itsModifiers;
		}

	protected:

		// (Required) default constructor.
		WCMotionEvent();

		// (Required) copy assignment.
		WCMotionEvent &operator=(const WCMotionEvent &other);
		friend class QtDisplayPanel;		/*** to access assignment operator to save "last event" ***/

	private:

		// Store the pixel position of the event here.
		casacore::Int itsPixX, itsPixY;

		// Store the linear position of the event here.
		casacore::Double itsLinX, itsLinY;

		// Store the world position of the event here.
		casacore::Vector<casacore::Double> itsWorld;

		// Store the button and keyboard modifier masks here.
		casacore::uInt itsModifiers;

	};


} //# NAMESPACE CASA - END

#endif


