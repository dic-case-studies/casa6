//# DisplayTool.cc: base class for event-based tools in the display classes
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

#include <casa/Exceptions/Error.h>
#include <display/DisplayEvents/DisplayTool.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

// Destructor.
	DisplayTool::~DisplayTool() {
	}

// Set/get the key to catch.
	void DisplayTool::setKey(const Display::KeySym &keysym) {
		itsKeySym = keysym;
		chooseKeyModifier();
	}
	Display::KeySym DisplayTool::getKey() const {
		return itsKeySym;
	}

// Constructor.
	DisplayTool::DisplayTool(const Display::KeySym &keysym) :
		itsKeySym(keysym) {
		chooseKeyModifier();
	}

// (Required) copy constructor.
	DisplayTool::DisplayTool(const DisplayTool &other) :
		itsKeySym(other.itsKeySym),
		itsKeyModifier(other.itsKeyModifier) {
	}

// (Required) copy assignment.
	DisplayTool &DisplayTool::operator=(const DisplayTool &other) {
		if (this != &other) {
			itsKeySym = other.itsKeySym;
			itsKeyModifier = other.itsKeyModifier;
		}
		return *this;
	}

	void DisplayTool::chooseKeyModifier() {
		try {
			itsKeyModifier = Display::keyModifierFromKeySym(itsKeySym);
		} catch (AipsError& x) {
			// use x to avoid compiler warning
			//if (&x) {
				itsKeyModifier = static_cast<Display::KeyModifier>(0);
			//}
		}
	}

} //# NAMESPACE CASA - END

