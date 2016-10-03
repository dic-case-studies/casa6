//# DSPixelRectangle.cc : Implementation of an absolute pixel DSRectangle
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

#include <display/DisplayShapes/DSPixelRectangle.h>
#include <display/DisplayShapes/DSWorldRectangle.h>
#include <display/DisplayShapes/DSScreenRectangle.h>

#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/QuantumHolder.h>
#include <casa/Quanta/UnitMap.h>

#include <casa/Containers/Record.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

	DSPixelRectangle::DSPixelRectangle() :
		DSRectangle() {

		UnitMap::putUser("pix", UnitVal(1.0), "absolute pixels");
	}

	DSPixelRectangle::DSPixelRectangle(const Record& settings) :
		DSRectangle() {
		UnitMap::putUser("pix", UnitVal(1.0), "absolute pixels");

		setOptions(settings);
	}

	DSPixelRectangle::DSPixelRectangle(DSScreenRectangle& other) :
		DSRectangle() {

		UnitMap::putUser("pix", UnitVal(1.0), "absolute pixels");

		DSRectangle::setOptions(other.getRawOptions());
	}

	DSPixelRectangle::DSPixelRectangle(DSWorldRectangle& other) :
		DSRectangle() {

		UnitMap::putUser("pix", UnitVal(1.0), "absolute pixels");

		DSRectangle::setOptions(other.getRawOptions());

	}


	DSPixelRectangle::~DSPixelRectangle() {

	}

	Record DSPixelRectangle::getOptions() {
		Record toReturn;

		toReturn = DSRectangle::getOptions();

		DisplayShapeWithCoords::floatPointToPix(toReturn, "center");
		DisplayShapeWithCoords::floatToPix(toReturn, "majoraxis");
		DisplayShapeWithCoords::floatToPix(toReturn, "minoraxis");


		// Shouldn't happen (should never be defined) .. but why not
		if (toReturn.isDefined("coords")) {
			toReturn.removeField("coords");
		}

		toReturn.define("coords", "pix");

		return toReturn;
	}

	Bool DSPixelRectangle::setOptions(const Record& settings) {

		Bool localChange = false;
		Record toSet = settings;

		if (settings.isDefined("coords")) {
			if (settings.asString("coords") != "pix") {
				throw(AipsError("I (DSPixelArrow) was expecting an option record which"
				                " had coords == \'pix\'. Please use a \'lock\' function"
				                " to change my co-ord system"));
			}
		}

		floatPointFromPix(toSet, "center");
		floatFromPix(toSet, "majoraxis");
		floatFromPix(toSet, "minoraxis");

		if (DSRectangle::setOptions(toSet)) {
			localChange = true;
		}

		return localChange;
	}




} //# NAMESPACE CASA - END

