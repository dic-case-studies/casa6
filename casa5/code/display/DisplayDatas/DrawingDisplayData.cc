//# DrawingDisplayData.cc: interactive drawing DisplayData
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

#include <casa/aips.h>
#include <casa/Exceptions.h>
#include <casa/Containers/Record.h>
#include <casa/Logging/LogIO.h>
#include <casa/System/AipsrcValue.h>
#include <display/DisplayDatas/DrawingDisplayMethod.h>
#include <display/DisplayDatas/DrawingDisplayData.h>
#include <display/DisplayDatas/DDDRectangle.h>
#include <display/DisplayDatas/DDDEllipse.h>
#include <display/DisplayDatas/DDDPolygon.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

	DrawingDisplayData::DrawingDisplayData(const Display::KeySym keysym) :
		PassiveCachingDD(),
		itsKeySym(keysym),
		itsObjectWhichIsShowingHandles(0) {
		setCaching(false);
		try {
			itsKeyModifier = Display::keyModifierFromKeySym(itsKeySym);
		} catch(AipsError x) {
			itsKeyModifier = (Display::KeyModifier)0;
		}
		AipsrcValue<Double>::find(itsDoubleClickInterval,
		                          "display.controls.doubleclickinterval",
		                          Double(0.5));
		installDefaultOptions();
	}

	DrawingDisplayData::~DrawingDisplayData() { }

	void DrawingDisplayData::setDefaultOptions() {
		PassiveCachingDD::setDefaultOptions();
		installDefaultOptions();
	}

	Bool DrawingDisplayData::setOptions(Record &rec, Record &recOut) {
		Bool ret = PassiveCachingDD::setOptions(rec, recOut);
		Bool localchange = false, error;

		localchange = (readOptionRecord(itsOptionsLabelPosition, error, rec,
		                                "labelposition") || localchange);

		ret = (ret || localchange);
		return ret;
	}

	Record DrawingDisplayData::getOptions( bool scrub ) const {
		Record rec = PassiveCachingDD::getOptions(scrub);

		Record labelposition;
		labelposition.define("dlformat", "labelposition");
		labelposition.define("listname", "Label position");
		labelposition.define("ptype", "choice");
		Vector<String> vlabelposition(2);
		vlabelposition(0) = "none";
		vlabelposition(1) = "centre";
		labelposition.define("popt", vlabelposition);
		labelposition.define("default", "none");
		labelposition.define("value", itsOptionsLabelPosition);
		labelposition.define("allowunset", false);
		rec.defineRecord("labelposition", labelposition);

		return rec;
	}

	CachingDisplayMethod *DrawingDisplayData::newDisplayMethod(
	    WorldCanvas *worldCanvas,
	    AttributeBuffer *wchAttributes,
	    AttributeBuffer *ddAttributes,
	    CachingDisplayData *dd) {
		return new DrawingDisplayMethod(worldCanvas, wchAttributes,
		                                ddAttributes, dd);
	}

	AttributeBuffer DrawingDisplayData::optionsAsAttributes() {
		AttributeBuffer buffer = PassiveCachingDD::optionsAsAttributes();
		return buffer;
	}
	void DrawingDisplayData::refreshEH(const WCRefreshEvent &ev) {
		for ( auto temp : itsDDDOList ) ((DDDObject *)temp)->operator()(ev);
		PassiveCachingDD::refreshEH(ev);
	}

	void DrawingDisplayData::addObject(const Record &description) {
		// YTBI: need to check that no existing object in the list has the
		// same id as this new one
		if (!description.isDefined("type")) {
			throw(AipsError("No 'type' field in object description Record"));
		}
		String type;
		description.get("type", type);
		DDDObject* dddObject= 0;
		if (type == "rectangle") {
			dddObject = new DDDRectangle(description, this);
		} else if (type == "ellipse") {
			dddObject = new DDDEllipse(description, this);
		} else if (type == "polygon") {
			dddObject = new DDDPolygon(description, this);
		} else {
			throw(AipsError("Unknown 'type' field in object description Record"));
		}

		// preferentially add to start of list for fast access to most
		// recently added object
		itsDDDOList.push_front((void *)dddObject);

		// install event handlers
		addPositionEventHandler(dddObject);

		// call refresh
		refresh();

	}

	Record DrawingDisplayData::description(const Int objectID) {
		Bool found = false;
		Record rec;
		for ( auto temp : itsDDDOList ) {
			if ( objectID == ((DDDObject *)temp)->objectID( ) ) {
				found = true;
				rec = ((DDDObject *)temp)->description();
				break;
			}
		}
//
		if (!found) {
			throw(AipsError("Couldn't find object with given id"));
		}
//
		return rec;
	}

	void DrawingDisplayData::setDescription(const Int objectID,
											const Record &rec) {
		Bool found = false;
		for ( auto temp : itsDDDOList ) {
			if ( objectID == ((DDDObject *)temp)->objectID( ) ) {
				found = true;
				((DDDObject *)temp)->setDescription(rec);
				break;
			}
		}

		if (!found) {
			LogIO os(LogOrigin("DrawingDisplayDatas", "setDescription(...)",
							   WHERE));
			os << LogIO::WARN << "Could not find object with given id" << LogIO::POST;
		}
	}

	void DrawingDisplayData::removeObject(const Int objectID) {
		std::list<void*> orig = itsDDDOList;
		std::list<void*> to_remove; 
		itsDDDOList.clear( );
		std::partition_copy( orig.begin( ), orig.end( ),
							 std::back_inserter(to_remove),
							 std::back_inserter(itsDDDOList),
							 [&](void *vp){return ((DDDObject *)vp)->objectID( ) == objectID;} );
		for ( void *vp : to_remove ) {
			auto temp = (DDDObject *) vp;
			if ( itsObjectWhichIsShowingHandles == temp )
				itsObjectWhichIsShowingHandles = 0;
			temp->showHandles(false, false);				// this removes motion EH
			removePositionEventHandler(*temp);
		}

		if ( to_remove.size( ) == 0 ) {
			LogIO os(LogOrigin("DrawingDisplayDatas", "removeObject(...)",
							   WHERE));
			os << LogIO::WARN << "Could not find object with given id" << LogIO::POST;
		} else {
			for ( void *vp : to_remove ) {
				auto temp = (DDDObject *) vp;
				delete temp;
			}
			// call refresh
			refresh();
		}
	}

	void DrawingDisplayData::setHandleState(DDDObject *item, const Bool state) {

		bool found = std::any_of( itsDDDOList.begin( ), itsDDDOList.end( ), [&](void *vp){return vp == item;} );

		if (!found) throw(AipsError("Cannot find object in list"));

		if (state) {
			if (itsObjectWhichIsShowingHandles) {
				// switch other off
				itsObjectWhichIsShowingHandles->showHandles(false, false);
			}
			itsObjectWhichIsShowingHandles = item;
		} else {
			if (itsObjectWhichIsShowingHandles == item) {
				itsObjectWhichIsShowingHandles = 0;
			}
		}
		refresh();
	}

// Set the key to catch.
	void DrawingDisplayData::setKey(const Display::KeySym &keysym) {
		itsKeySym = keysym;
		try {
			itsKeyModifier = Display::keyModifierFromKeySym(itsKeySym);
		} catch (AipsError x) {
			itsKeyModifier = (Display::KeyModifier)0;
		}
	}

	void DrawingDisplayData::doubleClick(const Int objectID) {
		cerr << "double click on object " << objectID << " detected" << endl;
	}

	DrawingDisplayData::DrawingDisplayData(const DrawingDisplayData &o) : PassiveCachingDD(o) {
	}

	void DrawingDisplayData::operator=(const DrawingDisplayData &/*other*/) {
	}

	void DrawingDisplayData::installDefaultOptions() {
		itsOptionsLabelPosition = "none";
	}

} //# NAMESPACE CASA - END

