//# QtOldMouseTools.qo.h: Qt versions of display library mouse tools.
//# Copyright (C) 2005
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


#ifndef QTOLDMOUSETOOLS_H
#define QTOLDMOUSETOOLS_H

#include <casa/aips.h>
#include <display/DisplayEvents/MWCRTRegion.h>
#include <display/DisplayEvents/MWCETRegion.h>
#include <display/DisplayEvents/MWCPTRegion.h>
#include <display/Display/PanelDisplay.h>
#include <casa/Containers/Record.h>

#include <graphics/X11/X_enter.h>
#  include <QObject>
#include <graphics/X11/X_exit.h>

namespace casa {



// <synopsis>
// Nothing yet: it may prove useful for Qt-signal-emitting mouse tools
// (which are MWCTools or possibly PCTools) to have a common base.
// </synopsis>
	class QtOldMouseTool: public QObject {

		Q_OBJECT	//# Allows slot/signal definition.  Must only occur in
		//# implement/.../*.h files; also, makefile must include
		//# name of this file in 'mocs' section.

	public:

		QtOldMouseTool() : QObject() {  }
		~QtOldMouseTool() {  }

	};




// <synopsis>
// QtOldRTRegion is the Rectangle Region mouse tool that sends a signal
// when a new rectangle is ready.
// </synopsis>
	class QtOldRTRegion: public QtOldMouseTool, public MWCRTRegion {

		Q_OBJECT	//# Allows slot/signal definition.  Must only occur in
		//# implement/.../*.h files; also, makefile must include
		//# name of this file in 'mocs' section.

	public:

		QtOldRTRegion(PanelDisplay* pd) : QtOldMouseTool(), MWCRTRegion(), pd_(pd) {  }

		~QtOldRTRegion() {  }

		// Retrieve the current rectangular mouse region record and WCH, if any.
		// (If nothing is ready, returns false -- be sure to check before using
		// return parameters.  See implementation for mouseRegion casacore::Record format).
		bool getMouseRegion(casacore::Record& mouseRegion, WorldCanvasHolder*& wch);

	signals:

		// See regionReady() implementation for format of the record.  (For some
		// uses, a connecting slot may be able to do without the WCH* parameter).
		void mouseRegionReady(casacore::Record mouseRegion, WorldCanvasHolder*);
		void echoClicked(casacore::Record);

	protected:

		// Signals mouseRegionReady with an appropriate casacore::Record, when
		// called by base class in response to user selection with the mouse.
		// See implementation for format of the record.
		virtual void regionReady();

		virtual void clicked(int x, int y);
		virtual void doubleClicked(int x, int y);
		//virtual void rectangleReady();
		//virtual void handleEvent(DisplayEvent& ev);
		//virtual void keyPressed(const WCPositionEvent &ev);

		PanelDisplay* pd_;	// (Kludge... zIndex inaccessible from WC...)

	};


// <synopsis>
// QtOldELRegion is the Ellipse Region mouse tool that sends a signal
// when a new circle is ready.
// </synopsis>
	class QtOldELRegion: public QtOldMouseTool, public MWCETRegion {

		Q_OBJECT	//# Allows slot/signal definition.  Must only occur in
		//# implement/.../*.h files; also, makefile must include
		//# name of this file in 'mocs' section.

	public:

		QtOldELRegion(PanelDisplay* pd) : QtOldMouseTool(), MWCETRegion(), pd_(pd) {  }

		~QtOldELRegion() {  }

		// Retrieve the current circular mouse region record and WCH, if any.
		// (If nothing is ready, returns false -- be sure to check before using
		// return parameters.  See implementation for mouseRegion casacore::Record format).
		bool getMouseRegion(casacore::Record& mouseRegion, WorldCanvasHolder*& wch);

	signals:

		// See regionReady() implementation for format of the record.  (For some
		// uses, a connecting slot may be able to do without the WCH* parameter).
		void mouseRegionReady(casacore::Record mouseRegion, WorldCanvasHolder*);
		void echoClicked(casacore::Record);

	protected:

		// Signals mouseRegionReady with an appropriate casacore::Record, when
		// called by base class in response to user selection with the mouse.
		// See implementation for format of the record.
		virtual void regionReady();

		virtual void clicked(int x, int y);
		virtual void doubleClicked(int x, int y);

		PanelDisplay* pd_;	// (Kludge... zIndex inaccessible from WC...)

	};


// <synopsis>
// QtOldPTRegion is the Polygon Region mouse tool that sends a signal
// when a new polygon is ready.
// </synopsis>
	class QtOldPTRegion: public QtOldMouseTool, public MWCPTRegion {

		Q_OBJECT	//# Allows slot/signal definition.  Must only occur in
		//# implement/.../*.h files; also, makefile must include
		//# name of this file in 'mocs' section.

	public:

		QtOldPTRegion(PanelDisplay* pd) : QtOldMouseTool(), MWCPTRegion(), pd_(pd) {  }

		~QtOldPTRegion() {  }

		// Retrieve the current polygon mouse region record and WCH, if any.
		// (If nothing is ready, returns false -- be sure to check before using
		// return parameters.  See implementation for mouseRegion casacore::Record format).
		bool getMouseRegion(casacore::Record& mouseRegion, WorldCanvasHolder*& wch);

	signals:

		// See regionReady() implementation for format of the record.  (For some
		// uses, a connecting slot may be able to do without the WCH* parameter).
		void mouseRegionReady(casacore::Record mouseRegion, WorldCanvasHolder*);
		void echoClicked(casacore::Record);

	protected:

		// This callback is invoked by the base when the user double-clicks
		// inside a polygon defined previously (but see also polygonReady(),
		// below).  This implementation emits the Qt signal mouseRegionReady()
		// with an appropriate casacore::Record defining the user's polygon mouse selection.
		// See implementation for format of the record.
		virtual void regionReady();
		virtual void clicked(int x, int y);
		virtual void doubleClicked(int x, int y);

		//virtual void handleEvent(DisplayEvent& ev);
		//virtual void keyPressed(const WCPositionEvent &ev);

		///////This is not correct !
		//
		// This callback is invoked by the base when the polygon is first
		// defined (by a double-click at last point) or when the mouse is
		// released after a move/resize.  It was unused in glish.  For Qt,
		// this will also signal that the polygon region has been fully
		// 'selected/made ready' (which will preclude the need for the
		// user to double-click [again] inside the polygon to select it).
		//virtual void polygonReady() { regionReady();  }

		// is this fix to 1393?
		virtual void polygonReady() { }

		PanelDisplay* pd_;	// (Kludge... zIndex inaccessible from WC...)

	};




} //# NAMESPACE CASA - END

#endif

