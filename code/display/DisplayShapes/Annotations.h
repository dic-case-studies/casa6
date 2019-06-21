//# Annotations.h : Control class for DisplayShape
//# Copyright (C) 1998,1999,2000,2001,2002,2003
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

#ifndef TRIALDISPLAY_ANNOTATIONS_H
#define TRIALDISPLAY_ANNOTATIONS_H

#include <casa/aips.h>

#include <display/DisplayEvents/PCPositionEH.h>
#include <display/DisplayEvents/PCMotionEH.h>
#include <display/DisplayEvents/PCRefreshEH.h>
#include <display/DisplayEvents/WCRefreshEH.h>
#include <casa/OS/Timer.h>

#include <display/Display/WorldCanvas.h>
#include <list>

#include <casa/Containers/Record.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>

namespace casa { //# NAMESPACE CASA - BEGIN

	class DisplayShape;
	class PanelDisplay;


// <summary>
// Control class for a collection of DisplayShapeInterfaces.
// </summary>
//
// <prerequisite>
// <li> <linkto class="DisplayShapeInterface">DisplayShapeInterface</linkto>
// </prerequisite>
//
// <etymology>
// Annotations is a method by which a collection of shapes can be managed
// on a pixel canvas.
// </etymology>
//
// <synopsis>
// Annotations registers itself as an event handler on whichever pixelcanvas
// it is supplied. It can then add, delete or modify shapes based on function
// calls and user mouse input.
//
// The class is essentially a casacore::PtrBlock of DisplayShapeInterfaces, with
// pixelcanvas event handling to control them.
// </synopsis>
//
// <motivation>
// To allow a collection of shapes to be displayed and cotrolled by one class.
// </motivation>
//
// <example>
// <srcblock>
// </srcblock>
// </example>


	class Annotations : public PCPositionEH,
		public PCMotionEH,
		public PCRefreshEH,
		public WCRefreshEH {

	public:
		// Used to determine the state of the class. The state may change due
		// to mouse events, or function calls.
		enum State {Nothing, Creation, WcCreation, Handle, OtherHandle, Move,
		            OtherMove
		           };

		// This is used to handle the creation of shapes by mouse events. Some
		// types of shapes require special treatment when constructed with
		// the mouse.
		enum CreationType {Normal, Arrow, Polygon, PolyLine, Marker};

		// Destructor, and constructor, taking the panelDisplay to which the
		// Annotator will attach itself to, and which mouse button to listen to.
		// <group>
		Annotations(PanelDisplay* panDisp,
		            const Display::KeySym& keysym = Display::K_Pointer_Button1,
		            const casacore::Bool useEH = true);
		virtual ~Annotations();
		// </group>

		// The event handlers.
		// <group>
		virtual void operator()(const PCPositionEvent& ev);
		virtual void operator()(const PCMotionEvent& ev);
		virtual void operator()(const PCRefreshEvent& ev);
		virtual void operator()(const WCRefreshEvent& ev);
		// </group>

		// Redraw all the shapes. Turns all handles off if noHandles is true
		virtual void draw(PixelCanvas* pc, const casacore::Bool noHandles = false);

		// Instantly creates a shape on screen, and hence requires
		// certain information in its record such as the center of the
		// new object on screen. Refresh == false can be used to surpress a
		// refresh.
		virtual void newShape(casacore::Record& settings, const casacore::Bool& refresh = true);

		// Set the key assigned to control shapes.
		virtual void setKey(const Display::KeySym& newKey);

		// Create a shape, based on the information contained in the record.
		// This function does not require the center to be specified, as once
		// createShape is called, the user can "drag" out the shape required.
		// If createShape is called, and the shape is no longer required,
		// the function 'cancelShapes' can be used to reset the state of the
		// annotator.
		virtual void createShape(casacore::Record& settings);

		// Disable event handling (with the exception of refresh calls). This stops
		// the user from moving shapes etc via the mouse motion / position event
		// handlers. These do not recognize nested calls (e.g. two disables,
		// followed by an enable will enable the handlers).
		// <group>
		virtual void enable();
		virtual void disable();
		// </group>

		// Cancels creation of any pending shapes, unselects all shapes and
		// turns off drawing of all handles.
		virtual void cancelShapes();

		// Delete the shape at the specified point in the collection of shapes
		// The current (selected) shape can be determined by using the
		// activeShape() function.
		virtual casacore::Bool deleteShape(const casacore::uInt& toDel);

		// These allow shapes to be "locked" together. Adding a locked shape
		// to the current shape means that whenever the current shape
		// moves, so too will the specified shape. removeLockedFromCurrent
		// can be used to remove this effect.
		// <group>
		virtual void removeLockedFromCurrent(const casacore::uInt& removeMe);
		virtual void addLockedToCurrent(const casacore::uInt& addMe);
		// </group>

		// Return the currently selected shapes' index
		virtual casacore::Int activeShape();

		// Return a brief summary of all current shapes
		virtual casacore::Record shapesSummary();

		// Set the options for the specified shape. The user selected shape
		// can be determined by call the 'activeShape' method. To see a list of
		// available options for different types of shapes / objects, see their
		// individual documentation.
		virtual void setShapeOptions(const casacore::uInt& whichShape,
		                             const casacore::Record& newOptions);

		// Returns a record containing the options for the specified shape.
		// The user selected shape can be determined by calling the 'activeShape'
		// method. To see a list of the options returned by shapes / objects,
		// see their individual documentation.
		virtual casacore::Record getShapeOptions(const casacore::uInt& whichShape);

		// Returns a record which contains a sub record (containing options) for
		// each shape. e.g. To obtain options for the first shape:
		// <srcblock>
		// casacore::Record all = myAnnotator.getAllOptions();
		// casacore::Record shape = all.subRecord(0);
		// </srcblock>
		virtual casacore::Record getAllOptions();

		// Deletes all existing shapes, and rebuilds new shapes based on records
		// in sub records of the supplied record. A record returned from
		// getAllOptions() can therefore be used to copy a set of shapes.
		virtual void setAllOptions(const casacore::Record& newSettings);

		// Called when a refresh is required. It copies back to front buffer,
		// then causes a refresh with reason = backCopiedToFront.
		virtual void update(PixelCanvas* pc);

		// Overwrite this to fire glish events.
		virtual void annotEvent(const casacore::String& /*event*/) {};

		// Functions for changing the co-ords system of a shape
		// <group>
		virtual casacore::Bool revertToPix(const casacore::Int& whichOne);
		virtual casacore::Bool revertToFrac(const casacore::Int& whichOne);
		virtual casacore::Bool lockToWC(const casacore::Int& whichOne);
		// </group>

		// Returns a list of shapes annotations knows about
		virtual casacore::Record availableShapes();


	protected:

	private:

		// My private variables...

		// The PanelDisplay I was given on construction and some stuff I got out
		// of it
		// <group>
		PanelDisplay* itsPanelDisplay;
		PixelCanvas* itsPC;
		// </group>

		// My copy of the list of Worldcanvases. I can check this against
		// itsPanelDisplay to see if anything has changed.
		std::list<WorldCanvas* > itsWCs;

		// casacore::List of shapes
		casacore::PtrBlock<DisplayShape*> itsShapes;

		// Some overall settings
		// <group>
		Display::KeySym itsKey;
		casacore::Bool itsUseEH;
		// </group>

		// Stuff to remember my current state:
		// <group>
		Annotations::State itsState;
		Annotations::CreationType itsCreation;
		casacore::Bool itsShapeMoved;
		casacore::Int itsActiveShape;
		casacore::Int itsActiveHandle;
		casacore::uInt itsX, itsY;
		casacore::Bool itsEnabled;
		casacore::Bool itsRefreshedYet;
		// </group>
		//


		// Internal functions
		// <group>
		casacore::Bool changeCoordSys(const casacore::Int shapeIndex, const casacore::String& type,
		                    const casacore::String& currentCoords, const casacore::String& toCoords);

		DisplayShape* constructShape(const casacore::Record& settings);
		void handleKeyUp();
		void handleCreation(const PCPositionEvent& ev);
		void handleNormalCreation(const casacore::Vector<casacore::Float>& createPix);
		void handleMarkerCreation(const casacore::Vector<casacore::Float>& createPix);
		void handlePolyLineCreation(const casacore::Vector<casacore::Float>& createPix);
		void handlePolygonCreation(const casacore::Vector<casacore::Float>& createPix);
		void handleArrowCreation(const casacore::Vector<casacore::Float>& createPix);
		void select(const casacore::Int i);
		casacore::Bool determineState(const PCPositionEvent& ev);
		void registerToWCs();
		casacore::Bool validateWCs();
		void changedWC();
		casacore::Bool validShape(const casacore::Record& shape);
		void polyLineToPolygon(const casacore::Int whichShape);
		// </group>
	};

// Ostream, mainly for debugging the state of the annotator.
	std::ostream& operator << (std::ostream& os, Annotations::State st);


} //# NAMESPACE CASA - END

#endif
