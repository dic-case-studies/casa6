//# MultiPVTool.h: Base class for MultiWorldCanvas event-based rectangle tools
//# Copyright (C) 2013
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

#ifndef DISPLAY_MULTIPVTOOL_H_
#define DISPLAY_MULTIPVTOOL_H_

#include <list>
#include <casa/aips.h>
#include <display/DisplayEvents/RegionTool.h>
#include <display/DisplayEvents/DTVisible.h>
#include <display/region/RegionCreator.h>
#include <display/region/RegionSourceFactory.h>


namespace casa { //# NAMESPACE CASA - BEGIN

	namespace viewer {
		class QtRegionDock;
	}

// <summary>
// Base class for MultiWorldCanvas event-based rectangle tools
// </summary>
//
// <use visibility=export>
//
// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>
//
// <prerequisites>
//   <li> WCTool
// </prerequisites>
//
// <etymology>
// MultiPVTool stands for MultiWorldCanvas PVLine Tool
// </etymology>
//
// <synopsis>
// This class adds to its base WCTool to provide a tool for drawing,
// resizing and moving rectangles on a WorldCanvas.  While MultiPVTool
// is not abstract, it performs no useful function.  The programmer
// should derive from this class, and implement the functions
// doubleInside and doubleOutside, which are called when the user
// double-clicks the key or mouse button inside or outside an existing
// rectangle respectively.  It is up to the programmer to decide what
// double clicks inside and outside the rectangle correspond to,
// although it is recommended that a double click inside correspond to
// the main action of the tool, and a double click outside correspond
// to a secondary action of the tool, if indeed a secondary action
// exists.
//
// The rectangle is drawn by dragging the mouse from one corner to
// the diagonally opposite corner.  Once constructed, the rectangle
// can be resized by dragging its corners, or relocated by dragging
// inside the rectangle.  The rectangle is removed from the display
// when the Esc key is pressed.
// </synopsis>
//
// <example>
// </example>
//
// <motivation>
// Many activities on the WorldCanvas will be based on the user drawing
// a rectangle, and then proceeding to some action with that rectangle.
// A nice example is zooming.
// </motivation>

	class MultiPVTool : public RegionTool, public DTVisible, public viewer::RegionCreator {

	public:

		// Constructor
		MultiPVTool( viewer::RegionSourceFactory *rsf, PanelDisplay* pd,
		             Display::KeySym keysym = Display::K_Pointer_Button1, const casacore::Bool persistent = false );

		// Destructor
		virtual ~MultiPVTool();

		// Switch the tool off - this erases the rectangle, if any,
		// and calls the base class disable.
		virtual void disable();

		// reset to non-existent, non-active rectangle.
		// Refreshes if necessary to erase (unless skipRefresh==true  In
		// that case, the caller should do the refresh itself).
		// (Does not unregister from WCs or disable future event handling).
		virtual void reset(casacore::Bool skipRefresh=false);

		// Is a rectangle currently defined?
		virtual casacore::Bool rectangleDefined() {
			return itsPVLineExists;
		}

		viewer::RegionSource *getRegionSource( ) {
			return rfactory;
		}

		void checkPoint( WorldCanvas *wc, State &state );

		// called when the user (read GUI user) indicates that a region should be deleted...
		void revokeRegion( viewer::Region * );

		// returns a set which indicates regions this creator creates...
		const std::set<viewer::region::RegionTypes> &regionsCreated( ) const;

		bool create( viewer::region::RegionTypes /*region_type*/, WorldCanvas */*wc*/, const std::vector<std::pair<double,double> > &/*pts*/,
		             const std::string &/*label*/, viewer::region::TextPosition /*label_pos*/, const std::vector<int> &/*label_off*/,
		             const std::string &/*font*/, int /*font_size*/, int /*font_style*/, const std::string &/*font_color*/,
		             const std::string &/*line_color*/, viewer::region::LineStyle /*line_style*/, unsigned int /*line_width*/,
		             bool /*annotation*/, VOID */*region_specific_state*/ );

		RegionToolTypes type( ) const {
			return PVLINETOOL;
		}

	protected:

		// Functions called by the base class event handling operators--and
		// normally only those.  This is the input that controls the rectangle's
		// appearance and action.  When the rectangle is ready and double-click
		// is received, the doubleInside/Outside routine will be invoked.
		// <group>
		virtual void keyPressed(const WCPositionEvent &ev);
		virtual void keyReleased(const WCPositionEvent &ev);
		virtual void otherKeyPressed(const WCPositionEvent &ev);
		virtual void moved(const WCMotionEvent &ev, const viewer::region::region_list_type & /*selected_regions*/);
		// </group>

		// draw the rectangle (if any) on the object's currently active WC.
		// Only to be called by the base class refresh event handler.  Derived
		// objects should use refresh() if they need to redraw, but even that
		// is normally handled automatically.
		virtual void draw(const WCRefreshEvent&/*ev*/, const viewer::region::region_list_type & /*selected_regions*/);

		// Output callback functions--to be overridden in derived class.
		// Called when there is a double click inside/outside the rectangle
		// <group>
		virtual void doubleInside() { };
		virtual void doubleOutside() { };
		// </group>

		// Called when a rectangle is ready and not being
		// edited.  (Unused so far on the glish level (12/01)).
		virtual void rectangleReady() { };


		// Retrieve the rectangle coordinates, in screen pixels.
		// Anchor (if applicable) is (x1,y1).  Valid during the output callbacks;
		// to be used by them, as well as internally.
		virtual void get(casacore::Int &x1, casacore::Int &y1, casacore::Int &x2, casacore::Int &y2) const ;

		virtual bool checkType( viewer::region::RegionTypes t ) {
			return t == viewer::region::PVLineRegion;
		}
		virtual std::shared_ptr<viewer::PVLine> allocate_region( WorldCanvas *wc, double x1, double y1, double x2, double y2, VOID *region_specific_state ) const;


		viewer::RegionSource *rfactory;

	private:
		typedef std::list<std::shared_ptr<viewer::PVLine> > pvlinelist;

		void update_stats(const WCMotionEvent &ev);
		void start_new_rectangle( WorldCanvas *, int x, int y );

		// set the pixel coordinates of the rectangle
		virtual void set(const casacore::Int &x1, const casacore::Int &y1, const casacore::Int &x2, const casacore::Int &y2);

		// get only the anchor point
		virtual void get(casacore::Int &x1, casacore::Int &y1) const;

		// do we have a rectangle yet? (if true, itsCurrentWC, itsEmitted, P1, and P2 are valid)
		casacore::Bool itsPVLineExists;

		// was the button pressed in the rectangle (or, if none, in an active WC)
		// and not yet released/reset?
		casacore::Bool itsActive;
		// itsActive is being replaced by resizing_region
		std::shared_ptr<viewer::PVLine> resizing_region;
		std::shared_ptr<viewer::PVLine> creating_region;
		int resizing_region_handle;

		// (valid only if itsActive==true):
		// true = being moved     false = being resized
		casacore::Bool itsMoving;
		// itsMoving is being replaced by moving_regions
		pvlinelist moving_regions;
		double moving_linx_;
		double moving_liny_;

		// (valid only if itsPVLineExists==true)
		// Has doubleInside/Outside been called for this rectangle?  If so, a
		// key press outside the rectangle will start a new rectangle, as if
		// itsPVLineExists were false.
		// However, a key press inside the rectangle will reset
		// itsEmitted to false, allowing the rectangle to be reused
		// (possibly moved or resized, and emitted again).
		casacore::Bool itsEmitted;

		// (Linear) coordinates of the rectangle (invariant over zooms, but not
		// coordinate system changes.  To do: support the WorldCoordinateChange
		// refresh reason, and reset this tool when it occurs).
		casacore::Vector<casacore::Double> itsP1, itsP2;

		// storage of the handle (pixel) coordinates
		casacore::Vector<casacore::Int> itsHX, itsHY;

		// position that move started from
		casacore::Int itsBaseMoveX, itsBaseMoveY;

		// store the times of the last two presses here:
		casacore::Double itsLastPressTime, its2ndLastPressTime;

		pvlinelist rectangles;
		PanelDisplay *pd_;

	};

} //# NAMESPACE CASA - END

#endif
