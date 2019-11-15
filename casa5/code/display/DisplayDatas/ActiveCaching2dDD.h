//# ActiveCaching2dDD.h: active 2d implementation of a CachingDisplayData
//# Copyright (C) 1999,2000,2001,2002,2003,2004
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

#ifndef TRIALDISPLAY_ACTIVECACHING2DDD_H
#define TRIALDISPLAY_ACTIVECACHING2DDD_H

#include <casa/aips.h>
#include <display/Display/DisplayCoordinateSystem.h>
#include <display/DisplayDatas/CachingDisplayData.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// Class providing active behaviour for 2d CachingDisplayDatas.
// </summary>

// <synopsis>
// An "active" DisplayData is one which is able to negotiate
// coordinates on a WorldCanvas, and subsequently provide coordinate
// conversion facilities.  This particular implementation is
// of a two-dimensional (2d) active DisplayData.
// </synopsis>

	class ActiveCaching2dDD : public CachingDisplayData {

	public:

		// Constructor.  Caller must provide a DisplayCoordinateSystem, and the
		// pixel range (via <src>pixelblc</src> and <src>pixeltrc</src>)
		// over which the DisplayData is expected to function.  The primary
		// use of <src>pixelblc</src> and <src>pixeltrc</src> is to provide
		// limits for the unzoomed WorldCanvas.  <src>coordsys</src> must
		// have only two world and pixel axes, otherwise an exception is
		// thrown, and likewise, <src>pixblc</src> and <src>pixtrc</src>
		// must both be of length two.
		ActiveCaching2dDD(const DisplayCoordinateSystem &coordsys,
		                  const casacore::Vector<casacore::Double> &pixblc,
		                  const casacore::Vector<casacore::Double> &pixtrc);

		// Destructor.
		virtual ~ActiveCaching2dDD();

		// casacore::Coordinate transformation handlers, called by WorldCanvasHolder.
		// <group>
		virtual casacore::Bool linToWorld(casacore::Vector<casacore::Double> &world, const casacore::Vector<casacore::Double> &lin);
		virtual casacore::Bool worldToLin(casacore::Vector<casacore::Double> &lin, const casacore::Vector<casacore::Double> &world);
		// </group>

		// casacore::Format a string containing coordinate information at the given
		// world coordinate.
		virtual casacore::String showPosition(const casacore::Vector<casacore::Double> &world,
		                            const casacore::Bool &displayAxesOnly = false);

		// World axis information suppliers.
		// <group>
		virtual casacore::Vector<casacore::String> worldAxisNames() const;
		virtual casacore::Vector<casacore::String> worldAxisUnits() const;
		// </group>

		// Return the number of display elements (ie. drawable images) in
		// this DisplayData.  Both return 1 because this is defined to be a
		// two-dimensional DisplayData, and can only have one view.  That
		// is, there is no third axis to iterate over for multiple views.
		// <group>
		virtual casacore::uInt nelements(const WorldCanvasHolder &/*wcHolder*/) const {
			return nelements();
		}
		virtual casacore::uInt nelements() const {
			return 1;
		}
		// </group>

		// Install the default options for this DisplayData.
		virtual void setDefaultOptions();

		// Apply options stored in <src>rec</src> to the DisplayData.  A
		// return value of <src>true</src> means a refresh is needed.
		// <src>recOut</src> contains any fields which were implicitly
		// changed as a result of the call to this function.
		virtual casacore::Bool setOptions(casacore::Record &rec, casacore::Record &recOut);

		// Retrieve the current and default options and parameter types.
		virtual casacore::Record getOptions( bool scrub=false ) const;

		// Negotiatiate WorldCanvas linear coordinate system when asked to
		// do so by the WorldCanvasHolder.
		virtual casacore::Bool sizeControl(WorldCanvasHolder &wcHolder,
		                         AttributeBuffer &holderBuf);

		// Determine whether DD can draw on the current coordinate system
		// of the given WC[H].
		virtual casacore::Bool conformsToCS(const WorldCanvas& wc);

		// Tidy up the elements of this DisplayData.
		virtual void cleanup();

	protected:

		// Return the current options of this DisplayData as an
		// AttributeBuffer.  The caller must delete the returned buffer.
		virtual AttributeBuffer optionsAsAttributes();

		// (Required) default constructor.
		ActiveCaching2dDD();

		// (Required) copy constructor.
		ActiveCaching2dDD(const ActiveCaching2dDD &other);

		// (Required) copy assignment.
		void operator=(const ActiveCaching2dDD &other);

		// Install a DisplayCoordinateSystem and limits.
		void setCoordinateSystem(const DisplayCoordinateSystem &coordsys,
		                         const casacore::Vector<casacore::Double> &pixblc,
		                         const casacore::Vector<casacore::Double> &pixtrc);

		// Identify a specified world axis by its casacore::Coordinate type, number,
		// and axis in that coordinate.
		void identifyWorldAxis(casacore::Coordinate::Type &type, casacore::Int &coordinate,
		                       casacore::Int &axisincoord, const casacore::uInt worldaxisnum);

		// Generate a casacore::String code for the specified world axis, describing
		// the coordinate type, number and axis in that coordinate.
		casacore::String codeWorldAxis(const casacore::uInt worldaxisnum);

		// Return the world axis increments.
		casacore::Vector<casacore::Double> worldAxisIncrements() const;

		// Return whether tracking is currently world (T) or pixel (F) coordinates.
		virtual casacore::Bool worldCoordTracking() const {
			return itsShowWorldCoordinate;
		}

		// Store the 2d DisplayCoordinateSystem here.
		DisplayCoordinateSystem itsCoordinateSystem;

		// Store the 2d pixel corners here.
		casacore::Vector<casacore::Double> itsPixelBlc, itsPixelTrc;

	private:

		// Option: aspect ratio for pixels.
		casacore::String itsOptionsAspect;

		// Position tracking
		casacore::String itsSpectralUnit;
		casacore::String itsVelocityType;
		casacore::Bool itsAbsolute;
		casacore::Bool itsShowWorldCoordinate;
		casacore::Bool itsFractionalPixels;

		// Set Spectral formatting
		void setSpectralFormatting (DisplayCoordinateSystem& cSys,
		                            const casacore::String& velTypeString,
		                            const casacore::String& unitString);


	};


} //# NAMESPACE CASA - END

#endif
