//# WCCSAxisLabeller.h: labelling axes using a DisplayCoordinateSystem on a WC
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

#ifndef TRIALDISPLAY_WCCSAXISLABELLER_H
#define TRIALDISPLAY_WCCSAXISLABELLER_H

#include <casa/aips.h>
#include <display/Display/DisplayCoordinateSystem.h>
#include <display/DisplayCanvas/WCAxisLabeller.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// Base class for WorldCanvas axis labelling using a DisplayCoordinateSystem.
// </summary>
//
// <synopsis>
// This (base) class adds to the interface of WCAxisLabeller functions
// which support the use/provision of a DisplayCoordinateSystem to assist with
// axis labelling.
// </synopsis>

	class WCCSAxisLabeller : public WCAxisLabeller {

	public:
		enum SpecAxisType { // taken from the casacore::FITS spectral coordinate type codes
		    FREQ,
		    VELO,
		    WAVE,
		    AWAV,
		};

		// Constructor
		WCCSAxisLabeller();

		// Destructor
		virtual ~WCCSAxisLabeller();

		// Install a DisplayCoordinateSystem.
		//#dk note: See casacore::Bool useWCCS, below; in some cases the WorldCanvas's
		//#   own CS is now used to draw labels, although the CS set here is
		//#   still needed for certain things (default user interface, e.g.)
		//#   when the WC is unknown.  (10/07)
		virtual void setCoordinateSystem(const DisplayCoordinateSystem& coordsys);

		// Get the DisplayCoordinateSystem.
		virtual DisplayCoordinateSystem coordinateSystem() const {
			return itsCoordinateSystem;
		}

		// Has a CS been set?
		casacore::Bool hasCoordinateSystem() const {
			return itsHasCoordinateSystem;
		};

		// Setting this true allows the labeller to use the WorldCanvas's
		// own CS to draw labels (although itsCoordinateSystem is still
		// needed for certain things at present).  Default: false.
		//#dk (See WCCSNLAxisLabeller::draw() and usage in PADD::setupElements()).
		//#
		//# (Kludge upon kludge, I know... all stemming from the original design
		//# flaw: not realizing that WC needs its _own CS_ (related to but _not
		//# the same_ as the Images'), and then trying to tack one on later...).
		casacore::Bool useWCCS;

		// install the default options for this labeller.
		virtual void setDefaultOptions();

		// apply options stored in rec to the labeller; return value
		// true means a refresh is needed.  Any fields added to the
		// updatedOptions argument are options which have changed in
		// some way due to the setting of other options - ie. they
		// are context sensitive.
		virtual casacore::Bool setOptions(const casacore::Record &rec, casacore::Record &updatedOptions);

		// retrieve the current and default options and parameter types.
		virtual casacore::Record getOptions() const;

		// return the X and Y label text - over-ridden from base class
		// <group>
		//# virtual casacore::String xAxisText(WorldCanvas* wc=0) const;
		//# virtual casacore::String yAxisText(WorldCanvas* wc=0) const;
		//# (Compiler whines unless you do it this way instead... grr...).
		virtual casacore::String xAxisText(WorldCanvas* wc) const;
		virtual casacore::String yAxisText(WorldCanvas* wc) const;
		virtual casacore::String xAxisText() const {
			return xAxisText(0);
		}
		virtual casacore::String yAxisText() const {
			return yAxisText(0);
		}
		// </group>

		virtual casacore::String zLabelType() const {
			return itsZLabelType;
		};

		virtual casacore::String zLabelPos() const {
			return itsZLabelPos;
		};

		virtual void setZIndex(casacore::Int zindex) {
			itsZIndex = zindex;
		};

		// DD 'Absolute Pixel Coordinates', e.g. channel numbers, are internally
		// 0-based (they begin numbering at 0), but 'Absolute Pixel coordinates'
		// have traditionally been displayed as 1-based in the glish viewer.
		// uiBase_, and related methods uiBase() and setUIBase(), allow newer
		// (python/Qt-based) code to cause such labelling to be produced with
		// 0-based values instead.  Unless setUIBase(0) is called, the
		// traditional 1-based labelling behavior is retained by default.
		//
		// If you are using 0-basing for 'Absolute Pixel casacore::Coordinate' labelling,
		// you should call setUIBase(0), before using draw().
		// <group>
		virtual casacore::Int uiBase() const {
			return uiBase_;
		}

		virtual void setUIBase(casacore::Int uibase) {
			if(uibase==0 || uibase==1) uiBase_ = uibase;
		}
		// </group>


		const casacore::String &spectralunitStr( ) const {
			return itsSpectralUnit;
		}

        casacore::Int spectralprec( ) const {
            return itsSpectralPrecision;
        }

	protected:
		casacore::Bool itsAbsolute;
		casacore::Bool itsWorldAxisLabels;
		mutable WCCSAxisLabeller::SpecAxisType itsSpecAxisType;
		casacore::Int itsZIndex;

		// Set spectral state onto given CS
		void setSpectralState(DisplayCoordinateSystem& cs) const;

		// Set direction state onto given CS
		void setDirectionState(DisplayCoordinateSystem& cs) const;

	private:

		DisplayCoordinateSystem itsCoordinateSystem;
		casacore::Bool itsHasCoordinateSystem;
		casacore::Int    itsSpectralPrecision;
		casacore::String itsSpectralUnit;
		casacore::String itsSpectralQuantity;
		casacore::String itsSpectralTypeUnit;
		casacore::String itsDirectionUnit;
		casacore::String itsDirectionSystem;
		casacore::String itsFrequencySystem;
		casacore::String itsZLabelType;
		casacore::String itsZLabelPos;
		casacore::String itsRestValue;      // rest frequency or wavelength (value and unit)
		casacore::String itsRestUnit;       // unit for rest frequency or wavelength

		casacore::Int uiBase_;		// (initialized to 1; see uiBase(), above).

		// Generate axis text for specified axis
		//# casacore::String axisText(casacore::Int worldAxis, WorldCanvas* wc=0) const;
		//# (Compiler whines unless you do it this way instead... grr...).
		casacore::String axisText(casacore::Int worldAxis, WorldCanvas* wc) const;
		casacore::String axisText(casacore::Int worldAxis) const {
			return axisText(worldAxis, 0);
		}

		// Set new spectral state in itsCoordinateSystem
		void setSpectralState() {
			setSpectralState(itsCoordinateSystem);
		}

		// Set new direction state in itsCoordinateSystem
		void setDirectionState() {
			setDirectionState(itsCoordinateSystem);
		}

		// Set absolute/relative state in itsCoordinateSystem
		void setAbsRelState();

		// "optical velocity [m/s]" --> "optical velocity" and "m/s"
		void distributeTypeUnit() ;

		static const casacore::String FRAME_REST;
	};


} //# NAMESPACE CASA - END

#endif
