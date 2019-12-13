//# DSWorldArrow.h : Implementation of a world coords DSArrow
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
#ifndef TRIALDISPLAY_DSWORLDARROW_H
#define TRIALDISPLAY_DSWORLDARROW_H

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>

#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/QuantumHolder.h>
#include <display/Display/PanelDisplay.h>

#include <display/DisplayShapes/DSArrow.h>
#include <display/DisplayShapes/DisplayShapeWithCoords.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// Implementation of a world coords DSArrow
// </summary>
	class PanelDisplay;
	class DSPixelArrow;
	class DSScreenArrow;

	class DSWorldArrow : public DSArrow , public DisplayShapeWithCoords {

	public:

		DSWorldArrow();
		DSWorldArrow(const casacore::Record& settings, PanelDisplay* pd);
		DSWorldArrow(DSScreenArrow& other, PanelDisplay* pd);
		DSWorldArrow(DSPixelArrow& other, PanelDisplay* pd);

		virtual ~DSWorldArrow();

		virtual void recalculateScreenPosition();

		// So we can update our WCs
		// <group>
		virtual void move(const casacore::Float& dX, const casacore::Float& dY);
		virtual void setCenter(const casacore::Float& xPos, const casacore::Float& yPos);
		virtual void rotate(const casacore::Float& angle);
		virtual void scale(const casacore::Float& scaleFactor);
		virtual void setStartPoint(const casacore::Vector<casacore::Float>& startPoint);
		virtual void setEndPoint(const casacore::Vector<casacore::Float>& endPoint);
		virtual void changePoint(const casacore::Vector<casacore::Float>&pos, const casacore::Int n);
		virtual void changePoint(const casacore::Vector<casacore::Float>& pos);
		virtual void draw(PixelCanvas* pc);
		// </group>

		virtual casacore::Bool setOptions(const casacore::Record& settings);
		virtual casacore::Record getOptions();



		virtual casacore::Record getRawOptions() {
			return DSArrow::getOptions();
		}

		virtual PanelDisplay* panelDisplay()  {
			return itsPD;
		}
	private:

		// Based on a WC option record
		//  virtual WorldCanvas* chooseWC(const casacore::Record& settings, PanelDisplay* pd);

		// Based on a pixel center
		virtual WorldCanvas* chooseWC(const casacore::Float& startXPos,
		                              const casacore::Float& startYPos,
		                              const casacore::Float& endXPos, const casacore::Float& endYPos,
		                              PanelDisplay* pd);

		// The paneldisplay from which I may choose an appropriate WC
		PanelDisplay* itsPD;

		// The WC of my choosing
		WorldCanvas* itsWC;

		// The center of the marker in world co-ords.
		casacore::Vector<casacore::Quantum<casacore::Double> > itsWorldStart;
		casacore::Vector<casacore::Quantum<casacore::Double> > itsWorldEnd;

		void updateWCoords();
	};

} //# NAMESPACE CASA - END

#endif
