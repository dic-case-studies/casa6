//# LatticeAsMarker.h: Class to display a lattice as markers
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002
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

#ifndef TRIALDISPLAY_LATTICEASMARKER_H
#define TRIALDISPLAY_LATTICEASMARKER_H

//# aips includes:
#include <casacore/casa/aips.h>
#include <casacore/casa/Quanta/Unit.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>

//# trial includes:
#include <images/Images/ImageInterface.h>

//# display library includes:
#include <display/DisplayDatas/LatticePADD.h>

namespace casa { //# NAMESPACE CASA - BEGIN

//# forwards:
	template <class T> class LatticePADMMarker;

// <summary>Class to manage the drawing of vector maps of slices from AIPS++ Lattices</summary>
//
// <use visibility=export>
//
// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>
//
// <prerequisite>
// <li> <linkto class="LatticePADisplayData">LatticePADisplayData</linkto>
// <li> <linkto class="DisplayData">DisplayData</linkto>
// </prerequisite>
//
// <etymology>
// From casacore::Lattice and Marker (shape)
// </etymology>
//
// <synopsis>
// The purpose of this class is to draw pixels as markers such as squares
// where the size of the marker reflects the pixel value.
// Presently only squares are available.  Positive values
// get solid squares, negative values open squares.
//
// At construction, any axes in the data can be mapped to the X and Y
// axes of the display device (see the <linkto
// class="PrincipalAxesDD">PrincipalAxesDD</linkto> class).  For data
// with more than two dimensions, a third axis in the dataset can be
// selected for generating a sequence of maps along: this is known
// as the "movie" axis.  Animation (see the <linkto
// class="Animator">Animator</linkto> class) will cause different
// slices of the data to be selected from along this axis.  After
// construction, the axis settings of a LatticeAsMarker object can be
// modified at a later stage.
//
//# </synopsis>
//
// <example>
// A LatticeAsMarker object could be constructed and used as follows:
// <srcblock>
// casacore::PagedImage<casacore::Complex> *pimage = new casacore::PagedImage<casacore::Complex>(casacore::String("test.im"));
// DisplayData *dd;
// casacore::uInt ndim = pimage->ndim();
// if (ndim < 2) {
//   throw(casacore::AipsError(casacore::String("Image has less than two dimensions")));
// } else if (ndim == 2) {
//   dd = (DisplayData *)(new LatticeAsMarker<casacore::Complex>(pimage, 0, 1));
// } else {
//   casacore::IPosition fixedPos(ndim);
//   fixedPos = 0;
//   dd = (DisplayData *)(new LatticeAsMarker<casacore::Complex>(pimage, 0, 1, 2,fixedPos));
// }
// // wcHolder is an existing WorldCanvasHolder *...
// wcHolder->addDisplayData(ddata);
// wcHolder->refresh();
// </srcblock>
// </example>
//
// <motivation>
// Displaying 2-dimensional slices of a lattice-based data volume is
// a standard display requirement for astronomical data visualization
// and presentation.
// </motivation>
//
// <templating arg=T>
// </templating>
//
// <thrown>
// </thrown>
//
// <todo asof="2001/11/15">
// We need a comprehensive set of marker shapes defined by
// an enum and implemented on the pixel canvas with native
// DL functions.
// </todo>

	template <class T> class LatticeAsMarker : public LatticePADisplayData<T> {

	public:
		// casacore::Array-based constructors: >2d and 2d.  xAxis and yAxis specify
		// which axis in the array (0-based) should be mapped to X and Y
		// on the display device: ie. 2-d slices of the data to be displayed
		// have these axes.  mAxis specifies the "movie" axis, which is the axis
		// along which different slices are taken.  fixedPos is an IPosition
		// having the same length as the number of dimensions in the array,
		// and indicate the fixed axis values for axes in the data that are
		// not specified as xAxis, yAxis or mAxis.
		// <group>
		LatticeAsMarker(casacore::Array<T>* array, const casacore::uInt xAxis,
		                const casacore::uInt yAxis, const casacore::uInt mAxis,
		                const casacore::IPosition fixedPos);
		LatticeAsMarker(casacore::Array<T>* array, const casacore::uInt xAxis,
		                const casacore::uInt yAxis);
		// </group>

		// Image-based constructors: >2d and 2d.  xAxis and yAxis specify
		// which axis in the image (0-based) should be mapped to X and Y
		// on the display device: ie. 2-d slices of the data to be displayed
		// have these axes.  mAxis specifies the "movie" axis, which is the axis
		// along which different slices are taken.  fixedPos is an IPosition
		// having the same length as the number of dimensions in the image,
		// and indicate the fixed axis values for axes in the data that are
		// not specified as xAxis, yAxis or mAxis.
		// <group>
		LatticeAsMarker(std::shared_ptr<casacore::ImageInterface<T> > image, const casacore::uInt xAxis,
		                const casacore::uInt yAxis, const casacore::uInt mAxis,
		                const casacore::IPosition fixedPos);
		LatticeAsMarker(std::shared_ptr<casacore::ImageInterface<T> > image, const casacore::uInt xAxis,
		                const casacore::uInt yAxis);
		// </group>

		// Destructor
		virtual ~LatticeAsMarker();

		// Create the various elements in the sequence of displayable
		// maps.  This is called upon construction as well as whenever
		// the display and/or movie axes are changed via a call to
		// PrincipalAxesDD::setAxes.
		//virtual void setupElements(casacore::IPosition fixedPos = casacore::IPosition(2));
		virtual void setupElements();

		// install the default options for display
		virtual void setDefaultOptions();

		// Apply the options stored in the provided casacore::Record to the
		// LatticeAsMarker object.  If the return value is true, then
		// some options have changed, and a refresh is needed to update
		// the display.
		virtual casacore::Bool setOptions(casacore::Record& rec, casacore::Record& recOut);

		// Retrieve the currently set options, and their types, default
		// values, and any help text associated with each parameter.  This
		// information can be used to generate form-type graphical user
		// interfaces or command-line interfaces to set the options with
		// prompts.
		virtual casacore::Record getOptions( bool scrub=false ) const;

		// Return the DisplayData type; used by the WorldCanvasHolder to
		// determine the order of drawing.
		virtual Display::DisplayDataType classType() {
			return Display::Vector;
		}

		// class name
		virtual casacore::String className() const {
			return casacore::String("LatticeAsMarker");
		}


	private:
		casacore::Float itsScale;          // Amplitude scale factor
		casacore::Float itsLineWidth;      // Line width of vectors
		casacore::Int itsIncX;             // Increment in X (pixels)
		casacore::Int itsIncY;             // Increment in Y (pixels)
		casacore::String itsColor;         // Color of vectors
		casacore::String itsMarker;        // Marker shape (should use an enum when we have some more)
//
		friend class LatticePADMMarker<T>;

		//# Make parent members known.
	public:
		using LatticePADisplayData<T>::nelements;
		using LatticePADisplayData<T>::nPixelAxes;
		using LatticePADisplayData<T>::fixedPosition;
		using LatticePADisplayData<T>::displayAxes;
		using LatticePADisplayData<T>::dataShape;
		using LatticePADisplayData<T>::getMinAndMax;
		using LatticePADisplayData<T>::readOptionRecord;
	protected:
		using LatticePADisplayData<T>::setNumImages;
		using LatticePADisplayData<T>::DDelement;
	};


} //# NAMESPACE CASA - END

#ifndef AIPS_NO_TEMPLATE_SRC
#include <display/DisplayDatas/LatticeAsMarker.tcc>
#endif //# AIPS_NO_TEMPLATE_SRC
#endif
