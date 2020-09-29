//# LatticeAsContour.h: Class to display lattice objects as contoured images
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

#ifndef TRIALDISPLAY_LATTICEASCONTOUR_H
#define TRIALDISPLAY_LATTICEASCONTOUR_H

//# aips includes:
#include <casacore/casa/aips.h>

//# trial includes:
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/casa/Arrays/Array.h>

//# display library includes:
#include <display/DisplayDatas/LatticePADD.h>
#include <display/Display/DParameterRange.h>

namespace casacore{

	class Record;
}

namespace casa { //# NAMESPACE CASA - BEGIN

//# forwards:
	template <class T> class LatticePADMContour;

// <summary>Class to manage the drawing of contour maps of slices from AIPS++ Lattices</summary>
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
// The purpose of this class is to draw "contour" maps of data that
// is "lattice"-based.
// </etymology>
//
// <synopsis>
// This class should be used to display contour maps---ie. line
// drawings connecting data points of equal intensity (or level)---of
// two-dimensional slices of data extracted from AIPS++ Lattices or
// Arrays having two or more dimensions.  Thus, this is the class to
// use to display standard zeroth-order moment map contour overlays
// for use on Digitized Sky Survey images, for example.
//
// At construction, any axes in the data can be mapped to the X and Y
// axes of the display device (see the <linkto
// class="PrincipalAxesDD">PrincipalAxesDD</linkto> class).  For data
// with more than two dimensions, a third axis in the dataset can be
// selected for generating a sequence of maps along: this is known
// as the "movie" axis.  Animation (see the <linkto
// class="Animator">Animator</linkto> class) will cause different
// slices of the data to be selected from along this axis.  After
// construction, the axis settings of a LatticeAsContour object can be
// modified at a later stage.
//
// The LatticeAsContour object supports a number of options which can
// be set or retrieved using the setOptions and getOptions functions.
// These functions simply accept a casacore::Record, which can be converted from
// a GlishRecord: this is done in the <linkto
// class="GTkDisplayData">GTkDisplayData</linkto> class.  The options
// for the LatticeAsContour class are:
// <li> levels: a casacore::Vector<casacore::Float> (or casacore::Vector<casacore::Double>) of one or more
// elements, being the data
// values at which the contours are drawn.  Depending on type,
// the values in the vector are interpreted as absolute or
// fractions between the data minimum and maximum.
// <li> type: a casacore::String, one of "frac" or "abs", indicating whether the
// contour values in "levels" (multiplied by "scale") are fractions of the
// data range between the minimum and maximum, or are instead absolute
// contour levels in the native units of the data.
// <li> scale: a casacore::Float or casacore::Double which provides an additional scale
// factor to apply to the contour levels in "levels."
// <li> line: a positive Integer specifying the line width of
// contours.
// <li> dashneg: a Boolean value, which if true, will force contours
// at negative data values to be drawn in dashed line style.
// <li> dashpos: a Boolean value, which if true, will force contours
// at positive data values to be drawn in dashed line style.
// <li> color: a casacore::String which is the color with which to draw the
// contours.  A valid X Color is required.
//# <li> complexmode: this is a casacore::String, and is only relevant for
//# LatticeAsContour<casacore::Complex> or LatticeAsContour<casacore::DComplex>
//# instantantiations.  One of "phase", "real", "imaginary" or
//# "magnitude" is appropriate for this option, and indicates how
//# complex data values should be translated to real pixel values.
//#
//# LatticeAsContour is templated, and can be used to draw casacore::Complex or Real
//# Images or Arrays.  For casacore::Complex data, the default complexmode is
//# "magnitude."
//# </synopsis>
//
// <example>
// A LatticeAsContour object could be constructed and used as follows:
// <srcblock>
// casacore::PagedImage<casacore::Float> *pimage = new casacore::PagedImage<casacore::Float>(casacore::String("test.im"));
// DisplayData *dd;
// casacore::uInt ndim = pimage->ndim();
// if (ndim < 2) {
//   throw(casacore::AipsError(casacore::String("Image has less than two dimensions")));
// } else if (ndim == 2) {
//   dd = (DisplayData *)(new LatticeAsContour<casacore::Float>(pimage, 0, 1));
// } else {
//   casacore::IPosition fixedPos(ndim);
//   fixedPos = 0;
//   dd = (DisplayData *)(new LatticeAsContour<casacore::Float>(pimage, 0, 1, 2,
//                                                   fixedPos));
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
// <todo asof="yyyy/mm/dd">
// </todo>

	template <class T> class LatticeAsContour : public LatticePADisplayData<T> {

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
		LatticeAsContour(casacore::Array<T> *array, const casacore::uInt xAxis,
		                 const casacore::uInt yAxis, const casacore::uInt mAxis,
		                 const casacore::IPosition fixedPos);
		LatticeAsContour(casacore::Array<T> *array, const casacore::uInt xAxis,
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
		LatticeAsContour(std::shared_ptr<casacore::ImageInterface<T> > image, const casacore::uInt xAxis, const casacore::uInt yAxis, const casacore::uInt mAxis, const casacore::IPosition fixedPos, viewer::StatusSink *sink=0 );
		LatticeAsContour(std::shared_ptr<casacore::ImageInterface<T> > image, const casacore::uInt xAxis,
		                 const casacore::uInt yAxis);
		// </group>

		// Destructor
		virtual ~LatticeAsContour();

		// Create the various elements in the sequence of displayable
		// maps.  This is called upon construction as well as whenever
		// the display and/or movie axes are changed via a call to
		// PrincipalAxesDD::setAxes.
		//virtual void setupElements(casacore::IPosition fixedPos = casacore::IPosition(2));
		virtual void setupElements();

		// install the default options for display
		virtual void setDefaultOptions();

		// Apply the options stored in the provided casacore::Record to the
		// LatticeAsContour object.  If the return value is true, then
		// some options have changed, and a refresh is needed to update
		// the display.
		virtual casacore::Bool setOptions(casacore::Record &rec, casacore::Record &recOut);

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
			return casacore::String("LatticeAsContour");
		}

		// Actual selected contour levels
		virtual casacore::Vector<casacore::Float> levels();

		// Actual levels, in casacore::String format.  If precision is unspecified,
		// one that is low enough not to be cluttered is used.
		virtual casacore::String levelString(casacore::Int prec=0);

		// Adds contour level information to the standard position tracking
		// string from PADD.
		virtual casacore::String showPosition(const casacore::Vector<casacore::Double> &wld, const casacore::Bool &abs,
		                            const casacore::Bool &dsp);

		using LatticePADisplayData<T>::dataUnit;


	protected:

		// Construct user option DisplayParameters (for min/max contour.)
		// (To be used by constructors only.)
		virtual void constructParameters_();

		// Set standard limits/values for contour sliders.  If recOut is provided,
		// they will be set onto it in a manner suitable for updating gui via
		// setOptions.
		virtual void setStdContourLimits_(casacore::Record* recOut=0);


	private:

		// Relative contour levels.
		casacore::Vector<casacore::Float> itsLevels;

		// Min and Max actual contours. linearly scaled from itsLevels to fit these.
		// <group>
		DParameterRange<casacore::Float> *itsBaseContour;
		DParameterRange<casacore::Float> *itsUnitContour;
		// </group>

		casacore::Float itsLine;
		casacore::Bool itsDashNeg;
		casacore::Bool itsDashPos;
		casacore::String itsColor;

		friend class LatticePADMContour<T>;

		//# Make parent members known.
	public:
		using LatticePADisplayData<T>::nelements;
		using LatticePADisplayData<T>::nPixelAxes;
		using LatticePADisplayData<T>::fixedPosition;
		using LatticePADisplayData<T>::displayAxes;
		using LatticePADisplayData<T>::dataShape;
		using LatticePADisplayData<T>::getMinAndMax;
		using LatticePADisplayData<T>::readOptionRecord;
		using LatticePADisplayData<T>::getDataMin;
		using LatticePADisplayData<T>::getDataMax;

	protected:
		using LatticePADisplayData<T>::setNumImages;
		using LatticePADisplayData<T>::DDelement;
		using LatticePADisplayData<T>::datamax;
	};


} //# NAMESPACE CASA - END

#ifndef AIPS_NO_TEMPLATE_SRC
#include <display/DisplayDatas/LatticeAsContour.tcc>
#endif //# AIPS_NO_TEMPLATE_SRC
#endif
