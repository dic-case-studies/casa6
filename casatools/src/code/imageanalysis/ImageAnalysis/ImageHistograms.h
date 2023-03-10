//# ImageHistograms.h: generate histograms from an image
//# Copyright (C) 1996,1997,1999,2000,2001
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
//# $Id: ImageHistograms.h 20229 2008-01-29 15:19:06Z gervandiepen $

#ifndef IMAGES_IMAGEHISTOGRAMS_H
#define IMAGES_IMAGEHISTOGRAMS_H


//# Includes
#include <casacore/casa/aips.h>
#include <casacore/lattices/LatticeMath/LatticeHistograms.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/BasicSL/String.h>


namespace casacore{

template <class T> class ImageInterface;
class IPosition;
}

namespace casa {

//# Forward Declarations


// <summary>
// Displays histograms of regions from an image.
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>

// <prerequisite>
//   <li> <linkto class=casacore::ImageInterface>ImageInterface</linkto>
//   <li> <linkto class=casacore::LatticeHistograms>LatticeHistograms</linkto>
// </prerequisite>

// <etymology>
// This is a class designed to display histograms from images
// </etymology>

// <synopsis>
// This class enable you to display and/or retrieve histograms evaluated over 
// specified regions from an image.  The dimension of the region is arbitrary, but 
// the size of each dimension is always the size of the corresponding image axis.
// The histograms are displayed as a function of location of the axes not
// used to evaluate the histograms over.  The axes which you evaluate the histograms 
// over are called the cursor axes, the others are called the display axwes.
//
// For example, consider an image cube (call the axes xyz or [0,1,2]).  You could 
// display histograms from xy planes (cursor axes [0,1]) as a function of z (display
// axes [2]).   Or  you could retrieve histograms from the z axis (cursor axes [2])
// for each [x,y] location (display axes [0,1]).
//
// The hard work is done by casacore::LatticeHistograms which this class (clumsily) inherits.
// It generates a "storage lattice" into which it writes the histograms.
// It is from this storage image that the plotting and retrieval
// arrays are drawn.  The storage image is either in core or on disk
// depending upon its size (if > 10% of memory given by .aipsrc system.resources.memory
// then it goes into a disk-based casacore::PagedArray).  If on disk,  the
// storage image is deleted when the <src>ImageHistograms</src> 
// object destructs.    
//
// See casacore::LatticeHistograms for most of the useful public interface.  ImageHistograms
// exists only so that it can write some world coordinate information to the plots
// and logger.
//
// <note role=tip>
// Note that for complex images, real and imaginary are treated independently.
// They are binned and plotted separately.
// </note>
//
// <note role=tip>
// If you ignore return error statuses from the functions that set the
// state of the class, the internal status of the class is set to bad.
// This means it will just  keep on returning error conditions until you
// explicitly recover the situation.   A message describing the last   
// error condition can be recovered with function errorMessage.

// </note>
// </synopsis>

// <example>
// <srcBlock>
//// Construct casacore::PagedImage from file name
//
//      casacore::PagedImage<casacore::Float> inImage(inName);
//   
//// Construct histogram object
//      
//      casacore::LogOrigin or("myClass", "myFunction(...)", WHERE);
//      casacore::LogIO os(or);
//      ImageHistograms<casacore::Float> histo(inImage, os);
//      
//// Set cursor axes to see statistics of yz planes (0 relative)
//
//      casacore::Vector<casacore::Int> cursorAxes(2)
//      cursorAxes(0) = 1;
//      cursorAxes(1) = 2;
//      if (!histo.setAxes(cursorAxes)) return 1;
//
//// Set to list and plot mean, sigma and rms
//
//      if (!histo.setList(true)) return 1;
//      casacore::String device = "/xs";
//      casacore::Vector<casacore::Int> nxy(2);
//      nxy(0) = 3;
//      nxy(1) = 3;
//      if (!histo.setPlotting(device, nxy)) return 1;
// 
//// Now activate actual listing and plotting
// 
//      if (!histo.display ()) return 1;
//
//// Retrieve histograms into array
//
//      casacore::Array<casacore::Float> values, counts;
//      if (!histo.getHistograms(values, counts)) return 1;
//
// </srcBlock>
// In this example, a <src>casacore::PagedImage</src> is constructed.  We set the cursor axes 
// to be the y and z axes so we make a histogram of each yz plane as a function 
// of x location on the  device "/xs" (no longer supported) with 9 subplots per page.
// After the plotting we also retrieve the histograms into an array.
// </example>

// <motivation>
// The generation of histograms from an image is a basic and necessary capability.
// </motivation>
//
// <todo asof="2000/04/04">
//   <li> Make ascii listing of histograms as well as plots if desired
// </todo>
//


template <class T> class ImageHistograms : public casacore::LatticeHistograms<T> {
public:

// Constructor takes the image and a <src>casacore::LogIO</src> object for logging.
// You can also specify whether you want to see progress meters or not.
// You can force the storage image to be disk based, otherwise
// the decision for core or disk is taken for you.
   ImageHistograms(const casacore::ImageInterface<T>& image, 
                   casacore::LogIO& os,
                   casacore::Bool showProgress=true,
                   casacore::Bool forceDisk=false);

// Constructor takes the image only. In the absence of a logger you get no messages.
// This includes error messages and potential listing of statistics.
// You can specify whether you want to see progress meters or not.
// You can force the storage image to be disk based, otherwise
// the decision for core or disk is taken for you.
   ImageHistograms(const casacore::ImageInterface<T>& image, 
                   casacore::Bool showProgress=true,
                   casacore::Bool forceDisk=false);

// Copy constructor (copy semantics)
   ImageHistograms(const ImageHistograms<T> &other);

// Destructor
  virtual ~ImageHistograms ();

// Assignment operator (copy semantics)
   ImageHistograms<T> &operator=(const ImageHistograms<T> &other);

// Set a new image.  A return value of <src>false</src> indicates the 
// image had an invalid type or that the internal status of the class is bad.
   casacore::Bool setNewImage (const casacore::ImageInterface<T>& image);

private:
   casacore::LogIO os_p;
   const casacore::ImageInterface<T>* pInImage_p;

   // Make a string with pixel and world coordinates of display axes
   virtual casacore::String writeCoordinates(const casacore::IPosition& histPos) const;

  //# Make members of parent class known.
protected:
  using casacore::LatticeHistograms<T>::locHistInLattice;
  using casacore::LatticeHistograms<T>::error_p;
  using casacore::LatticeHistograms<T>::goodParameterStatus_p;
  using casacore::LatticeHistograms<T>::displayAxes_p;
  using casacore::LatticeHistograms<T>::cursorAxes_p;
};

}

#ifndef CASACORE_NO_AUTO_TEMPLATES
#include <imageanalysis/ImageAnalysis/ImageHistograms.tcc>
#endif //# CASACORE_NO_AUTO_TEMPLATES
#endif
