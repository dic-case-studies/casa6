//# AspMatrixCleaner.h: Minor Cycle for Asp deconvolution
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003
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
//#
//# $Id: AspMatrixCleaner.h Genie H. 2020-04-06 <mhsieh@nrao.edu $

#ifndef SYNTHESIS_ASPMATRIXCLEANER_H
#define SYNTHESIS_ASPMATRIXCLEANER_H

//# Includes
#include <scimath/Mathematics/FFTServer.h>
#include <synthesis/MeasurementEquations/MatrixCleaner.h>

namespace casa { //# NAMESPACE CASA - BEGIN

class AspMatrixCleaner : public MatrixCleaner
{
public:
  // Create a cleaner : default constructor
  AspMatrixCleaner();

  // The destructor does nothing special.
  ~AspMatrixCleaner();

  casacore::Bool setaspcontrol(const casacore::Int niter,
      const casacore::Float gain, const casacore::Quantity& aThreshold,
      const casacore::Quantity& fThreshold);

  // Calculate the convolutions of the dirty image
  // this is the same as MatrixCleaner::makeDirtyScales
  // except this doesn't check itsCleanType
  void makedirtyscales();

  // Clean an image.
  //return value gives you a hint of what's happening
  //  1 = converged
  //  0 = not converged but behaving normally
  // -1 = not converged and stopped on cleaning consecutive smallest scale
  // -2 = not converged and either large scale hit negative or diverging
  // -3 = clean is diverging rather than converging
  casacore::Int aspclean(casacore::Matrix<casacore::Float> & model, casacore::Bool doPlotProgress=false);

  // helper functions for ASP
  float getPsfGaussianWidth(casacore::ImageInterface<casacore::Float>& psf);

  void setInitScaleXfrs(/*const casacore::Array<casacore::Float> arrpsf,*/ const casacore::Float width);

  // calculate the convolutions of the psf with the initial scales
  void setInitScalePsfs();

  casacore::Bool setInitScaleMasks(const casacore::Array<casacore::Float> arrmask, const casacore::Float& maskThreshold = 0.9);

  void maxDirtyConvInitScales(float& strengthOptimum, int& optimumScale, casacore::IPosition& positionOptimum);

  bool isGoodAspen(casacore::Float amp, casacore::Float scale, casacore::IPosition center, casacore::Float threshold);

  std::vector<casacore::Float> getActiveSetAspen();

  //void defineAspScales(const casacore::Vector<casacore::Float>& scales);
  void defineAspScales(std::vector<casacore::Float>& scaleSizes);


//protected:
private:

  using MatrixCleaner::validatePsf;
  using MatrixCleaner::makeScale;
  using MatrixCleaner::findMaxAbs;
  using MatrixCleaner::findMaxAbsMask;
  using MatrixCleaner::makeBoxesSameSize;
  using MatrixCleaner::itsGain;
  using MatrixCleaner::itsMaxNiter;
  using MatrixCleaner::itsThreshold;
  using MatrixCleaner::itsMask;
  using MatrixCleaner::itsPositionPeakPsf;
  using MatrixCleaner::itsSmallScaleBias;
  using MatrixCleaner::itsScaleMasks;
  using MatrixCleaner::itsScaleXfrs;
  using MatrixCleaner::itsScalesValid;
  using MatrixCleaner::itsNscales;
  using MatrixCleaner::itsMaskThreshold;
  using MatrixCleaner::itsDirty;
  using MatrixCleaner::itsXfr;
  using MatrixCleaner::itsScaleSizes;
  using MatrixCleaner::itsScales;
  casacore::Block<casacore::Matrix<casacore::Float> > itsInitScales;
  casacore::Block<casacore::Matrix<casacore::Complex> > itsInitScaleXfrs;
  using MatrixCleaner::itsPsfConvScales;
  using MatrixCleaner::itsDirtyConvScales;
  casacore::Block<casacore::Matrix<casacore::Float> > itsDirtyConvInitScales;
  casacore::Block<casacore::Matrix<casacore::Float> > itsInitScaleMasks;
  casacore::Block<casacore::Matrix<casacore::Float> > itsPsfConvInitScales;

  using MatrixCleaner::itsIteration;
  using MatrixCleaner::itsStartingIter;
  using MatrixCleaner::itsFracThreshold;
  using MatrixCleaner::itsMaximumResidual;
  using MatrixCleaner::itsStrengthOptimum;
  using MatrixCleaner::itsTotalFlux;
  using MatrixCleaner::index;

  using MatrixCleaner::destroyScales;
  casacore::Bool destroyAspScales();
  casacore::Bool destroyInitScales();
  using MatrixCleaner::destroyMasks;
  casacore::Bool destroyInitMasks();

  using MatrixCleaner::itsIgnoreCenterBox;
  using MatrixCleaner::itsStopAtLargeScaleNegative;
  using MatrixCleaner::itsStopPointMode;
  using MatrixCleaner::itsDidStopPointMode;
  using MatrixCleaner::psfShape_p;
  using MatrixCleaner::noClean_p;

  // set to 0, 1.5*,5*,10*width for initial scales in Asp
  std::vector<casacore::Float> itsInitScaleSizes;
  std::vector<casacore::Float> itsAspScaleSizes; // permanent list for making model image
  std::vector<casacore::Float> itsAspAmplitude;
  std::vector<casacore::IPosition> itsAspCenter;
  casacore::Int itsNInitScales;
  double itsPrevLBFGSGrad; // flor gradient clipping if LBFGS gradient explodes
  unsigned int itsNumIterNoGoodAspen;
};

} //# NAMESPACE CASA - END

#endif
