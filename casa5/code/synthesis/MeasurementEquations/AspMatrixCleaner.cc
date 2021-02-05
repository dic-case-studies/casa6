//# Copyright (C) 1997-2010
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
//# $Id:  $AspMatrixCleaner.cc

// Same include list as in MatrixCleaner.cc
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
//#include <casa/Arrays/ArrayIO.h>
#include <casa/BasicMath/Math.h>
#include <casa/BasicSL/Complex.h>
#include <casa/Logging/LogIO.h>
#include <casa/OS/File.h>
#include <casa/Containers/Record.h>

#include <lattices/LRegions/LCBox.h>
#include <casa/Arrays/Slicer.h>
#include <scimath/Mathematics/FFTServer.h>
#include <casa/OS/HostInfo.h>
#include <casa/Arrays/ArrayError.h>
#include <casa/Arrays/ArrayIter.h>
#include <casa/Arrays/VectorIter.h>

#include <casa/Utilities/GenSort.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>
#include <casa/Utilities/Fallible.h>

#include <casa/BasicSL/Constants.h>

#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogMessage.h>

#include <synthesis/MeasurementEquations/MatrixCleaner.h>
#include <synthesis/TransformMachines/StokesImageUtil.h>
#include <synthesis/TransformMachines2/Utils.h>
#include <coordinates/Coordinates/TabularCoordinate.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Additional include files
#include <algorithm>
#include <chrono>

#include<synthesis/MeasurementEquations/AspMatrixCleaner.h>

//for gsl bfgs
//#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_min.h>
//#include <gsl/gsl_blas.h>
#include <synthesis/MeasurementEquations/gslobjfunc.h>

// for lbfgs
/*#include <limits>
#include <LBFGS.h>
#include <LBFGSB.h>
//#include <synthesis/MeasurementEquations/lbfgsAsp.h>
#include <synthesis/MeasurementEquations/lbfgsGaussian.h>
//#include <synthesis/MeasurementEquations/lbfgsAspZhang.h>
using namespace LBFGSpp;*/

// for gsl
using Eigen::VectorXd;

//for CppNumericalSolvers
/*#include "./cppoptlib/solver/lbfgsbsolver.h"
#include <synthesis/MeasurementEquations/objfunc.h>
using namespace cppoptlib;*/

// for alglib
//#include "stdafx.h"
//#include "optimization.h"
//#include "objfunc.h"
//using namespace alglib;

using namespace casacore;
using namespace std::chrono;
namespace casa { //# NAMESPACE CASA - BEGIN
AspMatrixCleaner::AspMatrixCleaner():
  MatrixCleaner(),
  itsInitScaleSizes(0),
  itsAspScaleSizes(0),
  itsAspAmplitude(0),
  itsNInitScales(5),
  itsPrevLBFGSGrad(0.0),
  //itsNumIterNoGoodAspen(0),
  itsPsfWidth(0.0),
  itsUseZhang(false),
  itsSwitchedToHogbom(false),
  itsNumHogbomIter(0),
  itsNthHogbom(0),
  itsSwitchedToMS(false),
  //itsStrenThres(0.19), //m31_gaussian_trynorm8
  //itsStrenThres(0.4), //m31_gaussian_trynorm8
  //itsStrenThres(0.0001), //m31norm_hog, 5k->10k
  //itsStrenThres(0.003), //SNorm_hog
  //itsStrenThres(0.0001), //SNorm_nohog
  //itsStrenThres(0.00036), //SNorm_hog 5k->10k
  //itsStrenThres(0.00003), //m31Norm_nohog
  itsStrenThres(15.8), //GSL with der
  itsOptimumScale(0),
  itsOptimumScaleSize(0.0),
  itsPeakResidual(1000.0), // temp. should this be changed to MAX?
  itsOrigDirty( ),
  itsFusedThreshold(0.0)
{
  itsInitScales.resize(0);
  itsInitScaleXfrs.resize(0);
  itsDirtyConvInitScales.resize(0);
  itsInitScaleMasks.resize(0);
  itsPsfConvInitScales.resize(0);
  itsNumIterNoGoodAspen.resize(0);
  itsAspCenter.resize(0);
  itsAspGood.resize(0);
  itsGoodAspActiveSet.resize(0);
  itsGoodAspAmplitude.resize(0);
  itsGoodAspCenter.resize(0);
  itsPrevAspActiveSet.resize(0);
  itsPrevAspAmplitude.resize(0);
  itsUsedMemoryMB = double(HostInfo::memoryUsed()/2014);
  maxPsfConvInitScales.resize(0);
  itsNormMethod = casa::refim::SynthesisUtils::getenv("ASP_NORM", itsDefaultNorm);
}

AspMatrixCleaner::~AspMatrixCleaner()
{
  destroyAspScales();
  destroyInitMasks();
  //genie destroyInitScales here?
  if(!itsMask.null())
    itsMask=0;
}

/*void AspMatrixCleaner::makedirtyscales()
{
  LogIO os(LogOrigin("AspMatrixCleaner", "makedirtyscales()", WHERE));

  if(!itsScalesValid || itsNscales < 1 || itsDirty.null() || (itsNscales != Int(itsScaleXfrs.nelements())) )
    return;

  if( (psfShape_p) != (itsDirty->shape()))
    throw(AipsError("PSF and Dirty array are not of the same shape"));

  // "only" have 0 scale size -> switch to hogbom
  if (itsNscales == 1)
  {
    itsDirtyConvScales.resize(1, true);
    itsDirtyConvScales[0] = Matrix<Float>(itsDirty->shape());
    itsDirtyConvScales[0] = *itsDirty;
    return;
  }

  Matrix<Complex> dirtyFT;
  FFTServer<Float,Complex> fft(itsDirty->shape());
  fft.fft0(dirtyFT, *itsDirty);
  itsDirtyConvScales.resize(itsNscales, true);
  Int scale=0;

  //Having problem with fftw and omp
  // #pragma omp parallel default(shared) private(scale) firstprivate(fft)
  {
    //#pragma omp  for
    // Dirty*scale
    for (scale=0; scale<itsNscales; scale++)
    {
      os << "Calculating convolutions for scale size " << itsScaleSizes(scale) << LogIO::POST;
      Matrix<Complex> cWork;

      itsDirtyConvScales[scale]=Matrix<Float>(itsDirty->shape());
      cWork=((dirtyFT)*(itsScaleXfrs[scale]));
      fft.fft0((itsDirtyConvScales[scale]), cWork, false);
      fft.flip((itsDirtyConvScales[scale]), false, false);
      cout << "Calculating itsDirtyConvScales for scale size " << itsScaleSizes(scale) << endl;
    }
  }

}*/

void AspMatrixCleaner::makedirtyscale()
{
  LogIO os(LogOrigin("AspMatrixCleaner", "makedirtyscale()", WHERE));

  //if(!itsScalesValid || itsNscales < 1 || itsDirty.null() || (itsNscales != Int(itsScaleXfrs.nelements())) )
  //  return;

  //if( (psfShape_p) != (itsDirty->shape()))
  //  throw(AipsError("PSF and Dirty array are not of the same shape"));

  // 0 scale size -> switch to hogbom
  //if (itsOptimumScaleSize == 0)
  //{
    itsDirtyConvScales.resize(1, false);
    itsDirtyConvScales[0] = Matrix<Float>(itsDirty->shape());
    itsDirtyConvScales[0] = *itsDirty;
    return;
  //}

  /*Matrix<Complex> dirtyFT;
  FFTServer<Float,Complex> fft(itsDirty->shape());
  fft.fft0(dirtyFT, *itsDirty);
  itsDirtyConvScales.resize(1, false);

  // Dirty*scale

  os << "Calculating convolutions for scale size " << itsOptimumScaleSize << LogIO::POST;
  Matrix<Complex> cWork;

  itsDirtyConvScales[0]=Matrix<Float>(itsDirty->shape());
  cWork=((dirtyFT)*(itsScaleXfrs[0]));
  fft.fft0((itsDirtyConvScales[0]), cWork, false);
  fft.flip((itsDirtyConvScales[0]), false, false);
  cout << "Calculating itsDirtyConvScales for scale size " << itsOptimumScaleSize << endl;*/
}

Bool AspMatrixCleaner::setaspcontrol(const Int niter,
           const Float gain,
           const Quantity& aThreshold,
           const Quantity& fThreshold)
{
  itsMaxNiter=niter;
  itsGain=gain;
  itsThreshold=aThreshold;
  itsFracThreshold=fThreshold;
  return true;
}

// Do the clean as set up
Int AspMatrixCleaner::aspclean(Matrix<Float>& model,
                         Bool /*showProgress*/)
{
  AlwaysAssert(model.shape()==itsDirty->shape(), AipsError);

  LogIO os(LogOrigin("AspMatrixCleaner", "aspclean()", WHERE));

  cout << "Enter aspclean, itsNscales " << itsNscales << endl;
  os << LogIO::NORMAL1 << "Asp clean algorithm" << LogIO::POST;


  Int scale;
  Vector<Float> scaleBias(itsNscales);
  // scaleBias is 1.0 for AAsp for now
  for (scale=0; scale < itsNscales; scale++)
    scaleBias(scale) = 1.0;

  AlwaysAssert(itsScalesValid, AipsError);
  //no need to use all cores if possible
  Int nth = itsNscales;
#ifdef _OPENMP

    nth = min(nth, omp_get_max_threads());

#endif

  // Define a subregion for the inner quarter
  IPosition blcDirty(model.shape().nelements(), 0);
  IPosition trcDirty(model.shape()-1);

  if(!itsMask.null())
  {
    os << "Cleaning using given mask" << LogIO::POST;
    if (itsMaskThreshold < 0)
    {
        os << LogIO::NORMAL
           << "Mask thresholding is not used, values are interpreted as weights"
           <<LogIO::POST;
    }
    else
    {
      // a mask that does not allow for clean was sent
      if(noClean_p)
        return 0;

      os << LogIO::NORMAL
         << "Cleaning pixels with mask values above " << itsMaskThreshold
         << LogIO::POST;
    }

    Int nx=model.shape()(0);
    Int ny=model.shape()(1);

    AlwaysAssert(itsMask->shape()(0)==nx, AipsError);
    AlwaysAssert(itsMask->shape()(1)==ny, AipsError);
    Int xbeg=nx-1;
    Int ybeg=ny-1;
    Int xend=0;
    Int yend=0;
    for (Int iy=0;iy<ny;iy++)
    {
      for (Int ix=0;ix<nx;ix++)
      {
        if((*itsMask)(ix,iy)>0.000001)
        {
          xbeg=min(xbeg,ix);
          ybeg=min(ybeg,iy);
          xend=max(xend,ix);
          yend=max(yend,iy);
        }
      }
    }

    if (!itsIgnoreCenterBox)
    {
      if((xend - xbeg)>nx/2)
      {
        xbeg=nx/4-1; //if larger than quarter take inner of mask
        os << LogIO::WARN << "Mask span over more than half the x-axis: Considering inner half of the x-axis"  << LogIO::POST;
      }
      if((yend - ybeg)>ny/2)
      {
        ybeg=ny/4-1;
        os << LogIO::WARN << "Mask span over more than half the y-axis: Considering inner half of the y-axis" << LogIO::POST;
      }
      xend=min(xend,xbeg+nx/2-1);
      yend=min(yend,ybeg+ny/2-1);
    }

    blcDirty(0)=xbeg;
    blcDirty(1)=ybeg;
    trcDirty(0)=xend;
    trcDirty(1)=yend;
  }
  else
  {
    if (itsIgnoreCenterBox) {
      os << LogIO::NORMAL << "Cleaning entire image" << LogIO::POST;
      os << LogIO::NORMAL1 << "as per MF/WF" << LogIO::POST; // ???
    }
    else {
      os << "Cleaning inner quarter of the image" << LogIO::POST;
      for (Int i=0;i<Int(model.shape().nelements());i++)
      {
        blcDirty(i)=model.shape()(i)/4;
        trcDirty(i)=blcDirty(i)+model.shape()(i)/2-1;
        if(trcDirty(i)<0)
          trcDirty(i)=1;
      }
    }
  }
  LCBox centerBox(blcDirty, trcDirty, model.shape());



  // Start the iteration
  Float totalFlux=0.0;
  Int converged=0;
  Int stopPointModeCounter = 0;
  Float tmpMaximumResidual = 0.0;



  os << "Starting iteration"<< LogIO::POST;
  vector<Float> tempScaleSizes;
  itsIteration = itsStartingIter; // 0

  FFTServer<Float,Complex> fft(psfShape_p);
  Matrix<Float> itsScale0 = Matrix<Float>(psfShape_p);
  //makeScale(itsScale0, 0.0);
  //makeScaleImage(itsScale0, 0.0);
  Matrix<Complex>itsScaleXfr0 = Matrix<Complex> ();
  //fft.fft0(itsScaleXfr0, itsScale0);

  Matrix<Float> itsScale = Matrix<Float>(psfShape_p);
  Matrix<Complex>itsScaleXfr = Matrix<Complex> ();

  // Define a subregion so that the peak is centered
  IPosition support(model.shape());
  support(0) = max(Int(itsInitScaleSizes[itsNInitScales-1] + 0.5), support(0));
  support(1) = max(Int(itsInitScaleSizes[itsNInitScales-1] + 0.5), support(1));

  IPosition inc(model.shape().nelements(), 1);

  for (Int ii = itsStartingIter; ii < itsMaxNiter; ii++)
  {
    cout << "cur iter " << itsIteration << " max iter is "<<
            itsMaxNiter << endl;
    itsIteration++;

    // genie asp new ...
    // make single optimized scale
    cout << "aspclean: making scale " << itsOptimumScaleSize << endl;

    if (itsSwitchedToHogbom)
    {
      makeScaleImage(itsScale0, 0.0, itsStrengthOptimum, itsPositionOptimum);
      fft.fft0(itsScaleXfr0, itsScale0);
    	itsScale = 0.0;
    	itsScale = itsScale0;
    	itsScaleXfr.resize();
      itsScaleXfr = itsScaleXfr0;
    }
    else
    {
      makeScaleImage(itsScale, itsOptimumScaleSize, itsStrengthOptimum, itsPositionOptimum);
      //makeScaleImage(itsScale, itsOptimumScale, itsStrengthOptimum, itsPositionOptimum);
      itsScaleXfr.resize();
      fft.fft0(itsScaleXfr, itsScale);
    }

    //genie
    // trigger hogbom when itsStrengthOptimum is small enough
    // consider scale 5e-7 down every time this is triggered to see if imaging is improved
    //if (!itsSwitchedToHogbom && itsStrengthOptimum < 5e-7) // G55 value, no box
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 0/*itsStrenThres*/) //try
    //if (!itsSwitchedToHogbom && itsStrengthOptimum < 1e-7) // G55 value, with box
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 4) // old M31 value
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 0.55) // M31 value - new Asp: 5k->10k good
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 0.38) // M31 value - new Asp + new normalization: 5k->10k good
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 0.001) // M31 value - new Asp + gaussian
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 0) // M31 value - new Asp + gaussian
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 2.8) // M31 value-new asp: 1k->5k
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 0.002) // G55 value, new Asp, old normalization
    //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 0.0002) // G55 value, new Asp, sanjay's normalization
    //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 6) // GSL without der M31
    //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 2.2e-4 && abs(itsStrengthOptimum) < 1e-4) // GSL without der, G55 , with new norm this is not needed.
    //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 1.1e-4 && abs(itsStrengthOptimum) < 1e-4) // GSL with scaled der, G55 , with new norm this is not needed.
    //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 8e-5 && abs(itsStrengthOptimum) < 1e-4) // GSL with der, G55 , with new norm this is not needed.
    if (itsNormMethod == 1) // only Nomr Method 1 needs hogbom for speedup
    {
      //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 7e-3 && abs(itsStrengthOptimum) < 1e-7) // GSL with der, Points ,with new norm this is not needed.
      //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 0.138) //GSL M100 channel 22
      //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 0.15) // GSL M100 channel 22 & 23 & others
      //if(!itsSwitchedToHogbom && (abs(itsPeakResidual) < 2.5 || abs(itsStrengthOptimum) < 1e-3)) // GSL, CygA
      cout << "itsFusedThreshold " << itsFusedThreshold << endl;
      if(!itsSwitchedToHogbom && (abs(itsPeakResidual) < itsFusedThreshold || abs(itsStrengthOptimum) < 1e-3)) // GSL, CygA
      {
  	    cout << "Switch to hogbom b/c optimum strength is small enough: " << itsStrenThres << endl;
  	    //itsStrenThres = itsStrenThres/3.0; //box3
  	    //itsStrenThres = itsStrenThres/1.5; //box4, sNorm, SNorm2
  	    //itsStrenThres = itsStrenThres - 0.0005; //Snorm4, SNorm5

  	    switchedToHogbom();
      }
      /*if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 1.3e-3) // GSL with der, Points ,with new norm this is not needed.
      {
        cout << "Switch to hogbom b/c optimum strength is small enough: " << itsStrenThres << endl;
        switchedToHogbom();
      }*/
    }

    if (!itsSwitchedToHogbom)
    {
	    if (itsNumIterNoGoodAspen.size() >= 10)
	  	  itsNumIterNoGoodAspen.pop_front(); // only track the past 10 iters
	    if (itsOptimumScaleSize == 0)
	      itsNumIterNoGoodAspen.push_back(1); //genie Zhang 2018 fused-Asp approach
	    else
	      itsNumIterNoGoodAspen.push_back(0);
    }


    //itsStrengthOptimum /= scaleBias(optimumScale);
    //itsStrengthOptimum /= (vecWork_p[optimumScale])(posMaximum[optimumScale]);

    //AlwaysAssert(itsOptimumScale < itsNscales, AipsError);

    // Now add to the total flux
    totalFlux += (itsStrengthOptimum*itsGain);
    itsTotalFlux = totalFlux;

    if(ii == itsStartingIter)
    {
      itsMaximumResidual = abs(itsStrengthOptimum);
      tmpMaximumResidual = itsMaximumResidual;
      os << "Initial maximum residual is " << itsMaximumResidual;
      if( !itsMask.null() )
        os << " within the mask ";

      os << LogIO::POST;
    }

    // Various ways of stopping:
    //    1. stop if below threshold
    /*if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 1e-7/*threshold()* /)
    {
    	cout << "Reached stopping threshold " << threshold() << " at iteration "<< ii << endl;
      os << "Reached stopping threshold " << threshold() << " at iteration "<<
            ii << LogIO::POST;
      os << "Optimum flux is " << abs(itsStrengthOptimum) << LogIO::POST;
      converged = 1;
      break;
    }*/
    // with new norm this is not needed
    if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 2e-4 /*&& abs(itsStrengthOptimum) < 1e-7*//*threshold()*/)
    {
      cout << "Reached stopping threshold " << threshold() << " at iteration "<< ii << endl;
      os << "Reached stopping threshold " << threshold() << " at iteration "<<
            ii << LogIO::POST;
      os << "Optimum flux is " << abs(itsStrengthOptimum) << LogIO::POST;
      converged = 1;
      break;
    }
    //    2. negatives on largest scale?
    if ((itsNscales > 1) && itsStopAtLargeScaleNegative  &&
        //optimumScale == (itsNscales-1) &&
    	  itsOptimumScale == (itsNInitScales-1) &&
        itsStrengthOptimum < 0.0)
    {
      os << "Reached negative on largest scale" << LogIO::POST;
      converged = -2;
      break;
    }
    //  3. stop point mode at work
    if (itsStopPointMode > 0)
    {
      if (itsOptimumScale == 0)
        stopPointModeCounter++;
      else
        stopPointModeCounter = 0;

      if (stopPointModeCounter >= itsStopPointMode)
      {
        os << "Cleaned " << stopPointModeCounter <<
          " consecutive components from the smallest scale, stopping prematurely"
           << LogIO::POST;
        itsDidStopPointMode = true;
        converged = -1;
        break;
      }
    }
    //4. Diverging large scale
    //If actual value is 50% above the maximum residual. ..good chance it will not recover at this stage
    if(((abs(itsStrengthOptimum)-abs(tmpMaximumResidual)) > (abs(tmpMaximumResidual)/2.0))
       && !(itsStopAtLargeScaleNegative))
    {
      cout << "Diverging due to large scale?" << endl;
      os << "Diverging due to large scale?" << LogIO::POST;
      os << "itsStrengthOptimum " << itsStrengthOptimum << " tmp " << tmpMaximumResidual << LogIO::POST;
       //clean is diverging most probably due to the large scale
      converged=-2;
      break;
    }
    //5. Diverging for some other reason; may just need another CS-style reconciling
    if((abs(itsStrengthOptimum)-abs(tmpMaximumResidual)) > (abs(tmpMaximumResidual)/2.0))
    {
      os << "Diverging due to unknown reason" << LogIO::POST;
      converged=-3;
      break;
    }

      if (itsIteration == itsStartingIter + 1)
        os << "iteration    MaximumResidual   CleanedFlux" << LogIO::POST;
      if ((itsIteration % (itsMaxNiter/10 > 0 ? itsMaxNiter/10 : 1)) == 0)
      {
        //Good place to re-up the fiducial maximum residual
        //tmpMaximumResidual=abs(itsStrengthOptimum);
        os << itsIteration <<"      "<<itsStrengthOptimum<<"      "
           << totalFlux <<LogIO::POST;
      }


    // genie Continuing: subtract the peak that we found from the dirty image
    // Define a subregion so that the peak is centered
    /*IPosition support(model.shape());
    support(0) = max(Int(itsInitScaleSizes[itsNInitScales-1] + 0.5), support(0));
    support(1) = max(Int(itsInitScaleSizes[itsNInitScales-1] + 0.5), support(1));

    IPosition inc(model.shape().nelements(), 1);*/
    //cout << "support " << support.asVector()  << endl;

    // genie test
    IPosition blc(itsPositionOptimum - support/2);
    IPosition trc(itsPositionOptimum + support/2 - 1);
    LCBox::verify(blc, trc, inc, model.shape());

    /*IPosition blcPsf(blc + itsPositionPeakPsf - itsPositionOptimum);
    IPosition trcPsf(trc + itsPositionPeakPsf - itsPositionOptimum);
    LCBox::verify(blcPsf, trcPsf, inc, model.shape());
    makeBoxesSameSize(blc, trc, blcPsf, trcPsf);*/
    IPosition blcPsf(blc);
    IPosition trcPsf(trc);
    //cout << "itsPositionPeakPsf " << itsPositionPeakPsf << " blc " << blc;
    //cout << " trc " << trc << " blcPsf " << blcPsf << " trcPsf " << trcPsf << endl;
    blcDirty = blc;
    trcDirty = trc;

    Matrix<Float> modelSub = model(blc, trc);
    Float scaleFactor;
    scaleFactor = itsGain * itsStrengthOptimum;
    //scaleFactor = itsGain;
    Matrix<Float> scaleSub = (itsScale)(blcPsf,trcPsf);
    modelSub += scaleFactor * scaleSub; //MS'way
    //cout << "before dirtySub(262,291) = " << dirtySub(262,291) << " itsGain " << itsGain << " itsStrengthOptimum " << itsStrengthOptimum << endl;
    //cout << "before itsDirty(262,291) = " << (*itsDirty)(262,291) << endl;

    // Now update the residual image
    //PSF * model
    Matrix<Complex> cWork;
    //cout << "aspclean: Calculating Psf convolution for model " << itsOptimumScaleSize << endl;
    cWork = ((*itsXfr)*(itsScaleXfr)); //Asp's
    Matrix<Float> itsPsfConvScale = Matrix<Float>(psfShape_p);
    fft.fft0(itsPsfConvScale, cWork, false);
    fft.flip(itsPsfConvScale, false, false); //need this if conv with 1 scale; don't need this if conv with 2 scales
    Matrix<Float> psfSub = (itsPsfConvScale)(blcPsf, trcPsf);
    Matrix<Float> dirtySub=(*itsDirty)(blc,trc);

    float maxvalue;
    IPosition peakpos;
    findMaxAbs(psfSub, maxvalue, peakpos);
    cout << "psfSub pos " << peakpos << " maxval " << psfSub(peakpos) << endl;
    cout << "dirtySub pos " << peakpos << " val " << dirtySub(peakpos) << endl;
    findMaxAbs(itsPsfConvScale, maxvalue, peakpos);
    cout << "itsPsfConvScale pos " << peakpos << " maxval " << itsPsfConvScale(peakpos) << endl;
    cout << "itsDirty pos " << peakpos << " val " << (*itsDirty)(peakpos) << endl;
    /*findMaxAbs(dirtySub, maxvalue, peakpos);
    cout << "dirtySub pos " << peakpos << " maxval " << dirtySub(peakpos) << endl;
    findMaxAbs((*itsDirty), maxvalue, peakpos);
    cout << "itsDirty pos " << peakpos << " maxval " << (*itsDirty)(peakpos) << endl;*/
    cout << "itsPositionOptimum " << itsPositionOptimum << endl;
    cout << "itsStrengthOptimum " << itsStrengthOptimum << endl;
    cout << " maxPsfSub " << max(fabs(psfSub)) << " maxPsfConvScale " << max(fabs(itsPsfConvScale)) << " itsGain " << itsGain << endl;
    dirtySub -= scaleFactor * psfSub;

    // further update the model and residual with the remaining aspen of the active-set
    if (itsOptimumScaleSize != 0)
    {
    	bool aspenUnchanged = true;
    	if ((Int)itsGoodAspActiveSet.size() <= 1)
    		aspenUnchanged = false;

      for (scale = 0; scale < (Int)itsGoodAspActiveSet.size() - 1; ++scale)
      // -1 because we counted the latest aspen in the previous step already
      {
        if (itsPrevAspActiveSet[scale] == itsGoodAspActiveSet[scale])
          continue;

        cout << "I have active-set to adjust" << endl;
        aspenUnchanged = false;
        // "center" is unchanged for aspen
        // genie test
        IPosition blc1(itsGoodAspCenter[scale] - support/2);
        IPosition trc1(itsGoodAspCenter[scale] + support/2 - 1);
        LCBox::verify(blc1, trc1, inc, model.shape());

        /*IPosition blcPsf1(blc1 + itsPositionPeakPsf - itsGoodAspCenter[scale]);
        IPosition trcPsf1(trc1 + itsPositionPeakPsf - itsGoodAspCenter[scale]);
        LCBox::verify(blcPsf1, trcPsf1, inc, model.shape());
        makeBoxesSameSize(blc1, trc1, blcPsf1, trcPsf1);*/
        IPosition blcPsf1(blc1);
        IPosition trcPsf1(trc1);

        Matrix<Float> modelSub1 = model(blc1, trc1);
        Matrix<Float> dirtySub1 = (*itsDirty)(blc1,trc1);

        // First remove the previous values of aspen in the active-set
        cout << "aspclean: restore with previous scale " << itsPrevAspActiveSet[scale];
        cout << " amp " << itsPrevAspAmplitude[scale] << endl;
        //cout << " at blcPsf1 " << blcPsf1.asVector() << " trcPsf1 " << trcPsf1.asVector();
        //cout << "at blc1 " << blc1.asVector() << " trc1 " << trc1.asVector() << endl;

        //makeScale(itsScale, itsPrevAspActiveSet[scale]);
        makeScaleImage(itsScale, itsPrevAspActiveSet[scale], itsPrevAspAmplitude[scale], itsGoodAspCenter[scale]);
        itsScaleXfr.resize();
        fft.fft0(itsScaleXfr, itsScale);
        Matrix<Float> scaleSubPrev = (itsScale)(blcPsf1,trcPsf1);
        const float scaleFactorPrev = itsGain * itsPrevAspAmplitude[scale];
        //const float scaleFactorPrev = itsGain;
        // restore the model image...
        modelSub1 -= scaleFactorPrev * scaleSubPrev;
        // restore the residual image
        Matrix<Complex> cWorkPrev;
        cWorkPrev = ((*itsXfr)*(itsScaleXfr));
        Matrix<Float> itsPsfConvScalePrev = Matrix<Float>(psfShape_p);
        fft.fft0(itsPsfConvScalePrev, cWorkPrev, false);
        fft.flip(itsPsfConvScalePrev, false, false); //need this if conv with 1 scale; don't need this if conv with 2 scales
        Matrix<Float> psfSubPrev = (itsPsfConvScalePrev)(blcPsf1, trcPsf1);
        dirtySub1 += scaleFactorPrev * psfSubPrev;

        // Then update with the new values of aspen in the active-set
        cout << "aspclean: update with new scale " << itsGoodAspActiveSet[scale];
        cout << " amp " << itsGoodAspAmplitude[scale] << endl;
        //makeScale(itsScale, itsGoodAspActiveSet[scale]);
        makeScaleImage(itsScale, itsGoodAspActiveSet[scale], itsGoodAspAmplitude[scale], itsGoodAspCenter[scale]);
        itsScaleXfr.resize();
        fft.fft0(itsScaleXfr, itsScale);
        Matrix<Float> scaleSubNew = (itsScale)(blcPsf1,trcPsf1);
        const float scaleFactorNew = itsGain * itsGoodAspAmplitude[scale];
        //const float scaleFactorNew = itsGain;
        // Now do the addition of the active-set scales to the model image...
        modelSub1 += scaleFactorNew * scaleSubNew;
        // Now subtract the active-set scales from the residual image
        Matrix<Complex> cWorkNew;
        cWorkNew = ((*itsXfr)*(itsScaleXfr));
        Matrix<Float> itsPsfConvScaleNew = Matrix<Float>(psfShape_p);
        fft.fft0(itsPsfConvScaleNew, cWorkNew, false);
        fft.flip(itsPsfConvScaleNew, false, false); //need this if conv with 1 scale; don't need this if conv with 2 scales
        Matrix<Float> psfSubNew = (itsPsfConvScaleNew)(blcPsf1, trcPsf1);
        dirtySub1 -= scaleFactorNew * psfSubNew;
      } //update updating model/residual from aspen in active-set

      // switch to hogbom if aspen is unchaned?
	    /*if (!itsSwitchedToHogbom && aspenUnchanged)
	    {
	    	cout << "Switched to hogbom b/c aspen is unchanged" << endl;
	    	switchedToHogbom();
	    }*/
    }

    //cout << "after dirtySub(262,291) = " << dirtySub(262,291) << endl;
    //cout << "after itsDirty(262,291) = " << (*itsDirty)(262,291) << endl;
    Float maxVal=0;
    IPosition posmin((*itsDirty).shape().nelements(), 0);
    Float minVal=0;
    IPosition posmax((*itsDirty).shape().nelements(), 0);
    minMaxMasked(minVal, maxVal, posmin, posmax, (*itsDirty), itsInitScaleMasks[0]);
    itsPeakResidual = (fabs(maxVal) > fabs(minVal)) ? fabs(maxVal) : fabs(minVal);
    cout << "current peakres " << itsPeakResidual << " posmin " << posmin << " val " << (*itsDirty)(posmin);
    cout << " posmax " << posmax << " val " << (*itsDirty)(posmax) << endl;
    cout << "after: itsDirty optPos " << (*itsDirty)(itsPositionOptimum) << endl;

    //genie If we switch to hogbom (i.e. only have 0 scale size)
    // there is no need to do the following Aspen update

    //if (itsSwitchedToHogbom && itsNscales == 1)
    if (itsSwitchedToHogbom)
    {
      if (itsNumHogbomIter == 0)
      {
        itsSwitchedToHogbom = false;
        cout << "switched back to Asp." << endl;
      }
      else
        itsNumHogbomIter -= 1;

      //continue;
    }

    //genie Now update the actual residual image
    // At this point, itsDirty is not updated. Only itsDirtyConvScales is updated.
    //cout << "after dirtySub(100,100) = " << itsDirtyConvScales[0](100,100) << endl;
    //cout << "after itsdirty(100,100) = " << (*itsDirty)(100,100) << endl;

    tempScaleSizes.clear();
    tempScaleSizes = getActiveSetAspen();
    /*for (scale = 0; scale < int(tempScaleSizes.size()); scale++)
      cout << "2. getActiveSetAspen[" << scale << "] " << tempScaleSizes[scale] << endl;
    cout << "# tempScaleSizes " << tempScaleSizes.size() << endl;*/
    tempScaleSizes.push_back(0.0); // put 0 scale
    defineAspScales(tempScaleSizes);
    //makePsfScales();
    //makeScaleMasks();
    //makedirtyscale();
  }
  // End of iteration

  // memory used
  //itsUsedMemoryMB = double(HostInfo::memoryUsed()/1024);
  //cout << "Memory allocated in aspclean " << itsUsedMemoryMB << " MB." << endl;

  if(!converged)
    os << "Failed to reach stopping threshold" << LogIO::POST;

  //genie model should be updated here once?
  //Matrix<Float> modelSub = model(blc, trc);
  //Matrix<Float> dirtySub=(*itsDirty)(blc,trc);
  //Float scaleFactor;
  //scaleFactor = itsGain * itsStrengthOptimum; // MS's way
  //Matrix<Float> scaleSub = (itsScale)(blcPsf,trcPsf); //MS's way
  //modelSub += scaleFactor * scaleSub; //MS'way

  /*#pragma omp parallel default(shared) private(scale) num_threads(nth)
  {
    #pragma omp for
    for (scale = 0; scale < itsNscales; ++scale)
    {
    	Matrix<Float> scaleSub = (itsScales[scale])(blcPsf,trcPsf);
    	if (itsScaleSizes[scale] == 0)
    	  scaleFactor = itsGain * itsStrengthOptimum;
    	else
    	  scaleFactor = itsGain * itsGoodAspAmplitude[scale];

      // Now do the addition of the active-set scales to the model image...
      modelSub += scaleFactor * scaleSub;
    }
  }*/

  return converged;
}


Bool AspMatrixCleaner::destroyAspScales()
{
  destroyScales();

  for(uInt scale=0; scale < itsDirtyConvInitScales.nelements(); scale++)
    itsDirtyConvInitScales[scale].resize();

  itsDirtyConvInitScales.resize(0, true);

  return true;
}

Bool AspMatrixCleaner::destroyInitScales()
{
  for(uInt scale=0; scale < itsInitScales.nelements(); scale++)
    itsInitScales[scale].resize();
  for(uInt scale=0; scale < itsInitScaleXfrs.nelements(); scale++)
    itsInitScaleXfrs[scale].resize();
  for(uInt scale=0; scale < itsPsfConvInitScales.nelements(); scale++)
    itsPsfConvInitScales[scale].resize();

  itsInitScales.resize(0, true);
  itsInitScaleXfrs.resize(0, true);
  itsPsfConvInitScales.resize(0, true);
  maxPsfConvInitScales.resize(0);

  return true;
}

Bool AspMatrixCleaner::destroyInitMasks()
{
  for(uInt scale=0; scale<itsInitScaleMasks.nelements();scale++)
    itsInitScaleMasks[scale].resize();

  itsInitScaleMasks.resize(0);

  return true;
}


float AspMatrixCleaner::getPsfGaussianWidth(ImageInterface<Float>& psf)
{
  GaussianBeam beam;
  try
  {
      StokesImageUtil::FitGaussianPSF(psf, beam);
  }
  catch(AipsError &x)
  {
    LogIO os( LogOrigin("AspMatrixCleaner","getPsfGaussianWidth",WHERE) );
    os << "Error in fitting a Gaussian to the PSF : " << x.getMesg() << LogIO::POST;
    throw( AipsError("Error in fitting a Gaussian to the PSF" + x.getMesg()) );
  }

  CoordinateSystem cs = psf.coordinates();
  String dirunit = cs.worldAxisUnits()(0);
  Vector<String> unitas=cs.worldAxisUnits();
  unitas(0)="arcsec";
  unitas(1)="arcsec";
  cs.setWorldAxisUnits(unitas);

  cout << "major width " << beam.getMajor("arcsec") << " in " << cs.worldAxisUnits()(0) << endl;
  cout << "minor width " << beam.getMinor("arcsec") << endl;
  cout << " pixel sizes are " << abs(cs.increment()(0)) << " and ";
  cout << abs(cs.increment()(1)) << endl;
  const auto xpixels = beam.getMajor("arcsec") / abs(cs.increment()(0));
  const auto ypixels = beam.getMinor("arcsec") / abs(cs.increment()(1));
  cout << "xpixels " << xpixels << " ypixels " << ypixels << endl;
  cout << "init width " << float(ceil((xpixels + ypixels)/2)) << endl;

  itsPsfWidth = float(ceil((xpixels + ypixels)/2));

  return itsPsfWidth;
}

// Make a single scale size image by Gaussian
void AspMatrixCleaner::makeInitScaleImage(Matrix<Float>& iscale, const Float& scaleSize)
{

  const Int nx = iscale.shape()(0);
  const Int ny = iscale.shape()(1);
  iscale = 0.0;

  const Double refi = nx/2;
  const Double refj = ny/2;

  if(scaleSize == 0.0)
    iscale(Int(refi), Int(refj)) = 1.0;
  else
  {
    AlwaysAssert(scaleSize>0.0, AipsError);

    const Int mini = max(0, (Int)(refi - scaleSize));
    const Int maxi = min(nx-1, (Int)(refi + scaleSize));
    const Int minj = max(0, (Int)(refj - scaleSize));
    const Int maxj = min(ny-1, (Int)(refj + scaleSize));
    //cout << "makeScaleImage: scalesize " << scaleSize << " mini " << mini << " maxi " << maxi << " minj " << minj << " maxj " << maxj << endl;
    cout << "makeInitScaleImage: initscalesize " << scaleSize << endl;

    Gaussian2D<Float> gbeam(1.0/(sqrt(2*M_PI)*scaleSize), 0, 0, scaleSize, 1, 0);
    /*for (int j = minj; j <= maxj; j++)
    {
      for (int i = mini; i <= maxi; i++)
      {*/ //with this the max is 0
    float area = 0.0;
    float areag = 0.0;
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        const int px = i - refi;
        const int py = j - refj;
        ////iscale(i,j) = gbeam(px, py);
        iscale(i,j) = (1.0/(sqrt(2*M_PI)*scaleSize))*exp(-(pow(i-refi,2) + pow(j-refj,2))*0.5/pow(scaleSize,2)); //this is for 1D, but represents Sanjay's and gives good init scale
        //iscale(i,j) = (1.0/(2*M_PI*pow(scaleSize,2)))*exp(-(pow(i-refi,2) + pow(j-refj,2))*0.5/pow(scaleSize,2)); // this is for 2D, gives unit area but bad init scale (always 0)
        area += iscale(i,j);
        areag += gbeam(px, py);
      }
    }

    //debug
    /*float maxima;
    IPosition pos;
    findMaxAbs(iscale, maxima, pos);
    cout << "peak(iscale) " << max(fabs(iscale));
    cout << " maxima " << maxima << " pos " << pos << " iscale(pos) " << iscale(pos);
    cout << " area " << area << " areag " << areag << endl;*/
  }

}

// Make a single scale size image by Gaussian
void AspMatrixCleaner::makeScaleImage(Matrix<Float>& iscale, const Float& scaleSize, const Float& amp, const IPosition& center)
{

  const Int nx = iscale.shape()(0);
  const Int ny = iscale.shape()(1);
  iscale = 0.0;

  const Double refi = nx/2;
  const Double refj = ny/2;

  cout << "try to makeScaleImage: scalesize " << scaleSize << " center " << center << " amp " << amp << endl;
  if(scaleSize == 0.0)
    //iscale(Int(refi), Int(refj)) = 1.0;
    iscale(Int(center[0]), Int(center[1])) = 1.0;// genie check on me
  else
  {
    AlwaysAssert(scaleSize>0.0, AipsError);

    /*const Int mini = max(0, (Int)(refi - scaleSize));
    const Int maxi = min(nx-1, (Int)(refi + scaleSize));
    const Int minj = max(0, (Int)(refj - scaleSize));
    const Int maxj = min(ny-1, (Int)(refj + scaleSize));*/
    //cout << "makeScaleImage: scalesize " << scaleSize << " mini " << mini << " maxi " << maxi << " minj " << minj << " maxj " << maxj << endl;
    cout << "makeScaleImage: scalesize " << scaleSize << " center " << center << " amp " << amp << endl;

    const Float d = sqrt(pow(1.0/itsPsfWidth, 2) + pow(1.0/scaleSize, 2));
    //Gaussian2D<Float> gbeam(amp / (sqrt(2*M_PI)*scaleSize), center[0], center[1], scaleSize, 1, 0);
    Gaussian2D<Float> gbeam(1.0 / (sqrt(2*M_PI)*scaleSize), center[0], center[1], scaleSize, 1, 0);
    //Gaussian2D<Float> gbeam(amp * sqrt(2*M_PI)/d, center[0], center[1], scaleSize * d, 1, 0);
    //Gaussian2D<Float> gbeam(amp / (2*M_PI), center[0], center[1], scaleSize, 1, 0);
    //Gaussian2D<Float> gbeam(amp / pow(2,scaleSize-1), center[0], center[1], itsInitScaleSizes[scaleSize], 1, 0);
    /*for (int j = minj; j <= maxj; j++)
    {
      for (int i = mini; i <= maxi; i++)
      {*/ //with this the max is 0
    //double area = 0.0;
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        const int px = i;
        const int py = j;
        ////iscale(i,j) = gbeam(px, py);
        iscale(i,j) = (1.0/(sqrt(2*M_PI)*scaleSize))*exp(-(pow(i-center[0],2) + pow(j-center[1],2))*0.5/pow(scaleSize,2));
        //iscale(i,j) = (1.0/(2*M_PI*pow(scaleSize,2)))*exp(-(pow(i-center[0],2) + pow(j-center[1],2))*0.5/pow(scaleSize,2));
        //area += iscale(i,j);
      }
    }
    //debug
    /*float maxima;
    IPosition pos;
    findMaxAbs(iscale, maxima, pos);
    cout << "peak(iscale) " << max(fabs(iscale)) << " area " << area;
    cout << " maxima " << maxima << " pos " << pos << " iscale(pos) " << iscale(pos) << endl;*/
  }
}

/*void AspMatrixCleaner::makeAspScales()
{
  LogIO os(LogOrigin("AspMatrixCleaner", "makeAspScales()", WHERE));
  if(itsNscales < 1)
    throw(AipsError("Scales have to be set"));
  if(itsXfr.null())
    throw(AipsError("Psf is not defined"));
  //cout << "before destroy, size = " << itsNscales << " valid " << itsScalesValid << endl;
  destroyScales();
  //cout << "after destroy, size = " << itsNscales << " valid " << itsScalesValid << endl;
  itsScales.resize(itsNscales, true);
  itsScaleXfrs.resize(itsNscales, true);
  FFTServer<Float,Complex> fft(psfShape_p);

  for(int scale=0; scale < itsNscales; scale++)
  {
    itsScales[scale] = Matrix<Float>(psfShape_p);
    cout << "scale size[ " << scale << "] = " << itsScaleSizes(scale) << endl;
    makeScale(itsScales[scale], itsScaleSizes(scale));

    itsScaleXfrs[scale] = Matrix<Complex> ();
    fft.fft0(itsScaleXfrs[scale], itsScales[scale]);
  }

  itsScalesValid = true;
}*/


void AspMatrixCleaner::setInitScaleXfrs(/*const Array<Float> arrpsf, */const Float width)
{
  if(itsInitScales.nelements() > 0)
    destroyAspScales();

  /*Matrix<Float> tempMat;
  tempMat.reference(arrpsf);
  psfShape_p.resize(0, false);
  psfShape_p = tempMat.shape();*/

  // try 0, width, 2width, 4width and 8width
  itsInitScaleSizes.resize(itsNInitScales, false);
  itsInitScaleSizes = {0.0f, width, 2.0f*width, 4.0f*width, 8.0f*width};
  //itsInitScaleSizes = {0.0f, 2.0f, 4.0f, 8.0f, 16.0f};
  itsInitScales.resize(itsNInitScales, false);
  itsInitScaleXfrs.resize(itsNInitScales, false);
  FFTServer<Float,Complex> fft(psfShape_p);
  for (int scale = 0; scale < itsNInitScales; scale++)
  {
    itsInitScales[scale] = Matrix<Float>(psfShape_p);
    //makeScale(itsInitScales[scale], itsInitScaleSizes[scale]);
    makeInitScaleImage(itsInitScales[scale], itsInitScaleSizes[scale]);
    //cout << "made itsInitScales[" << scale << "] = " << itsInitScaleSizes[scale] << endl;
    //cout << "itsInitScales[" << scale << "](262,291) = " << (itsInitScales[scale])(262,291) << endl;
    //cout << "max itsInitScales[" << scale << "] = " << max(abs(itsInitScales[scale])) << endl;
    itsInitScaleXfrs[scale] = Matrix<Complex> ();
    fft.fft0(itsInitScaleXfrs[scale], itsInitScales[scale]);
    //cout << "itsInitScaleXfrs[" << scale << "](262,291) = " << (itsInitScaleXfrs[scale])(262,291) << endl;
    //cout << "max itsInitScaleXfrs[" << scale << "] = " << max(abs(itsInitScaleXfrs[scale])) << endl;
  }
  cout << "Norm Method is " << itsNormMethod << endl;
}

// calculate the convolutions of the psf with the initial scales
void AspMatrixCleaner::setInitScalePsfs()
{
  itsPsfConvInitScales.resize((itsNInitScales+1)*(itsNInitScales+1), false);
  itsNscales = itsNInitScales; // # initial scales. This will be updated in defineAspScales later.
  FFTServer<Float,Complex> fft(psfShape_p);

  Matrix<Complex> cWork;

  for (Int scale=0; scale < itsNInitScales; scale++)
  {
    //cout << "Calculating convolutions of psf for initial scale size " << itsInitScaleSizes[scale] << endl;
    //PSF * scale
    itsPsfConvInitScales[scale] = Matrix<Float>(psfShape_p);
    cWork=((*itsXfr)*(itsInitScaleXfrs[scale])*(itsInitScaleXfrs[scale]));
    fft.fft0((itsPsfConvInitScales[scale]), cWork, false);
    fft.flip(itsPsfConvInitScales[scale], false, false);

    for (Int otherscale = scale; otherscale < itsNInitScales; otherscale++)
    {
      AlwaysAssert(index(scale, otherscale) < Int(itsPsfConvInitScales.nelements()),
       AipsError);

      // PSF *  scale * otherscale
      itsPsfConvInitScales[index(scale,otherscale)] = Matrix<Float>(psfShape_p);
      cWork=((*itsXfr)*(itsInitScaleXfrs[scale])*(itsInitScaleXfrs[otherscale]));
      fft.fft0(itsPsfConvInitScales[index(scale,otherscale)], cWork, false);
    }
  }
}

// Set up the masks for the initial scales (i.e. 0, 1.5width, 5width and 10width)
Bool AspMatrixCleaner::setInitScaleMasks(const Array<Float> arrmask, const Float& maskThreshold)
{
  LogIO os(LogOrigin("AspMatrixCleaner", "setInitScaleMasks()", WHERE));

  destroyMasks();

  Matrix<Float> mask(arrmask);
  itsMask = new Matrix<Float>(mask.shape());
  itsMask->assign(mask);
  itsMaskThreshold = maskThreshold;
  noClean_p=(max(*itsMask) < itsMaskThreshold) ? true : false;

  if(itsMask.null() || noClean_p)
    return false;

  // make scale masks
  if(itsInitScaleSizes.size() <= 1)
  {
    os << "Initial scales are not yet set - cannot set initial scale masks"
       << LogIO::EXCEPTION;
  }

  AlwaysAssert((itsMask->shape() == psfShape_p), AipsError);
  FFTServer<Float,Complex> fft(itsMask->shape());

  Matrix<Complex> maskFT;
  fft.fft0(maskFT, *itsMask);
  itsInitScaleMasks.resize(itsNInitScales);
  // Now we can do all the convolutions
  Matrix<Complex> cWork;
  for (int scale=0; scale < itsNInitScales; scale++)
  {
    // Mask * scale
    // Allow only 10% overlap by default, hence 0.9 is a default mask threshold
    // if thresholding is not used, just extract the real part of the complex mask
    itsInitScaleMasks[scale] = Matrix<Float>(itsMask->shape());
    cWork=((maskFT)*(itsInitScaleXfrs[scale]));
    fft.fft0(itsInitScaleMasks[scale], cWork, false);
    fft.flip(itsInitScaleMasks[scale], false, false);
    for (Int j=0 ; j < (itsMask->shape())(1); ++j)
    {
      for (Int k =0 ; k < (itsMask->shape())(0); ++k)
      {
        if(itsMaskThreshold > 0)
          (itsInitScaleMasks[scale])(k,j) =  (itsInitScaleMasks[scale])(k,j) > itsMaskThreshold ? 1.0 : 0.0;
      }
    }
    Float mysum = sum(itsInitScaleMasks[scale]);
    if (mysum <= 0.1) {
      os << LogIO::WARN << "Ignoring initial scale " << itsInitScaleSizes[scale] <<
  " since it is too large to fit within the mask" << LogIO::POST;
    }

  }

   Int nx = itsInitScaleMasks[0].shape()(0);
   Int ny = itsInitScaleMasks[0].shape()(1);

   /* Set the edges of the masks according to the scale size */
   // Set the values OUTSIDE the box to zero....
  for(Int scale=0; scale < itsNInitScales; scale++)
  {
    Int border = (Int)(itsInitScaleSizes[scale]*1.5);
    // bottom
    IPosition blc1(2, 0 , 0 );
    IPosition trc1(2, nx-1, border);
    IPosition inc1(2, 1);
    LCBox::verify(blc1, trc1, inc1, itsInitScaleMasks[scale].shape());
    (itsInitScaleMasks[scale])(blc1, trc1) = 0.0;
    // top
    blc1[0]=0; blc1[1] = ny-border-1;
    trc1[0] = nx-1; trc1[1] = ny-1;
    LCBox::verify(blc1, trc1, inc1, itsInitScaleMasks[scale].shape());
    (itsInitScaleMasks[scale])(blc1, trc1) = 0.0;
    // left
    blc1[0]=0; blc1[1]=border;
    trc1[0]=border; trc1[1] = ny-border-1;
    LCBox::verify(blc1, trc1, inc1, itsInitScaleMasks[scale].shape());
    (itsInitScaleMasks[scale])(blc1, trc1) = 0.0;
    // right
    blc1[0] = nx-border-1; blc1[1]=border;
    trc1[0] = nx; trc1[1] = ny-border-1;
    LCBox::verify(blc1, trc1, inc1, itsInitScaleMasks[scale].shape());
    (itsInitScaleMasks[scale])(blc1,trc1) = 0.0;
  }

  return true;
}

void AspMatrixCleaner::maxDirtyConvInitScales(float& strengthOptimum, int& optimumScale, IPosition& positionOptimum)
{
  LogIO os(LogOrigin("AspMatrixCleaner", "maxDirtyConvInitScales()", WHERE));

  const int nx = itsDirty->shape()[0];
  const int ny = itsDirty->shape()[1];
  //cout << "nx " << nx << " ny " << ny << endl;
  IPosition blcDirty(itsDirty->shape().nelements(), 0);
  IPosition trcDirty(itsDirty->shape() - 1);

  if(!itsMask.null())
  {
    os << "Finding initial scales for Asp using given mask" << LogIO::POST;
    if (itsMaskThreshold < 0)
    {
        os << LogIO::NORMAL
           << "Mask thresholding is not used, values are interpreted as weights"
           <<LogIO::POST;
    }
    else
    {
      // a mask that does not allow for clean was sent
      if(noClean_p)
        return;

      os << LogIO::NORMAL
         << "Finding initial scales with mask values above " << itsMaskThreshold
         << LogIO::POST;
    }

    AlwaysAssert(itsMask->shape()(0) == nx, AipsError);
    AlwaysAssert(itsMask->shape()(1) == ny, AipsError);
    Int xbeg=nx-1;
    Int ybeg=ny-1;
    Int xend=0;
    Int yend=0;
    for (Int iy=0;iy<ny;iy++)
    {
      for (Int ix=0;ix<nx;ix++)
      {
        if((*itsMask)(ix,iy)>0.000001)
        {
          xbeg=min(xbeg,ix);
          ybeg=min(ybeg,iy);
          xend=max(xend,ix);
          yend=max(yend,iy);
        }
      }
    }
    blcDirty(0)=xbeg;
    blcDirty(1)=ybeg;
    trcDirty(0)=xend;
    trcDirty(1)=yend;
  }
  else
    os << LogIO::NORMAL << "Finding initial scales using the entire image" << LogIO::POST;

  // Find the peak residual
  /*IPosition gip;
  gip = IPosition(2, nx, ny);
  casacore::Block<casacore::Matrix<casacore::Float> > vecWork_p;
  vecWork_p.resize(0);
  vecWork_p.resize(itsNInitScales);

  for (int i = 0; i < itsNInitScales; i++)
    vecWork_p[i].resize(gip);*/

  Vector<Float> maxima(itsNInitScales);
  Block<IPosition> posMaximum(itsNInitScales);
  int scale;
  int nth = itsNInitScales;
  #ifdef _OPENMP
    nth = min(nth, omp_get_max_threads());
  #endif

  //genie
  // Find the peaks of the convolved Psfs
  /*Vector<Float> maxPsfConvInitScales(itsNInitScales);
  Int naxes = psfShape_p.nelements();

  #pragma omp parallel default(shared) private(scale) num_threads(nth)
  {
    #pragma omp for
    for (scale=0; scale < itsNInitScales; scale++)
    {
      IPosition positionPeakPsfConvInitScales(naxes, 0);

      findMaxAbs(itsPsfConvInitScales[scale], maxPsfConvInitScales(scale),
      positionPeakPsfConvInitScales);
     }
  } //End pragma parallel
  for (scale=0; scale < itsNInitScales; scale++)
  {
    if (maxPsfConvInitScales(scale) < 0.0)
    {
      os << "As Peak of PSF is negative, you should change the initial scales with a smaller scale size"
   << LogIO::SEVERE;
      return;
    }
  }*/
  //genie

  //#pragma omp parallel default(shared) private(scale) num_threads(nth)
  //{
  //  #pragma omp for // genie pragma seems to sometimes return wrong value to maxima on tesla
    /*float maxvalue;
    IPosition peakpos;
    findMaxAbs((*itsDirty), maxvalue, peakpos);
    cout << "orig itsDirty pos " << peakpos << " maxval " << (*itsDirty)(peakpos) << endl;*/
    Float maxVal=0;
    IPosition posmin((*itsDirty).shape().nelements(), 0);
    Float minVal=0;
    IPosition posmax((*itsDirty).shape().nelements(), 0);
    minMaxMasked(minVal, maxVal, posmin, posmax, (*itsDirty), itsInitScaleMasks[0]);
    cout << "orig itsDirty : min " << minVal << " max " << maxVal << endl;
    cout << "posmin " << posmin << " posmax " << posmax << endl;

    for (scale=0; scale < itsNInitScales; ++scale)
    {
      // Find absolute maximum for the dirty image
      //cout << "in omp loop for scale : " << scale << " : " << blcDirty << " : " << trcDirty << " :: " << itsDirty->shape().nelements() << endl;
      /*Matrix<Float> work = (vecWork_p[scale])(blcDirty, trcDirty);
      work = 0.0;

      //genie debug
      Float maxVal=0;
      IPosition posmin(vecWork_p[scale].shape().nelements(), 0);
      Float minVal=0;
      IPosition posmax(vecWork_p[scale].shape().nelements(), 0);
      minMaxMasked(minVal, maxVal, posmin, posmax, vecWork_p[scale], itsInitScaleMasks[scale]);
      cout << "before InitScaleVal " << scale << ": min " << minVal << " max " << maxVal << endl;
      cout << "posmin " << posmin << " posmax " << posmax << endl;
      //genie

      work = work + (itsDirtyConvInitScales[scale])(blcDirty, trcDirty);*/
      maxima(scale) = 0;
      posMaximum[scale] = IPosition(itsDirty->shape().nelements(), 0);
      //cout << "makedirtyinitscale mask shape: " << itsInitScaleMasks[scale].shape() << endl;

      //genie debug
      /*maxVal=0;
      posmin = IPosition(vecWork_p[scale].shape().nelements(), 0);
      minVal=0;
      posmax = IPosition(vecWork_p[scale].shape().nelements(), 0);
      minMaxMasked(minVal, maxVal, posmin, posmax, vecWork_p[scale], itsInitScaleMasks[scale]);
      cout << "after InitScaleVal " << scale << ": min " << minVal << " max " << maxVal << endl;
      cout << "posmin " << posmin << " posmax " << posmax << endl;*/

      Float maxVal=0;
      Float minVal=0;
      IPosition posmin(itsDirty->shape().nelements(), 0);
      IPosition posmax(itsDirty->shape().nelements(), 0);
      minMaxMasked(minVal, maxVal, posmin, posmax, itsDirtyConvInitScales[scale], itsInitScaleMasks[scale]);
      cout << "DirtyConvInitScale " << scale << ": min " << minVal << " max " << maxVal << endl;
      cout << "posmin " << posmin << " posmax " << posmax << endl;
      //genie

      if (!itsMask.null())
      {
        //findMaxAbsMask(vecWork_p[scale], itsInitScaleMasks[scale],
        //  maxima(scale), posMaximum[scale]);
        findMaxAbsMask(itsDirtyConvInitScales[scale], itsInitScaleMasks[scale],
          maxima(scale), posMaximum[scale]);
      }
      else
        //findMaxAbs(vecWork_p[scale], maxima(scale), posMaximum[scale]);
        findMaxAbs(itsDirtyConvInitScales[scale], maxima(scale), posMaximum[scale]);

      if (itsNormMethod == 2)
      {
        if (scale > 0)
        {
  	      float normalization;
  	      //normalization = 2 * M_PI / pow(itsInitScaleSizes[scale], 2); //sanjay's
  	      //normalization = 2 * M_PI / pow(itsInitScaleSizes[scale], 1); // this looks good on M31 but bad on G55
  	      //normalization = sqrt(2 * M_PI) / pow(itsInitScaleSizes[scale], 1); // M31 new asp gold
          normalization = sqrt(2 * M_PI *itsInitScaleSizes[scale]); //GSL trial. Need to recover and re-norm later
  	      maxima(scale) /= normalization;
  	      //cout << "normalization[" << scale << "] " << normalization << endl;
        } //sanjay's code doesn't normalize peak here, but this works well with GSL with derivatives
      }
      /*maxima(scale) /= maxPsfConvInitScales(scale); // mimic MS-clean. This seems wrong b/c Asp doesn't use scaleBias
      cout << "maxPsfconvinitscale[" << scale << "] = " << maxPsfConvInitScales(scale) << endl;
      cout << "after: maxima[" << scale << "] = " << maxima(scale) << endl;*/
    }
  //}//End parallel section

  // Find the peak residual among the 4 initial scales, which will be the next Aspen
  for (scale = 0; scale < itsNInitScales; scale++)
  {
    if(abs(maxima(scale)) > abs(strengthOptimum)) //genie, bug was from comparing to itsStrengthOptimum
    //if(abs(maxima(scale)) > abs(itsStrengthOptimum))
    {
      optimumScale = scale;
      strengthOptimum = maxima(scale);
      positionOptimum = posMaximum[scale];
    }
  }

  if (optimumScale > 0)
  {
    //const float normalization = 2 * M_PI / (pow(1.0/itsPsfWidth, 2) + pow(1.0/itsInitScaleSizes[optimumScale], 2)); // sanjay
    const float normalization = sqrt(2 * M_PI / (pow(1.0/itsPsfWidth, 2) + pow(1.0/itsInitScaleSizes[optimumScale], 2))); // this is good. Seems it can be in the loop as well.

    // norm method 2 recovers the optimal strength and then normalize it to get the init guess
    if (itsNormMethod == 2)
      strengthOptimum *= sqrt(2 * M_PI *itsInitScaleSizes[optimumScale]); // this is needed if we also first normalize and then compare.

    strengthOptimum /= normalization;
    cout << "normalization " << normalization << " strengthOptimum " << strengthOptimum << endl;
  } //genie's

  AlwaysAssert(optimumScale < itsNInitScales, AipsError);
}

// GSL
Float AspMatrixCleaner::isGoodAspen(Float amp, Float scale, IPosition center)
{
    const int nX = itsDirty->shape()(0);
    const int nY = itsDirty->shape()(1);
    const int refi = nX/2;
    const int refj = nY/2;
    const double sigma5 = 5 * scale / 2;
    const int minI = std::max(0, (int)(center[0] - sigma5));
    const int maxI = std::min(nX-1, (int)(center[0] + sigma5));
    const int minJ = std::max(0, (int)(center[1] - sigma5));
    const int maxJ = std::min(nY-1, (int)(center[1] + sigma5));

    casacore::Matrix<casacore::Float> Asp(nX, nY);
    Asp = 0.0;
    casacore::Matrix<casacore::Float> dAsp(nX, nY);
    dAsp = 0.0;

    for (int j = 0; j <= nY-1; j++)
    {
      for (int i = 0; i <= nX-1; i++)
      {
        const int px = i;
        const int py = j;

        Asp(i,j) = (1.0/(sqrt(2*M_PI)*fabs(scale)))*exp(-(pow(i-center[0],2) + pow(j-center[1],2))*0.5/pow(scale,2));
        dAsp(i,j)= Asp(i,j) * (((pow(i-center[0],2) + pow(j-center[1],2)) / pow(scale,2) - 1) / fabs(scale));
      }
    }

    casacore::Matrix<casacore::Complex> AspFT;
    casacore::FFTServer<casacore::Float,casacore::Complex> fft(itsDirty->shape());
    fft.fft0(AspFT, Asp);

    casacore::Matrix<casacore::Complex> cWork;
    cWork = AspFT * (*itsXfr);
    casacore::Matrix<casacore::Float> AspConvPsf(itsDirty->shape(), (casacore::Float)0.0);
    fft.fft0(AspConvPsf, cWork, false);
    fft.flip(AspConvPsf, false, false); //need this

    casacore::Matrix<casacore::Complex> dAspFT;
    fft.fft0(dAspFT, dAsp);
    casacore::Matrix<casacore::Complex> dcWork;
    dcWork = dAspFT * (*itsXfr);
    casacore::Matrix<casacore::Float> dAspConvPsf(itsDirty->shape(), (casacore::Float)0.0);
    fft.fft0(dAspConvPsf, dcWork, false);
    fft.flip(dAspConvPsf, false, false); //need this

    casacore::Matrix<casacore::Float> GradAmp(itsDirty->shape(), (casacore::Float)0.0);
    casacore::Matrix<casacore::Float> GradScale(itsDirty->shape(), (casacore::Float)0.0);

    // calculate the length of the direvative vector
    Float lenDirVec = 0.0;
    /*for (int j = minJ; j <= maxJ; j++)
    {
      for (int i = minI; i <= maxI; i++)
      {*/
    for (int j = 0; j <= nY-1; j++)
    {
      for (int i = 0; i <= nX-1; i++)
      {
        // generate derivatives of amplitude
        GradAmp(i,j) = (-2) * (*itsDirty)(i,j) * AspConvPsf(i,j);
        // generate derivative of scale
        GradScale(i,j) = (-2) * amp * (*itsDirty)(i,j) * dAspConvPsf(i,j);

        lenDirVec += sqrt(pow(GradAmp(i,j), 2) + pow(GradScale(i,j), 2));
      }
    }

    return lenDirVec;
}

vector<Float> AspMatrixCleaner::getActiveSetAspen()
{
  LogIO os(LogOrigin("AspMatrixCleaner", "getActiveSetAspen()", WHERE));

  if(int(itsInitScaleXfrs.nelements()) == 0)
    throw(AipsError("Initial scales for Asp are not defined"));

  /*if (itsSwitchedToMS)
  	return itsGoodAspActiveSet;*/

  /*if (itsSwitchedToHogbom)
    return {};

  if (!itsSwitchedToHogbom &&
  	  itsNumIterNoGoodAspen.size()>= 10 &&
  	  accumulate(itsNumIterNoGoodAspen.begin(), itsNumIterNoGoodAspen.end(), 0) >= 5)
  {
    switchedToHogbom();
    return {};
  } */

  if (!itsSwitchedToHogbom &&
  	  accumulate(itsNumIterNoGoodAspen.begin(), itsNumIterNoGoodAspen.end(), 0) >= 5)
  {
  	cout << "Switched to hogbom because of frequent small components." << endl;
    switchedToHogbom();
  }

  if (itsSwitchedToHogbom)
  	itsNInitScales = 1;
  else
  	itsNInitScales = itsInitScaleSizes.size();

  // Dirty * initial scales
  Matrix<Complex> dirtyFT;
  FFTServer<Float,Complex> fft(itsDirty->shape());
  fft.fft0(dirtyFT, *itsDirty);
  itsDirtyConvInitScales.resize(0);
  itsDirtyConvInitScales.resize(itsNInitScales); // 0, 1.5width, 5width and 10width

  for (int scale=0; scale < itsNInitScales; scale++)
  {
    Matrix<Complex> cWork;
    //cout << "scale " << scale << " itsInitScaleptr " << &(itsInitScaleXfrs[scale]) << endl;

    itsDirtyConvInitScales[scale] = Matrix<Float>(itsDirty->shape());
    cWork=((dirtyFT)*(itsInitScaleXfrs[scale]));
    fft.fft0((itsDirtyConvInitScales[scale]), cWork, false);
    fft.flip((itsDirtyConvInitScales[scale]), false, false);

    cout << "remake itsDirtyConvInitScales " << scale << " max itsInitScales[" << scale << "] = " << max(fabs(itsInitScales[scale])) << endl;
    cout << " max itsInitScaleXfrs[" << scale << "] = " << max(fabs(itsInitScaleXfrs[scale])) << endl;
  }

  float strengthOptimum = 0.0;
  int optimumScale = 0;
  IPosition positionOptimum(itsDirty->shape().nelements(), 0);
  itsGoodAspActiveSet.resize(0);
  itsGoodAspAmplitude.resize(0);
  itsGoodAspCenter.resize(0);
  itsPrevAspActiveSet.resize(0);
  itsPrevAspAmplitude.resize(0);

  maxDirtyConvInitScales(strengthOptimum, optimumScale, positionOptimum);
  cout << "Initial maximum residual is " << strengthOptimum;
  cout << " its itsDirty is " << (*itsDirty)(positionOptimum);
  cout << " at location " << positionOptimum[0] << " " << positionOptimum[1];
  cout << " " << positionOptimum[2] << " and scale: " << optimumScale << endl;

  // memory used
  //itsUsedMemoryMB = double(HostInfo::memoryUsed()/1024);
  //cout << "Memory allocated in getActiveSetAspen " << itsUsedMemoryMB << " MB." << endl;

  itsStrengthOptimum = strengthOptimum;
  itsPositionOptimum = positionOptimum;
  itsOptimumScale = optimumScale;
  itsOptimumScaleSize = itsInitScaleSizes[optimumScale];

  // initial scale size = 0 gives the peak res, so we don't
  // need to do the LBFGS optimization for it
  //if (itsNumIterNoGoodAspen.size() >= 10)
  //	itsNumIterNoGoodAspen.pop_front(); // only track the past 10 iters
  if (itsOptimumScale == 0)
  {
  	//if (!itsSwitchedToHogbom)
      //itsNumIterNoGoodAspen.push_back(1); //genie Zhang 2018 fused-Asp approach, handled by aspclean instead
    //return itsAspScaleSizes;
    return {};
  }
  else
  {

    //genie:
    // grab the existing Aspen from class variables, itsAspAmp and itsAspScale
    // and also add the new Aspen
    // abs(strengthOptimum), itsInitScaleSizes[optimumScale] to x
    // Also, push back the new center (positionOptimum) to itsAspCenter
    AlwaysAssert(itsAspScaleSizes.size() == itsAspAmplitude.size(), AipsError);
    AlwaysAssert(itsAspScaleSizes.size() == itsAspCenter.size(), AipsError);

    // avoid duplicated (i.e. same location, Amp and size) Aspen for speed up
    // i.e. only add the unique apsen to the permanent list
    /*if (!(find(itsAspAmplitude.begin(), itsAspAmplitude.end(), strengthOptimum) != itsAspAmplitude.end() &&
        find(itsAspScaleSizes.begin(), itsAspScaleSizes.end(), itsInitScaleSizes[optimumScale]) != itsAspScaleSizes.end() &&
        find(itsAspCenter.begin(), itsAspCenter.end(), positionOptimum) != itsAspCenter.end()))
    {
      itsAspAmplitude.push_back(strengthOptimum);
      itsAspScaleSizes.push_back(itsInitScaleSizes[optimumScale]);
      itsAspCenter.push_back(positionOptimum);
      itsAspGood.push_back(true); // true until compare with threshold
    }*/

    // test scale up genie
    /*for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
      itsAspAmplitude[i] = itsAspAmplitude[i] * 1e8;*/


    // heuristiclly determine active set for speed up
    Float resArea = 0.0;
    Int nX = itsDirty->shape()(0);
    Int nY = itsDirty->shape()(1);

    for (Int j = 0; j < nY; ++j)
    {
      for(Int i = 0; i < nX; ++i)
        resArea += abs((*itsDirty)(i,j));
    }
    //const Float lamda = 2e4; //M31
    //const Float lamda = 5e3; //M31 gold2
    //const Float lamda = 7e3; //M31 all fixed
    //const Float lamda = 200; //M31 new Asp
    const Float lamda = 318; //M31 new Asp - gold

    //const Float lamda = 13000.0; //G55
    //const Float lamda = 1120.0; //G55 tesla, 1e8, spw2
    //const Float lamda = 520.0; //G55 tesla, 1e8, spw3 used to be good
    //const Float lamda = 450.0; //G55 tesla, 1e8, spw3 nScales3
    //const Float lamda = 490.0; //G55 tesla, 1e8, spw3 nScales4/5/6/7/8/gold
    //const Float lamda = 485.8; // G55 tesla, 1e8, spw3 fixedScale, eigen
    //const Float lamda = 9600000; // G55 tesla, 1e8, spw3, new asp, old norm
    //const Float lamda = 5000000; // G55 tesla, 1e8, spw3, new asp, old norm
    //const Float lamda = 2400000; // G55 tesla, 1e8, spw3, new asp, SNorm
    //const Float lamda = 2600000; // G55 tesla, 1e8, spw3, new asp, SNorm2
    //const Float lamda = 130000; // G55 local, 1e8, spw3, new asp, permanent list
    //const Float lamda = 5000000; // G55 tesla, 1e8, spw3, new asp, sNorm3-5 m31 norm, gold

    const Float threshold = lamda * resArea;
    //const Float threshold = 30000;
    vector<Float> tempx;
    vector<IPosition> activeSetCenter;

    vector<pair<Float,int>> vp; //(LenDev, idx)
    // don't really need this since we use the latest aspen only
    /*
    for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
    {
    	//if (vp.size() >= 3) // limit the # active set to 3
    	//	break;

      if (!itsAspGood[i])
      	continue;

      //auto start = high_resolution_clock::now();

      const Float lenDirVec = isGoodAspen(itsAspAmplitude[i], itsAspScaleSizes[i], itsAspCenter[i]);
      //auto stop = high_resolution_clock::now();
      //auto duration = duration_cast<microseconds>(stop - start);
      //cout << "isGoodAspen runtime " << duration.count() << " ms" << endl;
      cout << "test: scale " << itsAspScaleSizes[i] << " amp " << itsAspAmplitude[i] << " center " << itsAspCenter[i] << " lenDirVec " << lenDirVec << " threshold " << threshold << endl;
    	//if (lenDirVec >= threshold)
      //if (lenDirVec >= 70000.0)
      //if (lenDirVec >= 200.0)
      if (lenDirVec >= 1e10)
    	{
    		//cout << "good: scale " << itsAspScaleSizes[i] << " amp " << itsAspAmplitude[i] << " center " << itsAspCenter[i] << " lenDirVec " << lenDirVec << " threshold " << threshold << endl;
    		vp.push_back({lenDirVec, i});
	    }
	    else
	    	itsAspGood[i] = false;
    }
    */

    sort(vp.begin(),vp.end(), [](const pair<Float,int> &l, const pair<Float,int> &r) {return l.first > r.first;});

    // select the top 5
    vector<int> goodvp;
    for (unsigned int i = 0; i < vp.size(); i++)
    {
      //if (i >= 5)
      if (i >= 20)
      //if (i >= 400)
        break;
      goodvp.push_back(vp[i].second);
    }
    sort(goodvp.begin(), goodvp.end(), [](const int &l, const int &r) {return l > r;});

    for (unsigned int i = 0; i < goodvp.size(); i++)
    {
      tempx.push_back(itsAspAmplitude[goodvp[i]]);
      tempx.push_back(itsAspScaleSizes[goodvp[i]]);
      activeSetCenter.push_back(itsAspCenter[goodvp[i]]);
      //cout << "temp erase aspen " << goodvp[i] << " AspScaleSize " << itsAspScaleSizes[goodvp[i]] << endl;
      itsAspAmplitude.erase(itsAspAmplitude.begin() + goodvp[i]);
      itsAspScaleSizes.erase(itsAspScaleSizes.begin() + goodvp[i]);
      itsAspCenter.erase(itsAspCenter.begin() + goodvp[i]);
      itsAspGood.erase(itsAspGood.begin() + goodvp[i]);
    }

    // the new aspen is always added to the active-set
    // remember to scale up the strength for G55
    //tempx.push_back(strengthOptimum * 1e8); //G55
    tempx.push_back(strengthOptimum); // M31
    tempx.push_back(itsInitScaleSizes[optimumScale]);
    activeSetCenter.push_back(positionOptimum);

    unsigned int length = tempx.size();

    // trial cppsolvers
    //AspObjFunc<float>::TVector x(length);
    //end trial

    // for lbfgs
    /*
    VectorXd x(length);
    // Bounds
    VectorXd lb(length);
    VectorXd ub(length);
    for (unsigned int i = 0; i < length; i+=2)
    {
      x[i] = tempx[i]; //Eigen::VectorXd needs to be assigned in this way
      x[i+1] = tempx[i+1];

      //lower bound for amp and scale
      lb[i] = -1000.0; // lower bound amp std::numeric_limits<double>::lowest();
      lb[i+1] = 1.0; // lower bound scale
      //upper bound for amp and scale
      ub[i] = 1000.0; // upper bound scale std::numeric_limits<double>::max();
      ub[i+1] = itsPsfWidth * 50; // upper bound scale

      // save aspen before optimization
      itsPrevAspAmplitude.push_back(tempx[i]); // M31 active-set amplitude before lbfgs
      //itsPrevAspAmplitude.push_back(tempx[i] / 1e8); // G55 active-set amplitude before lbfgs
      itsPrevAspActiveSet.push_back(tempx[i+1]); // prev active-set before lbfgs
    }

    cout << "Before: x = " << x.transpose() << endl;*/

    // GSL: set the initial guess
    gsl_vector *x = NULL;
    x = gsl_vector_alloc(length);
    gsl_vector_set_zero(x);

    for (unsigned int i = 0; i < length; i+=2)
    {
      gsl_vector_set(x, i, tempx[i]);
      gsl_vector_set(x, i+1, tempx[i+1]);

      // save aspen before optimization
      itsPrevAspAmplitude.push_back(tempx[i]); // M31 active-set amplitude before lbfgs
      //itsPrevAspAmplitude.push_back(tempx[i] / 1e8); // G55 active-set amplitude before lbfgs
      itsPrevAspActiveSet.push_back(tempx[i+1]); // prev active-set before lbfgs
    }

    // GSL optimization
    //fdf
    gsl_multimin_function_fdf my_func;
    gsl_multimin_fdfminimizer *s = NULL;
    // f only
    /*gsl_multimin_function my_func;
    gsl_multimin_fminimizer *s = NULL;*/

    // setupSolver
    ParamObj optParam(*itsDirty, *itsXfr, activeSetCenter);
    //ParamObj optParam(*itsOrigDirty, *itsXfr, activeSetCenter); // sanjay uses orig dirty
    ParamObj *ptrParam;
    ptrParam = &optParam;
    // fdf
    const gsl_multimin_fdfminimizer_type *T;
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T, length);
    my_func.n      = length;
    my_func.f      = my_f;
    my_func.df     = my_df;
    my_func.fdf    = my_fdf;
    my_func.params = (void *)ptrParam;
    // f only
    /*const gsl_multimin_fminimizer_type *T;
    T = gsl_multimin_fminimizer_nmsimplex2rand; //20.77, 31.63
    //T = gsl_multimin_fminimizer_nmsimplex2; //-1.72 313 -->bad
    s = gsl_multimin_fminimizer_alloc(T, length);
    my_func.n      = length;
    my_func.f      = my_f; //my_f
    my_func.params = (void *)ptrParam;*/

    // fdf
    const float InitStep = gsl_blas_dnrm2(x);
    gsl_multimin_fdfminimizer_set(s, &my_func, x, InitStep, 0.1/*1e-3*/);
    // f only
    /*gsl_vector *ss = NULL;
    ss = gsl_vector_alloc (length);
    gsl_vector_set_all (ss, gsl_blas_dnrm2(x));
    gsl_multimin_fminimizer_set(s, &my_func, x, ss);*/

    printf("\n---------- BFGS algorithm begin ----------\n");
    // fdf
    findComponent(5, s);

    // f only
    /*size_t iter = 0;
    int status = 0;
    double size;
    do
      {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
          break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);

        if (status == GSL_SUCCESS)
          {
            printf ("converged to minimum at\n");
          }

        printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
                iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                s->fval, size);
      }
    while (status == GSL_CONTINUE && iter < 20);*/
    printf("\n----------  BFGS algorithm end  ----------\n");
    // update x needed here?
    gsl_vector *optx = NULL;
    optx = gsl_multimin_fdfminimizer_x(s); //fdf
    //optx = gsl_multimin_fminimizer_x(s); // f only

    // end GSL optimization

    // trial cppsolvers
    /*AspObjFunc<float> fun(*itsOrigDirty, *itsXfr, activeSetCenter, length);
    bool probably_correct = fun.checkGradient(x);
    LbfgsbSolver<AspObjFunc<float>> solver;
    // change stopping criteria
    Criteria<float> crit = Criteria<float>::defaults(); // Create a Criteria class to set the solver's stop conditions
    crit.iterations = 10;
    crit.gradNorm = 5e-1;
    solver.setStopCriteria(crit);
    solver.minimize(fun, x);
    cout << "After: x = " << x.transpose() << endl;*/
    //end trial

    /*
    // lbfgs optimization
    ////LBFGSParam<double> param;
    LBFGSBParam<double> param;
    param.epsilon = 5e-3;
    //param.max_linesearch = 10;
    //param.min_step = 1e-30;
    param.max_iterations = 5; //M31
    //param.max_linesearch = 128;
    param.gclip = itsPrevLBFGSGrad;
    //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE; genie: doesn't work
    ////LBFGSSolver<double> solver(param);
    LBFGSBSolver<double> solver(param); //with bounds can cause amp/scale not updating
    // because scale cannot go negative but sometimes neg scale can later be corrected

    //if (itsUseZhang)
    //AspZhangObjFunc fun(*itsDirty, activeSetCenter); // Asp 2016
    //else
    GaussianObjFunc fun(*itsDirty, *itsXfr, activeSetCenter, itsGain); //Asp 2004
    ////AspObjFunc fun(*itsDirty, *itsXfr, activeSetCenter, itsGain); //Asp 2004

    //auto start = high_resolution_clock::now();

    double fx;
    double gclip;
    ////const int niter = solver.minimize(fun, x, fx, gclip);
    const int niter = solver.minimize(fun, x, fx, lb, ub, gclip); // with bounds

    //auto stop = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(stop - start);
    //cout << "lbfgs runtime " << duration.count() << " ms" << endl;

    // use the initial gradient as a roll back gradient if there is
    // gradient exploding in lbfgs
    if (itsPrevLBFGSGrad == 0.0)
    {
      cout << "lbfgs gradient exploded!" << endl;
      itsPrevLBFGSGrad = gclip;
    }
    //std::cout << "itsPrevLBFGSGrad " << itsPrevLBFGSGrad << std::endl;
    cout << "After: x = " << x.transpose() / *<< ", # iters " << niter * /<< endl;*/

    // x is not changing in LBFGS which causes convergence issue
    // so we use the previous x instead. This is a hack suggested by
    // the online community
    //bool becomesNegScale = false;
    /*if (fx == -999.0)
    { cout << "fx == -999; has convergence issue? Reset x to the previous one." << endl;
      becomesNegScale = true;
      for (unsigned int i = 0; i < length; i++)
        x[i] = tempx[i];
    }*/

    //prevent lbfgs overadjusting
    /*bool overflow = false;
    for (unsigned int i = 0; i < length; i+= 2)
    {
    	if ((fabs(x[i+1]) / fabs(tempx[i+1])) >= 2.0 ||
    		  (fabs(x[i]) / fabs(tempx[i])) >= 2.0     ||
    		  x[i+1] <= 1)
    	{
    		overflow = true;
    		cout << "overflow by: " << tempx[i+1] << " -> " << x[i+1] << " or " << tempx[i] << " -> " << x[i]  << endl;
    		break;
    	}
    }
    if (overflow)
    {
    	cout << "lbfgs overflow!!! Switched to hogbom" << endl;
    	switchedToHogbom();
    	for (unsigned int i = 0; i < length; i+= 2)
	    {
	    	x[i+1] = tempx[i+1];
	    	x[i] = tempx[i];
	    }
    }*/

    // put the updated x back to the class variables, itsAspAmp and itsAspScale
    for (unsigned int i = 0; i < length; i+= 2)
    {
      double amp = gsl_vector_get(optx, i);
      double scale = gsl_vector_get(optx, i+1);
      //scale = (scale = abs(scale)) < 0.4 ? 0.4 : scale;
      scale = (scale = fabs(scale)) < 0.4 ? 0 : scale;

      itsAspAmplitude.push_back(amp);
      itsAspScaleSizes.push_back(scale); //permanent list that doesn't get clear
      itsGoodAspAmplitude.push_back(amp); // M31 active-set amplitude
      //itsGoodAspAmplitude.push_back(x[i] / 1e8); // G55 active-set amplitude
      itsGoodAspActiveSet.push_back(scale); // active-set
      itsAspCenter.push_back(activeSetCenter[i/2]);
      if (scale < 0.4)
        itsAspGood.push_back(false);
      else
        itsAspGood.push_back(true);

      cout << "optx = " << amp << " " << scale << endl;
    }

    itsStrengthOptimum = itsAspAmplitude[itsAspAmplitude.size() -1]; //gsl_vector_get(optx, length - 2); // x[length - 2]; // the latest aspen is the last element of x
    // the above line affects residual image (brighter) (0 scales so far seems problemastic)
    itsOptimumScaleSize = itsAspScaleSizes[itsAspScaleSizes.size() -1]; //gsl_vector_get(optx, length - 1); //x[length - 1]; // the latest aspen is the last element of x
    itsGoodAspCenter = activeSetCenter;

    // debug
    cout << "opt strengthOptimum " << itsStrengthOptimum << " opt size " << itsOptimumScaleSize << endl;

    // free GSL stuff
    gsl_multimin_fdfminimizer_free(s); //fdf
    //gsl_multimin_fminimizer_free(s); // f only
    gsl_vector_free(x);
    //gsl_vector_free(optx); // free this causes seg fault!!!

  } // finish bfgs optimization

  AlwaysAssert(itsGoodAspCenter.size() == itsGoodAspActiveSet.size(), AipsError);

  // debug info
  /*for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
  {
    //cout << "After opt AspApm[" << i << "] = " << itsAspAmplitude[i] << endl;
    //cout << "After opt AspScale[" << i << "] = " << itsAspScaleSizes[i] << endl;
    //cout << "After opt AspCenter[" << i << "] = " << itsAspCenter[i] << endl;
    cout << "AspScale[ " << i << " ] = " << itsAspScaleSizes[i] << " center " << itsAspCenter[i] << endl;
  }*/
  // test scale back genie
  /*for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
    itsAspAmplitude[i] = itsAspAmplitude[i] / 1e8;*/

  //return itsAspScaleSizes; // this returns permanent list
  return itsGoodAspActiveSet; // return optimized scale
}

// Define the Asp scales without doing anything else
// user will call make makePsfScales and makeDirtyScales like an adult in the know
void AspMatrixCleaner::defineAspScales(vector<Float>& scaleSizes)
{
  /*if(itsScales.nelements()>0) {
    destroyAspScales();
  }

  destroyMasks();*/ //genie do I need this? prob not

  // try speed up here by removing duplicated scale sizes
  /*itsNscales = scaleSizes.nelements();
  itsScaleSizes.resize(itsNscales);
  itsScaleSizes = scaleSizes;  // make a copy that we can call our own
  GenSort<Float>::sort(itsScaleSizes);*/

  // this is good but needs to be reconsidered when we need to use itsGoodAspAmplitude
  sort(scaleSizes.begin(), scaleSizes.end());
  //scaleSizes.erase(unique(scaleSizes.begin(), scaleSizes.end()), scaleSizes.end()); // remove the duplicated scales
  scaleSizes.erase(unique(scaleSizes.begin(),scaleSizes.end(),[](Float l, Float r) { return abs(l - r) < 1e-3; }), scaleSizes.end());

  for (unsigned int i = 0; i < scaleSizes.size(); i++)
    cout << "defineAspScales scaleSizes[" << i << "] = " << scaleSizes[i] << endl;

  itsNscales = Int(scaleSizes.size());
  itsScaleSizes.resize(itsNscales);
  itsScaleSizes = Vector<Float>(scaleSizes);  // make a copy that we can call our own

  // analytically calculate component scale by Asp 2016
  if (itsUseZhang)
  {
    for (Int i = 0; i < itsNscales; i++)
    {
      if (itsScaleSizes[i] >= itsPsfWidth)
        itsScaleSizes[i] = sqrt(pow(itsScaleSizes[i], 2) - pow(Float(itsPsfWidth), 2));
    }
  }
  // end Asp 2016

  itsScalesValid = true;  //genie? It's false in MS clean
}

void AspMatrixCleaner::switchedToHogbom()
{
	itsSwitchedToHogbom = true;
	//itsSwitchedToHogbom = false;
  itsNthHogbom += 1;
  itsNumIterNoGoodAspen.resize(0);
  //itsNumHogbomIter = ceil(100 + 50 * (exp(0.05*itsNthHogbom) - 1)); //zhang's formula
  //itsNumHogbomIter = ceil(50 + 200 * (exp(0.05*itsNthHogbom) - 1)); //genie's formula, SNorm1-3
  itsNumHogbomIter = ceil(50 + 2 * (exp(0.05*itsNthHogbom) - 1)); //genie's formula, SNorm1-3
  //itsNumHogbomIter = ceil(65 + 200 * (exp(0.05*itsNthHogbom) - 1)); //genie's formula, SNorm4 5k
  //itsNumHogbomIter = ceil(100 + 50 * (exp(0.05*itsNthHogbom) - 1)); //SNorm5
  //itsNumHogbomIter = 1000000; //fake it for M31
  cout << "Run hogbom for " << itsNumHogbomIter << " iterations." << endl;
}

void AspMatrixCleaner::setOrigDirty(const Matrix<Float>& dirty){
  itsOrigDirty=new Matrix<Float>(dirty.shape());
  itsOrigDirty->assign(dirty);
}


void AspMatrixCleaner::setMaxPsfConvInitScales()
{
  itsPsfConvInitScales.resize(itsNInitScales, false);
  FFTServer<Float,Complex> fft(psfShape_p);
  Matrix<Complex> cWork;

  for (Int scale=0; scale < itsNInitScales; scale++)
  {
    //cout << "Calculating convolutions of psf for initial scale size " << itsInitScaleSizes[scale] << endl;
    //PSF * scale
    itsPsfConvInitScales[scale] = Matrix<Float>(psfShape_p);
    cWork=((*itsXfr)*(itsInitScaleXfrs[scale])*(itsInitScaleXfrs[scale]));
    fft.fft0((itsPsfConvInitScales[scale]), cWork, false);
    fft.flip(itsPsfConvInitScales[scale], false, false);
  }

  // Find the peaks of the convolved Psfs
  maxPsfConvInitScales.resize(0);
  maxPsfConvInitScales.resize(itsNInitScales);
  Int naxes = psfShape_p.nelements();

  for (int scale=0; scale < itsNInitScales; scale++)
  {
    IPosition positionPeakPsfConvInitScales(naxes, 0);

    findMaxAbs(itsPsfConvInitScales[scale], maxPsfConvInitScales(scale),
    positionPeakPsfConvInitScales);
  }

  for (int scale=0; scale < itsNInitScales; scale++)
  {
    if (maxPsfConvInitScales(scale) < 0.0)
    {
      cout << "As Peak of PSF is negative, you should change the initial scales with a smaller scale size"
   << endl;
      return;
    }
  }
}

void AspMatrixCleaner::testBFGS(const Matrix<Float>& psf)
{
  const int N = 512;
  helper();
  Matrix<Float> simpleDirty(N, N, 0.0);
  Matrix<Float> target(N, N, 0.0);
  Matrix<Float> simplePsf(N, N, 0.0);
  Matrix<Complex> simplePsfFT;
  const Float scaleSize = 20.0;
  //IPosition simplePos(2, N/2, N/2);
  IPosition simplePos(2, 261, 289);
  vector<IPosition> simCenter;
  Matrix<Float> msDirty(N, N, 0.0);
  msDirty = *itsOrigDirty;

  simplePsf(N/2, N/2) = 1.0;
  FFTServer<Float,Complex> fft(simplePsf.shape());
  fft.fft0(simplePsfFT, simplePsf);

  Gaussian2D<Float> gbeam(1.0/(sqrt(2*M_PI)*scaleSize), N/2, N/2, scaleSize, 1, 0);

  for (int j = 0; j < N; j++)
  {
    for (int i = 0; i < N; i++)
    {
      const int px = i;
      const int py = j;
      simpleDirty(i,j) = 3 * gbeam(px, py);
      target(i,j) = 3 * (1.0/(sqrt(2*M_PI)*scaleSize))*exp(-(pow(i-N/2,2) + pow(j-N/2,2))*0.5/pow(scaleSize,2));
      //target(i,j) = 3 * (1.0/(2*M_PI*pow(scaleSize,2)))*exp(-(pow(i-N/2,2) + pow(j-N/2,2))*0.5/pow(scaleSize,2));
    }
  }
  // create complicated dirty
  Matrix<Complex> DirtyFT;
  FFTServer<Float,Complex> fftd(target.shape());
  fftd.fft0(DirtyFT, target);
  Matrix<Complex> cWork;
  cWork = ((*itsXfr) * (DirtyFT));
  Matrix<Float> compDirty(target.shape(), (Float)0.0);
  fftd.fft0(compDirty, cWork, false);
  fftd.flip(compDirty, false, false); //need this

  /*gsl_multimin_function_fdf my_func;
  gsl_multimin_fdfminimizer *s = NULL;*/ // fdf
  gsl_multimin_function my_func;
  gsl_multimin_fminimizer *s = NULL;

  // this block should be equivalent of setupSolver
  simCenter.push_back(simplePos);
  //ParamObj optParam(compDirty, *itsXfr, simCenter);
  ParamObj optParam(msDirty, *itsXfr, simCenter);
  ParamObj *ptrParam;
  ptrParam = &optParam;
  /*const gsl_multimin_fdfminimizer_type *T;
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc(T, 2);
  my_func.n      = 2;
  my_func.f      = my_f; //my_f
  my_func.df     = my_df; //my_df
  my_func.fdf    = my_fdf; // my_fdf
  my_func.params = (void *)ptrParam;*/ //fdf
  const gsl_multimin_fminimizer_type *T;
  //T = gsl_multimin_fminimizer_nmsimplex2rand; //20.77, 31.63
  T = gsl_multimin_fminimizer_nmsimplex2; //-1.72 313
  s = gsl_multimin_fminimizer_alloc(T, 2);
  my_func.n      = 2;
  my_func.f      = my_f; //my_f
  my_func.params = (void *)ptrParam;
  cout << "msDirty " << msDirty(256,256) << endl;
  cout << "*itsXfr " << (*itsXfr)(256,256) << endl;


  /*setupSolver(&s,
      gsl_multimin_fdfminimizer_vector_bfgs2,
      &my_func,
      1, // # aspen
      par,
      &compDirty,
      &simplePsfFT);*/

  // Set the initial guess
  gsl_vector *x = NULL;
  x = gsl_vector_alloc(2);
  gsl_vector_set_zero(x);
  //gsl_vector_set(x, 0, 1.0);
  //gsl_vector_set(x, 1, 16.0);
  gsl_vector_set(x, 0, 115.677);
  gsl_vector_set(x, 1, 20.0);

  /*const float InitStep = gsl_blas_dnrm2(x);
  gsl_multimin_fdfminimizer_set(s, &my_func, x, InitStep, 1e-3);*/ //fdf
  gsl_vector *ss = NULL;
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, gsl_blas_dnrm2(x));
  gsl_multimin_fminimizer_set(s, &my_func, x, ss);

  printf("\n---------- BFGS algorithm begin ----------\n");
  //findComponent(10, s);
  size_t iter = 0;
  int status = 0;
  double size;
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
              iter,
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 20);
  printf("\n----------  BFGS algorithm end  ----------\n");

  //gsl_multimin_fdfminimizer_free(s);
  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
}

} //# NAMESPACE CASA - END
