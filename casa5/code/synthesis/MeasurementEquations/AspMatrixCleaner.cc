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
#include <synthesis/MeasurementEquations/gslobjfunc.h>

// for gsl
using Eigen::VectorXd;

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
  os << LogIO::NORMAL1 << "Asp clean algorithm" << LogIO::POST;


  Int scale;

  AlwaysAssert(itsScalesValid, AipsError);
  //no need to use all cores if possible
  Int nth = itsNscales;
#ifdef _OPENMP

    nth = min(nth, omp_get_max_threads());

#endif

  // Define a subregion for the inner quarter. No longer used
  /*IPosition blcDirty(model.shape().nelements(), 0);
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
  LCBox centerBox(blcDirty, trcDirty, model.shape());*/

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
  Matrix<Complex>itsScaleXfr0 = Matrix<Complex> ();

  Matrix<Float> itsScale = Matrix<Float>(psfShape_p);
  Matrix<Complex>itsScaleXfr = Matrix<Complex> ();

  // Define a subregion so that the peak is centered
  IPosition support(model.shape());
  support(0) = max(Int(itsInitScaleSizes[itsNInitScales-1] + 0.5), support(0));
  support(1) = max(Int(itsInitScaleSizes[itsNInitScales-1] + 0.5), support(1));
  IPosition inc(model.shape().nelements(), 1);

  for (Int ii = itsStartingIter; ii < itsMaxNiter; ii++)
  {
    //cout << "cur iter " << itsIteration << " max iter is "<< itsMaxNiter << endl;
    itsIteration++;

    // make single optimized scale image
    os << "Asp-Clean: making scale " << itsOptimumScaleSize << LogIO::POST;

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
      itsScaleXfr.resize();
      fft.fft0(itsScaleXfr, itsScale);
    }

    // trigger hogbom when itsStrengthOptimum or itsPeakResidual is small enough
    if (itsNormMethod == 1) // only Norm Method 1 needs hogbom for speedup
    {
    	//if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 0.001) // M31 value - new Asp + gaussian
      //if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < 2.8) // M31 value-new asp: 1k->5k
      //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 8e-5 && abs(itsStrengthOptimum) < 1e-4) // G55 
      //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 7e-3 && abs(itsStrengthOptimum) < 1e-7) // Points
      //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 0.138) //GSL M100 channel 22
      //if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 0.15) // GSL M100 channel 22 & 23 & others
      //if (!itsSwitchedToHogbom && (abs(itsPeakResidual) < 2.5 || abs(itsStrengthOptimum) < 1e-3)) // GSL, CygA
      
      if(!itsSwitchedToHogbom && (abs(itsPeakResidual) < itsFusedThreshold 
         || abs(itsStrengthOptimum) < (5e-4 * itsFusedThreshold))) // GSL, CygA. 
      	// 5e-4 is a experimental number here assuming under that threshold itsStrengthOptimum is too small to take affect.  
      {
  	    cout << "Switch to hogbom b/c optimum strength is small enough: " << itsFusedThreshold << endl;
  	    switchedToHogbom();
      }
      /*if (!itsSwitchedToHogbom && abs(itsPeakResidual) < 1.3e-3) // Points 
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
	      itsNumIterNoGoodAspen.push_back(1); // Zhang 2018 fused-Asp approach
	    else
	      itsNumIterNoGoodAspen.push_back(0);
    }

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
    //    1. stop if below threshold. 1e-6 is an experimental number
    if (!itsSwitchedToHogbom && abs(itsStrengthOptimum) < (1e-6 * itsFusedThreshold))
    {
    	//cout << "Reached stopping threshold " << 1e-6 * itsFusedThreshold << " at iteration "<< ii << endl;
      os << "Reached stopping threshold " << 1e-6 * itsFusedThreshold << " at iteration "<<
            ii << LogIO::POST;
      os << "Optimum flux is " << abs(itsStrengthOptimum) << LogIO::POST;
      converged = 1;
      break;
    }
    //    2. negatives on largest scale?
    if ((itsNscales > 1) && itsStopAtLargeScaleNegative &&
    	  itsOptimumScale == (itsNInitScales - 1) &&
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

    
    IPosition blc(itsPositionOptimum - support/2);
    IPosition trc(itsPositionOptimum + support/2 - 1);
    LCBox::verify(blc, trc, inc, model.shape());
    IPosition blcPsf(blc);
    IPosition trcPsf(trc);
    //blcDirty = blc;
    //trcDirty = trc;

    // Update the model image
    Matrix<Float> modelSub = model(blc, trc);
    Float scaleFactor;
    scaleFactor = itsGain * itsStrengthOptimum;   
    Matrix<Float> scaleSub = (itsScale)(blcPsf,trcPsf);
    modelSub += scaleFactor * scaleSub; 

    // Now update the residual image
    // PSF * model
    Matrix<Complex> cWork;
    cWork = ((*itsXfr)*(itsScaleXfr)); //Asp's
    Matrix<Float> itsPsfConvScale = Matrix<Float>(psfShape_p);
    fft.fft0(itsPsfConvScale, cWork, false);
    fft.flip(itsPsfConvScale, false, false); //need this if conv with 1 scale; don't need this if conv with 2 scales
    Matrix<Float> psfSub = (itsPsfConvScale)(blcPsf, trcPsf);
    Matrix<Float> dirtySub=(*itsDirty)(blc,trc);

    /* debug info
    float maxvalue;
    IPosition peakpos;
    findMaxAbs(psfSub, maxvalue, peakpos);
    cout << "psfSub pos " << peakpos << " maxval " << psfSub(peakpos) << endl;
    cout << "dirtySub pos " << peakpos << " val " << dirtySub(peakpos) << endl;
    findMaxAbs(itsPsfConvScale, maxvalue, peakpos);
    cout << "itsPsfConvScale pos " << peakpos << " maxval " << itsPsfConvScale(peakpos) << endl;
    cout << "itsDirty pos " << peakpos << " val " << (*itsDirty)(peakpos) << endl;
    findMaxAbs(dirtySub, maxvalue, peakpos);
    cout << "dirtySub pos " << peakpos << " maxval " << dirtySub(peakpos) << endl;
    findMaxAbs((*itsDirty), maxvalue, peakpos);
    cout << "itsDirty pos " << peakpos << " maxval " << (*itsDirty)(peakpos) << endl;
    cout << "itsPositionOptimum " << itsPositionOptimum << endl;    
    cout << " maxPsfSub " << max(fabs(psfSub)) << " maxPsfConvScale " << max(fabs(itsPsfConvScale)) << " itsGain " << itsGain << endl;*/
    os << "itsStrengthOptimum " << itsStrengthOptimum << LogIO::POST;

    // subtract the peak that we found from the dirty image
    dirtySub -= scaleFactor * psfSub;

    // further update the model and residual with the remaining aspen of the active-set
    // This is no longer needed since we found out using the last Aspen to update model/residual is good enough
    /*if (itsOptimumScaleSize != 0)
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
        IPosition blc1(itsGoodAspCenter[scale] - support/2);
        IPosition trc1(itsGoodAspCenter[scale] + support/2 - 1);
        LCBox::verify(blc1, trc1, inc, model.shape());

        IPosition blcPsf1(blc1);
        IPosition trcPsf1(trc1);

        Matrix<Float> modelSub1 = model(blc1, trc1);
        Matrix<Float> dirtySub1 = (*itsDirty)(blc1,trc1);

        // First remove the previous values of aspen in the active-set
        cout << "aspclean: restore with previous scale " << itsPrevAspActiveSet[scale];
        cout << " amp " << itsPrevAspAmplitude[scale] << endl;

        makeScaleImage(itsScale, itsPrevAspActiveSet[scale], itsPrevAspAmplitude[scale], itsGoodAspCenter[scale]);
        itsScaleXfr.resize();
        fft.fft0(itsScaleXfr, itsScale);
        Matrix<Float> scaleSubPrev = (itsScale)(blcPsf1,trcPsf1);
        const float scaleFactorPrev = itsGain * itsPrevAspAmplitude[scale];
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
        makeScaleImage(itsScale, itsGoodAspActiveSet[scale], itsGoodAspAmplitude[scale], itsGoodAspCenter[scale]);
        itsScaleXfr.resize();
        fft.fft0(itsScaleXfr, itsScale);
        Matrix<Float> scaleSubNew = (itsScale)(blcPsf1,trcPsf1);
        const float scaleFactorNew = itsGain * itsGoodAspAmplitude[scale];
        
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
	    / *if (!itsSwitchedToHogbom && aspenUnchanged)
	    {
	    	cout << "Switched to hogbom b/c aspen is unchanged" << endl;
	    	switchedToHogbom();
	    }* / 
    }*/

    Float maxVal=0;
    IPosition posmin((*itsDirty).shape().nelements(), 0);
    Float minVal=0;
    IPosition posmax((*itsDirty).shape().nelements(), 0);
    minMaxMasked(minVal, maxVal, posmin, posmax, (*itsDirty), itsInitScaleMasks[0]);
    itsPeakResidual = (fabs(maxVal) > fabs(minVal)) ? fabs(maxVal) : fabs(minVal);
    os << "current peakres " << itsPeakResidual << LogIO::POST;
    //cout << "after: itsDirty optPos " << (*itsDirty)(itsPositionOptimum) << endl;

    // If we switch to hogbom (i.e. only have 0 scale size),
    // we still need to do the following Aspen update to get the new optimumStrength
    if (itsSwitchedToHogbom)
    {
      if (itsNumHogbomIter == 0)
      {
        itsSwitchedToHogbom = false;
        os << "switched back to Asp." << LogIO::POST;
      }
      else
        itsNumHogbomIter -= 1; 
    }

    tempScaleSizes.clear();
    tempScaleSizes = getActiveSetAspen();
    tempScaleSizes.push_back(0.0); // put 0 scale
    defineAspScales(tempScaleSizes);
  }
  // End of iteration

  // memory used
  //itsUsedMemoryMB = double(HostInfo::memoryUsed()/1024);
  //cout << "Memory allocated in aspclean " << itsUsedMemoryMB << " MB." << endl;

  if(!converged)
    os << "Failed to reach stopping threshold" << LogIO::POST;

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

// Make a single initial scale size image by Gaussian
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

    /*const Int mini = max(0, (Int)(refi - scaleSize));
    const Int maxi = min(nx-1, (Int)(refi + scaleSize));
    const Int minj = max(0, (Int)(refj - scaleSize));
    const Int maxj = min(ny-1, (Int)(refj + scaleSize));*/
    cout << "makeInitScaleImage: initscalesize " << scaleSize << endl;

    //Gaussian2D<Float> gbeam(1.0/(sqrt(2*M_PI)*scaleSize), 0, 0, scaleSize, 1, 0);

    // has to make the whole scale image
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        const int px = i - refi;
        const int py = j - refj;
        //iscale(i,j) = gbeam(px, py); // gbeam with the above def is equivalent to the following
        iscale(i,j) = (1.0/(sqrt(2*M_PI)*scaleSize))*exp(-(pow(i-refi,2) + pow(j-refj,2))*0.5/pow(scaleSize,2)); //this is for 1D, but represents Sanjay's and gives good init scale
        //iscale(i,j) = (1.0/(2*M_PI*pow(scaleSize,2)))*exp(-(pow(i-refi,2) + pow(j-refj,2))*0.5/pow(scaleSize,2)); // this is for 2D, gives unit area but bad init scale (always picks 0)
      }
    }

  }
}

// Make a single scale size image by Gaussian
void AspMatrixCleaner::makeScaleImage(Matrix<Float>& iscale, const Float& scaleSize, const Float& amp, const IPosition& center)
{

  const Int nx = iscale.shape()(0);
  const Int ny = iscale.shape()(1);
  iscale = 0.0;

  if(scaleSize == 0.0)
    iscale(Int(center[0]), Int(center[1])) = 1.0;
  else
  {
    AlwaysAssert(scaleSize>0.0, AipsError);

    /* const Double refi = nx/2;
    const Double refj = ny/2;
    const Int mini = max(0, (Int)(refi - scaleSize));
    const Int maxi = min(nx-1, (Int)(refi + scaleSize));
    const Int minj = max(0, (Int)(refj - scaleSize));
    const Int maxj = min(ny-1, (Int)(refj + scaleSize));*/
    cout << "makeScaleImage: scalesize " << scaleSize << " center " << center << " amp " << amp << endl;

    ////Gaussian2D<Float> gbeam(1.0 / (sqrt(2*M_PI)*scaleSize), center[0], center[1], scaleSize, 1, 0);

    // all of the following was trying to replicate Sanjay's code but doesn't work
    //const Float d = sqrt(pow(1.0/itsPsfWidth, 2) + pow(1.0/scaleSize, 2));
    //Gaussian2D<Float> gbeam(amp / (sqrt(2*M_PI)*scaleSize), center[0], center[1], scaleSize, 1, 0);
    //Gaussian2D<Float> gbeam(amp * sqrt(2*M_PI)/d, center[0], center[1], scaleSize * d, 1, 0);
    //Gaussian2D<Float> gbeam(amp / (2*M_PI), center[0], center[1], scaleSize, 1, 0);
    //Gaussian2D<Float> gbeam(amp / pow(2,scaleSize-1), center[0], center[1], itsInitScaleSizes[scaleSize], 1, 0);

    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        const int px = i;
        const int py = j;
        //iscale(i,j) = gbeam(px, py); // this is equivalent to the following with the above gbeam definition 
        //this is for 1D, but represents Sanjay's and gives good init scale
        // Note that "amp" is not used in the expression
        iscale(i,j) = (1.0/(sqrt(2*M_PI)*scaleSize))*exp(-(pow(i-center[0],2) + pow(j-center[1],2))*0.5/pow(scaleSize,2));
        
        // this is for 2D, gives unit area but bad init scale (always picks 0)
        //iscale(i,j) = (1.0/(2*M_PI*pow(scaleSize,2)))*exp(-(pow(i-center[0],2) + pow(j-center[1],2))*0.5/pow(scaleSize,2));
      }
    }

  }
}


void AspMatrixCleaner::setInitScaleXfrs(const Float width)
{
  if(itsInitScales.nelements() > 0)
    destroyAspScales();

  // try 0, width, 2width, 4width and 8width
  itsInitScaleSizes.resize(itsNInitScales, false);
  itsInitScaleSizes = {0.0f, width, 2.0f*width, 4.0f*width, 8.0f*width};
  itsInitScales.resize(itsNInitScales, false);
  itsInitScaleXfrs.resize(itsNInitScales, false);
  FFTServer<Float,Complex> fft(psfShape_p);
  for (int scale = 0; scale < itsNInitScales; scale++)
  {
    itsInitScales[scale] = Matrix<Float>(psfShape_p);
    makeInitScaleImage(itsInitScales[scale], itsInitScaleSizes[scale]);
    //cout << "made itsInitScales[" << scale << "] = " << itsInitScaleSizes[scale] << endl;
    itsInitScaleXfrs[scale] = Matrix<Complex> ();
    fft.fft0(itsInitScaleXfrs[scale], itsInitScales[scale]);
  }
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

  /* It seems we don't need the following to define a region. Using minMaxMasked should be sufficient.

  const int nx = itsDirty->shape()[0];
  const int ny = itsDirty->shape()[1];
 
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
    os << LogIO::NORMAL << "Finding initial scales using the entire image" << LogIO::POST;*/

  Vector<Float> maxima(itsNInitScales);
  Block<IPosition> posMaximum(itsNInitScales);
  
  /*int nth = itsNInitScales;
  #ifdef _OPENMP
    nth = min(nth, omp_get_max_threads());
  #endif*/

  //#pragma omp parallel default(shared) private(scale) num_threads(nth)
  //{
  //  #pragma omp for // genie pragma seems to sometimes return wrong value to maxima on tesla
  
    /* debug info
    Float maxVal=0;
    IPosition posmin((*itsDirty).shape().nelements(), 0);
    Float minVal=0;
    IPosition posmax((*itsDirty).shape().nelements(), 0);
    minMaxMasked(minVal, maxVal, posmin, posmax, (*itsDirty), itsInitScaleMasks[0]);
    cout << "orig itsDirty : min " << minVal << " max " << maxVal << endl;
    cout << "posmin " << posmin << " posmax " << posmax << endl; */

    for (int scale = 0; scale < itsNInitScales; ++scale)
    {
      // Find absolute maximum for the dirty image   
      maxima(scale) = 0;
      posMaximum[scale] = IPosition(itsDirty->shape().nelements(), 0);

      /* debug info
      Float maxVal=0;
      Float minVal=0;
      IPosition posmin(itsDirty->shape().nelements(), 0);
      IPosition posmax(itsDirty->shape().nelements(), 0);
      minMaxMasked(minVal, maxVal, posmin, posmax, itsDirtyConvInitScales[scale], itsInitScaleMasks[scale]);
      cout << "DirtyConvInitScale " << scale << ": min " << minVal << " max " << maxVal << endl;
      cout << "posmin " << posmin << " posmax " << posmax << endl; */

      if (!itsMask.null())
      {
        findMaxAbsMask(itsDirtyConvInitScales[scale], itsInitScaleMasks[scale],
          maxima(scale), posMaximum[scale]);
      }
      else
        findMaxAbs(itsDirtyConvInitScales[scale], maxima(scale), posMaximum[scale]);

      if (itsNormMethod == 2)
      {
        if (scale > 0)
        {
  	      float normalization;
  	      //normalization = 2 * M_PI / pow(itsInitScaleSizes[scale], 2); //sanjay's
  	      //normalization = 2 * M_PI / pow(itsInitScaleSizes[scale], 1); // this looks good on M31 but bad on G55
          normalization = sqrt(2 * M_PI *itsInitScaleSizes[scale]); //GSL. Need to recover and re-norm later
  	      maxima(scale) /= normalization;
        } //sanjay's code doesn't normalize peak here. 
         // Norm Method 2 may work fine with GSL with derivatives, but Norm Method 1 is still the preferred approach.
      }
    }
  //}//End parallel section

  // Find the peak residual among the 4 initial scales, which will be the next Aspen
  for (int scale = 0; scale < itsNInitScales; scale++)
  {
    if(abs(maxima(scale)) > abs(strengthOptimum))
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
  } 

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

  if (!itsSwitchedToHogbom &&
  	  accumulate(itsNumIterNoGoodAspen.begin(), itsNumIterNoGoodAspen.end(), 0) >= 5)
  {
  	os << "Switched to hogbom because of frequent small components." << LogIO::POST;
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
  itsDirtyConvInitScales.resize(itsNInitScales); // 0, 1width, 2width, 4width and 8width

  for (int scale=0; scale < itsNInitScales; scale++)
  {
    Matrix<Complex> cWork;

    itsDirtyConvInitScales[scale] = Matrix<Float>(itsDirty->shape());
    cWork=((dirtyFT)*(itsInitScaleXfrs[scale]));
    fft.fft0((itsDirtyConvInitScales[scale]), cWork, false);
    fft.flip((itsDirtyConvInitScales[scale]), false, false);

    // cout << "remake itsDirtyConvInitScales " << scale << " max itsInitScales[" << scale << "] = " << max(fabs(itsInitScales[scale])) << endl;
    // cout << " max itsInitScaleXfrs[" << scale << "] = " << max(fabs(itsInitScaleXfrs[scale])) << endl;
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
  os << "Initial maximum residual is " << strengthOptimum  << " and scale: " << optimumScale << LogIO::POST;
  // cout << " its itsDirty is " << (*itsDirty)(positionOptimum);
  // cout << " at location " << positionOptimum[0] << " " << positionOptimum[1] << " " << positionOptimum[2];

  // memory used
  //itsUsedMemoryMB = double(HostInfo::memoryUsed()/1024);
  //cout << "Memory allocated in getActiveSetAspen " << itsUsedMemoryMB << " MB." << endl;

  itsStrengthOptimum = strengthOptimum;
  itsPositionOptimum = positionOptimum;
  itsOptimumScale = optimumScale;
  itsOptimumScaleSize = itsInitScaleSizes[optimumScale];

  // initial scale size = 0 gives the peak res, so we don't
  // need to do the LBFGS optimization for it
  if (itsOptimumScale == 0)
    return {};
  else
  {
    AlwaysAssert(itsAspScaleSizes.size() == itsAspAmplitude.size(), AipsError);
    AlwaysAssert(itsAspScaleSizes.size() == itsAspCenter.size(), AipsError);


    // heuristiclly determine active set for speed up
    Float resArea = 0.0;
    Int nX = itsDirty->shape()(0);
    Int nY = itsDirty->shape()(1);

    for (Int j = 0; j < nY; ++j)
    {
      for(Int i = 0; i < nX; ++i)
        resArea += abs((*itsDirty)(i,j));
    }
    const Float lamda = 318; //M31 new Asp - gold

    const Float threshold = lamda * resArea;
    
    vector<Float> tempx;
    vector<IPosition> activeSetCenter;

    vector<pair<Float,int>> vp; //(LenDev, idx)
   

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
    tempx.push_back(strengthOptimum);
    tempx.push_back(itsInitScaleSizes[optimumScale]);
    activeSetCenter.push_back(positionOptimum);

    unsigned int length = tempx.size();

    //cout << "Before: x = " << x.transpose() << endl;

    // GSL: set the initial guess
    gsl_vector *x = NULL;
    x = gsl_vector_alloc(length);
    gsl_vector_set_zero(x);

    for (unsigned int i = 0; i < length; i+=2)
    {
      gsl_vector_set(x, i, tempx[i]);
      gsl_vector_set(x, i+1, tempx[i+1]);

      // save aspen before optimization
      itsPrevAspAmplitude.push_back(tempx[i]); // active-set amplitude before lbfgs
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
    // update x is needed here.
    gsl_vector *optx = NULL;
    optx = gsl_multimin_fdfminimizer_x(s); //fdf
    //optx = gsl_multimin_fminimizer_x(s); // f only
    // end GSL optimization

    // put the updated x back to the class variables, itsAspAmp and itsAspScale
    for (unsigned int i = 0; i < length; i+= 2)
    {
      double amp = gsl_vector_get(optx, i);
      double scale = gsl_vector_get(optx, i+1);
      scale = (scale = fabs(scale)) < 0.4 ? 0 : scale;

      itsAspAmplitude.push_back(amp);
      itsAspScaleSizes.push_back(scale); // permanent list that doesn't get clear
      itsGoodAspAmplitude.push_back(amp); // active-set amplitude
      itsGoodAspActiveSet.push_back(scale); // active-set
      itsAspCenter.push_back(activeSetCenter[i/2]);
      if (scale < 0.4)
        itsAspGood.push_back(false);
      else
        itsAspGood.push_back(true);

      cout << "optx = " << amp << " " << scale << endl;
    }

    itsStrengthOptimum = itsAspAmplitude[itsAspAmplitude.size() -1]; // the latest aspen is the last element of x
    itsOptimumScaleSize = itsAspScaleSizes[itsAspScaleSizes.size() -1]; // the latest aspen is the last element of x
    itsGoodAspCenter = activeSetCenter;

    // debug
    cout << "opt strengthOptimum " << itsStrengthOptimum << " opt size " << itsOptimumScaleSize << endl;

    // free GSL stuff
    gsl_multimin_fdfminimizer_free(s); //fdf
    //gsl_multimin_fminimizer_free(s); // f only
    gsl_vector_free(x);
    //gsl_vector_free(optx); // Dont do it. Free this causes seg fault!!!

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
void AspMatrixCleaner::defineAspScales(vector<Float>& scaleSizes)
{
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

  itsScalesValid = true; 
}

void AspMatrixCleaner::switchedToHogbom()
{
	itsSwitchedToHogbom = true;
  itsNthHogbom += 1;
  itsNumIterNoGoodAspen.resize(0);
  //itsNumHogbomIter = ceil(100 + 50 * (exp(0.05*itsNthHogbom) - 1)); // zhang's formula
  itsNumHogbomIter = ceil(50 + 2 * (exp(0.05*itsNthHogbom) - 1)); // genie's formula
  cout << "Run hogbom for " << itsNumHogbomIter << " iterations." << endl;
}

void AspMatrixCleaner::setOrigDirty(const Matrix<Float>& dirty){
  itsOrigDirty=new Matrix<Float>(dirty.shape());
  itsOrigDirty->assign(dirty);
}


} //# NAMESPACE CASA - END
