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
#include <casa/Arrays/ArrayIO.h>
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
#include <coordinates/Coordinates/TabularCoordinate.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Additional include files
#include <algorithm>
#include<synthesis/MeasurementEquations/AspMatrixCleaner.h>

#include <LBFGS.h>
#include <synthesis/MeasurementEquations/lbfgsAsp.h>

using namespace casacore;
using Eigen::VectorXd;
using namespace LBFGSpp;

namespace casa { //# NAMESPACE CASA - BEGIN

AspMatrixCleaner::AspMatrixCleaner():
  MatrixCleaner(),
  itsInitScaleSizes(0),
  itsAspScaleSizes(0),
  itsAspAmplitude(0),
  itsNInitScales(4),
  itsPrevLBFGSGrad(0.0),
  itsNumIterNoGoodAspen(0),
  itsSwitchedToHogbom(false)
{
  itsInitScales.resize(0);
  itsInitScaleXfrs.resize(0);
  itsDirtyConvInitScales.resize(0);
  itsInitScaleMasks.resize(0);
  itsPsfConvInitScales.resize(0);
}

AspMatrixCleaner::~AspMatrixCleaner()
{
  destroyAspScales();
  destroyInitMasks();
  if(!itsMask.null())
    itsMask=0;
}

void AspMatrixCleaner::makedirtyscales()
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
      //os << "Calculating convolutions for scale size " << itsScaleSizes(scale) << LogIO::POST;
      Matrix<Complex> cWork;

      itsDirtyConvScales[scale]=Matrix<Float>(itsDirty->shape());
      cWork=((dirtyFT)*(itsScaleXfrs[scale]));
      fft.fft0((itsDirtyConvScales[scale]), cWork, false);
      fft.flip((itsDirtyConvScales[scale]), false, false);
    }
  }
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

  Int nScalesToClean = itsNscales;
  os << LogIO::NORMAL1 << "AAsp clean algorithm" << LogIO::POST;


  Int scale;
  Vector<Float> scaleBias(nScalesToClean);
  // scaleBias is 1.0 for AAsp for now
  for (scale=0; scale < nScalesToClean; scale++)
    scaleBias(scale) = 1.0;

  AlwaysAssert(itsScalesValid, AipsError);
  ////no need to use all cores if possible
  Int nth = nScalesToClean;
#ifdef _OPENMP

    nth = min(nth, omp_get_max_threads());

#endif
  // Find the peaks of the convolved Psfs
  Vector<Float> maxPsfConvScales(nScalesToClean);
  Int naxes = model.shape().nelements();

#pragma omp parallel default(shared) private(scale) num_threads(nth)
  {
    #pragma omp for
    for (scale=0; scale<nScalesToClean; scale++)
    {
      IPosition positionPeakPsfConvScales(naxes, 0);

      findMaxAbs(itsPsfConvScales[scale], maxPsfConvScales(scale),
      positionPeakPsfConvScales);

      //   cout  << "MAX  " << scale << "    " << positionPeakPsfConvScales
      //      << "  " << maxPsfConvScales(scale) << "\n"
      //      << endl;
     }
  } //End pragma parallel
  for (scale=0; scale<nScalesToClean; scale++)
  {
    if (maxPsfConvScales(scale) < 0.0)
    {
      os << "As Peak of PSF is negative, you should setscales again with a smaller scale size"
   << LogIO::SEVERE;
      return -1;
    }
  }

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

    //blcDirty(0)=xbeg> 0 ? xbeg-1 : 0;
    //blcDirty(1)=ybeg > 0 ? ybeg-1 : 0;
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

  Block<Matrix<Float> > scaleMaskSubs;
  if (!itsMask.null())
  {
    scaleMaskSubs.resize(itsNscales);
    for (Int is=0; is < itsNscales; is++)
      scaleMaskSubs[is] = ((itsScaleMasks[is]))(blcDirty, trcDirty);
  }

  // Start the iteration
  Vector<Float> maxima(nScalesToClean);
  Block<IPosition> posMaximum(nScalesToClean);
  Vector<Float> totalFluxScale(nScalesToClean);
  totalFluxScale=0.0;
  Float totalFlux=0.0;
  Int converged=0;
  Int stopPointModeCounter = 0;
  Int optimumScale=0;
  Float tmpMaximumResidual = 0.0;
  itsStrengthOptimum=0.0;
  IPosition positionOptimum(model.shape().nelements(), 0);

  Int nx=model.shape()(0);
  Int ny=model.shape()(1);
  IPosition gip;
  gip = IPosition(2, nx, ny);
  casacore::Block<casacore::Matrix<casacore::Float> > vecWork_p;
  vecWork_p.resize(nScalesToClean);
  for (Int i=0; i<nScalesToClean; i++)
    vecWork_p[i].resize(gip);

  os << "Starting iteration"<< LogIO::POST;
  vector<Float> tempScaleSizes;
  itsIteration = itsStartingIter; // 0

  for (Int ii = itsStartingIter; ii < itsMaxNiter; ii++)
  {
    cout << "cur iter " << itsIteration << " max iter is "<<
            itsMaxNiter << endl;
    itsIteration++;

    // Find the peak residual
    itsStrengthOptimum = 0.0;
    optimumScale = 0;

#pragma omp parallel default(shared) private(scale) num_threads(nth)
    {
#pragma omp  for
      for (scale=0; scale < min(nScalesToClean,itsNscales); ++scale)
      {
        // Find absolute maximum for the dirty image
        //cout << "clean: in omp loop for scale : " << scale << " : " << blcDirty << " : " << trcDirty << " :: " << itsDirtyConvScales.nelements() << endl;
        Matrix<Float> work = (vecWork_p[scale])(blcDirty,trcDirty);
        work = 0.0;
        work = work + (itsDirtyConvScales[scale])(blcDirty,trcDirty);
        maxima(scale)=0;
        posMaximum[scale]=IPosition(model.shape().nelements(), 0);

        //genie debug
        /*Float maxVal=0;
        IPosition posmin(vecWork_p[scale].shape().nelements(), 0);
        Float minVal=0;
        IPosition posmax(vecWork_p[scale].shape().nelements(), 0);
        minMaxMasked(minVal, maxVal, posmin, posmax, vecWork_p[scale], itsScaleMasks[scale]);
        cout << "clean: ScaleVal " << scale << ": min " << minVal << " max " << maxVal << endl;*/
        //genie

        if (!itsMask.null()) {
          findMaxAbsMask(vecWork_p[scale], itsScaleMasks[scale],
              maxima(scale), posMaximum[scale]);
          //cout << "before: maxima[" << scale << "] = " << maxima(scale) << endl;
        }
        else
          findMaxAbs(vecWork_p[scale], maxima(scale), posMaximum[scale]);

        // Remember to adjust the position for the window and for
        // the flux scale
        // cout << "scale " << scale << " maxPsfconvscale " << maxPsfConvScales(scale) << endl;
        // cout << "posmax " << posMaximum[scale] << " blcdir " << blcDirty << endl;
        // cout << "maxima " << maxima(scale) << " dirconvscale " << (itsDirtyConvScales[scale])(posMaximum[scale]) << endl;
        maxima(scale)/=maxPsfConvScales(scale);
        maxima(scale) *= scaleBias(scale);
        maxima(scale) *= (itsDirtyConvScales[scale])(posMaximum[scale]); //makes maxima(scale) positive to ensure correct scale is selected in itsStrengthOptimum for loop (next for loop).
        //cout << "maxPsfconvscale[" << scale << "] = " << maxPsfConvScales(scale) << endl;
        //cout << "after: maxima[" << scale << "] = " << maxima(scale) << endl;

        //posMaximum[scale]+=blcDirty;
      }
    } //End parallel section
    for (scale=0; scale < min(nScalesToClean,itsNscales); scale++)
    {
      if(abs(maxima(scale)) > abs(itsStrengthOptimum))
      {
        optimumScale=scale;
        itsStrengthOptimum = maxima(scale);
        positionOptimum = posMaximum[scale];
        cout << "clean: New optimum scale is " << itsScaleSizes(scale) << " itsStrengthOptimum " << itsStrengthOptimum << endl;
      }
    }

    itsStrengthOptimum /= scaleBias(optimumScale);
    itsStrengthOptimum /=  (itsDirtyConvScales[optimumScale])(posMaximum[optimumScale]);

    AlwaysAssert(optimumScale < nScalesToClean, AipsError);

    // Now add to the total flux
    totalFlux += (itsStrengthOptimum*itsGain);
    itsTotalFlux = totalFlux;
    totalFluxScale(optimumScale) += (itsStrengthOptimum*itsGain);

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
    if (abs(itsStrengthOptimum) < threshold())
    {
      os << "Reached stopping threshold " << threshold() << " at iteration "<<
            ii << LogIO::POST;
      os << "Optimum flux is " << abs(itsStrengthOptimum) << LogIO::POST;
      converged = 1;
      break;
    }
    //    2. negatives on largest scale?
    if ((nScalesToClean > 1) && itsStopAtLargeScaleNegative  &&
        optimumScale == (nScalesToClean-1) &&
        itsStrengthOptimum < 0.0)
    {
      os << "Reached negative on largest scale" << LogIO::POST;
      converged = -2;
      break;
    }
    //  3. stop point mode at work
    if (itsStopPointMode > 0)
    {
      if (optimumScale == 0)
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

      /*
      if(progress) {
        progress->info(false, itsIteration, itsMaxNiter, maxima,
        posMaximum, itsStrengthOptimum,
        optimumScale, positionOptimum,
        totalFlux, totalFluxScale,
        itsJustStarting );
        itsJustStarting = false;
      } else*/
    {
      if (itsIteration == itsStartingIter + 1)
        os << "iteration    MaximumResidual   CleanedFlux" << LogIO::POST;
      if ((itsIteration % (itsMaxNiter/10 > 0 ? itsMaxNiter/10 : 1)) == 0)
      {
        //Good place to re-up the fiducial maximum residual
        //tmpMaximumResidual=abs(itsStrengthOptimum);
        os << itsIteration <<"      "<<itsStrengthOptimum<<"      "
           << totalFlux <<LogIO::POST;
      }
    }

    // Continuing: subtract the peak that we found from all dirty images
    // Define a subregion so that the peak is centered
    IPosition support(model.shape());
    support(0) = max(Int(itsScaleSizes(itsNscales-1)+0.5), support(0));
    support(1) = max(Int(itsScaleSizes(itsNscales-1)+0.5), support(1));

    IPosition inc(model.shape().nelements(), 1);
    //cout << "support " << support.asVector()  << endl;
    //support(0)=1024;
    //support(1)=1024;
    //support(0)=min(Int(support(0)), Int(trcDirty(0)-blcDirty(0)));
    //support(1)=min(Int(support(1)), Int(trcDirty(1)-blcDirty(1)));
    // support(0)=min(Int(support(0)), (trcDirty(0)-blcDirty(0)+
    //        Int(2*abs(positionOptimum(0)-blcDirty(0)/2.0-trcDirty(0)/2.0))));
    //support(1)=min(Int(support(1)), (trcDirty(1)-blcDirty(1)+
    //        Int(2*abs(positionOptimum(1)-blcDirty(1)/2.0-trcDirty(1)/2.0))));

    IPosition blc(positionOptimum-support/2);
    IPosition trc(positionOptimum+support/2-1);
    LCBox::verify(blc, trc, inc, model.shape());

    // cout << "verify: blc " << blc.asVector() << " trc " << trc.asVector() << endl;

    IPosition blcPsf(blc+itsPositionPeakPsf-positionOptimum);
    IPosition trcPsf(trc+itsPositionPeakPsf-positionOptimum);
    LCBox::verify(blcPsf, trcPsf, inc, model.shape());
    makeBoxesSameSize(blc,trc,blcPsf,trcPsf);
    // cout << "samebox: blcPsf " << blcPsf.asVector() << " trcPsf " << trcPsf.asVector() << endl;
    // cout << "samebox: blc " << blc.asVector() << " trc " << trc.asVector() << endl;
    //    LCBox subRegion(blc, trc, model.shape());
    //  LCBox subRegionPsf(blcPsf, trcPsf, model.shape());

    Matrix<Float> modelSub = model(blc, trc);
    Matrix<Float> scaleSub = (itsScales[optimumScale])(blcPsf,trcPsf);

    // Now do the addition of this scale to the model image...
    Float scaleFactor;
    scaleFactor = itsGain * itsStrengthOptimum;
    modelSub += scaleFactor * scaleSub;

    // Now update the (residual image * scale)
    //cout << "before dirtySub(100,100) = " << itsDirtyConvScales[0](100,100) << endl;
    #pragma omp parallel default(shared) private(scale) num_threads(nth)
    {
      #pragma omp  for
      for (scale=0; scale < min(nScalesToClean,itsNscales); ++scale)
      {
        Matrix<Float> dirtySub=(itsDirtyConvScales[scale])(blc,trc);
        //AlwaysAssert(itsPsfConvScales[index(scale,optimumScale)], AipsError);
        Matrix<Float> psfSub=(itsPsfConvScales[index(scale,optimumScale)])(blcPsf, trcPsf);
        dirtySub -= scaleFactor*psfSub;
      }
    }//End parallel

    blcDirty = blc;
    trcDirty = trc;

    //genie If we switch to hogbom (i.e. only have 0 scale size)
    // there is no need to do the following Aspen update
    //if (nScalesToClean == 1)
    if (itsSwitchedToHogbom && itsNscales == 1)
      continue;

    //genie Now update the actual residual image
    // At this point, itsDirty is not updated. Only itsDirtyConvScales is updated.
    //cout << "after dirtySub(100,100) = " << itsDirtyConvScales[0](100,100) << endl;
    //cout << "after itsdirty(100,100) = " << (*itsDirty)(100,100) << endl;
    setDirty(residual()); // this updates itsDirty correctly
    //cout << "after2 itsdirty(100,100) = " << (*itsDirty)(100,100) << endl;
    // need to getActiveSetAspen, add 0 scale size
    // and then defineAspScales, makePsfScales, makeScaleMasks, makeDirtyScales?
    // need to re-think the logic v.s. what's done in ADAlgAAspClean.

    tempScaleSizes.clear();
    tempScaleSizes = getActiveSetAspen();
    /*for (scale = 0; scale < int(tempScaleSizes.size()); scale++)
      cout << "2. getActiveSetAspen[" << scale << "] " << tempScaleSizes[scale] << endl;
    cout << "# tempScaleSizes " << tempScaleSizes.size() << endl;*/
    tempScaleSizes.push_back(0.0); // put 0 scale
    //Vector<Float> scaleSizes(tempScaleSizes);
    //defineAspScales(scaleSizes);
    defineAspScales(tempScaleSizes);

    makePsfScales();
    makeScaleMasks();
    makedirtyscales();
    //genie
  }
  // End of iteration

  for (scale=0; scale < min(nScalesToClean,itsNscales); scale++) {
    os << LogIO::NORMAL
       << "  " << scale << "    " << totalFluxScale(scale)
       << LogIO::POST;
  }
  // Finish off the plot, etc.
  /*
  if(progress) {
    progress->info(true, itsIteration, itsMaxNiter, maxima, posMaximum,
       itsStrengthOptimum,
       optimumScale, positionOptimum,
       totalFlux, totalFluxScale);
  }
  */

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

  return true;
}

Bool AspMatrixCleaner::destroyInitMasks()
{
  for(uInt scale=0; scale<itsInitScaleMasks.nelements();scale++)
    itsInitScaleMasks[scale].resize();

  itsInitScaleMasks.resize(0);

  return true;
};


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

  return float(ceil((xpixels + ypixels)/2));
}

void AspMatrixCleaner::setInitScaleXfrs(/*const Array<Float> arrpsf, */const Float width)
{
  if(itsScales.nelements() > 0)
    destroyAspScales();

  /*Matrix<Float> tempMat;
  tempMat.reference(arrpsf);
  psfShape_p.resize(0, false);
  psfShape_p = tempMat.shape();*/

  // try 0, 1.5width, 5width and 10width
  itsInitScaleSizes.resize(itsNInitScales, false);
  itsInitScaleSizes = {0.0f, 1.5f*width, 5.0f*width, 10.0f*width};
  itsInitScales.resize(itsNInitScales, false);
  itsInitScaleXfrs.resize(itsNInitScales, false);
  FFTServer<Float,Complex> fft(psfShape_p);
  for (int scale = 0; scale < itsNInitScales; scale++)
  {
    itsInitScales[scale] = Matrix<Float>(psfShape_p);
    makeScale(itsInitScales[scale], itsInitScaleSizes[scale]);
    cout << "made itsInitScales[" << scale << "] = " << itsInitScaleSizes[scale] << endl;
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

  const int nx = itsDirty->shape()[0];
  const int ny = itsDirty->shape()[1];
  cout << "nx " << nx << " ny " << ny << endl;
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
  IPosition gip;
  gip = IPosition(2, nx, ny);
  casacore::Block<casacore::Matrix<casacore::Float> > vecWork_p;
  vecWork_p.resize(itsNInitScales);

  for (int i = 0; i < itsNInitScales; i++)
    vecWork_p[i].resize(gip);

  Vector<Float> maxima(itsNInitScales);
  Block<IPosition> posMaximum(itsNInitScales);
  int scale;
  int nth = itsNInitScales;
  #ifdef _OPENMP
    nth = min(nth, omp_get_max_threads());
  #endif

  //genie
  // Find the peaks of the convolved Psfs
  Vector<Float> maxPsfConvInitScales(itsNInitScales);
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
  }
  //genie

  #pragma omp parallel default(shared) private(scale) num_threads(nth)
  {
    #pragma omp for
    for (scale=0; scale < itsNInitScales; ++scale)
    {
      // Find absolute maximum for the dirty image
      //cout << "in omp loop for scale : " << scale << " : " << blcDirty << " : " << trcDirty << " :: " << itsDirty->shape().nelements() << endl;
      Matrix<Float> work = (vecWork_p[scale])(blcDirty, trcDirty);
      work = 0.0;
      work = work + (itsDirtyConvInitScales[scale])(blcDirty, trcDirty);
      maxima(scale) = 0;
      posMaximum[scale] = IPosition(itsDirty->shape().nelements(), 0);
      //cout << "makedirtyinitscale before: " << itsInitScaleMasks[scale].shape() << endl;

      //genie debug
      Float maxVal=0;
      IPosition posmin(vecWork_p[scale].shape().nelements(), 0);
      Float minVal=0;
      IPosition posmax(vecWork_p[scale].shape().nelements(), 0);
      minMaxMasked(minVal, maxVal, posmin, posmax, vecWork_p[scale], itsInitScaleMasks[scale]);
      //cout << "InitScaleVal " << scale << ": min " << minVal << " max " << maxVal << endl;
      //genie

      if (!itsMask.null())
      {
        findMaxAbsMask(vecWork_p[scale], itsInitScaleMasks[scale],
          maxima(scale), posMaximum[scale]);
      }
      else
        findMaxAbs(vecWork_p[scale], maxima(scale), posMaximum[scale]);

      // Remember to adjust the position for the window and for
      // the flux scale
      maxima(scale) /= maxPsfConvInitScales(scale);
      maxima(scale) *= (itsDirtyConvInitScales[scale])(posMaximum[scale]); //makes maxima(scale) positive to ensure correct scale is selected in strengthOptimum for loop (next for loop).
      //cout << "maxDirty: maxPsfconvinitscale[" << scale << "] = " << maxPsfConvInitScales(scale) << endl;
      //cout << "after: maxima[" << scale << "] = " << maxima(scale) << endl;
      //genie: it seems we need the above
    }
  }//End parallel section

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

  AlwaysAssert(optimumScale < itsNInitScales, AipsError);
}


bool AspMatrixCleaner::isGoodAspen(Float amp, Float scale, IPosition center, Float threshold)
{
  const int nX = itsDirty->shape()(0);
  const int nY = itsDirty->shape()(1);

  Matrix<Float> Asp(nX, nY);
  Gaussian2D<Float> gbeam(amp, center[0], center[1], scale, 1, 0);
  for (int j = 0; j < nY; ++j)
  {
    for(int i = 0; i < nX; ++i)
    {
      int px = i - nX/2;
      int py = j - nY/2;
      Asp(i,j) = gbeam(px, py);
    }
  }

  // gradient. 0: amplitude; 1: scale
  // generate derivative of amplitude
  Matrix<Float> GradAmp(nX, nY);
  Gaussian2D<Float> gbeamGradAmp(1, center[0], center[1], scale, 1, 0);
  for (int j = 0; j < nY; ++j)
  {
    for(int i = 0; i < nX; ++i)
    {
      int px = i - nX/2;
      int py = j - nY/2;
      GradAmp(i,j) = (-2) * gbeamGradAmp(px, py);
    }
  }
  Matrix<Float> Grad0 = product(transpose(*itsDirty), GradAmp);

  // generate derivative of scale
  Matrix<Float> GradScale(nX, nY);
  for (int j = 0; j < nY; ++j)
  {
    for(int i = 0; i < nX; ++i)
      GradScale(i, j) = (-2)*2*(pow(i-center[0],2) + pow(j-center[1],2))*Asp(i,j)/pow(scale,3);
  }
  Matrix<Float> Grad1 = product(transpose(*itsDirty), GradScale);

  // calculate the length of the direvative vector
  Float lenDirVec = 0.0;
  for (int j = 0; j < Grad0.shape()(1); ++j)
  {
    for(int i = 0; i < Grad0.shape()(0); ++i)
    {
      lenDirVec += sqrt(pow(Grad0(i,j), 2));
    }
  }

  for (int j = 0; j < Grad1.shape()(1); ++j)
  {
    for(int i = 0; i < Grad1.shape()(0); ++i)
    {
      lenDirVec += sqrt(pow(Grad1(i,j), 2));
    }
  }

  cout << "scale " << scale << " amp " << amp << " center " << center << " lenDirVec " << lenDirVec << " threshold " << threshold << endl;
  if (lenDirVec >= threshold)
  {
    //cout << "lenDirVec " << lenDirVec << " threshold " << threshold << endl;
    return true;
  }

  return false;

}


vector<Float> AspMatrixCleaner::getActiveSetAspen()
{
  LogIO os(LogOrigin("AspMatrixCleaner", "getActiveSetAspen()", WHERE));

  if(int(itsInitScaleXfrs.nelements()) == 0)
    throw(AipsError("Initial scales for Asp are not defined"));

  // if there is no "good" aspen put into the active-set
  // for more than 10 consecutive iterations, we switch
  // from Asp to hogbom
  if (itsNumIterNoGoodAspen > 10)
  {
    cout << "Switch to hogbom." << endl;
    itsSwitchedToHogbom = true;
    return {};
  }

  // Dirty * initial scales
  Matrix<Complex> dirtyFT;
  FFTServer<Float,Complex> fft(itsDirty->shape());
  fft.fft0(dirtyFT, *itsDirty);
  itsDirtyConvInitScales.resize(itsNInitScales); // 0, 1.5width, 5width and 10width

  for (int scale=0; scale < itsNInitScales; scale++)
  {
    Matrix<Complex> cWork;
    //cout << "scale " << scale << " itsInitScaleptr " << &(itsInitScaleXfrs[scale]) << endl;

    itsDirtyConvInitScales[scale] = Matrix<Float>(itsDirty->shape());
    cWork=((dirtyFT)*(itsInitScaleXfrs[scale]));
    fft.fft0((itsDirtyConvInitScales[scale]), cWork, false);
    fft.flip((itsDirtyConvInitScales[scale]), false, false);
  }

  float strengthOptimum = 0.0;
  int optimumScale = 0;
  IPosition positionOptimum(itsDirty->shape().nelements(), 0);

  maxDirtyConvInitScales(strengthOptimum, optimumScale, positionOptimum);
  cout << "Initial maximum residual is " << abs(strengthOptimum);
  cout << " at location " << positionOptimum[0] << " " << positionOptimum[1];
  cout << " " << positionOptimum[2] << " and scale: " << optimumScale << endl;

  // initial scale size = 0 gives the peak res, so we don't
  // need to do the LBFGS optimization for it
  if (optimumScale == 0)
    return itsAspScaleSizes;
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
    if (find(itsAspAmplitude.begin(), itsAspAmplitude.end(), strengthOptimum) != itsAspAmplitude.end() &&
        find(itsAspScaleSizes.begin(), itsAspScaleSizes.end(), itsInitScaleSizes[optimumScale]) != itsAspScaleSizes.end() &&
        find(itsAspCenter.begin(), itsAspCenter.end(), positionOptimum) != itsAspCenter.end())
    {
      return itsAspScaleSizes;
    }
    else //i.e. add the unique apsen to the permanent list
    {
      itsAspAmplitude.push_back(strengthOptimum);
      itsAspScaleSizes.push_back(itsInitScaleSizes[optimumScale]);
      itsAspCenter.push_back(positionOptimum);
    }

    // test scale up genie
    for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
      itsAspAmplitude[i] = itsAspAmplitude[i] * 1e8;


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
    const Float lamda = 7.65e-1; //G55
    //const Float lamda = 1.5; //G55
    const Float threshold = lamda * resArea;
    vector<Float> tempx;
    vector<IPosition> activeSetCenter;

    bool hasGoodAspen = false;
    for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
    {
      //cout << "scale " << i << " heuristic: " << endl;
      if (isGoodAspen(itsAspAmplitude[i], itsAspScaleSizes[i], itsAspCenter[i], threshold))
      {
        tempx.push_back(itsAspAmplitude[i]);
        tempx.push_back(itsAspScaleSizes[i]);
        activeSetCenter.push_back(itsAspCenter[i]);

        itsAspAmplitude.erase(itsAspAmplitude.begin() + i);
        itsAspScaleSizes.erase(itsAspScaleSizes.begin() + i);
        itsAspCenter.erase(itsAspCenter.begin() + i);

        hasGoodAspen = true;
      }
    }

    if (hasGoodAspen)
      itsNumIterNoGoodAspen = 0;
    else
      itsNumIterNoGoodAspen += 1;

    //
    //unsigned int length = (itsAspScaleSizes.size() + 1) * 2;  //genie recon
    //VectorXd x(length); //genie recon
    unsigned int length = tempx.size();
    VectorXd x(length);

    /*vector<Float> tempx;
    for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
    {
      tempx.push_back(itsAspAmplitude[i]);
      tempx.push_back(itsAspScaleSizes[i]);
    }*/
    //tempx.push_back(strengthOptimum);
    //tempx.push_back(itsInitScaleSizes[optimumScale]);
    //itsAspCenter.push_back(positionOptimum);

    for (unsigned int i = 0; i < length; i++)
      x[i] = tempx[i]; //Eigen::VectorXd needs to be assigned in this way

    cout << "Before: x = " << x.transpose() << endl;

    // lbfgs optimization
    LBFGSParam<double> param;
    param.epsilon = 1e-2;
    param.max_linesearch = 10;
    param.min_step = 1e-30;
    param.max_iterations = 2; //M31
    /*param.epsilon = 1e-5;
    param.max_linesearch = 10;
    param.min_step = 1e-30;
    param.max_iterations = 2;
    param.delta = 1e-15;*/ //G55
    param.gclip = itsPrevLBFGSGrad;
    //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE; genie: doesn't work
    LBFGSSolver<double> solver(param);

    //AspObjFunc fun(*itsDirty, *itsXfr, itsAspCenter);
    AspObjFunc fun(*itsDirty, *itsXfr, activeSetCenter);

    double fx;
    double gclip;
    /*int niter =*/ solver.minimize(fun, x, fx, gclip);

    //std::cout << niter << " iterations" << std::endl;
    // use the initial gradient as a roll back gradient if there is
    // gradient exploding in lbfgs
    if (itsPrevLBFGSGrad == 0.0)
      itsPrevLBFGSGrad = gclip;
    //std::cout << "itsPrevLBFGSGrad " << itsPrevLBFGSGrad << std::endl;
    std::cout << "After: x = " << x.transpose() << std::endl;
    //std::cout << "f(x) = " << fx << std::endl;
    //std::cout << "float is " << Float(x[1]) << endl;

    // x is not changing in LBFGS which causes convergence issue
    // so we use the previous x instead. This is a hack suggested by
    // the online community
    if (fx == 0)
    {
      for (unsigned int i = 0; i < length; i++)
        x[i] = tempx[i];
    }

    // put the updated x back to the class variables, itsAspAmp and itsAspScale
    //itsAspAmplitude.clear(); //cannot
    //itsAspScaleSizes.clear(); //cannot
    for (unsigned int i = 0; i < length; i+= 2)
    {
      itsAspAmplitude.push_back(x[i]);
      itsAspScaleSizes.push_back(x[i+1]);
      itsAspCenter.push_back(activeSetCenter[i/2]);
    }
  } // end of LBFGS optimization

  // debug info
  /*for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
  {
    //cout << "After opt AspApm[" << i << "] = " << itsAspAmplitude[i] << endl;
    //cout << "After opt AspScale[" << i << "] = " << itsAspScaleSizes[i] << endl;
    //cout << "After opt AspCenter[" << i << "] = " << itsAspCenter[i] << endl;
    cout << "AspScale[ " << i << " ] = " << itsAspScaleSizes[i] << " center " << itsAspCenter[i] << endl;
  }*/
  // test scale back genie
  for (unsigned int i = 0; i < itsAspAmplitude.size(); i++)
    itsAspAmplitude[i] = itsAspAmplitude[i] / 1e8;

  return itsAspScaleSizes; // return optimized scale
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

  sort(scaleSizes.begin(), scaleSizes.end());
  scaleSizes.erase(unique(scaleSizes.begin(), scaleSizes.end()), scaleSizes.end());
  //scaleSizes.erase(unique(scaleSizes.begin(),scaleSizes.end(),[](Float l, Float r) { return abs(l - r) < 1e-5; }), scaleSizes.end());

  itsNscales = Int(scaleSizes.size());
  itsScaleSizes.resize(itsNscales);
  itsScaleSizes = Vector<Float>(scaleSizes);  // make a copy that we can call our own

  itsScalesValid = true;  //genie? It's false in MS clean
}

} //# NAMESPACE CASA - END
