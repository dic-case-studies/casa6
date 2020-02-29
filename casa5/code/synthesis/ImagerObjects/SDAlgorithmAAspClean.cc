//# SDAlgorithmAAspClean.cc: Implementation of SDAlgorithmAAspClean classes
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
//# $Id$

#include <casa/Arrays/ArrayMath.h>
#include <casa/OS/HostInfo.h>
#include <synthesis/ImagerObjects/SDAlgorithmAAspClean.h>

#include <components/ComponentModels/SkyComponent.h>
#include <components/ComponentModels/ComponentList.h>
#include <images/Images/TempImage.h>
#include <images/Images/SubImage.h>
#include <images/Regions/ImageRegion.h>
#include <casa/OS/File.h>
#include <lattices/LEL/LatticeExpr.h>
#include <lattices/Lattices/TiledLineStepper.h>
#include <lattices/Lattices/LatticeStepper.h>
#include <lattices/Lattices/LatticeIterator.h>
#include <synthesis/TransformMachines/StokesImageUtil.h>
#include <coordinates/Coordinates/StokesCoordinate.h>
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>
#include <casa/OS/Directory.h>
#include <tables/Tables/TableLock.h>

#include<synthesis/ImagerObjects/SIMinorCycleController.h>

#include <casa/sstream.h>

#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogSink.h>

#include <casa/System/Choice.h>
#include <msvis/MSVis/StokesVector.h>


using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

  SDAlgorithmAAspClean::SDAlgorithmAAspClean( Vector<Float> scalesizes,
            Int stoppointmode ):
    SDAlgorithmBase(),
    itsMatPsf(), itsMatResidual(), itsMatModel(),
    itsCleaner(),
    itsScaleSizes(scalesizes),
    itsStopPointMode(stoppointmode),
    itsMCsetup(false)
  {
    itsAlgorithmName = String("aasp");   
  }

  SDAlgorithmAAspClean::~SDAlgorithmAAspClean()
  {
    
  }
  
  void SDAlgorithmAAspClean::initializeDeconvolver()
  {
    LogIO os(LogOrigin("SDAlgorithmAAspClean", "initializeDeconvolver", WHERE));
    AlwaysAssert((bool)itsImages, AipsError);

    itsImages->residual()->get( itsMatResidual, true );
    itsImages->model()->get( itsMatModel, true );
    itsImages->psf()->get( itsMatPsf, true );
    itsImages->mask()->get( itsMatMask, true );

    //// Initialize the MatrixCleaner.
    ///  ----------- do once ----------
    if( itsMCsetup == false)
    {
      const float width = itsCleaner.getPsfGaussianWidth(*(itsImages->psf()));
      itsCleaner.setInitScaleXfrs(itsMatPsf, width);
      itsCleaner.defineScales( itsScaleSizes ); // genie, this goes away

      // genie this is only done once 
      // FFT of 1R, 5R, 10R of psf is unchanged and only needs to be
      // computed once. This calls getPsfGaussianWidth 
      // and sets the new itsInitScaleXfrs
      //itsCleaner.setInitScaleXfrs(); 

      itsCleaner.stopPointMode( itsStopPointMode );
      itsCleaner.ignoreCenterBox( true ); // Clean full image

      // genie move the following 4 to be done at every minor cycle
      // it looks like they cannot be used for Asp directly
      Matrix<Float> tempMat;
      tempMat.reference( itsMatPsf );
      itsCleaner.setPsf(  tempMat );
      itsCleaner.makePsfScales();

      itsMCsetup=true;
    }

    // Parts to be repeated at each minor cycle start....
    itsCleaner.setcontrol(CleanEnums::MULTISCALE,0,0,0);/// Needs to come before makeDirtyScales
    // genie
    // find peak to determine a scale and put the optimized scale in the active set
    // by itsScaleSizes = getActiveSetAspen which convolve
    // cWork=((dirtyFT)*(itsInitScaleXfrs[scale]));
    // and then find the peak, scale, and optimes the obj function
    // Convolve psf with the active set (using the above 4) and do the following 

    Matrix<Float> tempmask(itsMatMask);
    itsCleaner.setMask( tempmask );

    Matrix<Float> tempMat1;
    tempMat1.reference( itsMatResidual );
    itsCleaner.setDirty( tempMat1 );
    itsCleaner.makeDirtyScales();
  }


  void SDAlgorithmAAspClean::takeOneStep( Float loopgain, 
					  Int cycleNiter, 
					  Float cycleThreshold, 
					  Float &peakresidual, 
					  Float &modelflux, 
					  Int &iterdone)
  {
    LogIO os( LogOrigin("SDAlgorithmAAspClean","takeOneStep", WHERE) );

    Quantity thresh(cycleThreshold, "Jy");
    itsCleaner.setcontrol(CleanEnums::MULTISCALE, cycleNiter, loopgain, thresh); 

    Matrix<Float> tempModel;
    tempModel.reference( itsMatModel );
    //save the previous model
    Matrix<Float> prevModel;
    prevModel=itsMatModel;

    //cout << "SDALMS,  matrix shape : " << tempModel.shape() << " array shape : " << itsMatModel.shape() << endl;

    // retval
    //  1 = converged
    //  0 = not converged but behaving normally
    // -1 = not converged and stopped on cleaning consecutive smallest scale
    // -2 = not converged and either large scale hit negative or diverging 
    // -3 = clean is diverging rather than converging
    itsCleaner.startingIteration( 0 );
    Int retval = itsCleaner.clean( tempModel );
    iterdone = itsCleaner.numberIterations();

    if( retval==-1 ) {os << LogIO::WARN << "AAspClean minor cycle stopped on cleaning consecutive smallest scale" << LogIO::POST; }
    if( retval==-2 ) {os << LogIO::WARN << "AAspClean minor cycle stopped at large scale negative or diverging" << LogIO::POST;}
    if( retval==-3 ) {os << LogIO::WARN << "AAspClean minor cycle stopped because it is diverging" << LogIO::POST; }

    ////This is going to be wrong if there is no 0 scale;
    Matrix<Float> residual(itsCleaner.residual(tempModel-prevModel));
    // account for mask as well
    peakresidual = max(abs(residual*itsMatMask));
    modelflux = sum( itsMatModel ); 
  }	    

  void SDAlgorithmAAspClean::finalizeDeconvolver()
  {
    (itsImages->residual())->put( itsMatResidual );
    (itsImages->model())->put( itsMatModel );
  }

} //# NAMESPACE CASA - END

