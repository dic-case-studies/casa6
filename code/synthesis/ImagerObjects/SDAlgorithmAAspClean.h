//# SDAlgorithmAAspClean.h: Definition for SDAlgorithmAAspClean
//# Copyright (C) 1996,1997,1998,1999,2000,2002
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
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#ifndef SYNTHESIS_SDALGORITHMAASPCLEAN_H
#define SYNTHESIS_SDALGORITHMAASPCLEAN_H

#include <ms/MeasurementSets/MeasurementSet.h>
//#include <synthesis/MeasurementComponents/SkyModel.h>
#include <casa/Arrays/Matrix.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/PagedImage.h>
#include <images/Images/TempImage.h>
#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>
#include <casa/System/PGPlotter.h>

#include<synthesis/ImagerObjects/SDAlgorithmBase.h>
#include<synthesis/ImagerObjects/SDAlgorithmHogbomClean.h>

namespace casa { //# NAMESPACE CASA - BEGIN

  /* Forware Declaration */
  class SIMinorCycleController;


  class SDAlgorithmAAspClean : public SDAlgorithmBase 
  {
  public:
    
    // Empty constructor
    SDAlgorithmAAspClean();
    virtual  ~SDAlgorithmAAspClean();
    
  protected:
    
    // Local functions to be overloaded by various algorithm deconvolvers.
    virtual void takeOneStep( casacore::Float loopgain, casacore::Int cycleNiter, casacore::Float cycleThreshold, casacore::Float &peakresidual, casacore::Float &modelflux, casacore::Int &iterdone );
    //    virtual void initializeDeconvolver( casacore::Float &peakresidual, casacore::Float &modelflux );
    virtual void initializeDeconvolver();
    virtual void finalizeDeconvolver();

    SDAlgorithmHogbomClean hogbom_p;

    ///    virtual void queryDesiredShape(casacore::Bool &onechan, casacore::Bool &onepol); // , nImageFacets.
    //    virtual void restorePlane();

    /*
    void findNextComponent( casacore::Float loopgain );
    void updateModel();
    void updateResidual();
    */

    /*
    casacore::SubImage<casacore::Float> itsResidual, itsPsf, itsModel, itsImage;
    casacore::Float itsComp;
    */
    //casacore::SubImage<casacore::Float> itsResidual, itsPsf, itsModel, itsImage;

    casacore::Array<casacore::Float> itsMatResidual, itsMatModel, itsMatPsf, itsMatMask;

    /*
    casacore::IPosition itsMaxPos;
    casacore::Float itsPeakResidual;
    casacore::Float itsModelFlux;

    casacore::Matrix<casacore::Float> itsMatMask;
    */

  };

} //# NAMESPACE CASA - END

#endif
