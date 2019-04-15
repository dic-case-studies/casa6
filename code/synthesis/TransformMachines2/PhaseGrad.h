// -*- C++ -*-
//# PhaseGrad.h: Definition of the PhaseGrad class
//# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
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
//
#ifndef SYNTHESIS_TRANSFORM2_PHASEGRAD_H
#define SYNTHESIS_TRANSFORM2_PHASEGRAD_H

#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <synthesis/TransformMachines2/CFBuffer.h>

// #include <coordinates/Coordinates/DirectionCoordinate.h>
// #include <images/Images/ImageInterface.h>
// #include <images/Images/TempImage.h>
// #include <msvis/MSVis/VisBuffer2.h>
// #include <synthesis/MeasurementComponents/SolvableVisCal.h>
// class SolvableVisJones;

namespace casa { //# NAMESPACE CASA - BEGIN
  namespace refim{
  //
  //-------------------------------------------------------------------------------------------
  //
  class PhaseGrad
  {
  public:
    PhaseGrad():field_phaseGrad_p(), antenna_phaseGrad_p(), cachedCFBPtr_p(NULL), maxCFShape_p(2)
    {cached_FieldOffset_p[0] = cached_FieldOffset_p[1] = cached_AntennaOffset_p[0] = cached_AntennaOffset_p[1] = 0;}

    ~PhaseGrad() {};

    PhaseGrad& operator=(const PhaseGrad& other);
    
    inline const casacore::Matrix<casacore::Complex>& getFieldPointingGrad() {return field_phaseGrad_p;}
    inline const casacore::Matrix<casacore::Complex>& getAntennaPointingGrad() {return antenna_phaseGrad_p;}
    //    void getPhaseGrad(casacore::Matrix<casacore::Complex>& fullPhaseGrad) {fullPhaseGrad = antenna_phaseGrad_p + sky_phaseGrad_p;}

    bool ComputeFieldPointingGrad(const casacore::Vector<double>& pointingOffset,
				  const casacore::CountedPtr<CFBuffer>& cfb,
				  const vi::VisBuffer2& vb);
    bool ComputeFieldPointingGrad(const casacore::Vector<casacore::Vector<double> >& pointingOffset,
				  const casacore::CountedPtr<CFBuffer>& cfb,
				  const vi::VisBuffer2& vb,
				  const int& row);
    // bool ComputeFieldPointingGrad(const casacore::Vector<double>& pointingOffset,
    // 				  const casacore::Vector<int>&cfShape,
    // 				  const casacore::Vector<int>& convOrigin,
    // 				  const double& /*cfRefFreq*/,
    // 				  const double& /*imRefFreq*/,
    // 				  const int& spwID, const int& fieldId);
    
    //  private:
    casacore::Matrix<casacore::Complex> field_phaseGrad_p;
    casacore::Matrix<casacore::Complex> antenna_phaseGrad_p;
    casacore::Double cached_FieldOffset_p[2], cached_AntennaOffset_p[2];
    CFBuffer* cachedCFBPtr_p;
    casacore::Vector<int> maxCFShape_p;
  };
  //
  //-------------------------------------------------------------------------------------------
  //
};
};
#endif
