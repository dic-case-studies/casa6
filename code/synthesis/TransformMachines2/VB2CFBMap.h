// -*- C++ -*-
// //# PhaseGrad.h: Definition of the PhaseGrad class
// //# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
// //# Associated Universities, Inc. Washington DC, USA.
// //#
// //# This library is free software; you can redistribute it and/or modify it
// //# under the terms of the GNU Library General Public License as published by
// //# the Free Software Foundation; either version 2 of the License, or (at your
// //# option) any later version.
// //#
// //# This library is distributed in the hope that it will be useful, but WITHOUT
// //# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// //# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
// //# License for more details.
// //#
// //# You should have received a copy of the GNU Library General Public License
// //# along with this library; if not, write to the Free Software Foundation,
// //# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
// //#
// //# Correspondence concerning AIPS++ should be addressed as follows:
// //#        Internet email: aips2-request@nrao.edu.
// //#        Postal address: AIPS++ Project Office
// //#                        National Radio Astronomy Observatory
// //#                        520 Edgemont Road
// //#                        Charlottesville, VA 22903-2475 USA
// //#
// //# $Id$
// //
#ifndef SYNTHESIS_TRANSFORM2_VB2CFBMAP_H
#define SYNTHESIS_TRANSFORM2_VB2CFBMAP_H

#include <casa/Arrays/Array.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <synthesis/TransformMachines2/CFBuffer.h>
#include <synthesis/TransformMachines2/CFStore2.h>
#include <synthesis/TransformMachines2/PhaseGrad.h>
#include <casa/Utilities/CountedPtr.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/TempImage.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <synthesis/MeasurementComponents/SolvableVisCal.h>
using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN
 namespace refim{
   using namespace CFDefs;
   class VB2CFBMap
   {
   public:
     VB2CFBMap(): vbRow2CFBMap_p(){vbRow2CFBMap_p(0) = 0;}
     
     ~VB2CFBMap() {};
     
     VB2CFBMap& operator=(const VB2CFBMap& other);
     casacore::CountedPtr<CFBuffer >& operator[](const int& i) {return vbRow2CFBMap_p[i];};
     
     inline casacore::Vector<casacore::CountedPtr<CFBuffer>>& getVBRow2CFBMap() {return vbRow2CFBMap_p;};	
     inline int nelements() {return vbRow2CFBMap_p.nelements();}

     virtual casacore::Int mapAntIDToAntType(const casacore::Int& /*ant*/) {return 0;};
     virtual casacore::Int makeVBRow2CFBMap(CFStore2& cfs,
					    const VisBuffer2& vb, const casacore::Quantity& dPA,
					    const casacore::Vector<casacore::Int>& dataChan2ImChanMap,
					    const casacore::Vector<casacore::Int>& dataPol2ImPolMap,
					    const casacore::Vector<casacore::Double>& pointingOffset);
     //   protected:
     casacore::Vector<casacore::CountedPtr<CFBuffer > > vbRow2CFBMap_p;
     
     //	VBRow2CFBMapType vbRow2CFBMap_p;
   };
 }
}
#endif	
