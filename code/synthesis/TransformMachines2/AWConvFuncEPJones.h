// -*- C++ -*-
//# AWConvFunc.h: Definition of the AWConvFunc class
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
#ifndef SYNTHESIS_TRANSFORM2_AWCONVFUNCEPJONES_H
#define SYNTHESIS_TRANSFORM2_AWCONVFUNCEPJONES_H

#include <synthesis/TransformMachines2/AWConvFunc.h>

namespace casa { //# NAMESPACE CASA - BEGIN
  namespace refim{
  //
  //-------------------------------------------------------------------------------------------
  //
  class AWConvFuncEPJones : public AWConvFunc
  {
  public:
    AWConvFuncEPJones(const casacore::CountedPtr<ATerm> ATerm,
		      const casacore::CountedPtr<PSTerm> psTerm,
		      const casacore::CountedPtr<WTerm> wTerm,
		      const casacore::Bool wbAWP=false,
		      const casacore::Bool conjPB=true):
      AWConvFunc(ATerm,psTerm,wTerm,wbAWP, conjPB)
    {};
    ~AWConvFuncEPJones() {};
    AWConvFuncEPJones& operator=(const AWConvFuncEPJones& other);

    virtual void makeConvFunction(const casacore::ImageInterface<casacore::Complex>& image,
				  const VisBuffer2& vb,
				  const casacore::Int wConvSize,
				  const casacore::CountedPtr<PolOuterProduct>& pop,
				  const casacore::Float pa,
				  const casacore::Float dpa,
				  const casacore::Vector<casacore::Double>& uvScale, const casacore::Vector<casacore::Double>& uvOffset,
				  const casacore::Matrix<casacore::Double>& vbFreqSelection,
				  CFStore2& cfs,
				  CFStore2& cfwts,
				  casacore::Bool fillCF=true);
    virtual casacore::Vector<casacore::Double> findPointingOffset(const casacore::ImageInterface<casacore::Complex>& image,
					      const VisBuffer2& vb);

  };
  //
  //-------------------------------------------------------------------------------------------
  //
};
};
#endif
