//# RFChunkStats.cc: this defines RFChunkStats
//# Copyright (C) 2000,2001,2002,2003
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
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Exceptions/Error.h>
#include <flagging/Flagging/Flagger.h>
#include <flagging/Flagging/RFChunkStats.h>
#include <msvis/MSVis/VisibilityIterator.h>
#include <msvis/MSVis/VisBuffer.h>
#include <casacore/casa/System/PGPlotter.h>
    

namespace casa { //# NAMESPACE CASA - BEGIN

template<class T> RFlagWord RFChunkStats::getCorrMask ( const casacore::Vector<T> &corrspec )
{
  RFlagWord mask=0;
  // loop over polzn spec
  for( casacore::uInt i=0; i<corrspec.nelements(); i++ )
  {
    // convert element of polspec to casacore::Stokes type
    casacore::Stokes::StokesTypes type = casacore::Stokes::type( corrspec(i) );
    if( type == casacore::Stokes::Undefined ){
      std::ostringstream oss;
      oss << corrspec(i);
      throw(casacore::AipsError( casacore::String("Unknown correlation type: ")+ casacore::String(oss)));
    }
    // find this type in current corrarizations
    casacore::Int icorr = findCorrType(type,corrtypes);
    if( icorr>=0 )
      mask |= (1<<icorr);
  }
  return mask;
}

} //# NAMESPACE CASA - END

