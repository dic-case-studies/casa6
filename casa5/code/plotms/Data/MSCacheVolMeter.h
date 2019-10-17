//# MSCacheVolMeter.h: Definition of MSCache Volume meter
//# Copyright (C) 2009
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
//# $Id: $
#ifndef MSCACHEVOLMETER_H_
#define MSCACHEVOLMETER_H_

#include <plotms/PlotMS/PlotMSConstants.h>
#include <plotms/PlotMS/PlotMSAveraging.h>

#include <casa/aips.h>
#include <casa/Arrays.h>
#include <casa/Containers/Block.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <ms/MeasurementSets/MeasurementSet.h>

namespace casa {

class MSCacheVolMeter {

public:

  // Constructor/Destructor
  MSCacheVolMeter();
  MSCacheVolMeter(const casacore::MeasurementSet& ms, const PlotMSAveraging ave,
		   const casacore::Vector<casacore::Vector<casacore::Slice> >& chansel,
		   const casacore::Vector<casacore::Vector<casacore::Slice> >& corrsel);
  ~MSCacheVolMeter();

  // reset (as if default ctor was run)
  void reset();

  // add in via a VisBuffer
  void add(const vi::VisBuffer2* vb);

  // add in via counts
  void add(casacore::Int DDID,casacore::Int nRows);

  // evaluate the volume for specified axes, and complain if 
  casacore::String evalVolume(std::map<PMS::Axis,casacore::Bool> axes,casacore::Vector<casacore::Bool> axesmask);
  casacore::String evalVolume(std::vector<casacore::IPosition> vbShapes, std::map<PMS::Axis,casacore::Bool> axes);

private:

  // The number of DATA_DESCRIPTIONs
  casacore::Int nDDID_;

  // Counters
  casacore::Vector<casacore::uInt64> nPerDDID_,nRowsPerDDID_,nChanPerDDID_,nCorrPerDDID_;

  // The number of antennas (max)
  casacore::Int nAnt_;

};

}

#endif /* MSCACHEVOLMETER_H_ */
