//# PlotMSCTAverager.h: class to average CT iterations for PlotMS
//# Copyright (C) 2020
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

#ifndef PLOTMSCTAVERAGER_H
#define PLOTMSCTAVERAGER_H

#include <plotms/PlotMS/PlotMSAveraging.h>
#include <synthesis/CalTables/CTIter.h>

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// A class to average NewCalTable chunks for PlotMS
// </summary>
//
// <use visibility=export>
//
// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>

// <prerequisite>
//   <li> CTIter
// </prerequisite>
//
// <etymology>
// </etymology>
//
// <synopsis>
// This class averages CalTable chunks together for PlotMS
// </synopsis>
//
// <example>
// </example>
//
// <motivation>
// </motivation>
//
// <thrown>
//    <li>
//    <li>
// </thrown>
//

class PlotMSCTAverager
{
public:
  // Construct from the averaging enabled, number of antennas and polarizations
  PlotMSCTAverager(
    PlotMSAveraging& averaging, casacore::Int nAnt, casacore::Int nPoln);

  // Null destructor
  ~PlotMSCTAverager();

  // Set debug messages to stdout
  inline void setDebug(casacore::Bool debug) {
    debug_ = debug;
  };

  // Cal table is baseline-based; assumes antenna-based unless this is set
  inline void setBaselineBased() {
    isAntennaBased_p = false;
  };

  // Accumulate a chunk
  inline void accumulate (ROCTIter& cti, std::vector<casacore::Slice>& chansel) {
    averaging_p.antenna() ? antennaAccumulate(cti, chansel) : simpleAccumulate(cti, chansel);
  };

  // Finalize averaging by writing values to CTMainRecord vector
  void finalizeAverage();

  // Return the result in NewCalTable filled with CTMainRecord vector
  void fillAvgCalTable(NewCalTable& tab);

  // CTMainRecord does not include chan or freq
  casacore::Int nchan();
  inline casacore::Vector<casacore::Int> chan() { return avgChan_; };
  inline casacore::Vector<casacore::Double> freq() { return avgFreq_; };

  // For channel averaging: for each output channel, row is vector of which channels were averaged together
  inline casacore::Array<casacore::Int> chansPerBin() { return chansPerBin_; }; // [nchan, nbin]
private:
  // Prohibit null constructor, copy constructor and assignment for now
  PlotMSCTAverager();
  PlotMSCTAverager& operator= (const PlotMSCTAverager&);
  PlotMSCTAverager (const PlotMSCTAverager&);

  // Initialize the next accumulation
  void initialize(ROCTIter& cti, std::vector<casacore::Slice>& chansel);

  // Different accumulate versions
  void simpleAccumulate (ROCTIter& cti, std::vector<casacore::Slice>& chansel);  // ordinary averaging
  void antennaAccumulate (ROCTIter& cti, std::vector<casacore::Slice>& chansel); // antenna-based averaging

  // Hash function to return a row index for a baseline-based cal table;
  // Returns 0 if baseline averaging, ant1 for antenna-based table
  casacore::Int baseline_index(const casacore::Int& ant1, const casacore::Int& ant2);

  // Convert r/i to a/p
  void convertToAP(casacore::Cube<casacore::Complex>& d);

  // Input averaging options
  PlotMSAveraging averaging_p;

  // Number of antennas, polarizations, and baselines
  casacore::Int nAnt_p, nPoln_p, nBlnMax_p;

  // Number of channels, selected channels, averaged channels
  casacore::Int nChan_p, nAvgChan_p, nSelChan_p;

  // For channel averaging and selection
  casacore::Bool avgChan_p, selChan_p;
  casacore::Int nChanPerBin_p;

  // Data is Complex or Float
  casacore::Bool isComplex_p;
  // Only ANTENNA1 varies
  casacore::Bool isAntennaBased_p;

  // Validation by baseline (if false, no attempt to accumulate this baseline)
  casacore::Vector<casacore::Bool> blnOK_p;

  // Accumulation helpers...
  // Used to set the time value for averaged chunk
  casacore::Double timeRef_p; // first time in averaged chunk
  casacore::Double minTimeOffset_p;
  casacore::Double maxTimeOffset_p;
  // Values for averaged chunk
  casacore::Int avgScan_p, avgField_p, avgSpw_p;

  // Keep track of initialization state
  casacore::Bool initialized_p;
 
  // Keep track of accumulation state; cannot finalize if no accumulation
  casacore::Bool isAccum_p;

  // Diagnostic print
  casacore::Bool debug_;

  // Mutable arrays, set in CTMainRecords when finalized
  // avg* vectors may be set to -1 if values are combined
  casacore::Vector<casacore::Int> avgScan_;
  casacore::Vector<casacore::Int> avgField_;
  casacore::Vector<casacore::Double> avgTime_;
  casacore::Vector<casacore::Int> avgSpw_;
  casacore::Vector<casacore::Int> avgChan_;
  casacore::Vector<casacore::Double> avgFreq_;
  casacore::Vector<casacore::Int> obsid_;
  casacore::Vector<casacore::Int> avgAntenna1_;
  casacore::Vector<casacore::Int> avgAntenna2_;
  casacore::Vector<casacore::Int> avgBaseline_;

  // Accumulated in accumulate()
  casacore::Cube<casacore::Complex> accumCParam_;
  casacore::Cube<casacore::Float> accumFParam_;
  casacore::Cube<casacore::Float> accumParamErr_;
  casacore::Cube<casacore::Float> accumSnr_;
  casacore::Cube<casacore::Float> accumWt_;
  casacore::Cube<casacore::Bool> avgFlag_;
  casacore::Vector<casacore::Double> accumFreq_; // if channel averaging
  casacore::Vector<casacore::Int> nAccumFreq_;   // if channel averaging

  // For channel averaging: nrow is number of averaged channels,
  // each row is vector of channels averaged into each averaged channel
  casacore::Array<casacore::Int> chansPerBin_;

  // Averaged results
  std::vector<CTMainRecord> main_rows_;
 
};


} //# NAMESPACE CASA - END

#endif


