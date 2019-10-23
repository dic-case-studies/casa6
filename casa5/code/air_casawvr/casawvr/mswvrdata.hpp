/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version January 2010.
   Maintained by ESO since 2013. 
   
   This file is part of LibAIR and is licensed under GNU Public
   License Version 2
   
   \file mswvrdata.hpp

   Tools for handling of WVR data in measurement sets
*/

#ifndef _LIBAIR_CASAWVR_MSWVRDATA_HPP__
#define _LIBAIR_CASAWVR_MSWVRDATA_HPP__

#include <set>
#include <vector>

#include "../src/apps/antennautils.hpp"
#include <ms/MeasurementSets/MeasurementSet.h>

namespace LibAIR2 {

  // Forward declarations
  class InterpArrayData;
  

  /** \brief Set of spectral windows
   */
  typedef std::set<size_t> SPWSet;

  /** Set of spectral windows associated with ALMA WVRs

      \bug At the moment I'm using an empirical algorithms that all
      SPWs with 4 channels are WVRs
   */
  SPWSet
  WVRSPWIDs(const casacore::MeasurementSet &ms);

  /** Set of DATA_DESC_IDs associated with WVRs
      among those given by spws
   */
  std::set<size_t>
  WVRDataDescIDs(const casacore::MeasurementSet &ms,
		 const std::vector<int> &spws);

  /** Number of spectral windows associated with ALMA WVRs
   */
  size_t nWVRSPWIDs(const casacore::MeasurementSet &ms);

  /** Set of antennas with data from WVRs available
   */
  AntSet
  WVRAntennas(const casacore::MeasurementSet &ms,
	      const std::vector<int> &wvrspws);

  /** Set of antennas with data from WVRs available
      
      This version uses the Feed table to figure out which antennas
      have WVRs
   */
  AntSet
  WVRAntennasFeedTab(const casacore::MeasurementSet &ms,
		     const std::vector<int> &wvrspws);

  /** Set of antennas with data from WVRs available
      
      This version obtains the set by iterating through the main
      table. This is slow, but allows correct result when the
      feedtable is corrupted
   */
  AntSet
  WVRAntennasMainTab(const casacore::MeasurementSet &ms,
		     const std::vector<int> &wvrspws);

  /** Add the antennas flagged in the ANTENNA table to the set
   */
  void WVRAddFlaggedAnts(const casacore::MeasurementSet &ms,
			 LibAIR2::AntSet &flaggedAnts);

  /** Time points, states, and field IDs at which WVR data have been
      recorded

   */
  void WVRTimeStatePoints(const casacore::MeasurementSet &ms,
			  std::vector<double> &times,
			  std::vector<size_t> &states,
			  std::vector<size_t> &field,
			  std::vector<size_t> &source,
			  const std::vector<int> &wvrspws,
			  const std::vector<size_t>& sortedI);

  /** Load all WVR data from a measurment set
      
   */
  InterpArrayData *loadWVRData(const casacore::MeasurementSet &ms,
			       const std::vector<int>& spws, 
			       std::vector<size_t> &sortedI,
			       std::set<int> &flaggedantsInMain,
			       double requiredUnflaggedFraction=0.8,
			       bool usepointing=true,
			       std::string offsetstable="");

}

#endif
