/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version January 2010.
   Maintained by ESO since 2013. 
   
   This file is part of LibAIR and is licensed under GNU Public
   License Version 2
   
   \file msutils.hpp

   Utilties for dealing with measurement sets and functions to extract
   non-WVR data from them. For handling of WVR data see mswvrdata.hpp
*/
#ifndef _LIBAIR_CASAWVR_MSUTILS_HPP__
#define _LIBAIR_CASAWVR_MSUTILS_HPP__

#include <vector>
#include <string>
#include <map>
#include <set>

#include <ms/MeasurementSets/MeasurementSet.h>

namespace LibAIR2 {

  /** Channel frequencies for spectral window spw     
   */
  void spwChannelFreq(const casacore::MeasurementSet &ms,
		      size_t spw,
		      std::vector<double> &fres);


  /** Retrieve row times, field IDs, and Source IDs
   */
  void fieldIDs(const casacore::MeasurementSet &ms,
		std::vector<double> &time,
		std::vector<int> &fieldID,
		std::vector<int> &sourceID,
		const std::vector<size_t> &sortedI
		);

  /** \brief Connection between state_id and the ScanIntent (or
      ObsMode) 
   */
  struct StateIntentMap:
    public std::map<size_t, std::string>
  {

  public:

  };

  std::ostream & operator<<(std::ostream &o, 
			    const StateIntentMap &t);

  /** \brief Load the map between state ids and scan intents from
      measurement set*/
  void scanIntents(const casacore::MeasurementSet &ms,
		   StateIntentMap &mi);

  /** Returns the set of State IDs that are known to have the WVR
      looking at the sky
   */
  std::set<size_t> skyStateIDs(const casacore::MeasurementSet &ms);

  /// Association between field number and field names
  typedef std::map<size_t, std::string > field_t;
  field_t getFieldNames(const casacore::MeasurementSet &ms);


  /// Association between source IDs and source names
  field_t getSourceNames(const casacore::MeasurementSet &ms);

  /** Get all fields associated with a particular source name
   */
  std::set<size_t> getSrcFields(const casacore::MeasurementSet &ms,
				const std::string &source);

  /** Map from field IDs to Source Ids. One way only as multiple
      fields can correspond to one source (thats the way CASA does
      it).
   */
  std::map<size_t, size_t> getFieldSrcMap(const casacore::MeasurementSet &ms);
      


}

#endif
