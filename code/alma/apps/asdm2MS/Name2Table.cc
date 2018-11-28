
/*
 * ALMA - Atacama Large Millimeter Array
 * (c) European Southern Observatory, 2002
 * (c) Associated Universities Inc., 2002
 * Copyright by ESO (in the framework of the ALMA collaboration),
 * Copyright by AUI (in the framework of the ALMA collaboration),
 * All rights reserved.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY, without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307  USA
 *
 * Warning!
 *  -------------------------------------------------------------------- 
 * | This is generated code!  Do not modify this file.                  |
 * | If you do, all changes will be lost when the file is re-generated. |
 *  --------------------------------------------------------------------
 *
 * File Name2Table.cpp
 */
#include <alma/apps/asdm2MS/ASDMTables.h>
#include <alma/apps/asdm2MS/Name2Table.h>

#include <string>
#include <regex>
#include <sstream>
#include <algorithm>
#include <iterator>

#include <casa/Logging.h>

using namespace std;
using namespace casacore;

map<string, ASDM_TABLE_BASE*> Name2Table::name2Table_;
bool Name2Table::init_ = Name2Table::init();
set<ASDM_TABLE_BASE*> Name2Table::table_;

bool Name2Table::init() {
	
	  name2Table_["Main"] = ASDM_MAIN::instance();	
	
	  name2Table_["AlmaRadiometer"] = ASDM_ALMARADIOMETER::instance();	
	
	  name2Table_["Annotation"] = ASDM_ANNOTATION::instance();	
	
	  name2Table_["Antenna"] = ASDM_ANTENNA::instance();	
	
	  name2Table_["CalAmpli"] = ASDM_CALAMPLI::instance();	
	
	  name2Table_["CalAntennaSolutions"] = ASDM_CALANTENNASOLUTIONS::instance();	
	
	  name2Table_["CalAppPhase"] = ASDM_CALAPPPHASE::instance();	
	
	  name2Table_["CalAtmosphere"] = ASDM_CALATMOSPHERE::instance();	
	
	  name2Table_["CalBandpass"] = ASDM_CALBANDPASS::instance();	
	
	  name2Table_["CalCurve"] = ASDM_CALCURVE::instance();	
	
	  name2Table_["CalData"] = ASDM_CALDATA::instance();	
	
	  name2Table_["CalDelay"] = ASDM_CALDELAY::instance();	
	
	  name2Table_["CalDevice"] = ASDM_CALDEVICE::instance();	
	
	  name2Table_["CalFlux"] = ASDM_CALFLUX::instance();	
	
	  name2Table_["CalFocus"] = ASDM_CALFOCUS::instance();	
	
	  name2Table_["CalFocusModel"] = ASDM_CALFOCUSMODEL::instance();	
	
	  name2Table_["CalGain"] = ASDM_CALGAIN::instance();	
	
	  name2Table_["CalHolography"] = ASDM_CALHOLOGRAPHY::instance();	
	
	  name2Table_["CalPhase"] = ASDM_CALPHASE::instance();	
	
	  name2Table_["CalPointing"] = ASDM_CALPOINTING::instance();	
	
	  name2Table_["CalPointingModel"] = ASDM_CALPOINTINGMODEL::instance();	
	
	  name2Table_["CalPosition"] = ASDM_CALPOSITION::instance();	
	
	  name2Table_["CalPrimaryBeam"] = ASDM_CALPRIMARYBEAM::instance();	
	
	  name2Table_["CalReduction"] = ASDM_CALREDUCTION::instance();	
	
	  name2Table_["CalSeeing"] = ASDM_CALSEEING::instance();	
	
	  name2Table_["CalWVR"] = ASDM_CALWVR::instance();	
	
	  name2Table_["ConfigDescription"] = ASDM_CONFIGDESCRIPTION::instance();	
	
	  name2Table_["CorrelatorMode"] = ASDM_CORRELATORMODE::instance();	
	
	  name2Table_["DataDescription"] = ASDM_DATADESCRIPTION::instance();	
	
	  name2Table_["DelayModel"] = ASDM_DELAYMODEL::instance();	
	
	  name2Table_["DelayModelFixedParameters"] = ASDM_DELAYMODELFIXEDPARAMETERS::instance();	
	
	  name2Table_["DelayModelVariableParameters"] = ASDM_DELAYMODELVARIABLEPARAMETERS::instance();	
	
	  name2Table_["Doppler"] = ASDM_DOPPLER::instance();	
	
	  name2Table_["Ephemeris"] = ASDM_EPHEMERIS::instance();	
	
	  name2Table_["ExecBlock"] = ASDM_EXECBLOCK::instance();	
	
	  name2Table_["Feed"] = ASDM_FEED::instance();	
	
	  name2Table_["Field"] = ASDM_FIELD::instance();	
	
	  name2Table_["Flag"] = ASDM_FLAG::instance();	
	
	  name2Table_["FlagCmd"] = ASDM_FLAGCMD::instance();	
	
	  name2Table_["Focus"] = ASDM_FOCUS::instance();	
	
	  name2Table_["FocusModel"] = ASDM_FOCUSMODEL::instance();	
	
	  name2Table_["FreqOffset"] = ASDM_FREQOFFSET::instance();	
	
	  name2Table_["GainTracking"] = ASDM_GAINTRACKING::instance();	
	
	  name2Table_["History"] = ASDM_HISTORY::instance();	
	
	  name2Table_["Holography"] = ASDM_HOLOGRAPHY::instance();	
	
	  name2Table_["Observation"] = ASDM_OBSERVATION::instance();	
	
	  name2Table_["Pointing"] = ASDM_POINTING::instance();	
	
	  name2Table_["PointingModel"] = ASDM_POINTINGMODEL::instance();	
	
	  name2Table_["Polarization"] = ASDM_POLARIZATION::instance();	
	
	  name2Table_["Processor"] = ASDM_PROCESSOR::instance();	
	
	  name2Table_["Pulsar"] = ASDM_PULSAR::instance();	
	
	  name2Table_["Receiver"] = ASDM_RECEIVER::instance();	
	
	  name2Table_["SBSummary"] = ASDM_SBSUMMARY::instance();	
	
	  name2Table_["Scale"] = ASDM_SCALE::instance();	
	
	  name2Table_["Scan"] = ASDM_SCAN::instance();	
	
	  name2Table_["Seeing"] = ASDM_SEEING::instance();	
	
	  name2Table_["Source"] = ASDM_SOURCE::instance();	
	
	  name2Table_["SpectralWindow"] = ASDM_SPECTRALWINDOW::instance();	
	
	  name2Table_["SquareLawDetector"] = ASDM_SQUARELAWDETECTOR::instance();	
	
	  name2Table_["State"] = ASDM_STATE::instance();	
	
	  name2Table_["Station"] = ASDM_STATION::instance();	
	
	  name2Table_["Subscan"] = ASDM_SUBSCAN::instance();	
	
	  name2Table_["SwitchCycle"] = ASDM_SWITCHCYCLE::instance();	
	
	  name2Table_["SysCal"] = ASDM_SYSCAL::instance();	
	
	  name2Table_["SysPower"] = ASDM_SYSPOWER::instance();	
	
	  name2Table_["TotalPower"] = ASDM_TOTALPOWER::instance();	
	
	  name2Table_["WVMCal"] = ASDM_WVMCAL::instance();	
	
	  name2Table_["Weather"] = ASDM_WEATHER::instance();	
	
  return true;
}

const set<ASDM_TABLE_BASE*>& Name2Table::find (const vector<string>& name, bool verbose) {

  // LogSinkInterface& lsif = LogSink::globalSink();
  ostringstream infostream; 

  // Empty the set just to be sure.
  table_.clear();
  
  // Process each string given as input
  for (unsigned int i = 0; i < name.size(); i++) {
    // Replace any possible wildcard (*) by a sequences '(.)*'
    stringstream ss(name.at(i));
    string token;
    bool first = true;
    ostringstream oss;
    while (getline(ss,token,'*')) {
      if (!first) {
        oss << "(.)*";
      } else {
        first = false;
      }
      oss << token;
    }
    // watch for a trailing wildcard, will be missed by the above
    if (name.back()=='*') oss << "(.)*";

    // Build a regexp out of this string written onto oss.
    std::regex expression(oss.str().c_str()); 

    // For each table name
    for (map<string, ASDM_TABLE_BASE*>::const_iterator iter = name2Table_.begin();
	 iter!= name2Table_.end(); iter++) {
      // Is there a match between the string and the name of the table ?
      std::cmatch what; 
      if(std::regex_match(iter->first.c_str(), what, expression)) {
      	string uppername = iter->first;
        std::transform(uppername.begin(),uppername.end(), uppername.begin(), [](unsigned char c) { return std::toupper(c);});
      	if (verbose) {
      		infostream.str("");
      		infostream << "An ASDM_" << uppername << " table will be added to the MS" << endl;
      		LogSink::postGlobally(LogMessage(infostream.str(),
					              			 LogOrigin("Name2Table", 
					                                   "find (const vector<string>& name, bool verbose=false)", 
					                                   WHERE)));
      	}
		table_.insert(iter->second);
      }
    }
  }
  return table_;
}
                                    
