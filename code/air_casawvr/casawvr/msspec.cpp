/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version June 2010. 
   Maintained by ESO since 2013.

   This file is part of LibAIR and is licensed under GNU Public
   License Version 2
   
   \file msspec.cpp
   
*/

#include "msspec.hpp"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>

namespace LibAIR2 {

  void loadSpec(const casacore::MeasurementSet &ms,
		const std::vector<int> &spws,
		MSSpec &s)
  {
    casacore::MSSpectralWindow specTable(ms.spectralWindow());
    // Number of spectral windows in this measurement set

    const size_t nspw=specTable.nrow();

    // Number of channels 
    casacore::ScalarColumn<casacore::Int> nc
      (specTable,
       casacore::MSSpectralWindow::columnName(casacore::MSSpectralWindow::NUM_CHAN));

    // Frequencies of the channels
    casacore::ArrayColumn<casacore::Double>
      chfreq(specTable,
	     casacore::MSSpectralWindow::columnName(casacore::MSSpectralWindow::CHAN_FREQ));

    std::vector<int> spwsout = spws;
    size_t nspwout = spwsout.size();
    if(nspw==0){
      for(size_t i=0; i<nspw; i++){
	if(nc(i)!=4){
	  spwsout.push_back(i);
	}
      }
      nspwout=spwsout.size();
    }

    s.spws.resize(nspwout);

    for (size_t ispw=0; ispw<nspwout; ++ispw)
    {
      size_t spwid=spwsout[ispw];

      s.spws[ispw].spwid=spwid;
      casacore::Array<casacore::Double> freq;
      chfreq.get(spwid, freq, casacore::True);
      s.spws[ispw].chf.resize(nc(spwid));
      for(size_t i=0; i< static_cast<size_t>(nc(spwid)); ++i)
      {
	s.spws[ispw].chf[i]=freq(casacore::IPosition(1,i));
      }
    }
  }

  std::ostream & 
  operator<<(std::ostream &os,
	     const MSSpec &s)
  {
    os<<"This MS has "<<s.spws.size()<<" spectral windows"
      <<std::endl;
    for (size_t spw=0; spw<s.spws.size(); ++spw)
    {
      os<<"SPW "<<spw<<" has "<<s.spws[spw].chf.size()<<" channels ";
      if (s.spws[spw].chf.size()==1)
      {
	os<<"at frequency "<<s.spws[spw].chf[0];
      }
      else
      {
	os<<"starting at frequency "<<s.spws[spw].chf[0]
	  <<" and with last at frequency "<<s.spws[spw].chf[s.spws[spw].chf.size()-1];
      }
      os<<std::endl;
    }
    return os;
  }

  std::map<size_t, size_t>
  SPWDataDescMap(const casacore::MeasurementSet &ms)
  {
    std::map<size_t, size_t> res;
    const casacore::MSDataDescription dd(ms.dataDescription());
    const size_t n=dd.nrow();

    casacore::ScalarColumn<casacore::Int>
      spwid(dd,
	    casacore::MSDataDescription::columnName(casacore::MSDataDescriptionEnums::SPECTRAL_WINDOW_ID));


    for(size_t i=0; i<n; ++i)
    {
      int spw;
      spwid.get(i, spw);
      res.insert(std::pair<size_t, size_t>(spw, i));
    }
    return res;
  }

  std::map<size_t, size_t>
  DataDescSPWMap(const casacore::MeasurementSet &ms)
  {
    std::map<size_t, size_t> res;
    const casacore::MSDataDescription dd(ms.dataDescription());
    const size_t n=dd.nrow();

    casacore::ScalarColumn<casacore::Int>
      spwid(dd,
	    casacore::MSDataDescription::columnName(casacore::MSDataDescriptionEnums::SPECTRAL_WINDOW_ID));


    for(size_t i=0; i<n; ++i)
    {
      int spw;
      spwid.get(i, spw);
      res.insert(std::pair<size_t, size_t>(i, spw));
    }
    return res;
  }

  void dataSPWs(const casacore::MeasurementSet &ms,
		std::vector<size_t> &spw,
		const std::vector<size_t> &sortedI)
  {
    std::map<size_t, size_t> map=DataDescSPWMap(ms);
    const casacore::ROMSMainColumns cols(ms);
    const casacore::ScalarColumn<casacore::Int> &dd=cols.dataDescId();
    const size_t nrows=dd.nrow();          
    spw.resize(nrows);
    for(size_t ii=0; ii<nrows; ++ii)
    {
      size_t i = sortedI[ii];

      spw[ii]=map[dd(i)];
    }
  }

  size_t numSPWs(const casacore::MeasurementSet &ms){
    casacore::MSSpectralWindow specTable(ms.spectralWindow());

    return specTable.nrow();
  }


}



