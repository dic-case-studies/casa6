/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version December 2010. 
   Maintained by ESO since 2013.

   This file is part of LibAIR and is licensed under GNU Public
   License Version 2
   
   \file wvrgcalfeedback.hpp

   Feedback to the user 

*/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "../src/libair_main.hpp"
#include "wvrgcalfeedback.hpp"


namespace LibAIR2 {

  void fatalMsg(const std::string &m)
  {
    std::cout<<"FATAL: "<<m
	     <<std::endl;
  }

  void errorMsg(const std::string &m)
  {
    std::cout<<"ERR: "<<m
	     <<std::endl;
  }

  void warnMsg(const std::string &m)
  {
    std::cout<<"WARN: "<<m
	     <<std::endl;
  }

  void printBanner(std::ostream& /*os*/)
  {
    std::cout<<std::endl
	     <<"WVRGCAL  -- Version "
	     <<LibAIR2::version()
	     <<std::endl
	     <<std::endl
	     <<"Developed by B. Nikolic at the University of Cambridge as part of EU FP6 ALMA Enhancement"
	     <<std::endl
	     <<"Maintained and extended since 2013 by the European Southern Observatory as part of the ALMA project"
	     <<std::endl
	     <<"GPLv2 License -- you have a right to the source code"
	     <<std::endl
	     <<std::endl;
    
  }

  void printUsedStates(const std::set<size_t> &useID)
  {
    std::cout<<std::endl
	     <<"Good state IDs: ";
    BOOST_FOREACH(const size_t &x, useID)
      std::cout<<x<<", ";
    std::cout<<std::endl;
    
  }


  std::ostream & boolYesNo(std::ostream &os, bool b)
  {
    if (b)
      os<<"Yes";
    else
      os<<"No";
    return os;
  }

  std::ostream & operator<<(std::ostream &os, const AntennaInfo &ai)
  {
    os<<ai.no<<"\t"<<ai.name<<"\t";
    boolYesNo(os, ai.haswvr);
    os<<"\t";
    boolYesNo(os,ai.flag);
    os<<boost::format("\t%4.3g") % (ai.pathRMS/1e-6)
      <<"\t"
      <<boost::format("\t%4.3g") % (ai.pathDisc/1e-6);
    return os;
  }

  AntITable::AntITable(const aname_t &names,
		       const LibAIR2::AntSet &flag,
                       const LibAIR2::AntSet &nowvr,                       
		       const std::vector<double> &rms,
		       const std::vector<double> &disc,
		       const LibAIR2::AntSet &interpolImpossibleAnts)
  {
    for(aname_t::const_iterator i=names.begin(); i !=names.end(); ++i)
    {
      AntennaInfo x;
      x.no= i->first;
      x.name= i->second;
      x.haswvr=(not (nowvr.count(x.no)));
      x.flag=flag.count(x.no);
      x.pathRMS=rms[x.no];
      x.pathDisc=disc[x.no];
      if((x.flag>0) && (interpolImpossibleAnts.count(x.no)>0))
      {
	x.pathRMS = 0.;
	x.pathDisc = 0.;
      }
      push_back(x);
    }
  }

  std::ostream & operator<<(std::ostream &os, const AntITable &at)
  {
    os<<"     Antenna/WVR information:                     "<<std::endl
      <<"-----------------------------------------------------------------------"<<std::endl;
    os<<"#"<<"\t"<<"Name"<<"\t"<<"WVR?"<<"\t"<<"Flag?"<<"\t"<<"RMS (um)" << "\t"<<"Disc (um)" 
      <<std::endl;
    BOOST_FOREACH(const AntennaInfo &x, at)
      os<<x<<std::endl;
    return os;
  }

  void printStatTimes(std::ostream &os,
		      const std::vector<double> &time,
		      const std::vector<std::pair<double, double> > &tmask)
  {
    double tbase=time[0];
    os<<"Times used for the statistics calculation (in seconds from first astro datum) "<<std::endl
      <<"-----------------------------------------------------------------------------"<<std::endl;
    for(size_t i=0; i<tmask.size(); ++i)
    {
      os<<"("<<tmask[i].first-tbase<<", "<<tmask[i].second-tbase<<")"
	<<std::endl;
    }

  }

}


