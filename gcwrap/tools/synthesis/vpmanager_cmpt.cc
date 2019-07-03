
/***
 * Framework independent implementation file for vpmanager...
 *
 * Implement the vpmanager component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 *
 * @author
 * @version 
 ***/

#include <iostream>
#include <vpmanager_cmpt.h>
#include <casa/Logging/LogIO.h>
#include <casa/BasicSL/String.h>
#include <casa/OS/Path.h>
#include <casa/OS/Directory.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/Assert.h>
#include <casa/Utilities/GenSort.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>
#include <casa/Quanta.h>
#include <casa/Quanta/MVTime.h>
#include <casa/Quanta/MVFrequency.h>
#include <casa/Quanta/Quantum.h>
#include <synthesis/MeasurementEquations/VPManager.h>
#include <synthesis/MeasurementEquations/AntennaResponses.h>

using namespace std;
using namespace casacore;
using namespace casa;

using namespace casacore;
namespace casac {

vpmanager::vpmanager()
{

  VPManager::Instance();
  itsLog = new LogIO();

}

vpmanager::~vpmanager()
{
  if(itsLog) delete itsLog;
  itsLog=0;
}

bool
vpmanager::saveastable(const std::string& tablename)
{
  bool rstat=false;
  try{
    rstat=VPManager::Instance()->saveastable(String(tablename));
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
vpmanager::loadfromtable(const std::string& tablename)
{
  bool rstat=false;
  try{
    rstat=VPManager::Instance()->loadfromtable(String(tablename));

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
vpmanager::summarizevps(const bool verbose)
{
  bool rstat=false;

  try {
    rstat = VPManager::Instance()->summarizevps(verbose);
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}


//TODO!
// list
// vpmanager::list_known_telescopes()
// {
//   return me.   // Use pbmath instead.
// }

::casac::record*
vpmanager::setcannedpb(const std::string& telescope, const std::string& othertelescope, const bool dopb, const std::string& commonpb, const bool dosquint, const ::casac::variant& paincrement, const bool usesymmetricbeam)
{
  ::casac::record* r=0;
  try{
    casacore::Quantity pai;
    if(toCasaString(paincrement)==casacore::String(""))
      pai=casacore::Quantity(720,"deg");
    else
      pai=casaQuantity(paincrement);

    casacore::Record rec;
    VPManager::Instance()->setcannedpb(String(telescope), String(othertelescope), dopb, String(commonpb),
			dosquint, pai, usesymmetricbeam, rec);
    r=fromRecord(rec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return r;
}

::casac::record*
vpmanager::setpbairy(const std::string& telescope, const std::string& othertelescope, const bool dopb, const ::casac::variant& dishdiam, const ::casac::variant& blockagediam, const  ::casac::variant& maxrad, const ::casac::variant& reffreq, const ::casac::variant& squintdir, const ::casac::variant& squintreffreq, const bool dosquint, const ::casac::variant& paincrement, const bool usesymmetricbeam)
{

  ::casac::record* r=0;
  try{
    casacore::Quantity rf;
    if(toCasaString(reffreq)==casacore::String(""))
      rf=casacore::Quantity(1.0,"GHz");
    else
      rf=casaQuantity(reffreq);
    casacore::Quantity srf;
    if(toCasaString(squintreffreq)==casacore::String(""))
      srf=casacore::Quantity(1.0,"GHz");
    else
      srf=casaQuantity(squintreffreq);
    casacore::Quantity pai;
    if(toCasaString(paincrement)==casacore::String(""))
      pai=casacore::Quantity(720,"deg");
    else
      pai=casaQuantity(paincrement);
    casacore::MDirection sdir;
    if(toCasaString(squintdir)==casacore::String(""))
      sdir=casacore::MDirection(casacore::Quantity(0,"deg"),casacore::Quantity(0,"deg")) ;
    else
      casaMDirection(squintdir, sdir);
    casacore::Quantity dd;
    if(toCasaString(dishdiam)==casacore::String(""))
      dd=casacore::Quantity(25.0,"m");
    else
      dd=casaQuantity(dishdiam);
    casacore::Quantity bd;
    if(toCasaString(blockagediam)==casacore::String(""))
      bd=casacore::Quantity(2.5,"m");
    else
      bd=casaQuantity(blockagediam);
    casacore::Quantity mr;
    if(toCasaString(maxrad)==casacore::String(""))
      mr=casacore::Quantity(0.8,"deg");
    else
      mr=casaQuantity(maxrad);
    casacore::Record rec;
    VPManager::Instance()->setpbairy(String(telescope), String(othertelescope), dopb, dd, bd, 
		      mr, rf, sdir, srf, dosquint, pai, usesymmetricbeam, rec);
    r=fromRecord(rec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return r;
}

::casac::record*
vpmanager::setpbcospoly(const std::string& telescope, const std::string& othertelescope, const bool dopb, const std::vector<double>& coeff, const std::vector<double>& scale, const ::casac::variant& maxrad, const ::casac::variant& reffreq, const std::string& isthispb, const ::casac::variant& squintdir, const ::casac::variant& squintreffreq, const bool dosquint, const ::casac::variant& paincrement, const bool usesymmetricbeam)
{
  ::casac::record* r=0;
  try{
    casacore::Quantity mr;
    if(toCasaString(maxrad)==casacore::String(""))
      mr=casacore::Quantity(0.8,"deg");
    else
      mr=casaQuantity(maxrad);
    casacore::Quantity rf;
    if(toCasaString(reffreq)==casacore::String(""))
      rf=casacore::Quantity(1.0,"GHz");
    else
      rf=casaQuantity(reffreq);
    casacore::MDirection sdir;
    if(toCasaString(squintdir)==casacore::String(""))
      sdir=casacore::MDirection(casacore::Quantity(0,"deg"),casacore::Quantity(0,"deg")) ;
    else
      casaMDirection(squintdir, sdir);
    casacore::Quantity srf;
    if(toCasaString(squintreffreq)==casacore::String(""))
      srf=casacore::Quantity(1.0,"GHz");
    else
      srf=casaQuantity(squintreffreq);
    casacore::Quantity pai;
    if(toCasaString(paincrement)==casacore::String(""))
      pai=casacore::Quantity(720,"deg");
    else
      pai=casaQuantity(paincrement);
    casacore::Record rec;
    VPManager::Instance()->setpbcospoly(String(telescope), String(othertelescope), dopb,
			 coeff, scale, mr, rf, isthispb, sdir, srf,
			 dosquint, pai, usesymmetricbeam, rec);
    r=fromRecord(rec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return r;
}

::casac::record*
vpmanager::setpbgauss(const std::string& telescope, const std::string& othertelescope, const bool dopb, const ::casac::variant& halfwidth, const ::casac::variant& maxrad, const ::casac::variant& reffreq, const std::string& isthispb, const ::casac::variant& squintdir, const ::casac::variant& squintreffreq, const bool dosquint, const ::casac::variant& paincrement, const bool usesymmetricbeam)
{

  ::casac::record* r=0;

  try{
    casacore::Quantity rf;
    if(toCasaString(reffreq)==casacore::String(""))
      rf=casacore::Quantity(1.0,"GHz");
    else
      rf=casaQuantity(reffreq);
    casacore::Quantity srf;
    if(toCasaString(squintreffreq)==casacore::String(""))
      srf=casacore::Quantity(1.0,"GHz");
    else
      srf=casaQuantity(squintreffreq);
    casacore::Quantity pai;
    if(toCasaString(paincrement)==casacore::String(""))
      pai=casacore::Quantity(720,"deg");
    else
      pai=casaQuantity(paincrement);
    casacore::MDirection sdir;
    if(toCasaString(squintdir)==casacore::String(""))
      sdir=casacore::MDirection(casacore::Quantity(0,"deg"),casacore::Quantity(0,"deg")) ;
    else
      casaMDirection(squintdir, sdir);
    casacore::Quantity hw;
    if(toCasaString(halfwidth)==casacore::String(""))
      hw=casacore::Quantity(0.5,"deg");
    else
      hw=casaQuantity(halfwidth);
    casacore::Quantity mr;
    if(toCasaString(maxrad)==casacore::String(""))
      mr=casacore::Quantity(0.8,"deg");
    else
      mr=casaQuantity(maxrad);
    casacore::Record rec;
    VPManager::Instance()->setpbgauss(String(telescope), String(othertelescope), dopb, hw, 
		      mr, rf, String(isthispb), sdir, srf, dosquint, pai, 
		      usesymmetricbeam, rec);
    r=fromRecord(rec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return r;
}

::casac::record*
vpmanager::setpbinvpoly(const std::string& telescope, const std::string& othertelescope, const bool dopb, const std::vector<double>& coeff, const ::casac::variant& maxrad, const ::casac::variant& reffreq, const std::string& isthispb, const ::casac::variant& squintdir, const ::casac::variant& squintreffreq, const bool dosquint, const ::casac::variant& paincrement, const bool usesymmetricbeam)
{
  ::casac::record* r=0;
  try{
    casacore::Quantity mr;
    if(toCasaString(maxrad)==casacore::String(""))
      mr=casacore::Quantity(0.8,"deg");
    else
      mr=casaQuantity(maxrad);
    casacore::Quantity rf;
    if(toCasaString(reffreq)==casacore::String(""))
      rf=casacore::Quantity(1.0,"GHz");
    else
      rf=casaQuantity(reffreq);
    casacore::MDirection sdir;
    if(toCasaString(squintdir)==casacore::String(""))
      sdir=casacore::MDirection(casacore::Quantity(0,"deg"),casacore::Quantity(0,"deg")) ;
    else
      casaMDirection(squintdir, sdir);
    casacore::Quantity srf;
    if(toCasaString(squintreffreq)==casacore::String(""))
      srf=casacore::Quantity(1.0,"GHz");
    else
      srf=casaQuantity(squintreffreq);
    casacore::Quantity pai;
    if(toCasaString(paincrement)==casacore::String(""))
      pai=casacore::Quantity(720,"deg");
    else
      pai=casaQuantity(paincrement);
    casacore::Record rec;
    VPManager::Instance()->setpbinvpoly(String(telescope), String(othertelescope), dopb,
			 coeff, mr, rf, isthispb, sdir, srf,
			 dosquint, pai, usesymmetricbeam, rec);
    r=fromRecord(rec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return r;
}

::casac::record*
vpmanager::setpbnumeric(const std::string& telescope, const std::string& othertelescope, const bool dopb, const std::vector<double>& vect, const ::casac::variant& maxrad, const ::casac::variant& reffreq, const std::string& isthispb, const ::casac::variant& squintdir, const ::casac::variant& squintreffreq, const bool dosquint, const ::casac::variant& paincrement, const bool usesymmetricbeam)
{
  ::casac::record* r=0;
  try{
    casacore::Quantity mr;
    if(toCasaString(maxrad)==casacore::String(""))
      mr=casacore::Quantity(0.8,"deg");
    else
      mr=casaQuantity(maxrad);
    casacore::Quantity rf;
    if(toCasaString(reffreq)==casacore::String(""))
      rf=casacore::Quantity(1.0,"GHz");
    else
      rf=casaQuantity(reffreq);
    casacore::MDirection sdir;
    if(toCasaString(squintdir)==casacore::String(""))
      sdir=casacore::MDirection(casacore::Quantity(0,"deg"),casacore::Quantity(0,"deg")) ;
    else
      casaMDirection(squintdir, sdir);
    casacore::Quantity srf;
    if(toCasaString(squintreffreq)==casacore::String(""))
      srf=casacore::Quantity(1.0,"GHz");
    else
      srf=casaQuantity(squintreffreq);
    casacore::Quantity pai;
    if(toCasaString(paincrement)==casacore::String(""))
      pai=casacore::Quantity(720,"deg");
    else
      pai=casaQuantity(paincrement);
    casacore::Record rec;
    VPManager::Instance()->setpbnumeric(String(telescope), String(othertelescope), dopb,
			 vect, mr, rf, isthispb, sdir, srf,
			 dosquint, pai, usesymmetricbeam, rec);
    r=fromRecord(rec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return r;
}

::casac::record*
vpmanager::setpbimage(const std::string& telescope, const std::string& othertelescope, const bool dopb, const std::string& realimage, const std::string& imagimage, const std::string& compleximage, const std::vector<std::string>& antnames)
{

  ::casac::record* r=0;
  try{
    casacore::Record rec;
    Vector <String> anames(toVectorString(antnames));
    ////CAS-5666 go round
    if(anames.nelements()==1 && anames[0]==casacore::String(""))
      anames[0]=casacore::String("*");
    VPManager::Instance()->setpbimage(String(telescope), String(othertelescope), dopb, 
				      String(realimage), String(imagimage), String(compleximage), anames, rec);
    r=fromRecord(rec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return r;
}

::casac::record*
vpmanager::setpbpoly(const std::string& telescope, const std::string& othertelescope, const bool dopb, const std::vector<double>& coeff, const ::casac::variant& maxrad, const ::casac::variant& reffreq, const std::string& isthispb, const ::casac::variant& squintdir, const ::casac::variant& squintreffreq, const bool dosquint, const ::casac::variant& paincrement, const bool usesymmetricbeam)
{
  ::casac::record* r=0;
  try{
    casacore::Quantity mr;
    if(toCasaString(maxrad)==casacore::String(""))
      mr=casacore::Quantity(0.8,"deg");
    else
      mr=casaQuantity(maxrad);
    casacore::Quantity rf;
    if(toCasaString(reffreq)==casacore::String(""))
      rf=casacore::Quantity(1.0,"GHz");
    else
      rf=casaQuantity(reffreq);
    casacore::MDirection sdir;
    if(toCasaString(squintdir)==casacore::String(""))
      sdir=casacore::MDirection(casacore::Quantity(0,"deg"),casacore::Quantity(0,"deg")) ;
    else
      casaMDirection(squintdir, sdir);
    casacore::Quantity srf;
    if(toCasaString(squintreffreq)==casacore::String(""))
      srf=casacore::Quantity(1.0,"GHz");
    else
      srf=casaQuantity(squintreffreq);
    casacore::Quantity pai;
    if(toCasaString(paincrement)==casacore::String(""))
      pai=casacore::Quantity(720,"deg");
    else
      pai=casaQuantity(paincrement);
    casacore::Record rec;
    VPManager::Instance()->setpbpoly(String(telescope), String(othertelescope), dopb,
		      coeff, mr, rf, isthispb, sdir, srf,
		      dosquint, pai, usesymmetricbeam, rec);
    r=fromRecord(rec);

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return r;
}

bool
vpmanager::setpbantresptable(const std::string& telescope, 
			     const std::string& othertelescope, 
			     const bool dopb, 
			     const std::string& antresppath)
{
  *itsLog << LogOrigin("vp", "setantresptable");

  Bool rval = VPManager::Instance()->setpbantresptable(telescope, othertelescope, dopb, antresppath);

  return rval;

}


bool vpmanager::reset(){

  *itsLog << LogOrigin("vp", "reset");

  VPManager::Instance()->reset();

  return true;
}

bool vpmanager::setuserdefault(const int vplistnum,
			       const std::string& telescope,
			       const std::string& anttype)
{

  *itsLog << LogOrigin("vp", "setuserdefault");

  bool rval = VPManager::Instance()->setuserdefault(vplistnum, telescope, anttype);

  return rval;

}

int vpmanager::getuserdefault(const std::string& telescope,
			      const std::string& anttype)
{

  int rval = -2;

  *itsLog << LogOrigin("vp", "getuserdefault");

  if(!VPManager::Instance()->getuserdefault(rval, telescope, anttype)){
    *itsLog << LogIO::WARN << "No default response for telescope \""
	    << telescope << "\", antenna type \"" << anttype << "\"" 
	    << LogIO::POST;
    rval = -2;
  }

  return rval;

}


std::vector<std::string> vpmanager::getanttypes(const std::string& telescope,
						const casac::variant& obstime, 
						const casac::variant& freq, 
						const casac::variant& obsdirection)
{
  std::vector<string> rval(0);
  Vector<String> antTypes;

  *itsLog << LogOrigin("vp", "getanttypes");

  try{    
	  casacore::MEpoch mObsTime;
	  casacore::MFrequency mFreq;
	  casacore::MDirection mObsDir;
    
    if(!casaMEpoch(obstime, mObsTime)){
      *itsLog << LogIO::SEVERE << "Could not interprete obstime parameter "
	      << toCasaString(obstime) << LogIO::POST;
      return rval;
    }
    if(!casaMFrequency(freq, mFreq)){
      *itsLog << LogIO::SEVERE << "Could not interprete freq parameter "
	      << toCasaString(freq) << LogIO::POST;
      return rval;
    }
    if(!casaMDirection(obsdirection, mObsDir)){
      if(toCasaString(obsdirection).empty()){
	casaMDirection("AZEL 0deg 0deg", mObsDir);
      }
      else{
	*itsLog << LogIO::SEVERE << "Could not interprete obsdirection parameter "
		<< toCasaString(obsdirection) << LogIO::POST;
	return rval;
      }
    }

    if(!VPManager::Instance()->getanttypes(antTypes, telescope, mObsTime, mFreq, mObsDir)){
      *itsLog << LogIO::SEVERE << "Error determining antenna types for telescope "
	      << telescope << LogIO::POST;
      return rval;
    }      

  } catch(AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  
  for(uInt i=0; i<antTypes.size(); i++){
    rval.push_back(antTypes(i));
  }

  return rval;
  
}
  

int vpmanager::numvps(const std::string& telescope,
		      const casac::variant& obstime, 
		      const casac::variant& freq, 
		      const casac::variant& obsdirection){

  int rval(-1);

  *itsLog << LogOrigin("vp", "numvps");

  try{

	  casacore::MEpoch mObsTime;
	  casacore::MFrequency mFreq;
	  casacore::MDirection mObsDir;
    
    if(!casaMEpoch(obstime, mObsTime)){
      *itsLog << LogIO::SEVERE << "Could not interprete obstime parameter "
	      << toCasaString(obstime) << LogIO::POST;
      return rval;
    }
    if(!casaMFrequency(freq, mFreq)){
      *itsLog << LogIO::SEVERE << "Could not interprete freq parameter "
	      << toCasaString(freq) << LogIO::POST;
      return rval;
    }
    if(!casaMDirection(obsdirection, mObsDir)){
      if(toCasaString(obsdirection).empty()){
	casaMDirection("AZEL 0deg 0deg", mObsDir);
      }
      else{
	*itsLog << LogIO::SEVERE << "Could not interprete obsdirection parameter "
		<< toCasaString(obsdirection) << LogIO::POST;
	return rval;
      }
    }

    rval = VPManager::Instance()->numvps(telescope, mObsTime, mFreq, mObsDir);

  } catch(AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return rval;

}


::casac::record* 
vpmanager::getvp(const std::string& telescope,
		 const std::string& antennatype, 
		 const casac::variant& obstime, 
		 const casac::variant& freq, 
		 const casac::variant& obsdirection){

  ::casac::record* r=0;

  *itsLog << LogOrigin("vp", "getvp");

  try{
    casacore::MEpoch mObsTime;
    casacore::MFrequency mFreq;
    casacore::MDirection mObsDir;
    Record rec;
    
    int nRefs = 0;

    if(casaMEpoch(obstime, mObsTime)){
      nRefs +=1;
    }
    else{
      if(toCasaString(obstime).empty()){
	casaMEpoch("2000/01/01T00:00:00", mObsTime);
      }
      else{
	*itsLog << LogIO::SEVERE << "Could not interprete obstime parameter "
		<< toCasaString(obstime) << LogIO::POST;
	return r;
      }
    }
    if(casaMFrequency(freq, mFreq)){
      nRefs +=1;
    }
    else{
      if(toCasaString(freq).empty()){
	casaMFrequency("TOPO 0Hz", mFreq);
      }
      else{
	*itsLog << LogIO::SEVERE << "Could not interprete freq parameter "
		<< toCasaString(freq) << LogIO::POST;
	return r;
      }
    }

    if(!casaMDirection(obsdirection, mObsDir)){
      if(toCasaString(obsdirection).empty()){
	casaMDirection("AZEL 0deg 90deg", mObsDir);
      }
      else{
	*itsLog << LogIO::SEVERE << "Could not interprete obsdirection parameter "
		<< toCasaString(obsdirection) << LogIO::POST;
	return r;
      }
    }

    if(nRefs==0){
      if(VPManager::Instance()->getvp(rec,telescope, antennatype)){
	r = fromRecord(rec);
      }
    }
    else{
      if(VPManager::Instance()->getvp(rec,telescope, mObsTime, mFreq, antennatype, mObsDir)){
	r = fromRecord(rec);
      }
    }

  } catch(AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return r;

}

::casac::record* 
vpmanager::getvps(const std::string& telescope,
		  const std::vector<std::string>& antennas, 
		  const casac::variant& obstimestart, 
		  const casac::variant& obstimeend, 
		  const casac::variant& minfreq, 
		  const casac::variant& maxfreq, 
		  const casac::variant& obsdirection){

  ::casac::record* r=0;

  *itsLog << LogOrigin("vp", "getvp");

  try{
    casacore::Vector<casacore::MEpoch> mObsTimeV;
    casacore::Vector<casacore::MFrequency> mFreqV;
    casacore::MDirection mObsDir;
    Record rec;
    
    int nRefs = 0;


    if(antennas.size()==0 || (antennas[0]=="" && antennas.size()==1)){
      *itsLog << LogIO::SEVERE << "Antenna name list is empty." << LogIO::POST;
      return r;
    }

    std::vector<casac::variant> obst;
    obst.push_back(obstimestart);
    obst.push_back(obstimeend);
    mObsTimeV.resize(2);

    for(uInt i=0; i<2; i++){
      if(casaMEpoch(obst[i], mObsTimeV[i])){
	nRefs +=1;
      }
      else{
	if(toCasaString(obst[i]).empty()){
	  casaMEpoch("2000/01/01T00:00:00", mObsTimeV[i]);
	}
	else{
	  *itsLog << LogIO::SEVERE << "Could not interprete obstime parameter "
		  << toCasaString(obst[i]) << LogIO::POST;
	  return r;
	}
      }
    }

    std::vector<casac::variant> obsf;
    obsf.push_back(minfreq);
    obsf.push_back(maxfreq);
    mFreqV.resize(2);

    for(uInt i=0; i<2; i++){
      if(casaMFrequency(obsf[i], mFreqV[i])){
	nRefs +=1;
      }
      else{
	if(toCasaString(obsf[i]).empty()){
	  casaMFrequency("TOPO 0Hz", mFreqV[i]);
	}
	else{
	  *itsLog << LogIO::SEVERE << "Could not interprete freq parameter "
		  << toCasaString(obsf[i]) << LogIO::POST;
	  return r;
	}
      }
    }

    if(!casaMDirection(obsdirection, mObsDir)){
      if(toCasaString(obsdirection).empty()){
	casaMDirection("AZEL 0deg 90deg", mObsDir);
      }
      else{
	*itsLog << LogIO::SEVERE << "Could not interprete obsdirection parameter "
		<< toCasaString(obsdirection) << LogIO::POST;
	return r;
      }
    }


    if(nRefs==0){
	*itsLog << LogIO::SEVERE << "None of the necessary input parameters was provided."
	        << LogIO::POST;
	return r;
    }
    else{
      casacore::Vector<Record> out_rec_list;
      casacore::Vector<casacore::Vector<uInt> > beam_index;

      if(VPManager::Instance()->getvps(out_rec_list, beam_index,
				       telescope, mObsTimeV, mFreqV, antennas, mObsDir)){
	if(beam_index.size()>0){
          // translate index into matrix
	  Matrix<uInt> beamIndex(beam_index.size(), beam_index[0].size(), -1);
	  for(uInt i=0; i<beam_index.size(); i++){
	    for(uInt j=0; j<beam_index[i].size(); j++){
	      beamIndex(i,j) = beam_index[i][j];
	    }
	  }
	  rec.define("beamindex", beamIndex); 

	  // store the records in the output record
	  for(uInt i=0; i<out_rec_list.size(); i++){
	    String bnum = String::toString(i);
	    bnum.trim();
	    rec.defineRecord(String("uniquebeam_")+bnum, out_rec_list[i]);
	  }

	  r = fromRecord(rec);
	}
	else{
	  *itsLog << LogIO::SEVERE << "No appropriate beams found." << LogIO::POST;
	}
      }
    }

  } catch(AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }

  return r;

}


bool 
vpmanager::createantresp(const std::string& imdir, 
			 const std::string& starttime, 
			 const std::vector<std::string>& bandnames, 
			 const std::vector<std::string>& bandminfreq, 
			 const std::vector<std::string>& bandmaxfreq)
{
  bool rstat(false);

  *itsLog << LogOrigin("vp", "createantresp");


  try{

//     cout << " imdir " << imdir << " starttime " << starttime << endl;
//     for(uInt i=0; i<bandnames.size(); i++){
//       cout << i << " bandnames " << bandnames[i] << " bandminfreq " << bandminfreq[i] << " bandmaxfreq " << bandmaxfreq[i] << endl;
//     }

    String imDir = toCasaString(imdir);
    Path imPath(imDir);
    String absPathName = imPath.absoluteName();
    Path absImPath(absPathName);
    if(!absImPath.isValid()){
      *itsLog << LogIO::SEVERE << "Invalid path: \"" << imDir  << "\"." << LogIO::POST;
      return rstat;
    }
      
    String startTime = toCasaString(starttime);
    casacore::MEpoch theStartTime;
    Quantum<Double> qt;
    if (MVTime::read(qt,startTime)) {
      MVEpoch mv(qt);
      theStartTime = casacore::MEpoch(mv, casacore::MEpoch::UTC);
    } else {
      *itsLog << LogIO::SEVERE << "Invalid time format: " 
	      << startTime << LogIO::POST;
      return rstat;
    }

    uInt nBands = bandnames.size();

    if((nBands != bandminfreq.size()) || (nBands != bandmaxfreq.size())){
      *itsLog << LogIO::SEVERE << "bandnames, bandminfreq, and bandmaxfreq need to have the same number of elements." 
	      << LogIO::POST;
      return rstat;
    }
      
    Vector<MVFrequency> bandMinFreqV(nBands), bandMaxFreqV(nBands);

    for(uInt i=0; i<nBands; i++){
      Quantum<Double> freqmin, freqmax;
      if(!Quantum<Double>::read(freqmin, toCasaString(bandminfreq[i]))){
	*itsLog << LogIO::SEVERE << "Invalid quantity in bandminfreq " << i << ": " << bandminfreq[i] 
		<< LogIO::POST;
	return rstat;
      }
      if(!Quantum<Double>::read(freqmax, toCasaString(bandmaxfreq[i]))){
	*itsLog << LogIO::SEVERE << "Invalid quantity in bandmaxfreq " << i << ": " << bandmaxfreq[i] 
		<< LogIO::POST;
	return rstat;
      }

      try{
	MVFrequency fMin(freqmin);
	bandMinFreqV(i) = fMin;
      }
      catch(AipsError x) {
	*itsLog << LogIO::SEVERE << "Error interpreting bandminfreq " << i << " as frequency: " << bandminfreq[i] << endl << x.getMesg() << LogIO::POST;
	return rstat;
      }
      try{
	MVFrequency fMax(freqmax);
	bandMaxFreqV(i) = fMax;
      }
      catch(AipsError x) {
	*itsLog << LogIO::SEVERE << "Error interpreting bandmaxfreq " << i << " as frequency: " << bandmaxfreq[i] << endl << x.getMesg() << LogIO::POST;
	return rstat;
      }
    }

    // all input parameters checked, now find all the images
    Directory imDirD(imDir);
    String imPattern = "*?_*?_*_*_*?_*?_*?_*?_*?_*?_*?_*?_*?_*?_*_*?.im";
    Vector<String> imNamesV = imDirD.find(Regex::fromPattern(imPattern), true, false); // follow symlinks, non-recursive
    
    uInt nIm = imNamesV.size();

    *itsLog << LogIO::NORMAL << "Found " << nIm << " response images in directory " << absPathName << LogIO::POST;

    {
      imPattern = "*.im";
      Vector<String> imNamesV2 = imDirD.find(Regex::fromPattern(imPattern), true, false);

      if(nIm==0){
	if(imNamesV2.nelements()==0){
	  *itsLog << LogIO::SEVERE << "No images found in directory " << imDir << LogIO::POST;
	}
	else{ // no images found with the correctly formatted name
	  *itsLog << LogIO::SEVERE << "No images with correctly formatted file name found in directory " << imDir << LogIO::POST;
	  *itsLog << LogIO::SEVERE << "Expected image file name format is:" << endl
		  << " obsname_beamnum_anttype_rectype_azmin_aznom_azmax_elmin_elnom_elmax_freqmin_freqnom_freqmax_frequnit_comment_functype.im" 
		  << LogIO::POST;
	}
	return rstat;
      }
      else if(imNamesV2.nelements()!=nIm){
	*itsLog << LogIO::WARN << "Not all the images found in directory " << imDir 
		<< endl << " have correctly formatted filenames. Continuing anyway ..." << LogIO::POST;
      }
    }

    // sort images by name

    GenSort<String>::sort(imNamesV);

    // disect the image file names

    Vector<String> obsnameV(nIm);
    Vector<uInt> beamnumV(nIm);
    Vector<String> anttypeV(nIm);
    Vector<String> rectypeV(nIm);
    Vector<Double> azminV(nIm),aznomV(nIm),azmaxV(nIm),elminV(nIm),elnomV(nIm),elmaxV(nIm),
      freqminV(nIm),freqnomV(nIm),freqmaxV(nIm);
    Vector<String> frequnitV(nIm);
    Vector<String> functypeV(nIm);

    for(uInt i=0; i<nIm; i++){
      
      String::size_type pos = 0;
      String::size_type pos2;
      
      pos = imNamesV(i).find("_", pos);
      obsnameV(i) = imNamesV(i).substr(0, pos); // use pos as length in this case

      casacore::MPosition Xpos;
      if (obsnameV(i).length() == 0 || 
	  !MeasTable::Observatory(Xpos,obsnameV(i))) {
	// unknown observatory
	*itsLog << LogIO::SEVERE << "Unknown observatory: \"" << obsnameV(i) << "\"." << LogIO::POST;
	return rstat;
      }

      pos2 = imNamesV(i).find("_", ++pos); // start one after the "_"
      beamnumV(i) = atoi(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      anttypeV(i) = imNamesV(i).substr(pos, pos2-pos);

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      rectypeV(i) = imNamesV(i).substr(pos, pos2-pos);

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      azminV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      aznomV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      azmaxV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      elminV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      elnomV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      elmaxV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      freqminV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      freqnomV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      freqmaxV(i) = atof(imNamesV(i).substr(pos, pos2-pos).c_str());

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); 
      frequnitV(i) = imNamesV(i).substr(pos, pos2-pos);

      pos = pos2;
      pos2 = imNamesV(i).find("_", ++pos); // skip the comment 
      pos = pos2;
      pos2 = imNamesV(i).find(".", ++pos); 
      functypeV(i) = imNamesV(i).substr(pos, pos2-pos);

//       cout << obsnameV(i) << ", " << beamnumV(i) << ", " << anttypeV(i) << ", "  << rectypeV(i) 
// 	   << ", AZ " << azminV(i) << ", " << aznomV(i) << ", " << azmaxV(i) 
// 	   << ", EL " << elminV(i) << ", " << elnomV(i) << ", " << elmaxV(i) 
// 	   << ", freq " << freqminV(i) << ", " << freqnomV(i) << ", " << freqmaxV(i) << ", " 
// 	   << frequnitV(i) << ", " << functypeV(i) << endl;
     
    }

    AntennaResponses aR(""); // empty table in memory
    uInt iRow = 0;
    Int currBeamId = 0;

    for(uInt iIm=0; iIm<nIm; iIm++){

      String currObsName = obsnameV(iIm);
      uInt currBeamnum = beamnumV(iIm);
      String currAntType = anttypeV(iIm);
      String currRecType = rectypeV(iIm);
      Double currAzMin = azminV(iIm);
      Double currAzNom = aznomV(iIm);
      Double currAzMax = azmaxV(iIm);
      Double currElMin = elminV(iIm);
      Double currElNom = elnomV(iIm);
      Double currElMax = elmaxV(iIm);

      if(iIm>0 && currObsName==obsnameV(iIm)){ // beam id should be unique for a given obsname
	currBeamId++;
      }
      else{
	currBeamId = 0;
      }

      // accumulate vectors

      vector<String> bNameAV;
      vector<MVFrequency> minFreqAV;
      vector<MVFrequency> nomFreqAV;
      vector<MVFrequency> maxFreqAV;
      vector<AntennaResponses::FuncTypes> fTypAV;
      vector<String> funcNameAV;
      vector<uInt> funcChannelAV;
      
      uInt ii=iIm;
      while(ii<nIm &&
	    currObsName == obsnameV(ii) &&
	    currBeamnum == beamnumV(ii) &&
	    currAntType == anttypeV(ii) &&
	    currRecType == rectypeV(ii) &&
	    currAzMin == azminV(ii) &&
	    currAzNom == aznomV(ii) &&
	    currAzMax == azmaxV(ii) &&
	    currElMin == elminV(ii) &&
	    currElNom == elnomV(ii) &&
	    currElMax == elmaxV(ii)){
	
	// determine band name
	MVFrequency currFreqMin(Quantum<Double>(freqminV(ii), Unit(frequnitV(ii))));
	MVFrequency currFreqNom(Quantum<Double>(freqnomV(ii), Unit(frequnitV(ii))));
	MVFrequency currFreqMax(Quantum<Double>(freqmaxV(ii), Unit(frequnitV(ii))));

	Int iFound = -1;
	for(uInt i=0; i<nBands; i++){
	  if(bandMinFreqV(i).get().getBaseValue() <= currFreqMin.get().getBaseValue() &&
	     currFreqMin.get().getBaseValue() <= bandMaxFreqV(i).get().getBaseValue() &&
	     bandMinFreqV(i).get().getBaseValue() <= currFreqNom.get().getBaseValue() &&
	     currFreqNom.get().getBaseValue() <= bandMaxFreqV(i).get().getBaseValue() &&
	     bandMinFreqV(i).get().getBaseValue() <= currFreqMax.get().getBaseValue() &&
	     currFreqMax.get().getBaseValue() <= bandMaxFreqV(i).get().getBaseValue()){
	    // found the band
	    iFound = i;
	    break;
	  }
	}
	if(iFound == -1){
	  *itsLog << LogIO::SEVERE << "Image " << imNamesV(ii) 
		  << " has band definition inconsistent with the available bands." << LogIO::POST;
	  return rstat;
	}
	String currBandName = bandnames[iFound];

	bNameAV.push_back(currBandName);
	minFreqAV.push_back(currFreqMin);
	nomFreqAV.push_back(currFreqNom);
	maxFreqAV.push_back(currFreqMax);
	fTypAV.push_back(AntennaResponses::FuncType(functypeV(ii)));
	if(AntennaResponses::FuncType(functypeV(ii))==AntennaResponses::INVALID){
	  *itsLog << LogIO::WARN << "Function type " << functypeV(ii) << " of image " << imNamesV(ii) 
		  << " is not recognised. Continuing anyway ..." << LogIO::POST;
	}			 
	funcNameAV.push_back(imNamesV(ii));
	funcChannelAV.push_back(0); // multi-channel images not yet supported

	ii++;
      } // end while

      iIm = ii-1; // update global image counter

      // sort accumulated vectors by nominal frequency
      Vector<uInt> sortIndex;
      uInt nSubBands = bNameAV.size();
      Vector<String> sortedBName(nSubBands);
      Vector<MVFrequency> sortedMinFreq(nSubBands);
      Vector<MVFrequency> sortedMaxFreq(nSubBands); 
      Vector<AntennaResponses::FuncTypes> sortedFTyp(nSubBands);
      Vector<String> sortedFuncName(nSubBands); 
      Vector<uInt> sortedFuncChannel(nSubBands);

      Vector<MVAngle> rotAngOffset(nSubBands,MVAngle(0.));

      Vector<MVFrequency> sortedNomFreq(nomFreqAV); // convert this std::vector to a Vector at the same time

      GenSortIndirect<MVFrequency>::sort(sortIndex, sortedNomFreq);

      for(uInt i=0; i<nSubBands; i++){
	sortedBName(i) = bNameAV[sortIndex(i)];
	sortedMinFreq(i) = minFreqAV[sortIndex(i)];
	sortedMaxFreq(i) = maxFreqAV[sortIndex(i)]; 
	sortedFTyp(i) = fTypAV[sortIndex(i)];
	sortedFuncName(i) = funcNameAV[sortIndex(i)]; 
	sortedFuncChannel(i) = funcChannelAV[sortIndex(i)];
	sortedNomFreq(i) = nomFreqAV[sortIndex(i)];
      }

      // fill table
      if(!aR.putRow(iRow, currObsName, currBeamId,
		    sortedBName, sortedMinFreq, sortedMaxFreq, 
		    sortedFTyp, sortedFuncName, 
		    sortedFuncChannel, sortedNomFreq,
		    rotAngOffset, // all zero, no sorting necessary
		    currAntType, theStartTime,
		    casacore::MDirection(casacore::Quantity(currAzNom, "deg"),
			       casacore::Quantity(currElNom, "deg"),
			       casacore::MDirection::AZEL),
		    casacore::MDirection(casacore::Quantity(currAzMin, "deg"),
			       casacore::Quantity(currElMin, "deg"),
			       casacore::MDirection::AZEL),
		    casacore::MDirection(casacore::Quantity(currAzMax, "deg"),
			       casacore::Quantity(currElMax, "deg"),
			       casacore::MDirection::AZEL),
		    currRecType,
		    currBeamnum
		    )
	 ){
	*itsLog << LogIO::SEVERE << "Error entering row " << iRow << " into antenna responses table." 
		<< LogIO::POST;
	return rstat;
      } // end if

      iRow++;

    } // end for

    // write table to disk

    aR.create(absPathName+"/AntennaResponses");

    *itsLog << LogIO::NORMAL << "Created antenna responses table " << absPathName+"/AntennaResponses" 
	    << endl << " with " << iRow << " rows." << LogIO::POST;
    
    rstat = true;

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

std::string 
vpmanager::getrespimagename(const std::string& telescope, 
			    const std::string& starttime, 
			    const std::string& frequency, 
			    const std::string& functype, 
			    const std::string& anttype, 
			    const std::string& azimuth,
			    const std::string& elevation, 
			    const std::string& rectype,
			    const int beamnumber)
{
  std::string rval("");

  *itsLog << LogOrigin("vp", "getrespimagename");
    
  try{

    String obsName(telescope);
    String antRespPath;
    if (!MeasTable::AntennaResponsesPath(antRespPath, obsName)) {
      // unknown observatory
      *itsLog << LogIO::SEVERE << "No antenna responses path set for observatory: \"" << obsName 
	      << "\" in the AntennaResponses column of the Observatories table." << LogIO::POST;
      return rval;
    }

    AntennaResponses aR(antRespPath);

    String functionImageName;
    uInt funcChannel;
    MFrequency nomFreq;
    AntennaResponses::FuncTypes fType;
    MVAngle rotAngOffset;

    String startTime = toCasaString(starttime);
    casacore::MEpoch theStartTime;
    Quantum<Double> qt;
    if (MVTime::read(qt,startTime)) {
      MVEpoch mv(qt);
      theStartTime = casacore::MEpoch(mv, casacore::MEpoch::UTC);
    } else {
      *itsLog << LogIO::SEVERE << "Invalid time format: " 
	      << startTime << LogIO::POST;
      return rval;
    }

    Quantum<Double> freqQ;
    if(!Quantum<Double>::read(freqQ, toCasaString(frequency))){
      *itsLog << LogIO::SEVERE << "Invalid quantity in parameter frequency " << frequency 
	      << LogIO::POST;
      return rval;
    }
    MFrequency freqM(freqQ, MFrequency::TOPO);
    AntennaResponses::FuncTypes requFType = AntennaResponses::FuncType(toCasaString(functype));
    String antennaType = String(anttype);
    Quantum<Double> az, el;
    if(!Quantum<Double>::read(az, toCasaString(azimuth))){
      *itsLog << LogIO::SEVERE << "Invalid quantity in parameter azimuth " << azimuth 
	      << LogIO::POST;
      return rval;
    }
    if(!Quantum<Double>::read(el, toCasaString(elevation))){
      *itsLog << LogIO::SEVERE << "Invalid quantity in parameter elevation " << elevation 
	      << LogIO::POST;
      return rval;
    }

    casacore::MDirection center = casacore::MDirection(az, el, casacore::MDirection::AZEL);
    String receiverType = toCasaString(rectype);

    if(!aR.getImageName(functionImageName, 
			funcChannel, 
			nomFreq, 
			fType,
			rotAngOffset,
			obsName,
			theStartTime,
			freqM,
			requFType,
			antennaType,
			center,
			receiverType,
			beamnumber)){
      *itsLog << LogIO::SEVERE << "Could not find appropriate response image." << LogIO::POST;
      return rval;
    }
      
    *itsLog << LogIO::NORMAL << "Found image " << functionImageName << endl
	    << " of type " << (Int)fType << ". Nominal frequency is " << nomFreq.getValue() 
	    << " Hz. Use channel " << funcChannel << "." << LogIO::POST;
    rval = functionImageName;

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
    RETHROW(x);
  }
  return rval;
}


} // casac namespace

