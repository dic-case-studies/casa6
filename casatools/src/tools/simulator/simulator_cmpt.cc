
/***
 * Framework independent implementation file for simulator...
 *
 * Implement the simulator component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 *
 * @author
 * @version 
 ***/

#include <iostream>
#include <simulator_cmpt.h>
#include <casa/Exceptions/Error.h>
#include <casa/Containers/Record.h>
#include <synthesis/MeasurementEquations/Simulator.h>
#include<casa/BasicSL/String.h>
#include<casa/Utilities/Assert.h>
#include<measures/Measures/MDirection.h>
#include<measures/Measures/MPosition.h>
#include<measures/Measures/MEpoch.h>
#include <measures/Measures/MeasureHolder.h>
#include<casa/Quanta/QuantumHolder.h>
#include<ms/MeasurementSets.h>
#include <casa/Logging/LogIO.h>
#include <imageanalysis/ImageAnalysis/ImageFactory.h>

using namespace std;
using namespace casacore;
using namespace casa;

using namespace casacore;
namespace casac {

simulator::simulator()
{

  itsSim=0;
  itsLog = new LogIO();

}

simulator::~simulator()
{

  if(itsSim !=0)
    delete itsSim;
  
  delete itsLog;

}

bool
simulator::open(const std::string& ms)
{
  Bool rstat(false);
  try {

    // In case already open, close it!
    close();

    String lems(ms);
    itsSim=new Simulator(lems);
    if(itsSim !=0)
      rstat=true;


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

bool
simulator::openfromms(const std::string& ms)
{

  Bool rstat(false);

  try {

    // In case already open, close it!
    close();

    MeasurementSet thems(ms,TableLock(TableLock::AutoLocking), Table::Update);
    itsSim=new Simulator(thems);
    if(itsSim !=0)
      rstat=true;


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

bool
simulator::close()
{
  Bool rstat(false);
  
  try {
  
    if(itsSim !=0)
      delete itsSim;
    itsSim=0;
    rstat=true;


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

  
}

bool
simulator::done()
{

Bool rstat(false);
  
  try {
  
    if(itsSim !=0)
      delete itsSim;
    itsSim=0;
    rstat=true;


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;


}

std::string
simulator::name()
{

  std::string nameString;
try {
  
    if(itsSim !=0)
      nameString=itsSim->name();


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

 return nameString;
    
}

bool
simulator::summary()
{

 Bool rstat(false);
  
  try {
  
    if(itsSim !=0)
      rstat=itsSim->summary();
   
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

std::string
simulator::type()
{

  return std::string("Simulator");
}

bool
simulator::settimes(const ::casac::variant& integrationtime, const bool usehourangle, const ::casac::variant& referencetime)
{
 Bool rstat(false);
  
  try {
  
    if(itsSim !=0){
	    casacore::MEpoch lepoch;
      if (!casaMEpoch(referencetime, lepoch)){
	*itsLog << LogIO::SEVERE 
		<< "Could not convert referencetime to an Epoch Measures"
		<< LogIO::POST;
	
      }
      casacore::Quantity qIntTime(casaQuantity(integrationtime));
      // Negative integration time crashes casa
      AlwaysAssert(qIntTime.getValue("s")>=0, AipsError);
      rstat=itsSim->settimes(qIntTime, usehourangle, lepoch);
    }


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

 
}

bool
simulator::observe(const std::string& sourcename, const std::string& spwname, const ::casac::variant& starttime, const ::casac::variant& stoptime,
		   const bool add_observation,
		   const bool state_sig,
		   const bool state_ref,
		   const double state_cal,
		   const double state_load,
		   const long state_sub_scan,
		   const std::string& state_obs_mode,
		   const std::string& observername,
		   const std::string& projectname)
// defaults are inserted to .h file from xml
//		   const bool add_observation=true,
//		   const bool state_sig=true,
//		   const bool state_ref=true,
//		   const double state_cal=0.,
//		   const double state_load=0.,
//		   const int state_sub_scan=0,
//		   const std::string& state_obs_mode="OBSERVE_TARGET.ON_SOURCE")
{
Bool rstat(false);
  
  try {
  
    if(itsSim !=0){
      casacore::Quantity qstarttime(casaQuantity(starttime));
      casacore::Quantity qstoptime(casaQuantity(stoptime));
      rstat=itsSim->observe(sourcename, spwname, qstarttime, qstoptime,
			    add_observation,state_sig,state_ref,state_cal,state_load,state_sub_scan,state_obs_mode,observername,projectname);
    }


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}


bool
//simulator::observemany(const std::vector<string>& sourcenames, const std::string& spwname, const ::casac::variant& starttimes, const ::casac::variant& stoptimes, const ::casac::variant& directions)
simulator::observemany(const std::vector<string>& sourcenames, const std::string& spwname, const std::vector<string>& starttimes, const std::vector<string>& stoptimes, const std::vector<string>& directions,
		       const bool add_observation,
		       const bool state_sig,
		       const bool state_ref,
		       const double state_cal,
		       const double state_load,
		       const long state_sub_scan,
		       const std::string& state_obs_mode,
		       const std::string& observername,
		       const std::string& projectname)
{
Bool rstat(false);
  
  try {
  
    if(itsSim !=0){
      casacore::String sspwname(spwname);
      Vector<String> ssourcenames(toVectorString(sourcenames));
//      Vector<String> ssourcenames();
//     Vector<String> inFiles;
//     if (infiles.type() == ::casac::variant::BOOLVEC) {
//	inFiles.resize(0);      // unset
//     } else if (infiles.type() == ::casac::variant::STRING) {
//	sepCommaEmptyToVectorStrings(inFiles, infiles.toString());
//     } else if (infiles.type() == ::casac::variant::STRINGVEC) {
//	inFiles = toVectorString(infiles.toStringVec());
//     } else {
//	*itsLog << LogIO::WARN << "Unrecognized infiles datatype"
//		<< LogIO::POST;
//     }
      Vector<String> sstarttimes(toVectorString(starttimes));
      Vector<String> sstoptimes(toVectorString(stoptimes));
      Vector<String> sdirections(toVectorString(directions));

      uInt nptg=sstarttimes.nelements();
      AlwaysAssert(sstoptimes.nelements() == nptg, AipsError);
      if (ssourcenames.nelements() == 1 ) {
	ssourcenames.resize(nptg,true);
      }
      AlwaysAssert(ssourcenames.nelements() == nptg, AipsError);

      if (sdirections.nelements() == 1) {
	sdirections.resize(nptg,"");
      }
      AlwaysAssert(sdirections.nelements() == nptg, AipsError);

      casacore::Vector<casacore::Quantity> qstarttimes(nptg);
      casacore::Vector<casacore::Quantity> qstoptimes(nptg);
      casacore::Vector<casacore::MDirection> mdirections(nptg);
      casacore::MDirection mdir,northPole;
      for (uInt i=0; i<nptg; i++) {
	qstarttimes[i]=casaQuantity(sstarttimes[i]);
	qstoptimes[i]=casaQuantity(sstoptimes[i]);
	if (sdirections[i].length() >0) {
	  if (!casaMDirection(sdirections[i], mdir)){
	    *itsLog << LogIO::SEVERE 
		    << "Could not convert direction to a Direction Measure."
		    << LogIO::POST;
	  }
	} else {
	  mdir=northPole;
	}
	mdirections[i]=mdir;
      }
      rstat=itsSim->observemany(ssourcenames, sspwname, qstarttimes, qstoptimes, mdirections,
				add_observation,state_sig,state_ref,state_cal,state_load,state_sub_scan,state_obs_mode,observername,projectname);
    }


  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;

}

bool
simulator::setlimits(const double shadowlimit, const ::casac::variant& elevationlimit)
{

  Bool rstat(false);
  try {
  
    if(itsSim !=0){
      casacore::Quantity qelevationlimit(casaQuantity(elevationlimit));
      rstat=itsSim->setlimits(shadowlimit, qelevationlimit);
    }

    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;

}

bool
simulator::setauto(const double autocorrwt)
{

  Bool rstat(false);
  try {
  
    if(itsSim !=0){
      rstat=itsSim->setauto(autocorrwt);
    }
    

  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }

  return rstat;
  
}


bool
simulator::setconfig(const std::string& telescopename, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& dishdiameter, const std::vector<double>& offset, const std::vector<std::string>& mount, const std::vector<std::string>& antname, const std::vector<std::string>& padname, const std::string& coordsystem, const ::casac::variant& referencelocation)
{
  Bool rstat(false);
  try {
  
    if(itsSim !=0){
	    casacore::MPosition mpos;
      if(!casaMPosition(referencelocation, mpos) && referencelocation.toString() != "[]"){
	*itsLog << LogIO::SEVERE 
		<< "Could not convert referencelocation "
		<< referencelocation.toString() << " to a Position Measures"
		<< LogIO::POST;

	
      }
      rstat=itsSim->setconfig(telescopename, x, y, z, dishdiameter, offset,  
			      toVectorString(mount), toVectorString(antname), 
			      toVectorString(padname), 
			      coordsystem, mpos );
    }
    else
      *itsLog << LogIO::SEVERE
              << "The sm tool must be opened first."
              << LogIO::POST;
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;
  
}

bool
simulator::setfeed(const std::string& mode, const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::string>& pol)
{

Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
      rstat=itsSim->setfeed(mode,x,y,toVectorString(pol));
    }
    
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;


}

bool
simulator::setfield(const std::string& sourcename, const ::casac::variant& sourcedirection, const std::string& calcode, const ::casac::variant& distance)
{

  Bool rstat(false);
  try {
    
    if(itsSim !=0){
    
	    casacore::MDirection mdir;
      if (!casaMDirection(sourcedirection, mdir)){
	*itsLog << LogIO::SEVERE 
		<< "Could not convert source direction to a Direction Measure."
		<< LogIO::POST;
	
      }
      casacore::Quantity qdist(casaQuantity(distance));
      rstat=itsSim->setfield(sourcename, mdir, calcode, qdist);
    }
    
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;



}

bool
simulator::setmosaicfield(const std::string& sourcename, 
			  const std::string& calcode, 
			  const ::casac::variant& fieldcenter, const long xmosp,
			  const long ymosp, const ::casac::variant& mosspacing, 
			  const ::casac::variant& distance)
{


Bool rstat(false);
 try {
    
   if(itsSim !=0){
     
	   casacore::MDirection mdir;
     if (!casaMDirection(fieldcenter, mdir)){
	*itsLog << LogIO::SEVERE 
		<< "Could not convert field center to a Direction Measures"
		<< LogIO::POST;
     }
     casacore::Quantity qdist(casaQuantity(distance));
     casacore::Quantity qmosspacing(casaQuantity(mosspacing));
     rstat=itsSim->setmosaicfield(sourcename, calcode, mdir, xmosp, ymosp, 
				  qmosspacing, qdist);
   }
   
   
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
    	   << LogIO::POST;
    RETHROW(x);
 }
  
 return rstat;


}

bool
simulator::setspwindow(const std::string& spwname, 
		       const ::casac::variant& freq,
		       const ::casac::variant& deltafreq, 
		       const ::casac::variant& freqresolution, 
		       const std::string& refcode,
		       const long nchannels, const std::string& stokes)
{

Bool rstat(false);
 try {
    
   if(itsSim !=0){
     
     casacore::Quantity qfreq(casaQuantity(freq));
     casacore::Quantity qdeltafreq(casaQuantity(deltafreq));
     casacore::Quantity qfreqres(casaQuantity(freqresolution));
     MFrequency::Types freqType;
     if (!refcode.empty()) {
       String code = refcode;
       code.upcase();
       if (!MFrequency::getType(freqType, code)) {
	 *itsLog << "Invalid frequency reference '" << code
		 << "'" << LogIO::EXCEPTION;
       }
     } else {
       freqType=MFrequency::TOPO;
     }
     rstat=itsSim->setspwindow(spwname, qfreq, qdeltafreq, qfreqres, freqType,
			       nchannels,stokes);
   }
   
   
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
		   << LogIO::POST;
    RETHROW(x);
 }
  
 return rstat;

}

bool
simulator::setdata(const std::vector<long>& spwid, const std::vector<long>& fieldid, const std::string& msselect)
{

 Bool rstat(false);
 try {
   
   if(itsSim !=0){
     rstat=itsSim->setdata(spwid, fieldid, msselect);
   }
    
    
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }
  
 return rstat;
 

}

bool
simulator::predict(const std::vector<std::string>& modelImage, const std::string& complist, const bool incremental)
{

  Bool rstat(false);
  try {
    
    for (uInt i=0; i<modelImage.size(); i++) { 
      if (!modelImage[i].empty()) {
	_imageF.reset();
	_imageC.reset();
    // NOTE. fromFile() now returns more than two pointers in the tuple	
    std::tie(_imageF, _imageC, std::ignore, std::ignore) = ImageFactory::fromFile(modelImage[i]);
    ThrowIf(
        ! (_imageF || _imageC),
        "Unsupported image data type"
    );

	
	auto unit =_imageF
	  ? _imageF->units().getName()
	  : _imageC->units().getName();
	String jyPix("Jy/pixel");
	
	if (unit!=jyPix) {
	  LogIO os(LogOrigin("simulator", "predict"));	
	  os << LogIO::WARN 
	     << "sm.predict will overwrite the brightness unit " << unit 
	     << " in your image " << modelImage[i]
	     << LogIO::POST;
	}
      }
    }
    _imageF.reset();
    _imageC.reset();
    
    if(itsSim !=0){
      rstat=itsSim->predict(toVectorString(modelImage), complist, incremental);
    }
    
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;

}

bool
simulator::setoptions(const std::string& ftmachine, const long cache, 
		      const long tile, const std::string& gridfunction, 
		      const ::casac::variant& location, const double padding,
		      const long facets, const double maxdata, 
		      const long wprojplanes)
{

  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
	    casacore::MPosition mpos;
      if (!casaMPosition(location, mpos) && (location.toString()!="[]")){
	*itsLog << LogIO::SEVERE 
		<< "Could not convert location "
		<< location.toString() << " to a Position Measures"
		<< LogIO::POST;
      }
      rstat=itsSim->setoptions(ftmachine, cache, tile, gridfunction, mpos, padding, facets, maxdata, wprojplanes);
    }
    
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;

}

bool
simulator::setvp(const bool dovp, const bool usedefaultvp, 
		 const std::string& vptable, const bool dosquint, 
		 const ::casac::variant& parangleinc, 
		 const ::casac::variant& skyposthreshold, const double pblimit)
{

  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      casacore::Quantity skyposthr(180, "deg");
      casacore::Quantity parang(360, "deg");
      if ((String(parangleinc.toString()) != casacore::String("")) &&
	  (String(parangleinc.toString()) != casacore::String("[]")) )
	parang=casaQuantity(parangleinc);
      if ((String(skyposthreshold.toString()) != casacore::String("")) &&
	  (String(skyposthreshold.toString()) != casacore::String("[]")))
	skyposthr=casaQuantity(skyposthreshold);
      rstat=itsSim->setvp(dovp, usedefaultvp, vptable, dosquint,parang, 
			  skyposthr, pblimit);
    }
    
    
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }
  
 return rstat;
   



}

bool
// simulator::corrupt()
simulator::corrupt()
{

  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      rstat=itsSim->corrupt();
    }
    
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;
}

bool
simulator::reset()
{
  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      rstat=itsSim->resetviscal();
      rstat=itsSim->resetimcal();
    }
    
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;

}

bool
simulator::setbandpass(const std::string& mode, const std::string& table, 
		       const ::casac::variant& interval, 
		       const std::vector<double>& amplitude)
{

  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
      casacore::Quantity qinter(casaQuantity(interval));
      rstat=itsSim->setbandpass(mode, table, qinter, amplitude);
    }
    
    
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }
  
 return rstat;
   

}

bool 
simulator::setapply(const std::string& table, 
		    const std::string& type, 
		    const double t, 		    
		    const ::casac::variant& field,
		    const std::string& interp,
		    const bool calwt,
		    const std::vector<long>& spwmap,
		    const double opacity)
{
  Bool rstat(false);
  try {
    
    LogIO os(LogOrigin("simulator", "setapply"));
    os << "Beginning setapply--------------------------" << LogIO::POST;

    // Forward to the Simulator object
    if (itsSim)
      rstat = itsSim->setapply(type,t,table,
			       "",toCasaString(field),
			       interp,calwt,spwmap,opacity);
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;

}



// RI TODO make pwv and windspeed variants here and quantities in Simulator.cc

bool
simulator::settrop(const std::string& mode, const std::string& table, 
		   const double pwv, const double deltapwv, 
		   const double beta, const double windspeed)
{
  Bool rstat(false);
  try {
    
    if(itsSim !=0){  
      //      casacore::Quantity qinter(casaQuantity(interval));
      rstat=itsSim->settrop(mode, table, pwv, deltapwv, beta, windspeed);
    }
    // RI TODO interpolation params have to get to SolvableVisCal::setApply.
    // RI TODO do we make the user call sm.setapply to deal with that, 
    // RI TODO or do we have it pass through here to VC::setApply? 
      
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }
 return rstat;  
}






bool
simulator::setgain(const std::string& mode, const std::string& table, 
		   const ::casac::variant& interval, 
		   const std::vector<double>& amplitude)
{
  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
      casacore::Quantity qint(casaQuantity(interval));
      rstat=itsSim->setgain(mode, table, qint, amplitude);
    }
    
    
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }
  
 return rstat;
   
}




bool 
simulator::setpointingerror(const std::string& epjtablename, 
			    const bool applypointingoffsets, 
			    const bool dopbcorrection){

  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
      rstat=itsSim->setpointingerror(epjtablename, applypointingoffsets, dopbcorrection);
    }
    
    
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }
  
 return rstat;
   

}



bool
simulator::setleakage(const std::string& mode, const std::string& table, 
		      //const ::casac::variant& interval, 
		      const std::vector<double>& amplitude,
		      const std::vector<double>& offset)
{

  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
      //casacore::Quantity qinter(casaQuantity(interval));
      //rstat=itsSim->setleakage(mode, table, qinter, amplitude);
      rstat=itsSim->setleakage(mode, table, amplitude, offset);
    }
    
    
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }
  
 return rstat;
   
   
}

bool
simulator::oldsetnoise(const std::string& mode, const std::string& table, 
		    const ::casac::variant& simplenoise, 
		    const double antefficiency, const double correfficiency, 
		    const double spillefficiency, const double tau, 
		    const double trx, const double tatmos, const double tcmb)
{
  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
      LogIO os(LogOrigin("simulator", "oldsetnoise"));
      os << LogIO::WARN << "Using deprecated ACoh Noise - this will dissapear in the future - please switch to sm.setnoise" << LogIO::POST;
      
      casacore::Quantity qnoise(casaQuantity(simplenoise));
      rstat=itsSim->oldsetnoise(mode, table, qnoise, antefficiency, correfficiency, spillefficiency, tau, trx, tatmos, tcmb);
    }
    
    
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }  
 return rstat;
}


bool
simulator::setnoise(const std::string& mode, 
		    const std::string& table, 
		    const ::casac::variant& simplenoise, 
		    // atm parameters
		    const ::casac::variant& pground,
		    const double relhum,
		    const ::casac::variant& altitude,
		    const ::casac::variant& waterheight,
		    const ::casac::variant& pwv,
		    // OR specify tau and tatmos 
		    const double tatmos, 
		    const double tau,
		    // antenna parameters
		    const double antefficiency, 
		    const double spillefficiency, 
		    const double correfficiency,    
		    const double trx, 
		    const double tground, 
		    const double tcmb,
		    const bool OTF,
		    const double senscoeff, 
		    const long rxtype
		    ) {
  Bool rstat(false);
  try {
    if(itsSim !=0){      
      casacore::Quantity qnoise(casaQuantity(simplenoise));
      casacore::Quantity qpress(casaQuantity(pground));
      casacore::Quantity qalt(casaQuantity(altitude));
      casacore::Quantity qwaterht(casaQuantity(waterheight));
      casacore::Quantity qpwv(casaQuantity(pwv));
#ifdef RI_DEBUG
      cout<<qnoise<<" "<<qpress<<" "<<qalt<<" "<<qwaterht<<" "<<qpwv<<endl;
#endif
      rstat=itsSim->setnoise(mode, table, qnoise, 
			     qpress,relhum,qalt,qwaterht,qpwv,
			     tatmos,tau,
			     antefficiency, spillefficiency, correfficiency, 
			     trx, tground, tcmb, OTF, senscoeff, rxtype);
    }    
 } catch  (AipsError x) {
   *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	   << LogIO::POST;
    RETHROW(x);
 }  
 return rstat;
}



bool
simulator::setpa(const std::string& mode, const std::string& table, 
		 const ::casac::variant& interval)
{


  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
      casacore::Quantity qinter(casaQuantity(interval));
      rstat=itsSim->setpa(mode, table, qinter);
    }
    
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;
   

}

bool
simulator::setseed(const long seed)
{

  Bool rstat(false);
  try {
    
    if(itsSim !=0){
      
      rstat=itsSim->setseed(seed);
    }
    
    
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() 
	    << LogIO::POST;
    RETHROW(x);
  }
  
  return rstat;
   


}



} // casac namespace

