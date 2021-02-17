/***
 * Framework independent implementation file for imager...
 *
 * Implement the imager component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 
 ***/

#include <iostream>
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/Assert.h>
#include <casa/Logging/LogIO.h>

#include <casa/OS/Directory.h>
#include <images/Images/PagedImage.h>

#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>

#include <synthesisutils_cmpt.h>

using namespace std;
using namespace casacore;
using namespace casa;

     
using namespace casacore;
namespace casac {

synthesisutils::synthesisutils() 
{
  itsUtils = new SynthesisUtilMethods();
}

synthesisutils::~synthesisutils()
{
  done();
}

  ::casac::record* synthesisutils::advisechansel(const ::casac::variant& freqstart, const ::casac::variant& freqend, const ::casac::variant& freqstep, const std::string& freqframe, const std::string& ephemtab, const std::string& msname, const long fieldid, const bool getfreqrange, const std::string& spwselection){

  casac::record* retval=0;
  LogIO theLog;
      try {

	if( msname.length() > 0){

	  Vector<Int>  spw;
	  Vector<Int>  start;
	  Vector<Int>  nchan;
	  MFrequency::Types tp;
	  String specframe(upcase(String(freqframe)));
	  
	  if(freqframe=="SOURCE"){
	    tp=MFrequency::REST;
	  }
	  else if(!MFrequency::getType(tp, freqframe))
	    throw(AipsError("Invalid frequency frame"));
	  Double fstart, fend, fstep;
	  /***if(freqstart.type()==::casac::variant::INT || freqend.type()==::casac::variant::INT || freqstep.type()==::casac::variant::INT)
	    throw(AipsError("Due to a known bug (CAS-13149) in casa parameter passing \nDo not pass a python integer as parameter here... use a period '.' at the end of the number to make sure it is interpreted as double")); 
          ***/
	  if(freqstart.type() > ::casac::variant::BOOL && freqstart.type() < ::casac::variant::COMPLEX){
	    fstart=freqstart.toDouble();
	  }
	  else
	    fstart=casaQuantity(freqstart).get("Hz").getValue();
	
	  //cerr << "FSTART " << fstart << "  quant " << casaQuantity(freqstart).get("Hz") << endl;
	  if(freqend.type() > ::casac::variant::BOOL && freqend.type() <  ::casac::variant::COMPLEX){
	    fend=freqend.toDouble();
	  }
	  else
	    fend=casaQuantity(freqend).get("Hz").getValue();
	  if(freqstep.type()> ::casac::variant::BOOL && freqstep.type() < ::casac::variant::COMPLEX){
	    fstep=freqstep.toDouble();
	  }
	  else
	    fstep=casaQuantity(freqstep).get("Hz").getValue();
	  if(itsUtils->adviseChanSel(fstart, fend, fstep, tp, spw, start, nchan, msname, ephemtab, fieldid, getfreqrange, spwselection)){
	    Record outRec;
	    if(!getfreqrange){
		
		outRec.define("spw", spw);
		outRec.define("start", start);
		outRec.define("nchan", nchan);
	      
	    }
	    else{
	      //outRec.define("freqstart", fstart);
	      //outRec.define("freqend", fend);
	      casacore::QuantumHolder qh(casacore::Quantity(fstart, "Hz"));
	      outRec.defineRecord("freqstart", qh.toRecord());
	      qh=casacore::QuantumHolder(casacore::Quantity(fend, "Hz"));
	      outRec.defineRecord("freqend", qh.toRecord());
	      
	    }
	    retval=fromRecord(outRec);
	  }
	} else {
	  theLog << LogIO::SEVERE << "No MeasurementSet has been assigned, please give a valid ms in msname" << LogIO::POST;
	}

	
      } catch  (AipsError x) {
          //*itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg() << LogIO::POST;
	  RETHROW(x);
      }


      return retval;



  }
  
  casac::record* synthesisutils::contdatapartition(const casac::record& selpars, const long npart)
{
  casac::record* rstat(0);

  try 
    {
      std::unique_ptr<casacore::Record> recpars(toRecord( selpars ));
      rstat = fromRecord(  itsUtils->continuumDataPartition( *recpars , npart ) );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

  casac::record* synthesisutils::cubedatapartition(const casac::record& selpars, const long npart, const ::casac::variant& fstart, const ::casac::variant&  fend, const string& frame)
{
  casac::record* rstat(0);

  try 
    {
      std::unique_ptr<casacore::Record> recpars(toRecord( selpars ));
      casacore::Quantity qstart(1.0, "GHz");
      casacore::Quantity qend(1.5,"GHz");
      if( (fstart.toString().size() != 0) && String(fstart.toString()) != String("[]"))
	qstart=casaQuantity(fstart);
      if( (fend.toString().size() != 0) && String(fend.toString()) != String("[]"))
	qend=casaQuantity(fend);
      
      casacore::MFrequency::Types eltype;
      casacore::MFrequency::getType(eltype, frame);
      rstat = fromRecord(  itsUtils->cubeDataPartition( *recpars , npart, qstart.getValue("Hz"), qend.getValue("Hz"), eltype ) );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

  casac::record* synthesisutils::cubeimagepartition(const casac::record& selpars, const long npart)
{
  casac::record* rstat(0);

  try 
    {
      std::unique_ptr<casacore::Record> recpars(toRecord( selpars ));
      rstat = fromRecord(  itsUtils->cubeImagePartition( *recpars , npart ) );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

  casac::record* synthesisutils::cubedataimagepartition(const casac::record& selpars, const casac::record& incsysrec,
                                                      const long npart, const long nchannel)
{
   casac::record* rstat(0);

   try
     {
        std::unique_ptr<casacore::Record> recselpars(toRecord( selpars ));
        std::unique_ptr<casacore::Record> recincsys(toRecord( incsysrec ));
        if (recincsys->nfields() != 0 ) {
          CoordinateSystem *incsys;
          incsys = CoordinateSystem::restore(*recincsys,"coordsys");
          Vector<CoordinateSystem> ocsysvec;
          Vector<Int> outnchanvec;
          rstat = fromRecord( itsUtils->cubeDataImagePartition( *recselpars, *incsys, npart, nchannel, ocsysvec, outnchanvec) );
        }
     }
   catch (AipsError x)
     {
       RETHROW(x);
     }

   return rstat;
}

  casac::record* synthesisutils::checkselectionparams(const casac::record& selpars)
{
  casac::record* rstat(0);

  try 
    {
      std::unique_ptr<casacore::Record> recpars(toRecord( selpars ));
      SynthesisParamsSelect pars;
      pars.fromRecord( *recpars );
      rstat = fromRecord(  pars.toRecord()  );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

  casac::record* synthesisutils::checkimageparams(const casac::record& impars)
{
  casac::record* rstat(0);

  try 
    {
      std::unique_ptr<casacore::Record> recpars(toRecord( impars ));
      SynthesisParamsImage pars;
      pars.fromRecord( *recpars );
      rstat = fromRecord(  pars.toRecord()  );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}


  casac::record* synthesisutils::checkgridparams(const casac::record& gridpars)
{
  casac::record* rstat(0);

  try 
    {
      std::unique_ptr<casacore::Record> recpars(toRecord( gridpars ));
      SynthesisParamsGrid pars;
      pars.fromRecord( *recpars );
      rstat = fromRecord(  pars.toRecord()  );
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}

  casac::record*  synthesisutils::updateimpars(const casac::record& impars)
{
   casac::record* rstat(0);

  try
    {
      SynthesisParamsImage pars;
      std::unique_ptr<casacore::Record> recpars(toRecord( impars ));
      rstat = fromRecord( pars.updateParams( *recpars ) );
    }
  catch  (AipsError x)
    {
      RETHROW(x);
    }

  return rstat;
}

/***
 bool 
 synthesisutils::makeimage(const casac::record& impars, const casac::record& selpars, const string& msname)
{
  bool rstat(0);

  try 
    {
      // Construct parameter object, and verify params.
      casacore::Record recpars = *toRecord( impars );
      SynthesisParamsImage pars;
      pars.fromRecord( recpars );  // will throw exception if parameters are invalid

      // Construct Coordinate system and make image. 
      MeasurementSet ms;
      //if( msname.length() > 0 && (Directory(msname)).exists() ) { ms = MeasurementSet(msname); }
      ms = MeasurementSet(msname);

      cout << "Making image : " << pars.imageName << " of shape : " << pars.shp() << endl;
      PagedImage<Float> diskimage( pars.shp(), pars.buildCoordinateSystem(ms), pars.imageName );

    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }

  return rstat;
}
***/


  long  synthesisutils::getOptimumSize(const long size)
{
   int rstat(size);

  try
    {
      rstat = SynthesisUtilMethods::getOptimumSize(size);
    }
  catch  (AipsError x)
    {
      RETHROW(x);
    }

  return rstat;
}
  
  bool  synthesisutils::fitPsfBeam(const string& imagename, const long nterms)
{
   int rstat(false);

  try
    {
      rstat = SynthesisUtilMethods::fitPsfBeam(imagename, nterms);
    }
  catch  (AipsError x)
    {
      RETHROW(x);
    }

  return rstat;
}
  


bool
synthesisutils::done()
{
  Bool rstat(false);

  try 
    {
      if (itsUtils)
	{
	  delete itsUtils;
	  itsUtils=NULL;
	}
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return rstat;
}

} // casac namespace
