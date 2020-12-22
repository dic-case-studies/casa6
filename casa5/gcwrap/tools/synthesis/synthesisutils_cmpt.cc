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
