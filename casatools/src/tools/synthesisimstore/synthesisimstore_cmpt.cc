/***
 * Framework independent implementation file for imager...
 *
 * Implement the imager component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 
 * @author Wes Young
 * @version 
 ***/

#include <iostream>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/ms/MeasurementSets.h>
#include <casacore/ms/MeasurementSets/MSHistoryHandler.h>
#include <casacore/casa/Logging/LogIO.h>

#include <synthesis/ImagerObjects/SIImageStore.h>

#include <synthesisimstore_cmpt.h>

using namespace std;
using namespace casacore;
using namespace casa;

     
using namespace casacore;
namespace casac {

 synthesisimstore::synthesisimstore()
{
  itsImStore = new SIImageStore();
  containsimage=false;
}

  synthesisimstore::synthesisimstore(casa::SIImageStore* imstore)
{
  itsImStore = imstore;  /// Or use the constructor from the internal pointers.
  containsimage=true;
}

casa::SIImageStore* synthesisimstore::getImStore()
{
  containsimage=false; // why??? get* methods should not have side effects....
  return itsImStore;
}


synthesisimstore::~synthesisimstore()
{
  done();
}

bool
synthesisimstore::done()
{
  Bool rstat(false);

  try 
    {
      if (itsImStore )
	{
	  if(containsimage)
	    {
	      delete itsImStore;
	    }
	  itsImStore=NULL;
	}
    } 
  catch  (AipsError x) 
    {
      RETHROW(x);
    }
  
  return rstat;
}



} // casac namespace
