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
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/Assert.h>
#include <ms/MeasurementSets.h>
#include <ms/MeasurementSets/MSHistoryHandler.h>
#include <casa/Logging/LogIO.h>

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

std::string synthesisimstore::getImageName(const std::string& imageId, long taylorTerm)
{
    if (imageId == "MASK") {
        return itsImStore->mask((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "PSF") {
        return itsImStore->psf((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "MODEL") {
        return itsImStore->model((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "RESIDUAL") {
        return itsImStore->residual((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "WEIGHT") {
        return itsImStore->weight((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "IMAGE") {
        return itsImStore->image((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "SUMWT") {
        return itsImStore->sumwt((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "GRIDWT") {
        throw AipsError("Retrieval of gridwt image name not supported at this time.");
    } else if (imageId == "PB") {
        return itsImStore->pb((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "FORWARDGRID") {
        return itsImStore->forwardGrid((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "BACKWARDGRID") {
        return itsImStore->backwardGrid((casacore::uInt)taylorTerm)->name();
    } else if (imageId == "IMAGEPBCOR") {
        return itsImStore->imagepbcor((casacore::uInt)taylorTerm)->name();
    } else {
        throw AipsError("Image id \""+imageId+"\" not recognized.");
    }
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
