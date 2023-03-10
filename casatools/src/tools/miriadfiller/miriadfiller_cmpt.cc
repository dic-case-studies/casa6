
/***
 * Framework independent implementation file for miriadfiller...
 *
 * Implement the miriadfiller component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 *
 * @author
 * @version 
 ***/

#include <iostream>
#include <miriadfiller_cmpt.h>

#include <miriad/Filling/Importmiriad.h>
#include <casacore/casa/Containers/Block.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/OS/Path.h>
#include <casacore/casa/OS/SymLink.h>
#include <casacore/casa/Quanta/MVEpoch.h>
#include <casacore/casa/Quanta/MVFrequency.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/BasicSL/String.h>
#include <stdcasa/StdCasa/CasacSupport.h>


using namespace std;

using namespace casacore;
using namespace casa;

namespace casac {

miriadfiller::miriadfiller()
{

}

miriadfiller::~miriadfiller()
{

}

void
miriadfiller::fill(const std::string& vis, 
		    const std::string& inputfile, 
		    const bool tsys, 
		    const std::vector<long>& spw, 
		    const std::string& vel, 
		    const bool linecal, 
		    const std::vector<long>& wide, 
		    const long debug) 
{
    if (inputfile=="")
      throw(AipsError("No miriad file given"));
    
    String visfile=vis, infile=inputfile;
    if(vis.size()==0) visfile=infile.before('.') + ".ms";
    // convert vector to Block to make Importmiriad happy
    Block<Int> spws(spw.size()), wides(wide.size());
    for (uint i=0; i<spw.size(); i++) spws[i]=spw[i];
    for (uint i=0; i<wide.size(); i++) wides[i]=wide[i];

    File t(inputfile);                      // only used for sanity checks
    if (!t.isDirectory())
      throw(AipsError("Input file does not appear to be a miriad dataset"));
  try {
    Importmiriad bf(infile,debug,tsys,false,linecal);

    bf.checkInput(spws,wides);
    bf.setupMeasurementSet(visfile);
    bf.fillAntennaTable();    // put first array in place


    // fill the main table
    bf.fillMSMainTable();


    // fill remaining tables
    bf.fillSyscalTable();      
    bf.fillSpectralWindowTable(vel);   
    bf.fillFieldTable();
    bf.fillSourceTable();
    bf.fillFeedTable();
    bf.fixEpochReferences();
    bf.fillObsTables();
    bf.close();
    //cout << "Importmiriad::close()  Created MeasurementSet " << vis << endl;
    cout << "Created MeasurementSet " << vis << endl;
  } 
  catch (AipsError x) {
    cerr << "Error: " << x.getMesg() << endl;
  } 
  return;
}

} // casac namespace

