/***
 * Framework independent implementation file for atcafiller...
 *
 * Implement the atcafiller component here.
 * 
 * // The filler for ATCA data in the RPFITS format
 * // Resurrected from the pre-CABB, pre-svn version.
 * // Updated for CABB, Jan 2015.
 *
 * @author
 * @version 
 ***/

//#include <iostream>

#include <atnf/atca/ATCAFiller.h>
#include <atcafiller_cmpt.h>
#include <stdcasa/StdCasa/CasacSupport.h>

using namespace std;

using namespace casacore;
using namespace casa;

namespace casac {

atcafiller::atcafiller()
{
  itsATCAFiller = new casa::ATCAFiller();
}


atcafiller::~atcafiller() 
{
  if (itsATCAFiller) {
    delete itsATCAFiller;
    itsATCAFiller=0;
  }
}

bool atcafiller::open(const std::string& msName,
		 const std::vector<std::string>& rpfitsFiles,
		 const std::vector<std::string>& options)
{
    Vector<String> const rpfitsFilesV(rpfitsFiles);
    Vector<String> const optionsV(options);
    itsATCAFiller->open(msName,rpfitsFilesV,optionsV,0);
    return true;
}

bool atcafiller::fill()
{
    return itsATCAFiller->fill();
}

bool atcafiller::select(long firstScan, long lastScan, 
       const std::vector<long>& spw, 
		   const Quantity& lowFreq, const Quantity& highFreq,
		   const std::vector<std::string>& fieldList,float edge)
{
    itsATCAFiller->scanRange(firstScan,lastScan);
    Vector<Int> const spwV(spw);
    itsATCAFiller->freqSel(spwV);
    itsATCAFiller->freqRange(casaQuantity(lowFreq).getValue("Hz"),
			    casaQuantity(highFreq).getValue("Hz"));
    Vector<String> const fieldListV(fieldList);
    itsATCAFiller->fields(fieldListV);
    itsATCAFiller->edge(edge);
    return true;
}

std::vector<double> atcafiller::opacities(const std::vector<double>& fGHz,float tempK,
  float humi, float press, float height)
{
  Vector<Double> const fGHzV(fGHz);
  Vector<Double> v=itsATCAFiller->opacities(fGHzV,tempK,humi,press,height);
  std::vector<double> out(0);
  v.tovector(out);
  return out;  
}

} // casac namespace

