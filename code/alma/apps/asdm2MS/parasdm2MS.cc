// NOTE: this code is not currently being used or developed. It _may_ have been checked to see if
// it still compiles but it should not be trusted as being up to date or remotely correct. It is
// being kept around as an an example of a previous attempt to do something like this.
// Use with extreme caution.

// #ifdef _OPENMP
// #include <omp.h>
// #endif 
//#else
//  #define omp_get_num_threads() 0
//  #define omp_get_thread_num() 0
//#endif

#define DDPRIORITY 1
#include <omp.h>
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <regex>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <ASDMAll.h>
#include <Misc.h>

#include "SDMBinData.h"
using namespace sdmbin;

#include <exception>
using namespace asdm;
#include "IllegalAccessException.h"

#include "UvwCoords.h"
#include "ASDM2MSFiller.h"

#include "measures/Measures/Stokes.h"
#include "measures/Measures/MFrequency.h"
using namespace casacore;

#include <tables/Tables/Table.h>
#include <tables/Tables/PlainTable.h>
#include <tables/Tables/TableCopy.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/OS/Path.h>
#include <casa/OS/Directory.h>
#include <casa/OS/DirectoryIterator.h>
#include <casa/OS/File.h>
#include <casa/OS/RegularFile.h>
#include <casa/OS/SymLink.h>


#include "CBasebandName.h"
#include "CCalibrationDevice.h"
using namespace CalibrationDeviceMod;
#include "CFrequencyReferenceCode.h"
#include "CPolarizationType.h"
#include "CPrimitiveDataType.h"
#include "CProcessorSubType.h"
#include "CProcessorType.h"
using namespace ProcessorTypeMod;
#include "CScanIntent.h"
#include "CSubscanIntent.h"
using namespace SubscanIntentMod;
#include "CSpectralResolutionType.h"
using namespace SpectralResolutionTypeMod;
#include "CStokesParameter.h"

#include "Name2Table.h"
#include "ASDMVerbatimFiller.h"

#include "SDMDataObjectStreamReader.h"
#include "SDMDataObjectReader.h"
#include "SDMDataObject.h"
using namespace asdmbinaries;

#include "TableStreamReader.h"
#include "asdm2MSGeneric.h"

#include "AnyValueMap.hpp"

#include "asdmstman/AsdmStMan.h"
#include "BDF2AsdmStManIndex.h"

#include "ScansParser.h"

#include	"time.h"			/* <time.h> */
#if	defined(__sysv__)
#include	<sys/time.h>
#elif	defined(__bsd__)
#define		ftime	FTIME
#include	<sys/timeb.h>			/* from system */
#undef		ftime
extern	void	ftime( struct timeb * );	/* this is the funtion */
#else
#endif


void	myTimer( double *cpu_time ,		/* cpu timer */
		 double *real_time ,		/* real timer */
		 int  *mode )			/* the mode */
{
  clock_t	tc;				/* clock time */
  double	ct;				/* cpu time in seconds */
  double	rt = 0.0 ;           		/* real time in seconds */

#if	defined(__sysv__)
  struct timeval 	Tp;
  struct timezone	Tzp;
#elif	defined(__bsd__)
  struct timeb tr;				/* struct from ftime */
#else
#endif
  tc = clock( );				/* get clock time */
  ct = (double)(tc) / (double)CLOCKS_PER_SEC;	/* to seconds */
#if	defined(__sysv__)
  gettimeofday( &Tp, &Tzp );			/* get timeofday */
  rt = (double) Tp.tv_sec + 0.000001 * (double) Tp.tv_usec;
#elif	defined(__bsd__)
  ftime( &tr );				/* get real time */
  rt = (double) tr.time + 0.001 * (double) tr.millitm;	/* in seconds */
#else
#endif
  if (*mode) {					/* calculate difference */
    (*cpu_time)  = ct - (*cpu_time);		/* cpu time */
    (*real_time) = rt - (*real_time);		/* real time */
  } else {
    (*cpu_time)  = ct;			/* set cpu time */
    (*real_time) = rt;			/* set real time */
  }
}


using namespace casacore;
//using namespace casa;

ASDM2MSFiller* msFiller;
string appName;
bool verbose = true;
bool isEVLA = false;
bool lazy = false;

//LogIO os;

#include <casa/Logging/StreamLogSink.h>
#include <casa/Logging/LogSink.h>

void info (const string& message) {
  
  if (!verbose){
    return;
  }

  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
}

void warning (const string& message) {
  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
}

void error(const string& message) {
  LogSink::postGlobally(LogMessage(message, LogOrigin(appName,WHERE), LogMessage::NORMAL));
  //os << LogIO::POST;
  exit(1);
}

#include <iostream>
#include <sstream>

/*
** A simplistic tracing toolbox.
*/
bool debug = (getenv("ASDM_DEBUG") != NULL);
vector<char> logIndent;
// #define LOGENTER(name) if (debug) {for_each(logIndent.begin(), logIndent.end(), cout << _1); logIndent.push_back('\t'); cout << #name ": entering" << endl;}
// #define LOGEXIT(name)  if (debug) {logIndent.pop_back(); for_each(logIndent.begin(), logIndent.end(), cout << _1); cout << #name ": exiting" << endl;}
// #define LOG(msg) if (debug) {for_each(logIndent.begin(), logIndent.end(), cout << _1); cout << msg << endl;}


ostringstream errstream;
ostringstream infostream;

using namespace std;
template<typename T> string TO_STRING(const T &v) {
    const T shift = 100000000;
    const T rnd = .000000005;
    char buffer[128];
    int whole = int(v);
    int fraction = int((v-whole+rnd)*shift);
    sprintf(buffer,"%d.%08d",whole,fraction);
    return string(buffer);
}

string TO_STRING(const int &v) {
    char buffer[128];
    sprintf(buffer,"%d",v);
    return string(buffer);
}


//
// A class to describe Exceptions.
//
class ASDM2MSException {  
public:
  /**
   * An empty contructor.
   */
  ASDM2MSException();
  
  /**
   * A constructor with a message associated with the exception.
   * @param m a string containing the message.
   */
  ASDM2MSException(string m);
  
  /**
   * The destructor.
   */
  virtual ~ASDM2MSException();
  
  /**
   * Returns the message associated to this exception.
   * @return a string.
   */
  string getMessage() const;
  
protected:
  string message;
  
};

inline ASDM2MSException::ASDM2MSException() : message ("ASDM2MSException") {}
inline ASDM2MSException::ASDM2MSException(string m) : message(m) {}
inline ASDM2MSException::~ASDM2MSException() {}
inline string ASDM2MSException::getMessage() const {
  return "ASDM2MSException : " + message;
}

// A facility to get rid of blanks at start and end of a string.
// 
string lrtrim(std::string& s,const std::string& drop = " ")
{
  std::string r=s.erase(s.find_last_not_of(drop)+1);
  return r.erase(0,r.find_first_not_of(drop));
}


// These classes provide mappings from some ALMA Enumerations to their CASA counterparts.
class StokesMapper {
private :
  Stokes::StokesTypes* sa;

public :
  StokesMapper();
  ~StokesMapper();

  static Stokes::StokesTypes value(StokesParameterMod::StokesParameter s);
  Stokes::StokesTypes* to1DArray(const vector<StokesParameterMod::StokesParameter>& v);
  static vector<Stokes::StokesTypes> toVectorST(const vector<StokesParameterMod::StokesParameter>& v);
  static vector<int> toVectorI(const vector<StokesParameterMod::StokesParameter>& v);
}
  ;

StokesMapper::StokesMapper() {
  sa = 0;
}

StokesMapper::~StokesMapper() {
  if (sa)  delete[] sa;
}

Stokes::StokesTypes StokesMapper::value(StokesParameterMod::StokesParameter s) {
  switch (s) {
  case I : return Stokes::I; 
  case Q : return Stokes::Q; 
  case U : return Stokes::U; 
  case V : return Stokes::V; 
  case RR : return Stokes::RR; 
  case RL : return Stokes::RL; 
  case LR : return Stokes::LR; 
  case LL : return Stokes::LL; 
  case XX : return Stokes::XX; 
  case XY : return Stokes::XY; 
  case YX : return Stokes::YX; 
  case YY : return Stokes::YY; 
  case RX : return Stokes::RX; 
  case RY : return Stokes::RY; 
  case LX : return Stokes::LX; 
  case LY : return Stokes::LY; 
  case XR : return Stokes::XR; 
  case XL : return Stokes::XL; 
  case YR : return Stokes::YR; 
  case YL : return Stokes::YL; 
  case PP : return Stokes::PP; 
  case PQ : return Stokes::PQ; 
  case QP : return Stokes::QP; 
  case QQ : return Stokes::QQ; 
  case RCIRCULAR : return Stokes::RCircular; 
  case LCIRCULAR : return Stokes::LCircular; 
  case LINEAR : return Stokes::Linear; 
  case PTOTAL : return Stokes::Ptotal; 
  case PLINEAR : return Stokes::Plinear; 
  case PFTOTAL : return Stokes::PFtotal; 
  case PFLINEAR : return Stokes::PFlinear; 
  case PANGLE : return Stokes::Pangle;  
  }
  return Stokes::Undefined;
}

Stokes::StokesTypes* StokesMapper::to1DArray(const vector<StokesParameterMod::StokesParameter>& v) {
  if (v.size() == 0) return 0;

  if (sa) {
    delete[] sa;
    sa = 0;
  }

  sa = new Stokes::StokesTypes[v.size()];
  for (unsigned int i = 0; i < v.size(); i++) 
    sa[i] = value(v.at(i));

  return sa;
}

vector<int> StokesMapper::toVectorI(const vector<StokesParameterMod::StokesParameter>& v) {
  vector<int> result;

  for (unsigned int i = 0; i < v.size(); i++)
    result.push_back(value(v[i]));

  return result;
}

vector<Stokes::StokesTypes> StokesMapper::toVectorST(const vector<StokesParameterMod::StokesParameter>& v) {
  vector<Stokes::StokesTypes> result;

  for (unsigned int i = 0; i < v.size(); i++)
    result.push_back(value(v[i]));

  return result;
}

class FrequencyReferenceMapper {

public :
  FrequencyReferenceMapper();
  ~FrequencyReferenceMapper();

  static MFrequency::Types value(FrequencyReferenceCodeMod::FrequencyReferenceCode frc);
}
  ;

FrequencyReferenceMapper::FrequencyReferenceMapper() {
  ;
}

FrequencyReferenceMapper::~FrequencyReferenceMapper() {
  ;
}

MFrequency::Types FrequencyReferenceMapper::value(FrequencyReferenceCodeMod::FrequencyReferenceCode frc) {
  switch (frc) {
  case FrequencyReferenceCodeMod::LABREST : 
    errstream.str("");
    errstream << " Can't map FrequencyReferenceCode::LABREST to an MFrequency::Types" << endl;
    error(errstream.str());
    break;
    
  case FrequencyReferenceCodeMod::LSRD : return MFrequency::LSRD; 
  case FrequencyReferenceCodeMod::LSRK : return MFrequency::LSRK; 
  case FrequencyReferenceCodeMod::BARY : return MFrequency::BARY;
  case FrequencyReferenceCodeMod::REST : return MFrequency::REST;
  case FrequencyReferenceCodeMod::GEO  : return MFrequency::GEO; 
  case FrequencyReferenceCodeMod::GALACTO : return MFrequency::GALACTO;
  case FrequencyReferenceCodeMod::TOPO : return MFrequency::TOPO; 
  }
  return MFrequency::TOPO; // Never happens.
}

class PolTypeMapper {
private :
  vector<string> polType;

public :
  PolTypeMapper();
  ~PolTypeMapper();

  static char value(PolarizationTypeMod::PolarizationType p);
  vector<string> toStringVector(const vector<PolarizationTypeMod::PolarizationType>& v);
};


PolTypeMapper::PolTypeMapper() {
  ;
}

PolTypeMapper::~PolTypeMapper() {
  ;
}

char PolTypeMapper::value(PolarizationTypeMod::PolarizationType p) {
  return (CPolarizationType::name(p)).at(0);
}

vector<string> PolTypeMapper::toStringVector(const vector<PolarizationTypeMod::PolarizationType>& v) {
  polType.clear();
  for (unsigned int i = 0; i < v.size(); i++) {
    polType.push_back((CPolarizationType::name(v.at(i))));
  }
  return polType;
}


// These classes provide methods to convert from vectors of "anything" into vectors of basic types.
class FConverter {
public:

  static vector<float> toVectorF(const vector<vector<float> >& vvF, bool tranpose=false);

  template<class T>
  static vector<float> toVectorF(const vector<vector<T> >& vvT, bool transpose=false) {
    vector<float> result;

    if (transpose == false) {
      //
      // Simply linearize the vector of vectors.
      //
      for (unsigned int i = 0; i < vvT.size(); i++)
	for (unsigned int j = 0; j < vvT.at(i).size(); j++)
	  result.push_back(vvT.at(i).at(j).get());
    }
    else {
      //
      // We want to transpose.
      // Let's consider the very general case where our vector of vectors does not represent
      // a rectangular matrix, i.e. all the elements of this vector are vectors with possibly
      // different size. 
      //
      unsigned int maxsize = 0;
      unsigned int cursize = 0;
      for (unsigned int i = 0; i < vvT.size(); i++) {
	cursize = vvT.at(i).size();
	if (cursize > maxsize) maxsize = cursize;
      }
      
      for (unsigned int i = 0; i < maxsize; i++)
	for (unsigned int j = 0; j < vvT.size(); j++)
	  if (i < vvT.at(j).size())
	    result.push_back(vvT.at(j).at(i).get());     
    }
    return result;
  }
};

float d2f(double d) { return (float) d; }

vector<float> FConverter::toVectorF(const vector< vector <float> > & vvF, bool transpose) {
  vector<float> result;
  if (transpose == false) {
    //
    // Simply linearize the vector of vectors.
    //
    for (unsigned int i = 0; i < vvF.size(); i++)
      for (unsigned int j = 0; j < vvF.at(i).size(); j++)
	result.push_back(vvF.at(i).at(j));
  }
  else {
    //
    // We want to transpose.
    // Let's consider the very general case where our vector of vectors does not represent
    // a rectangular matrix, i.e. all the elements of this vector are vectors with possibly
    // different size. 
    //
    unsigned int maxsize = 0;
    unsigned int cursize = 0;
    for (unsigned int i = 0; i < vvF.size(); i++) {
      cursize = vvF.at(i).size();
      if (cursize > maxsize) maxsize = cursize;
    }
    
    for (unsigned int i = 0; i < maxsize; i++)
      for (unsigned int j = 0; j < vvF.size(); j++)
	if (i < vvF.at(j).size())
	  result.push_back(vvF.at(j).at(i));     
  }
  
  return result;
}   

class DConverter {
  
public :

  static vector<double> toVectorD(const vector<vector<double> >& vv);
  template<class T>

  static vector<double> toVectorD(const vector<T>& v) {
    vector<double> result;
    for (typename vector<T>::const_iterator iter = v.begin(); iter != v.end(); ++iter)
      result.push_back(iter->get());
    return result;
  }

  template<class T>
  static vector<double> toVectorD(const vector<vector<T> >& vv) {
    vector<double> result;
    for (typename vector<vector<T> >::const_iterator iter = vv.begin(); iter != vv.end(); ++iter)
      for (typename vector<T>::const_iterator iiter = iter->begin(); iiter != iter->end(); ++iiter)
	result.push_back(iiter->get());
    return result;
  }

  template<class T>
  static vector<vector<double> > toMatrixD(const vector<vector<T> >& vv) {
    vector<vector<double> > result;
    vector<double> vD;
    for(vector<T> v: vv) {
      vD.clear();
      for(T x: v) {
	vD.push_back(x.get());
      }
      result.push_back(vD);
    }
    return result;
  }
};

vector<double> DConverter::toVectorD(const vector<vector<double> >& vv) {
  vector<double> result;
  
  for (vector<vector<double> >::const_iterator iter = vv.begin(); iter != vv.end(); ++iter)
    result.insert(result.end(), iter->begin(), iter->end());

  return result;
}

class IConverter {
  
public :
  vector<int> static toVectorI(vector<Tag>& v);
};

vector<int> IConverter::toVectorI(vector<Tag>& v) {
  vector<int> result(v.size());
  vector<Tag>::const_iterator iiter = v.begin();
  vector<int>::iterator oiter = result.begin();
  for (; iiter != v.end(); ++iiter, ++oiter)
    *oiter = iiter->getTagValue();

  return result;
}

class SConverter {
public :
  template<class Enum, class EnumHelper>
  static vector<string> toVectorS(const vector<Enum>& vEnum) {
    vector<string> result(vEnum.size());
    typename vector<Enum>::const_iterator iiter = vEnum.begin();
    vector<string>::iterator     oiter = result.begin();

    for (; iiter != vEnum.begin(); ++iiter, ++oiter)
      *oiter = EnumHelper::name(*iiter);

    return result;
  }
};

class  CConverter {
public :
  static vector<std::complex<float> > toVectorCF(const vector<vector<asdm::Complex> >& vv);
};

vector<std::complex<float> > CConverter::toVectorCF(const vector<vector<asdm::Complex> >& vv) {
  vector<std::complex<float> > result;

  for (vector<vector<asdm::Complex> >::const_iterator iter = vv.begin(); iter != vv.end(); ++iter) 
    result.insert(result.end(), iter->begin(), iter->end());

  return result;
}

class ComplexDataFilter {

public:
  ComplexDataFilter();
  virtual ~ComplexDataFilter();
  virtual float* to4Pol(int numPol, int numChan, float* cdata);

private:
  vector<float *> storage; 
};

ComplexDataFilter::ComplexDataFilter() {
  ;
}

ComplexDataFilter::~ComplexDataFilter() {
  for (unsigned int i = 0; i < storage.size(); i++) 
    delete[] storage.at(i);
}

float *ComplexDataFilter::to4Pol(int numCorr, int numChan, float* cdata) {
  // Do nothing if numCorr != 3
  if (numCorr != 3) return cdata;

  
  // Increase the storage size by appending nul chars.
  float* filtered = new float[ 2 * 4 * numChan ];

  storage.push_back(filtered);
  
  for (int i = 0; i < numChan; i++) {
    // The 1st row goes to the first row.
    filtered[ 8 * i ]     = cdata[ 6 * i ] ;
    filtered[ 8 * i + 1 ] = cdata[ 6 * i + 1 ] ;

    // The second row goes the second row.
    filtered [ 8 * i + 2 ] = cdata[ 6 * i + 2 ] ;
    filtered [ 8 * i + 3 ] = cdata[ 6 * i + 3 ] ;

    // The second row's conjugate goes to the third row.
    filtered [ 8 * i + 4 ] = cdata[ 6 * i + 2 ] ;
    filtered [ 8 * i + 5 ] = - cdata[ 6 * i + 3 ] ;

    // The third row goes to the third row.
    filtered [ 8 * i + 6 ] = cdata[ 6 * i + 4 ] ;
    filtered [ 8 * i + 7 ] = cdata[ 6 * i + 5 ] ;    
  }

  return filtered;
}

//
//  A collection of functions to combine target and offset
//  into a pointing direction in the Pointing table.
//
/**
 *
 * From a Fortran 90 subroutine, written by Didier Despois (1980)
 * and revised by Michel Perault (1984).
 *
 * Computes a coodinates conversion matrix for a rotation defined
 * by Euler angles psi, the and phi (values in radians).
 *
 * The product of the euclidean coordinates of a vector in the original
 * base by this matrix gives the coordinates in the rotated base.
 * 
 */
void eulmat(double psi, double the, double phi, vector<vector<double> >& mat) {
  double cpsi, spsi, cthe, sthe, cphi, sphi;
  double x1, x2;

  cpsi = cos(psi);
  spsi = sin(psi);
  cthe = cos(the);
  sthe = sin(the);
  cphi = cos(phi);
  sphi = sin(phi);

  x1 = spsi*cthe;
  x2 = cpsi * cthe; 
  mat.at(0).at(0) =  cpsi * cphi - x1 * sphi;
  mat.at(0).at(1) =  spsi * cphi + x2 * sphi;
  mat.at(0).at(2) =  sthe * sphi;

  mat.at(1).at(0) = -cpsi * sphi - x1 * cphi;
  mat.at(1).at(1) = -spsi * sphi + x2 * cphi;
  mat.at(1).at(2) = sthe * cphi;

  mat.at(2).at(0) =  spsi * sthe;
  mat.at(2).at(1) =  -cpsi * sthe;
  mat.at(2).at(2) =  cthe;
}

/**
 * Performs a product Matrix x Vector.
 *
 * Attention ! vFactor , mFactor and vResult must be correctly
 * sized ! And also vFactor and vResult must be different objects.
 */
void matvec (const vector<vector<double> >& mFactor, const vector<double>& vFactor, vector<double>& vResult) {
  vResult.at(0) = vFactor.at(0)*mFactor.at(0).at(0) + vFactor.at(1)*mFactor.at(0).at(1) + vFactor.at(2)*mFactor.at(0).at(2);
  vResult.at(1) = vFactor.at(0)*mFactor.at(1).at(0) + vFactor.at(1)*mFactor.at(1).at(1) + vFactor.at(2)*mFactor.at(1).at(2);
  vResult.at(2) = vFactor.at(0)*mFactor.at(2).at(0) + vFactor.at(1)*mFactor.at(2).at(1) + vFactor.at(2)*mFactor.at(2).at(2);
}

/**
 * Spherical to cartesian conversion.
 * Radius assumed to be equal to 1.
 *
 * Attention ! a and x must be correctly sized !
 */
void rect(const vector<double>& s, vector<double>& x) {
  x.at(0) = cos(s.at(0)) * cos(s.at(1));
  x.at(1) = sin(s.at(0)) * cos(s.at(1));
  x.at(2) = sin(s.at(1));
}

/**
 * Cartesian to spherical conversion.
 * Radius assumed to be equal to 1.
 *
 * Attention ! a and x must be correctly sized !
 */
void spher(const vector<double>& x, vector<double>& s) {
  s.at(0) = atan2(x.at(1) , x.at(0));
  s.at(1) = atan2(x.at(2), sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1)));
}


/** 
 * a rotation matrix from local (topocentric) coordinates to ITRF geocentric coordinates
 * lambda - longitude, phi - latitude 
 */
void topo2geomat(double lambda, double phi, vector<vector<double> >& mat) {
 
  double clam, slam, cphi, sphi;

  clam = cos(lambda);
  slam = sin(lambda);
  cphi = cos(phi);
  sphi = sin(phi);

  mat.at(0).at(0) =  -slam;
  mat.at(0).at(1) =  -sphi * clam;
  mat.at(0).at(2) = cphi * clam; 

  mat.at(1).at(0) = clam;
  mat.at(1).at(1) = -sphi * slam;
  mat.at(1).at(2) = cphi * slam;

  mat.at(2).at(0) =  0; 
  mat.at(2).at(1) =  cphi;
  mat.at(2).at(2) =  sphi;
}

/**
 * Reorder the collection of spectral windows ids.
 *
 * Given a dataset 'ds', the statement 'ds.spectralWindowTable.get()' returns
 * a vector of pointers on instances of 'SpectralWindowRow' and 
 * defines implicitely an ordered collection of Tag (of type TagType::SpectralWindow) instances
 * obtained by using 'getSpectralWindowId()' on each instance.
 *
 * This method partitions this ordered collection into two collections
 * based on the criterium : does the Tag appear in at least one row
 * of the DataDescription table in the SpectralWindowId attribute or not.
 *
 * The methods returns a vector of (SpactralWindow) Tag  obtained by appending the collections
 * of (SpectralWindow) Tag which appear in the DataDescription table
 * followed by the collection of (SpectralWindow) Tag which does not appear
 * in the DataDescription table. The int value which associated with a Tag
 * is the order number of the insertion of the pair (Tag, int) in the map.
 * 
 */

struct TagCmp {
  bool operator() (const Tag& t1, const Tag& t2) const {
    return t1.getTagValue() < t2.getTagValue();
  }
};

vector<Tag> reorderSwIds(const ASDM& ds) {
  vector<SpectralWindowRow *> swRs = ds.getSpectralWindow().get();
  vector<DataDescriptionRow *> ddRs = ds.getDataDescription().get();
  map<Tag, bool, TagCmp> isInDD;

  for (vector<SpectralWindowRow *>::size_type  i = 0; i < swRs.size(); i++) isInDD[swRs[i]->getSpectralWindowId()] = false;
  for (vector<DataDescriptionRow *>::size_type i = 0; i < ddRs.size(); i++) isInDD[ddRs[i]->getSpectralWindowId()] = true;

  vector<Tag> swIdsDD, swIdsNoDD;
  for (map<Tag, bool, TagCmp>::iterator iter = isInDD.begin(); iter != isInDD.end(); ++iter)
    if (iter->second) swIdsDD.push_back(iter->first);
    else swIdsNoDD.push_back(iter->first);

  vector<Tag> result (swIdsDD.begin(), swIdsDD.end());
  for (vector<Tag>::size_type i = 0; i < swIdsNoDD.size(); i++)
    result.push_back(swIdsNoDD[i]);

  return result;
}

// mapping SpwId to DDId
vector<int> getDDIdsFromSwId(const ASDM& ds, const int& spwid) {
  vector<int> DDIds;
  vector<DataDescriptionRow *> ddRs = ds.getDataDescription().get();
  //cerr<<"ddRs nrow="<<ddRs.size()<<endl;
  for (vector<DataDescriptionRow *>::size_type i = 0; i < ddRs.size(); i++) {
    int ddSwId = ddRs[i]->getSpectralWindowId().getTagValue();
    if (spwid == ddSwId) {
      DDIds.push_back(i);
    }
  }
  return DDIds;
}

void usage(char* command) {
  cout << "Usage : " << command << " DataSetName" << endl;
}

template<class T>
void checkVectorSize(const string& vectorAttrName,
		     const vector<T>& vectorAttr,
		     const string& sizeAttrName,
		     unsigned int sizeAttr,
		     const string& tableName,
		     unsigned int rowNumber) {
  if (vectorAttr.size() != sizeAttr) {
    errstream.str("");
    errstream << "In the '"
	      << tableName 
	      << " table, at row #"
	      << rowNumber
	      << ", I found '"
	      << vectorAttrName
	      << "' with a size of '"
	      << vectorAttr.size() 
	      << "', I was expecting it to be equal to the size of '"
	      << sizeAttrName
	      <<"' which is '"
	      <<"'"
	      <<sizeAttr
	      <<"'. I can't go further."
	      << endl;
    error(errstream.str());
  }
}

EnumSet<AtmPhaseCorrection> apcLiterals(const ASDM& ds) {
  EnumSet<AtmPhaseCorrection> result;

  vector<MainRow *> mRs = ds.getMain().get();
  
  for (unsigned int i = 0; i < mRs.size(); i++) {
    ConfigDescriptionRow * configDescriptionRow = mRs.at(i)->getConfigDescriptionUsingConfigDescriptionId();
    vector<AtmPhaseCorrection> apc = configDescriptionRow -> getAtmPhaseCorrection();
    for (unsigned int i = 0; i < apc.size(); i++)
      result.set(apc.at(i));
  }
  return result;
}

bool hasCorrectedData(const EnumSet<AtmPhaseCorrection>& es) {
  return es[AP_CORRECTED];
}

bool hasUncorrectedData(const EnumSet<AtmPhaseCorrection>& es) {
  return es[AP_UNCORRECTED];
}
		    
//
// A number of EnumSet to encode the different selection criteria.
//
EnumSet<CorrelationMode>         es_cm;
EnumSet<SpectralResolutionType>  es_srt;
EnumSet<TimeSampling>            es_ts;
Enum<CorrelationMode>            e_query_cm; 
EnumSet<AtmPhaseCorrection>      es_query_apc;    

// An EnumSet to store the different values of AtmPhaseCorrection present
// in the binary data (apc in datastruct).
//
EnumSet<AtmPhaseCorrection>      es_apc;

//
// By default the resulting MS will not contain compressed columns
// unless the 'compress' option has been given.
// 
bool                             withCompression = false;

//
// A function to determine if overTheTop is present in a given row of the Pointing table.
//
bool overTheTopExists(PointingRow* row) { return row->isOverTheTopExists(); }

map<int, int> swIdx2Idx ;                       // A map which associates old and new index of Spectral Windows before/after reordering.

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Some functions defined for the processing of the SysPower table 
// with 'functional programming' techniques.
//
int sysPowerAntennaId(SysPowerRow* row) {
  return row->getAntennaId().getTagValue();
}

int sysPowerSpectralWindowId(const SysPowerRow* row) {
  return swIdx2Idx[row->getSpectralWindowId().getTagValue()];
}

int sysPowerFeedId(const SysPowerRow* row) {
  return row->getFeedId();
}

double sysPowerMidTimeInSeconds(const SysPowerRow* row) {
  //if (isEVLA) {
  //  return row->getTimeInterval().getStartInMJD()*86400 ; 
  //}
  return row->getTimeInterval().getStartInMJD()*86400 + ((double) row->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond / 2.;
}

double sysPowerIntervalInSeconds(const SysPowerRow* row) {
  return ((double) row->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
}

int sysPowerNumReceptor(const SysPowerRow* row) {
  return row->getNumReceptor();
}

struct sysPowerCheckConstantNumReceptor {
private:
  unsigned int numReceptor;

public:
  sysPowerCheckConstantNumReceptor(unsigned int numReceptor) : numReceptor(numReceptor) {}
  bool operator()(SysPowerRow* row) {
    return numReceptor != (unsigned int) row->getNumReceptor();
  }
};

struct sysPowerCheckSwitchedPowerDifference {
private:
  unsigned int	numReceptor;
  bool		exists;

public:
  sysPowerCheckSwitchedPowerDifference(unsigned int numReceptor, bool exists) : numReceptor(numReceptor), exists(exists) {}
  bool operator()(SysPowerRow* row) {
    return (exists != row->isSwitchedPowerDifferenceExists()) || (row->getSwitchedPowerDifference().size() != numReceptor);
  }
};

struct sysPowerCheckSwitchedPowerSum {
private:
  unsigned int	numReceptor;
  bool		exists;

public:
  sysPowerCheckSwitchedPowerSum(unsigned int numReceptor, bool exists) : numReceptor(numReceptor), exists(exists) {}
  bool operator()(SysPowerRow* row) {
    return (exists != row->isSwitchedPowerSumExists()) || (row->getSwitchedPowerSum().size() != numReceptor);
  }
};

struct sysPowerCheckRequantizerGain {
private:
  unsigned int	numReceptor;
  bool		exists;

public:
  sysPowerCheckRequantizerGain(unsigned int numReceptor, bool exists) : numReceptor(numReceptor), exists(exists) {}
  bool operator()(SysPowerRow* row) {
    return (exists != row->isRequantizerGainExists()) || (row->getRequantizerGain().size() != numReceptor);
  }
};

struct sysPowerSwitchedPowerDifference {
private:
  vector<float>::iterator iter;
  
public:
  sysPowerSwitchedPowerDifference(vector<float>::iterator iter): iter(iter) {}
  void operator()(SysPowerRow* row) { 
    vector<float> tmp = row->getSwitchedPowerDifference();
    copy(tmp.begin(), tmp.end(), iter);
    iter += tmp.size();
  }
};

struct sysPowerSwitchedPowerSum {
private:
  vector<float>::iterator iter;
  
public:
  sysPowerSwitchedPowerSum(vector<float>::iterator iter): iter(iter) {}
  void operator()(SysPowerRow* row) { 
    vector<float> tmp = row->getSwitchedPowerSum();
    copy(tmp.begin(), tmp.end(), iter);
    iter += tmp.size();
  }
};


struct sysPowerRequantizerGain {
private:
  vector<float>::iterator iter;

public:
  sysPowerRequantizerGain(vector<float>::iterator iter): iter(iter) {}
  void operator()(SysPowerRow* row) { 
    vector<float> tmp = row->getRequantizerGain();
    copy(tmp.begin(), tmp.end(), iter);
    iter += tmp.size();
  }
};


/**
 * This function tries to calculate the sizes of the successives slices of the BDF 
 * to be processed given a) the overall size of the BDF and b) the approximative maximum
 * size that one wants for one slice.
 * @paramater BDFsize the overall size of the BDF expressed in bytes,
 * @parameter approxSizeInMemory the approximative size that one wants for one slice, expressed in byte.
 *
 * \par Note:
 * It tries to avoid a last slice, corresponding to 
 * the remainder in the euclidean division BDFsize / approxSizeInMemory.
 */
vector<uint64_t> sizeInMemory(uint64_t BDFsize, uint64_t approxSizeInMemory) {
  if (debug) cout << "sizeInMemory: entering" << endl;
  vector<uint64_t> result;
  uint64_t Q = BDFsize / approxSizeInMemory;
  if (Q == 0) { 
    result.push_back(BDFsize);
  }
  else {
    result.resize(Q, approxSizeInMemory);
    unsigned int R = BDFsize % approxSizeInMemory;
    if ( R > (Q * approxSizeInMemory / 5) )  {
      result.push_back(R); 
    }
    else {
      while (R > 0) 
	for (unsigned int i = 0; R >0 && i < result.size(); i++) {
	  result[i]++ ; R--;
	}
    }
  }
  if (debug) cout << "sizeInMemory: exiting" << endl;
  return result;
} 

vector<map<AtmPhaseCorrection,ASDM2MSFiller*> >  msFillers_v;

map<AtmPhaseCorrection, ASDM2MSFiller*> msFillers; // There will be one filler per value of the axis APC.

vector<int>	dataDescriptionIdx2Idx;
int		ddIdx;
map<MainRow*, int>     stateIdx2Idx;

set<int> SwIdUsed;


/**
 * This function creates a table dedicated to the storage of Ephemeris informations and attaches it to an existing measurement set.
 *
 * The format of the table is inherited from the format of the tables created by the service HORIZONS provided by the JPL (see http://ssd.jpl.nasa.gov/?horizons)
 * 
 * By construction the table will be located under the FIELD directory of the MS to which the table will be attached.
 */
Table *  buildAndAttachEphemeris(const string & name, vector<double> observerLocation, MS* attachMS_p) {
  TableDesc tableDesc;

  // The table keywords firstly.
  tableDesc.comment() = "An ephemeris table.";
  tableDesc.rwKeywordSet().define("MJD0", casacore::Double(0.0));
  tableDesc.rwKeywordSet().define("dMJD", casacore::Double(0.0));
  tableDesc.rwKeywordSet().define("NAME", "T.B.D");
  tableDesc.rwKeywordSet().define("GeoLong", casacore::Double(observerLocation[0] / 3.14159265 * 180.0));
  tableDesc.rwKeywordSet().define("GeoLat", casacore::Double(observerLocation[1] / 3.14159265 * 180.0));
  tableDesc.rwKeywordSet().define("GeoDist", casacore::Double(observerLocation[2]));
  
  // Then the fields definitions and keywords.
  ScalarColumnDesc<casacore::Double> mjdColumn("MJD");
  mjdColumn.rwKeywordSet().define("UNIT", "d");
  tableDesc.addColumn(mjdColumn);

  ScalarColumnDesc<casacore::Double> raColumn("RA");
  raColumn.rwKeywordSet().define("UNIT", "deg");
  tableDesc.addColumn(raColumn);

  ScalarColumnDesc<casacore::Double> decColumn("DEC");
  decColumn.rwKeywordSet().define("UNIT", "deg");
  tableDesc.addColumn(decColumn);
  
  ScalarColumnDesc<casacore::Double> rhoColumn("Rho");
  rhoColumn.rwKeywordSet().define("UNIT", "AU");
  tableDesc.addColumn(rhoColumn);

  ScalarColumnDesc<casacore::Double> radVelColumn("RadVel");
  radVelColumn.rwKeywordSet().define("UNIT", "AU/d");
  tableDesc.addColumn(radVelColumn);

  ScalarColumnDesc<casacore::Double> diskLongColumn("diskLong");
  diskLongColumn.rwKeywordSet().define("UNIT", "deg");
  tableDesc.addColumn(diskLongColumn);

  ScalarColumnDesc<casacore::Double> diskLatColumn("diskLat");
  diskLatColumn.rwKeywordSet().define("UNIT", "deg");
  tableDesc.addColumn(diskLatColumn);

  SetupNewTable tableSetup(attachMS_p->tableName() + "/FIELD/" + String(name),
			   tableDesc,
			   Table::New);

  Table* table_p = new Table(tableSetup, TableLock(TableLock::PermanentLockingWait));
  AlwaysAssert(table_p, AipsError);
  attachMS_p->rwKeywordSet().defineTable(name, *table_p);
  table_p->flush();

  return table_p;
}

void solveTridiagonalSystem(unsigned int		n,
			    const vector<double>&	a,
			    const vector<double>&	b,
			    const vector<double>&	c,
			    const vector<double>&	d,
			    vector<double>&             x) {
  LOGENTER("solveTridiagonalSystem");
  vector<double>   cprime(n-1);
  vector<double>   dprime(n);

  cprime[0] = c[0] / b[0];
  for (unsigned int i = 1 ; i < n-1; i++ ) 
    cprime[i] = c[i] / (b[i] - cprime[i-1] * a[i]);
  LOG("cprime computed");
  
  dprime[0] = d[0] / b[0];
  for (unsigned int i = 1; i < n; i++)
    dprime[i] = (d[i] - dprime[i-1] * a[i]) / (b[i] - cprime[i-1] * a[i]);
  LOG("dprime computed");

  x.clear(); x.resize(n);
  x[n-1] = dprime[n-1];
  for (int i = n-2; i >=0; i--) {
    x[i] = dprime[i] - cprime[i] * x[i+1];
  }
  cout << endl;
  LOGEXIT("solveTridiagonalSystem");
}

void linearInterpCoeff(uint32_t                   npoints,
		       const vector<double>&      time_v,
		       const vector<double>&      k_v,
		       vector<vector<double> >&   coeff_vv) {
  LOGENTER("linearInterpCoeff");
  coeff_vv.clear();
  coeff_vv.resize(npoints-1);
  
  for (uint32_t i = 0; i < npoints-1; i++) {
    vector<double> coeff_v (2);
    coeff_v[0] = k_v[i];
    coeff_v[1] = (k_v[i+1] - k_v[i]) / (time_v[i+1] - time_v[i]);
    coeff_vv[i] = coeff_v;
  }
  LOGEXIT("linearInterpCoeff");
}

void cubicSplineCoeff(unsigned int		npoints,
		      const vector<double>&     time_v,
		      const vector<double>&     k_v,
		      const vector<double>&	a_v,
		      const vector<double>&	b_v,
		      const vector<double>&	c_v,
		      const vector<double>&	d_v,
		      vector<vector<double> >&	coeff_vv) {

  LOGENTER("cubicSplineCoeff");
  
  vector<double> m_v(npoints);
  vector<double> x_v;

  coeff_vv.clear();
  coeff_vv.resize(npoints-1);

  solveTridiagonalSystem(npoints-1, a_v, b_v, c_v, d_v, x_v);
  m_v.clear();
  m_v.resize(npoints);
  m_v[0] = 0.0;
  m_v[npoints-1] = 0.0;
  for (unsigned int i = 1; i < npoints - 1; i++)
    m_v[i] = x_v[i-1];

  for (unsigned int i = 0 ; i < npoints - 1; i++) {
    vector<double> coeff_v (4);
    LOG(" delta k " + TO_STRING(k_v[i+1]-k_v[i]) + " delta t " + TO_STRING(time_v[i+1] - time_v[i]));
    LOG(" m_i = " + TO_STRING(m_v[i]) + ", m_i+1 = " + TO_STRING(m_v[i+1]));
    coeff_v[0] = k_v[i];
    coeff_v[1] = (k_v[i+1] - k_v[i]) / (time_v[i+1] - time_v[i]) - (time_v[i+1] - time_v[i]) * (2 * m_v[i] + m_v[i+1]) / 6.0;
    coeff_v[2] = m_v[i] / 2.0;
    coeff_v[3] = (m_v[i+1] - m_v[i]) / 6.0 / (time_v[i+1] - time_v[i]);
    coeff_vv[i] = coeff_v;
  }
  LOGEXIT("cubicSplineCoeff");
  return;
}

double evalPoly (unsigned int		numCoeff,
		 const vector<double>&	coeff,
		 double 		timeOrigin,
		 double 		time) {
  LOGENTER("evalPoly");
  LOG( "numCoeff=" + TO_STRING(numCoeff) + ", size of coeff=" + TO_STRING(coeff.size())) ;
  //
  // Let's use the Horner schema to evaluate the polynomial.
  double result = coeff[numCoeff-1];
  for (int i = numCoeff - 2; i >= 0; i--) 
    result = coeff[i] + result*(time-timeOrigin);

  LOGEXIT("evalPoly");
  return result;
}

typedef struct NumChanSelector {
public:
  NumChanSelector( ASDM * ds_p ) { this->ds_p = ds_p; } 
  ~NumChanSelector() {;}
  int operator()(Tag& dataDescriptionId) { return ds_p->getSpectralWindow().getRowByKey((ds_p->getDataDescription().getRowByKey(dataDescriptionId))->getSpectralWindowId())->getNumChan(); }
private:
  ASDM * ds_p;  
} NumChanSelector;

typedef struct NumCorrSelector {
public:
  NumCorrSelector( ASDM * ds_p ) { this->ds_p = ds_p; }
  ~NumCorrSelector() {;} 
  int operator()(Tag& dataDescriptionId) { return ds_p->getPolarization().getRowByKey((ds_p->getDataDescription().getRowByKey(dataDescriptionId))->getPolOrHoloId())->getNumCorr(); }
private:
  ASDM * ds_p;  
} NumCorrSelector;

template<class T, class S>
vector<T> dd2v(vector<Tag> dataDescriptionId_v, S s) {
  vector<T> result_v;
  std::transform(dataDescriptionId_v.begin(), dataDescriptionId_v.end(), std::back_inserter(result_v), s);
  return result_v;
}

template<class T> 
vector<T> emptyVector() {
  vector<T> result;
  return result;
}

map<Tag, AnyValueMap<string> > msGroup(ASDM * ds_p) {
  LOGENTER("msGroup");
  
  set<ConfigDescriptionRow * > cfgR_s;
  map<Tag, vector<string> > bdfNamesPerCfgId_m;
  vector<string> empty_v;

  const vector<MainRow*> & mR_v = ds_p->getMain().get();
  for (unsigned int i = 0; i < mR_v.size(); i++) {
    Tag cfgId = mR_v[i]->getConfigDescriptionId();
    cfgR_s.insert(ds_p->getConfigDescription().getRowByKey(cfgId));
    if (bdfNamesPerCfgId_m.find(cfgId) == bdfNamesPerCfgId_m.end() )
      bdfNamesPerCfgId_m[cfgId] = empty_v;
    string dataUID = mR_v[i]->getDataUID().getEntityId().toString();
    replace(dataUID.begin(),dataUID.end(),':','_');
    replace(dataUID.begin(),dataUID.end(),'/','_');
    bdfNamesPerCfgId_m[cfgId].push_back(dataUID);    
  }

  map<Tag, AnyValueMap<string> > msGroupPerConfigDesc_m;
  for (ConfigDescriptionRow * cfgR_p: cfgR_s) {
    AnyValueMap<string> msGroup_avm;
    msGroup_avm.setValue("spectralResolution", cfgR_p->getSpectralType());
    msGroup_avm.setValue("processorType", cfgR_p->getProcessorType());
    msGroup_avm.setValue("bdfNames", bdfNamesPerCfgId_m[cfgR_p->getConfigDescriptionId()]);
    msGroup_avm.setValue("antennas", cfgR_p->getAntennaId());
    msGroup_avm.setValue("ddId", cfgR_p->getDataDescriptionId());
    msGroup_avm.setValue("numChan", dd2v<uint32_t>(cfgR_p->getDataDescriptionId(), NumChanSelector(ds_p)));
    msGroup_avm.setValue("numCorr", dd2v<uint32_t>(cfgR_p->getDataDescriptionId(), NumCorrSelector(ds_p)));
    msGroupPerConfigDesc_m[cfgR_p->getConfigDescriptionId()] = msGroup_avm;
  }

  LOGEXIT("msGroup");
  return msGroupPerConfigDesc_m;

}

/**
 * Fill the Antenna table.
 * @parameter ds_p a pointer to the ASDM dataset.
 */
void fillAntenna(ASDM * ds_p) {
  LOGENTER("fillAntenna");
  try { 
    AntennaTable& antennaT = ds_p->getAntenna();
    AntennaRow*   r   = 0;

    int nAntenna = antennaT.size();
    infostream.str("");
    infostream << "The dataset has " << nAntenna << " antenna(s)...";
    info(infostream.str());
    
    for (int i = 0; i < nAntenna; i++) {
      if ((r = antennaT.getRowByKey(Tag(i, TagType::Antenna))) == 0){
	errstream.str("");
	errstream << "Problem while reading the Antenna table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl;
	error(errstream.str());
      }

      // The MS Antenna position is defined as the sum of the ASDM station position and
      // of the ASDM Antenna position after applying to it a coordinate system transformation.
      // Since the ASDM Antenna position is 0,0,0 for now, we only use the ASDM station position.
      // Update - 2012-03-22
      // Now the ASDM Antenna position contains non-zeros  so need to take account
      // for this now as shown below. For EVLA, this is still 0,0,0.       
      vector<Length> position = r->getStationUsingStationId()->getPosition();
      double xStation = position.at(0).get();
      double yStation = position.at(1).get();
      double zStation = position.at(2).get();
      
      // ---- transform antenna position to geocentric coordinates ----
      // Method 1 - assume z axis of antenna position lines up with station vector
      //            and transformation is done by a rotation matrix based on 
      //            geocntric longitude and latitude. Good enough for current 
      //            antenna position measurement accuracy.
      //
      // geocentric longitude  and latitude
      double glat = atan2(zStation,sqrt(xStation*xStation + yStation*yStation)); 
      double glon = atan2(yStation,xStation);

      // get ASDM Antenna position vector
      vector<Length> antPosition = r->getPosition();
      
      vector<double> cartesianAnt1(3, 0.0);
      vector<double> cartesianAnt2(3, 0.0);
      vector<vector<double> > matrixAnt3x3;
      for (unsigned int ii = 0; ii < 3; ii++) {
        matrixAnt3x3.push_back(cartesianAnt1);
      }
      cartesianAnt1[0] = antPosition.at(0).get();
      cartesianAnt1[1] = antPosition.at(1).get();
      cartesianAnt1[2] = antPosition.at(2).get();

      if (cartesianAnt1[0]!=0.0 || cartesianAnt1[1]!=0 || cartesianAnt1[2]!=0.0) {
        topo2geomat(glon,glat,matrixAnt3x3);
        matvec(matrixAnt3x3,cartesianAnt1,cartesianAnt2); 
      }

      //add antenna position: for now use ones obtained by the topo2geomat() 
      double xPosition = position.at(0).get() + cartesianAnt2.at(0);
      double yPosition = position.at(1).get() + cartesianAnt2.at(1);
      double zPosition = position.at(2).get() + cartesianAnt2.at(2);

      vector<Length> offset = r->getOffset();
      double xOffset = offset.at(0).get();
      double yOffset = offset.at(1).get();
      double zOffset = offset.at(2).get();

      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	// the return value here is not used
	/* int antenna_id = */ iter->second->addAntenna(r->getName(),
						  r->getStationUsingStationId()->getName(),
						  xPosition,
						  yPosition,
						  zPosition,
						  xOffset,
						  yOffset,
						  zOffset,
						  (float)r->getDishDiameter().get());	
      }
    }

    int numTrueAntenna = msFillers.begin()->second->ms()->antenna().nrow();
    if (numTrueAntenna) {
      infostream.str("");
      infostream << "converted in " << numTrueAntenna << " antenna(s)  in the measurement set(s)." ;
      info(infostream.str());
    }
  }
  catch (IllegalAccessException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }

  LOGEXIT("fillAntenna");
}

/**
 * This function fills the CalDevice table.
 *
 * @parameter ds_p a pointer to the ASDM dataset.
 * @param ignoreTime a boolean value to indicate if the selected scans are taken into account or if all the table is going to be processed.
 * @param selectedScanRow_v a vector of pointers on ScanRow used to determine which rows of SysPower are going to be processed.
 */
void fillCalDevice(ASDM* ds_p, bool ignoreTime, vector<ScanRow*>& selectedScanRow_v) {
  LOGENTER("fillCalDevice");
  try {
    const CalDeviceTable& calDeviceT = ds_p->getCalDevice();
    infostream.str("");
    infostream << "The dataset has " << calDeviceT.size() << " calDevice(s)...";

    rowsInAScanbyTimeIntervalFunctor<CalDeviceRow> selector(selectedScanRow_v);    
    const vector<CalDeviceRow *>& calDevices = selector(calDeviceT.get(), ignoreTime);
    if (!ignoreTime) 
      infostream << calDevices.size() << " of them in the selected exec blocks / scans ... ";
    
    info(infostream.str());
    
    for (vector<CalDeviceRow*>::const_iterator iter = calDevices.begin(); iter != calDevices.end(); iter++) {
      bool ignoreThisRow = false;
      unsigned int numCalload = 0;
      unsigned int numReceptor = 0;
      
      //
      // Let's make some checks on the attributes.
      errstream.str("");
      infostream.str("");

      //
      // Is numCalload > 0 ?
      if ((numCalload = (*iter)->getNumCalload()) <= 0) { 
	errstream << "In the table CalDevice, the attribute 'numCalload' in row #"
		  << (unsigned int) (iter - calDevices.begin())
		  << " has an invalid value '("
		  << numCalload << "'), a strictly positive value is expected."
		  << endl; 
	ignoreThisRow = true;
      }
      
      //
      // Do we have enough elements in calLoadNames ?
      vector<CalibrationDevice> temp = (*iter)->getCalLoadNames();
      vector<string> calLoadNames;
      if (temp.size() < numCalload) { 
	errstream  << "In the table CalDevice, the size of the attribute 'calLoadNames' in row #"
		   << (unsigned int) (iter - calDevices.begin())
		   << " is too small. It should be greater than or equal to the value of the atttribute 'numCalload' ("
		   << numCalload
		   <<")."
		   << endl;
	ignoreThisRow = true;
      }
      else {
	calLoadNames.resize(temp.size());
	std::transform(temp.begin(), temp.end(), calLoadNames.begin(), stringValue<CalibrationDevice, CCalibrationDevice>);	  
      }

      //
      // Do we have numReceptor ?
      if ((*iter)->isNumReceptorExists()) {
	numReceptor = (*iter)->getNumReceptor();
	if (numReceptor == 0) {
	  errstream << "In the table CalDevice, the value of the attribute 'numReceptor' in row #"
		    << (unsigned int) (iter - calDevices.begin())
		    << " is invalid (" 
		    << numReceptor 
		    << "). It is expected to be strictly positive."
		    << endl;
	  ignoreThisRow = true;
	}
      }

      //
      // Do we have calEff ?
      vector<vector<float> > calEff;
      if ((*iter)->isCalEffExists()) {
	//
	// Do we take it into account ?
	if (numReceptor == 0) {
	  infostream << "In the table CalDevice, the attribute 'calEff' is present in row #"
		     << (unsigned int) (iter - calDevices.begin())
		     << " but it will be ignored due to the fact that the attribute 'numReceptor' is null."
		     << endl;
	}
	else {
	  calEff = (*iter)->getCalEff();
	  
	  //
	  // Let's check the sizes of its content.
	  if (calEff.size() < numReceptor) {
	    errstream << "In the table CalDevice, the size of the attribute 'calEff' in row #"
		      << (unsigned int) (iter - calDevices.begin())
		      << " is too small. It should be greater than or equal to the value of the attribute 'numReceptor' ("
		      << numReceptor
		      <<")."
		      << endl;
	    ignoreThisRow = true;
	  }
	  else {
	    if (find_if(calEff.begin(), calEff.end(), size_lt<float>(numCalload)) != calEff.end()) {
	      errstream << "In the table CalDevice, the attribute 'calEff' in row #"
			<< (unsigned int) (iter - calDevices.begin())
			<< " has at least one element whose size is too small. All its elements should have their size"
			<< " greater then  or equal to the value of the attribute 'numCalload' ("
			<< numCalload
			<< ")."
			<< endl;
	      ignoreThisRow = true;
	    }
	  }
	}	
      }
      
      //
      // In priority let's see if we have coupledNoiseCal ?
      vector<vector<float> > coupledNoiseCal;
      if ((*iter)->isCoupledNoiseCalExists()) {
	//
	// Do we take it into account ?
	if (numReceptor == 0) {
	  infostream << "In the table CalDevice, the attribute 'coupledNoiseCal' is present in row #"
		     << (unsigned int) (iter - calDevices.begin())
		     << " but it will be ignored due to the fact that the attribute 'numReceptor' is null."
		     << endl;
	}
	else {
	  coupledNoiseCal = (*iter)->getCoupledNoiseCal();
	  
	  //
	  // Let's check the sizes of its content.
	  if (coupledNoiseCal.size() < numReceptor) {
	    errstream << "In the table CalDevice, the size of the attribute 'coupledNoiseCal' in row #"
		      << (unsigned int) (iter - calDevices.begin())
		      << " is too small. It should be greater than or equal to the value of the attribute 'numReceptor' ("
		      << numReceptor
		      <<")."
		      << endl;
	    ignoreThisRow = true;
	  }
	  else {
	    if (find_if(coupledNoiseCal.begin(), coupledNoiseCal.end(), size_lt<float>(numCalload)) != coupledNoiseCal.end()) {
	      errstream << "In the table CalDevice, the attribute 'coupledNoiseCal' in row #"
			<< (unsigned int) (iter - calDevices.begin())
			<< " has at least one element whose size is too small. All its elements should have their size"
			<< " greater than or equal to the value of the attribute 'numCalload' (=="
			<< numCalload
			<< ")."
			<< endl;
	      ignoreThisRow = true;
	    }
	  }
	}	
      }
      // Ok we don't have coupledNoiseCal , but maybe we have noiseCal ?
      else if ((*iter)->isNoiseCalExists()) {
	//
	// Do we take it into account ?
	vector<double> noiseCal = (*iter)->getNoiseCal();
	
	if (noiseCal.size() < numCalload) {
	  infostream << "In the table CalDevice, the size of the attribute 'noiseCal' in row #"
		     << (unsigned int) (iter - calDevices.begin())
		     << " is too small. It should be greater than or equal to the value of the attribute 'numCalload' ("
		     << numCalload
		     << ")."
		     << endl;
	  ignoreThisRow = true;
	}
	else {
	  // So yes we have a noiseCal attribute, then pretend we have coupledNoiseCal. 
	  // Artificially force numReceptor to 2 and fill coupledNoiseCal by replicating what we have in noiseCal :
	  // coupledNoiseCal[0] = noiseCal
	  // coupledNoiseCal[1] = noiseCal
	  //
	  // infostream << "In the table CalDevice  there is no attribute 'coupledNoiseCal' but there an attribute 'noiseCal' in row #"
	  // 	       << (unsigned int) (iter - calDevices.begin())
	  // 	       << " which we are going to use to fill the MS NOISE_CAL by replicating its values."
	  // 	       << endl;
	  
	  numReceptor = 2;
	  coupledNoiseCal.resize(numReceptor);
	  for (unsigned int iReceptor = 0; iReceptor < numReceptor; iReceptor++) {
	    coupledNoiseCal[iReceptor].resize(numCalload);
	    std::transform(noiseCal.begin(), noiseCal.begin()+numCalload, coupledNoiseCal[iReceptor].begin(), d2f);
	  } 
	}
      }
      
      //
      // Do we have temperatureLoad ?
      vector<double> temperatureLoad;
      if ((*iter)->isTemperatureLoadExists()) {
	vector<Temperature> temp = (*iter)->getTemperatureLoad();
	if (temp.size() < numCalload) {
	  errstream  << "In the table CalDevice, the size of the attribute 'temperatureLoad' in row #"
		     << (unsigned int) (iter - calDevices.begin())
		     << " is too small. It should be greater than or equal to the value of the atttribute 'numCalload' ("
		     << numCalload
		     <<")."
		     << endl;
	  ignoreThisRow = true;
	}
	else {
	  temperatureLoad.resize(temp.size());
	  std::transform(temp.begin(), temp.end(), temperatureLoad.begin(), basicTypeValue<Temperature, float>);	  
	}
      }
      
      if (errstream.str().size() > 0) 
	error(errstream.str());
      
      if (infostream.str().size() > 0)
	info(infostream.str());

      if (ignoreThisRow) {
	infostream.str("");
	infostream << "This row will be ignored." << endl;
	info(infostream.str());
	continue;
      }

      //
      // And finally we can add a new row to the MS CALDEVICE table.
      double interval = ((double) (*iter)->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  (*iter)->getTimeInterval().getStartInMJD()*86400 ;
      // }
      // else {
      time =  (*iter)->getTimeInterval().getStartInMJD()*86400 + interval / 2.0 ;
      //}
      
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator msIter = msFillers.begin();
	   msIter != msFillers.end();
	   ++msIter) {
	msIter->second->addCalDevice((*iter)->getAntennaId().getTagValue(),
				     (*iter)->getFeedId(),
				     swIdx2Idx[(*iter)->getSpectralWindowId().getTagValue()],
				     time,
				     interval,
				     numCalload,
				     calLoadNames,
				     numReceptor,
				     calEff,
				     coupledNoiseCal,
				     temperatureLoad);
      }      
    }
    unsigned int numMSCalDevices = (const_cast<casacore::MeasurementSet*>(msFillers.begin()->second->ms()))->rwKeywordSet().asTable("CALDEVICE").nrow();
    if (numMSCalDevices > 0) {
      infostream.str("");
      infostream << "converted in " << numMSCalDevices << " caldevice(s) in the measurement set.";
      info(infostream.str());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillCalDevice");
}

/**
 * This function fills the DataDescription table.
 * @parameter ds_p a pointer to the ASDM dataset.
 */
void fillDataDescription(ASDM* ds_p, vector<int>& polarizationIdx2Idx_v) {
  LOGENTER("fillDataDescription");
  try {
    DataDescriptionTable& ddT = ds_p->getDataDescription();
    DataDescriptionRow* r = 0;
    int nDataDescription = ddT.size();
    infostream.str("");
    infostream << "The dataset has " << nDataDescription << " data description(s)...";
    info(infostream.str());
    
    for (int i = 0; i < nDataDescription; i++) {
      if ((r=ddT.getRowByKey(Tag(i, TagType::DataDescription))) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the DataDescription table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	ddIdx = iter->second->addUniqueDataDescription(swIdx2Idx[r->getSpectralWindowId().getTagValue()],
						       polarizationIdx2Idx_v.at(r->getPolOrHoloId().getTagValue()));
      }
      dataDescriptionIdx2Idx.push_back(ddIdx); 
    }
    if (nDataDescription) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->dataDescription().nrow() << " data description(s)  in the measurement set(s)." ;
      info(infostream.str());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillDataDescription");
}

#define LOG_EPHEM(message) if (getenv("FILLER_LOG_EPHEM")) cout << message;

void fillEphemeris(ASDM* ds_p, uint64_t timeStepInNanoSecond) {
  LOGENTER("fillEphemeris");
  
  try {
    // Retrieve the Ephemeris table's content.
    EphemerisTable & eT = ds_p->getEphemeris();
    const vector<EphemerisRow *>&  eR_v = eT.get();
  
    infostream.str("");
    infostream << "The dataset has " << eT.size() << " ephemeris row(s)...";
    info(infostream.str());

    // Let's partition the vector of Ephemeris rows based on the value of the field ephemerisId.
    vector<EphemerisRow *> empty_v;
    map<int, vector<EphemerisRow *> > i2e_m; // A map which associates a value of ephemerisId to the vector of Ephemeris rows 
    // having this value in their field ephemerisId.
    set<int> ephemerisId_s;
    for(EphemerisRow * eR_p: eR_v) {
      int ephemerisId = eR_p -> getEphemerisId();
      if (i2e_m.find(ephemerisId) == i2e_m.end()) {
	ephemerisId_s.insert(ephemerisId);
	i2e_m[ephemerisId] = empty_v;
      }
      i2e_m[ephemerisId].push_back(eR_p);
    }

    // Let's create and fill the MS ephemeris tables.

    for(int ephemerisId: ephemerisId_s) {
      /**
       * Check if there is at least one ASDM::Field row refering to this ephemerisId.
       */
      const vector<FieldRow *> fR_v = ds_p->getField().get();
      vector<FieldRow *> relatedField_v;    // This vector elements will contain pointers to all the Fields refering to this ephemerisId.
      for(const FieldRow* fR_p: fR_v) {
	if (fR_p->isEphemerisIdExists() && (fR_p->getEphemerisId() == ephemerisId))
	  relatedField_v.push_back(const_cast<FieldRow *>(fR_p));
      }
    
      if (relatedField_v.size() == 0) {
	infostream.str("");
	infostream << "No ASDM Field found with 'ephemerisId=='"<< ephemerisId << "'. An ephemeris table will be created in the MS though, but it will not be connected a Field.";
	info(infostream.str());
      }

      // Let's use the value of 'fieldName' in first element (if any) of relatedField_v to build the name of the MS Ephemeris table to
      // be created.
      //
      string fieldName;
      if (relatedField_v.size() > 0)
	fieldName = relatedField_v[0]->getFieldName();

      /**
       * For each MS ephemeris table we need :
       *
       * 1) MJD0
       * 2) DMJD0
       * 3) NAME
       * 4) GeoLong in deg.
       * 5) GeoLat  in deg.
       * 6) GeoDist in km.
       * 7) origin (a textual summary of the origin of this ephemeris)
       *
       * We derive these values from the informations found in the first element of i2e_m[ephemerisId]
       */
      vector<EphemerisRow *>&	v		 = i2e_m[ephemerisId];    
      vector<double>		observerLocation = v[0]->getObserverLocation();
      double			geoLong		 = observerLocation[0] / 3.14159265 * 180.0; // in order to get degrees.
      double			geoLat		 = observerLocation[1] / 3.14159265 * 180.0; // in order to get degrees.
      double                      geoDist          = observerLocation[2] / 1000.0;             // in order to get km (supposedly above the reference ellipsoid)

      int64_t	t0ASDM = v[0]->getTimeInterval().getStart().get();	// The first time recorded for this ephemerisId.
      int64_t	q      = t0ASDM / timeStepInNanoSecond;
      int64_t	r      = t0ASDM % timeStepInNanoSecond;
      int64_t	t0MS   = t0ASDM;
      if ( r != 0 ) {  
	q		       = q + 1;
	t0MS	       = q * timeStepInNanoSecond;
      }

      double mjd0 = ArrayTime(t0MS).getMJD();
      double dmjd = 0.001;

      // Prepare the table keywords with the values computed above.
      TableDesc tableDesc;
    
      tableDesc.comment() = v[0]->getOrigin();
      tableDesc.rwKeywordSet().define("MJD0", casacore::Double(mjd0));
      tableDesc.rwKeywordSet().define("dMJD", casacore::Double(dmjd));
      tableDesc.rwKeywordSet().define("NAME", fieldName);
      tableDesc.rwKeywordSet().define("GeoLong", casacore::Double(geoLong));
      tableDesc.rwKeywordSet().define("GeoLat", casacore::Double(geoLat));
      tableDesc.rwKeywordSet().define("GeoDist", casacore::Double(geoDist));
    
      // Then the fields definitions and keywords.
      ScalarColumnDesc<casacore::Double> mjdColumn("MJD");
      mjdColumn.rwKeywordSet().define("UNIT", "d");
      tableDesc.addColumn(mjdColumn);
    
      ScalarColumnDesc<casacore::Double> raColumn("RA");
      raColumn.rwKeywordSet().define("UNIT", "deg");
      tableDesc.addColumn(raColumn);
    
      ScalarColumnDesc<casacore::Double> decColumn("DEC");
      decColumn.rwKeywordSet().define("UNIT", "deg");
      tableDesc.addColumn(decColumn);
    
      ScalarColumnDesc<casacore::Double> rhoColumn("Rho");
      rhoColumn.rwKeywordSet().define("UNIT", "AU");
      tableDesc.addColumn(rhoColumn);
    
      ScalarColumnDesc<casacore::Double> radVelColumn("RadVel");
      radVelColumn.rwKeywordSet().define("UNIT", "km/s");
      tableDesc.addColumn(radVelColumn);
    
      ScalarColumnDesc<casacore::Double> diskLongColumn("diskLong");
      diskLongColumn.rwKeywordSet().define("UNIT", "deg");
      tableDesc.addColumn(diskLongColumn);
    
      ScalarColumnDesc<casacore::Double> diskLatColumn("diskLat");
      diskLatColumn.rwKeywordSet().define("UNIT", "deg");
      tableDesc.addColumn(diskLatColumn);
    
  
      string tableName = "EPHEM"
	+ TO_STRING(ephemerisId)
	+ "_"
	+ fieldName
	+ "_"
	+ TO_STRING(mjd0)
	+ ".tab";

      map<AtmPhaseCorrection, Table*> apc2EphemTable_m;
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	string tablePath = (string) iter->second->ms()->tableName() + string("/FIELD/") + tableName;  
	SetupNewTable tableSetup(tablePath,
				 tableDesc,
				 Table::New);
    
	Table* table_p = new Table(tableSetup, TableLock(TableLock::PermanentLockingWait));
	AlwaysAssert(table_p, AipsError);
	(const_cast<casacore::MeasurementSet*>(iter->second->ms()))->rwKeywordSet().defineTable(tableName, *table_p);
	table_p->flush();
	apc2EphemTable_m[iter->first] = table_p;
	LOG("Empty ephemeris table '" + tablePath + "' created");
      }

      // 
      // Now it's time to fill the EPHEM table(s).
      //
    
      //
      // Check that for each polynomial column the degree is always null or never null and that the optional fields are always present or always absent.
      // And also verify that there is no "hole" in the time range covered by the sequence of ArrayTime intervals when the degree is == 0.
      //
      //
      bool	numPolyDirIsOne	   = v[0]->getNumPolyDir() == 1;
      bool	numPolyDistIsOne   = v[0]->getNumPolyDist() == 1;
      bool	radVelExists	   = v[0]->isRadVelExists() && v[0]->isNumPolyRadVelExists();
      bool	numPolyRadVelIsOne = radVelExists ? v[0]->getNumPolyRadVel() == 1 : false;

      vector<double> duration_v;  // In seconds
      vector<double> time_v;      //  "      "

      time_v.push_back(1.0e-09*v[0]->getTimeInterval().getStart().get());

      errstream.str("");
      for (unsigned int i = 1; i < v.size(); i++) {
	if (numPolyDirIsOne != (v[i]->getNumPolyDir() == 1)) {
	  errstream << "In the table Ephemeris the value of the field 'numPolyDir' is expected to be whether always equal to 1 or always greater than 1. This rule is violated at line #" << i <<"."; 
	  error(errstream.str());
	}

	if (numPolyDistIsOne != (v[i]->getNumPolyDist() == 1)) {
	  errstream << "In the table Ephemeris the value of the field 'numPolyDist' is expected to be whether always equal to 1 or always greater than 1. This rule is violated at line #" << i <<"."; 
	  error(errstream.str());
	}
      
	if (radVelExists != (v[i]->isRadVelExists() && v[i]->isNumPolyRadVelExists())) {
	  errstream << "In the table Ephemeris the fields 'radVel' and 'numPolyRadVel' are expected to be whether always absent or always present. This rule is violated at line #" << i <<".";
	  error(errstream.str());
	}

	if (radVelExists) {
	  if (numPolyRadVelIsOne != (v[i]->getNumPolyRadVel() == 1)) {
	    errstream << "In the table Ephemeris the value of the field 'numPolyRadVel' is expected to be whether always equal to 1 or always greater than 1. This rule is violated at line #" << i <<"."; 
	    error(errstream.str());
	  }	 
	}

	if (numPolyDirIsOne || numPolyDistIsOne || (radVelExists && numPolyRadVelIsOne)) {
	  int64_t start_i = v[i]->getTimeInterval().getStart().get() ;
	  int64_t start_i_1 = v[i-1]->getTimeInterval().getStart().get();
	  int64_t duration_i_1 = v[i-1]->getTimeInterval().getDuration().get();
	  if (start_i != (start_i_1 + duration_i_1)) {
	    infostream.str("");
	    infostream << "The value of 'timeInterval' at row #" << i-1 << " does not cover the time range up to the start time of the next row. The polynomial will be evaluated despite the presence of this 'hole'";
	    info(infostream.str());
	  }
	  duration_v.push_back(1.0e-09*(start_i - start_i_1));
	  time_v.push_back(1.0e-09*start_i); 
	}
      }
    
      LOG ("numPolyDirIsOne = " + TO_STRING(numPolyDirIsOne));
      LOG ("numPolyDistIsOne = " + TO_STRING(numPolyDistIsOne));
      LOG ("radVelExists = " + TO_STRING(radVelExists));
      LOG ("numPolyRadVelIsOne = " + TO_STRING(numPolyRadVelIsOne));
    

      //
      // The number of tabulated values (i.e. the number of rows in the MS Ephemerides) table depend on the 
      // degrees of each polynomial column. If there is at least one such degree which is equal to 1 then
      // we exclude the last element of the vector v, i.e. the last ArrayTimeInterval, since on this time interval
      // we would miss one end value for the interpolation (so far extrapolation is excluded).
      //
      if (numPolyDirIsOne || numPolyDistIsOne || (radVelExists && numPolyRadVelIsOne)) {
	// Then just "forget" the last element.
	LOG("Erasing the last element of v (size before = '" + TO_STRING(v.size()) + "')");
	v.erase(v.begin() + v.size() - 1);
	LOG("Erasing the last element of v (size after = '" + TO_STRING(v.size()) + "')");

	LOG("Erasing the last element of duration_v (size before = '" + TO_STRING(duration_v.size()) + "')");
	duration_v.erase(duration_v.begin() + duration_v.size() - 1);
	LOG("Erasing the last element of duration_v (size after = '" + TO_STRING(duration_v.size()) + "')");
      }

      // 
      // Determine the timely ordered sequence of indexes in v which will be used to tabulate the ephemeris data to be put into the MS table.
      //
      LOG("Prepare the time ordered sequence of indexes used to tabulate the ephemeris data to be written in the MS table.");
      typedef pair<uint32_t, int64_t> atiIdxMStime_pair;
      vector<atiIdxMStime_pair>  atiIdxMStime_v;

      uint32_t index = 0;  
      int64_t tMS = t0MS;
      atiIdxMStime_v.push_back(atiIdxMStime_pair(index, tMS));
      LOG ("size of atiIdxMStime_v="+TO_STRING(atiIdxMStime_v.size())+", index = "+TO_STRING(index)+", tMS = "+TO_STRING(tMS));
      tMS += timeStepInNanoSecond;

      int64_t  start =  v[index]->getTimeInterval().getStart().get();
      int64_t  end   =  start + v[index]->getTimeInterval().getDuration().get();
      do {
	if (tMS < end) {
	  atiIdxMStime_v.push_back(atiIdxMStime_pair(index, tMS));
	  tMS += timeStepInNanoSecond;
	  LOG ("size of atiIdxMStime_v="+TO_STRING(atiIdxMStime_v.size())+", index = "+TO_STRING(index)+", tMS = "+TO_STRING(tMS));
	}
	else {
	  index++;
	  end   =  v[index]->getTimeInterval().getStart().get() + v[index]->getTimeInterval().getDuration().get();
	}

      } while (index < v.size()-1);
    
      LOG("atiIdxMStime_v has " + TO_STRING(atiIdxMStime_v.size()) + " elements.");

      //
      // Prepare the coefficients which will be used for the tabulation.
 
      LOG("Prepare the coefficients which will be used for the tabulations.");
      vector<vector<double> >  raASDM_vv;
      vector<double>           raASDM_v;
      vector<vector<double> >  decASDM_vv;
      vector<double>           decASDM_v;
      vector<vector<double> >  distanceASDM_vv;
      vector<double>           distanceASDM_v;
      vector<vector<double> >  radVelASDM_vv;
      vector<double>           radVelASDM_v;
      vector<double>           empty_v;
      vector<double>           temp_v;
    
      cout.precision(10);
      for (unsigned int i = 0; i < v.size(); i++) {
	LOG_EPHEM("original " + TO_STRING (ArrayTime(v[i]->getTimeInterval().getStart().get()).getMJD()));
	vector<vector<double> > temp_vv = v[i]->getDir();
	if (numPolyDistIsOne) {
	  raASDM_v.push_back(temp_vv[0][0]/3.14159265*180.0);
	  decASDM_v.push_back(temp_vv[0][1]/3.14159265*180.0);
	  LOG_EPHEM (" " + TO_STRING(raASDM_v.back()) + " " + TO_STRING(decASDM_v.back()));
	}
	else {
	  raASDM_vv.push_back(empty_v);
	  decASDM_vv.push_back(empty_v);
	  for (int j = 0; j < v[i]->getNumPolyDir(); j++) {
	    raASDM_vv.back().push_back(temp_vv[j][0]/3.14159265*180.0);
	    decASDM_vv.back().push_back(temp_vv[j][1]/3.14159265*180.0);
	  }
	}

	temp_v = v[i]->getDistance();      
	if (numPolyDistIsOne) {
	  distanceASDM_v.push_back(temp_v[0] / 1.4959787066e11);           // AU
	  LOG_EPHEM (" " + TO_STRING(distanceASDM_v.back()));
	}
	else {
	  distanceASDM_vv.push_back(empty_v);
	  for (int j = 0; j < v[i]->getNumPolyDist(); j++)
	    distanceASDM_vv.back().push_back(temp_v[j] / 1.4959787066e11); // AU
	}

	if (radVelExists) {
	  temp_v = v[i]->getRadVel();
	  if (numPolyRadVelIsOne) { 
	    radVelASDM_v.push_back(temp_v[0] / 1000.0);                    // km/s
	    LOG_EPHEM(" " + TO_STRING(radVelASDM_v.back()));
	  }
	  else {
	    radVelASDM_vv.push_back(empty_v);
	    for (int j = 0; j < v[i]->getNumPolyRadVel(); j++)
	      radVelASDM_vv.back().push_back(temp_v[j]/1000);              // km/s
	  }	
	}
	LOG_EPHEM("\n");
      }

#if 0
      // Preparing the coefficients of cubic splines.
      vector<double> a_v(v.size()-2);
      vector<double> b_v(v.size()-1);
      vector<double> c_v(v.size()-2);
      vector<double> dRA_v(v.size()-2);
      vector<double> dDEC_v(v.size()-2);
      vector<double> dDist_v(v.size()-2);
      vector<double> dRadV_v(v.size()-2);

      //
      // The calculations below are done only once and will be useful for all the columns
      // requiring the cubic spline interpolation.
      //
      if (numPolyDirIsOne || numPolyDistIsOne || (radVelExists && numPolyRadVelIsOne)) {
	LOG("Prepare the tridiagonal system of equations.");
	LOG("a");
	for (unsigned int  i = 0; i < a_v.size(); i++) 
	  a_v[i] = duration_v[i];       

	LOG("b");
	for (unsigned int i = 0; i < b_v.size(); i++)
	  b_v[i] = 2.0 * (duration_v[i] + duration_v[i+1]);
      
	LOG("c");
	for (unsigned int i = 0; i < c_v.size(); i++) 
	  c_v[i] = duration_v[i+1];

	if (numPolyDirIsOne) {
	  LOG("dRA");
	  for (unsigned int i = 0; i < dRA_v.size(); i++)
	    dRA_v[i] = 6.0 * (raASDM_v[i+1] - raASDM_v[i]) / duration_v[i];

	  LOG("Compute the cubic spline coefficients for RAD");
	  cubicSplineCoeff(v.size(),
			   time_v,
			   raASDM_v,
			   a_v,
			   b_v,
			   c_v,
			   dRA_v,
			   raASDM_vv);
	  LOG("size of raASDM_vv="+TO_STRING(raASDM_vv.size()));
	  for(vector<double> temp_v: raASDM_vv) {
	    LOG("raASDM_v = [");
	    for(double temp: temp_v){
	      LOG(TO_STRING(temp)+" ");
	    }
	    LOG("]");
	  }
	
	  for (unsigned int i = 0; i < dDEC_v.size(); i++)
	    dDEC_v[i] = 6.0 * (decASDM_v[i+1] - decASDM_v[i]) / duration_v[i];

	  LOG("Compute the cubic spline coefficients for DEC");
	  cubicSplineCoeff(v.size(),
			   time_v,
			   decASDM_v,
			   a_v,
			   b_v,
			   c_v,
			   dDEC_v,
			   decASDM_vv);
	  LOG("size of decASDM_vv="+TO_STRING(decASDM_vv.size()));	
	}

	if (numPolyDistIsOne)  {
	  for (unsigned int i = 0; i < dRA_v.size(); i++)
	    dDist_v[i] = 6.0 * (distanceASDM_v[i+1] - distanceASDM_v[i]) / duration_v[i];

	  LOG("Compute the cubic spline coefficients for Dist");
	  cubicSplineCoeff(v.size(),
			   time_v,
			   distanceASDM_v,
			   a_v,
			   b_v,
			   c_v,
			   dDist_v,
			   distanceASDM_vv);
	  LOG("size of distanceASDM_vv="+TO_STRING(distanceASDM_vv.size()));		
	}
      
	if (radVelExists && numPolyRadVelIsOne) {
	  LOG("Compute the cubic spline coefficients for RadVel");
	  for (unsigned int i = 0; i < dRA_v.size(); i++) 
	    dRadV_v[i] = 6.0 * (radVelASDM_v[i+1] - radVelASDM_v[i]) / duration_v[i];

	  cubicSplineCoeff(v.size(),
			   time_v,
			   radVelASDM_v,
			   a_v,
			   b_v,
			   c_v,
			   dRadV_v,
			   radVelASDM_vv);
	  LOG("size of radVelASDM_vv="+TO_STRING(radVelASDM_vv.size()));
	}     
      }
      // End of preparing the coefficients of cubic splines.
#else
      // Preparing the coefficients of piecewise polynomial of degree 1.
      //
      // The calculations below are done only once and will be useful for all the columns
      // requiring the cubic spline interpolation.
      //
      if (numPolyDirIsOne || numPolyDistIsOne || (radVelExists && numPolyRadVelIsOne)) {
	if (numPolyDirIsOne) {
	  LOG("Compute the linear interpolation coefficients for RAD");
	  linearInterpCoeff(v.size(),
			    time_v,
			    raASDM_v,
			    raASDM_vv);
	
	  LOG("Compute the linear interpolation coefficients for DEC");
	  linearInterpCoeff(v.size(),
			    time_v,
			    decASDM_v,
			    decASDM_vv);
	}

	if (numPolyDistIsOne)  {
	  LOG("Compute the linear interolation coefficients for Dist");
	  linearInterpCoeff(v.size(),
			    time_v,
			    distanceASDM_v,
			    distanceASDM_vv);
	}
      
	if (radVelExists && numPolyRadVelIsOne) {
	  LOG("Compute the linear interpolation coefficients for RadVel");
	  linearInterpCoeff(v.size(),
			    time_v,
			    radVelASDM_v,
			    radVelASDM_vv);
	}     
      }
      // End of interpolating with piecewise polynomial of degree 1.
#endif
    
      vector<double> mjdMS_v;
      vector<double> raMS_v;
      vector<double> decMS_v;
      vector<double> distanceMS_v;
      vector<double> radVelMS_v;

      for (atiIdxMStime_pair atiIdxMStime: atiIdxMStime_v) {
	//
	// MJD
	mjdMS_v.push_back(ArrayTime(atiIdxMStime.second).getMJD());
	LOG_EPHEM( "resampled " + TO_STRING(mjdMS_v.back()));
	LOG("mjdMS_v -> "+TO_STRING(mjdMS_v.back()));
      
	double timeOrigin = 1.0e-09 * v[atiIdxMStime.first]->getTimeOrigin().get(); 
	double time       = 1.0e-09 * atiIdxMStime.second;
      
	LOG("timeOrigin="+TO_STRING(timeOrigin)+", time="+TO_STRING(time));

	//
	// RA / DEC
	LOG("Eval poly for RA");
	LOG("atiIdxMStime.first = " + TO_STRING(atiIdxMStime.first));
	raMS_v.push_back(evalPoly(raASDM_vv[atiIdxMStime.first].size(),
				  raASDM_vv[atiIdxMStime.first],
				  timeOrigin,
				  time));
	LOG_EPHEM(" " + TO_STRING(raMS_v.back()));
	LOG("raMS_v -> "+TO_STRING(raMS_v.back()));

	LOG("Eval poly for DEC");
	decMS_v.push_back(evalPoly(decASDM_vv[atiIdxMStime.first].size(),
				   decASDM_vv[atiIdxMStime.first],
				   timeOrigin,
				   time));
	LOG_EPHEM(" " + TO_STRING(decMS_v.back()));
      
	//
	// Distance
	LOG("Eval poly for distance");
	distanceMS_v.push_back(evalPoly(distanceASDM_vv[atiIdxMStime.first].size(),
					distanceASDM_vv[atiIdxMStime.first],
					timeOrigin,
					time));
	LOG_EPHEM(" " + TO_STRING(distanceMS_v.back()));
	//
	// Radvel
	if (radVelExists) { 
	  LOG("Eval poly for radvel");
	  radVelMS_v.push_back(evalPoly(radVelASDM_vv[atiIdxMStime.first].size(),
					radVelASDM_vv[atiIdxMStime.first],
					timeOrigin,
					time));
	  LOG_EPHEM(" " + TO_STRING(radVelMS_v.back()));
	}
	LOG_EPHEM("\n");
      }

      // Now the data are ready to be written to the MS Ephemeris table.
      // Let's proceed, using Slicers.
      //
    
      unsigned int numRows = raMS_v.size();
      Slicer slicer(IPosition(1, 0),
		    IPosition(1, numRows-1),
		    Slicer::endIsLast);

      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	Table * table_p = apc2EphemTable_m[iter->first];
	table_p->addRow(numRows);
	LOG ("Added "+TO_STRING(numRows)+" rows to table "+((string)table_p->tableName()));
      
	LOG("Filling column MJD");
	Vector<casacore::Double> MJD_V(IPosition(1, numRows), &mjdMS_v[0], SHARE);
	ScalarColumn<casacore::Double> MJD(*table_p, "MJD");
	MJD.putColumnRange(slicer, MJD_V);
      
	LOG("Filling column RA");
	Vector<casacore::Double> RA_V(IPosition(1, numRows), &raMS_v[0], SHARE);
	ScalarColumn<casacore::Double> RA(*table_p,  "RA");
	RA.putColumnRange(slicer, RA_V);
      
	LOG("Filling column DEC");
	Vector<casacore::Double> DEC_V(IPosition(1, numRows), &decMS_v[0], SHARE);
	ScalarColumn<casacore::Double> DEC(*table_p, "DEC");
	DEC.putColumnRange(slicer, DEC_V);
      
	LOG ("Filling column Rho");
	Vector<casacore::Double> Rho_V(IPosition(1, numRows), &distanceMS_v[0], SHARE);
	ScalarColumn<casacore::Double> Rho(*table_p, "Rho");
	Rho.putColumnRange(slicer, Rho_V);
      
	if (radVelExists) {
	  LOG ("Filling column RadVel");
	  Vector<casacore::Double> RadVel_V(IPosition(1, numRows), &radVelMS_v[0], SHARE);
	  ScalarColumn<casacore::Double> RadVel(*table_p, "RadVel");
	  RadVel.putColumnRange(slicer, RadVel_V);
	}
      
	table_p->flush();
	infostream.str("");
	infostream << "converted in " << table_p->nrow() << " ephemeris rows in the table '" << table_p->tableName() << "'.";
	info(infostream.str());
      }
    }
  }
  catch (IllegalAccessException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillEphemeris");
  
}

/** 
 * This function fills the MS Feed table.
 * @parameter ds_p a pointer to the ASDM dataset.
 * @ignoreTime a bool to indicate if a time range must be considered in the selection of the rows to process.
 * @selectedScanRow_v a vector of pointers to the ScanRow s  considered.
 */
void fillFeed(ASDM* ds_p, bool ignoreTime, vector<ScanRow *>& selectedScanRow_v) {
  LOGENTER("fillFeed");
  // Process the Feed table
  // Issues :
  //    - time (epoch) : at the moment it takes directly the time as it is stored in the ASDM.
  //    - focusLength (in AIPS++) is no defined.
  try {
    const FeedTable& feedT = ds_p->getFeed();
    FeedRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << feedT.size() << " feed(s)...";
    rowsInAScanbyTimeIntervalFunctor<FeedRow> selector(selectedScanRow_v);
    
    const vector<FeedRow *>& v = selector(feedT.get(), ignoreTime);
    if (!ignoreTime)
      infostream << v.size() << " of them in the exec blocks / selected scans ... ";
    
    info(infostream.str());
    int nFeed = v.size();
    for (int i = 0; i < nFeed; i++) {
      r = v.at(i);
      // For now we just adapt the types of the time related informations and compute a mid-time.
      //
      double interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond;
      // }
      // else {
      time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval/2.0;
      //}

      vector<double> beam_offset_ =  DConverter::toVectorD(r->getBeamOffset());
      vector<std::string> polarization_type_ = PolTypeMapper().toStringVector(r->getPolarizationTypes());
      vector<complex<float> > polarization_response_ = CConverter::toVectorCF(r->getPolResponse());
      vector<double> xyzPosition (3, 0.0);
      if (r->isPositionExists()) {
	vector<Length> position = r->getPosition();
	if (position.size() != 3) {
	  errstream.str("");
	  errstream << "The size of attribute position ('" 
		    << position.size()
		    << "') is not equal to 3. Can't go further."
		    << endl;
	  error(errstream.str());
	}
	
	xyzPosition = DConverter::toVectorD<Length>(position);
      }
      vector<double> receptor_angle_ = DConverter::toVectorD<Angle>(r->getReceptorAngle());
      
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addFeed((int) r->getAntennaId().getTagValue(),
			      r->getFeedId(),
			      swIdx2Idx[r->getSpectralWindowId().getTagValue()],
			      time,
			      interval,
			      r->getNumReceptor(), 
			      -1,             // We ignore the beamId array
			      beam_offset_,
			      polarization_type_,
			      polarization_response_,
			      xyzPosition,
			      receptor_angle_);
      } 
    }
    if (nFeed) {
      infostream.str("");
      infostream <<  "converted in " << msFillers.begin()->second->ms()->feed().nrow() << " feed(s) in the measurement set." ;
      info(infostream.str());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillFeed");
}

/** 
 * This function fills the MS Field table.
 * given :
 * @parameter ds_p a pointer to the ASDM dataset.
 * @parameter considerEphemeris take into account the reference to Ephemeris table(s).
 */
void fillField(ASDM* ds_p, bool considerEphemeris) {
  LOGENTER("fillField");
  vector<pair<int, int> > idxEphemerisId_v;

  try {
    FieldTable& fieldT = ds_p->getField();
    FieldRow* r = 0;
    int nField = fieldT.size();
    infostream.str("");
    infostream << "The dataset has " << nField << " field(s)...";
    info(infostream.str());
    
    for (int i = 0; i < nField; i++) {
      if ((r=fieldT.getRowByKey(Tag(i, TagType::Field))) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the Field table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
      
      string fieldName = r->getFieldName();
      string code = r->getCode();
      DirectionReferenceCodeMod::DirectionReferenceCode dirRefCode = DirectionReferenceCodeMod::J2000;
      if(r->isDirectionCodeExists()){
	dirRefCode = r->getDirectionCode();
	//cout << "found directionCode for field " << fieldName << ": ";
      }
      //       else{
      // 	cout << "No directionCode in input table. Assuming ";
      //       }
      string directionCode = CDirectionReferenceCode::name(dirRefCode);
      //cout << directionCode << endl;
      
      vector<vector<double> > delayDir     = DConverter::toMatrixD<Angle>(r->getDelayDir());
      vector<vector<double> > phaseDir     = DConverter::toMatrixD<Angle>(r->getPhaseDir());
      vector<vector<double> > referenceDir = DConverter::toMatrixD<Angle>(r->getReferenceDir());
      
      if (r->isEphemerisIdExists()) idxEphemerisId_v.push_back(pair<int, int>(i, r->getEphemerisId()));

      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addField( fieldName,
				code,
				r->isTimeExists() ? ((double) r->getTime().get()) / ArrayTime::unitsInASecond : 0.0,
				r->getNumPoly(),
				delayDir,
				phaseDir,
				referenceDir,
				directionCode,
				r->isSourceIdExists()?r->getSourceId():0 );
      }
    }
    
    if (considerEphemeris && (idxEphemerisId_v.size() > 0)) 
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->updateEphemerisIdInField(idxEphemerisId_v);
      }

    if (nField) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->field().nrow() << "  field(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  }
  catch (IllegalAccessException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }

  LOGEXIT("fillField");
}

/**
 * This function fills the FlagCmd table.
 * @parameter ds_p a pointer to the ASDM dataset.
 * @ignoreTime a bool to indicate if a time range must be considered in the selection of the rows to process.
 * @selectedScanRow_v a vector of pointers to the ScanRow s  considered.
 */
void fillFlagCmd(ASDM* ds_p, bool ignoreTime, vector<ScanRow *>& selectedScanRow_v) {
  LOGENTER("fillFlagCmd");
  try {
    const FlagCmdTable& flagCmdT  = ds_p->getFlagCmd();
    FlagCmdRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << flagCmdT.size() << " FlagCmd(s)...";
    rowsInAScanbyTimeIntervalFunctor<FlagCmdRow> selector(selectedScanRow_v);

    const vector<FlagCmdRow *>& v = selector(flagCmdT.get(), ignoreTime);
    if (!ignoreTime)
      infostream << v.size() << " of them in the exec blocks / selected scans ... ";

    info(infostream.str());
    int nFlagCmd = v.size();
    for (int i = 0; i < nFlagCmd; i++) {
      r = v.at(i);
      // For now we just adapt the types of the time related informations and compute a mid-time.
      //
      double interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond ;
      // }
      // else {
      time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval/2.0;
      //}
      string type = r->getType();
      string reason = r->getReason();
      string command = r->getCommand();
   
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addFlagCmd(time,
				 interval,
				 type,
				 reason,
				 r->getLevel(),
				 r->getSeverity(),
				 r->getApplied() ? 1 : 0,
				 command);
      }
    }
    if (nFlagCmd) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->flagCmd().nrow() << " in the measurement set." ;
      info(infostream.str());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillFlagCmd");
}

/**
 * This function fills the History table.
 * @parameter ds_p a pointer to the ASDM dataset.
 * @ignoreTime a bool to indicate if a time range must be considered in the selection of the rows to process.
 * @selectedScanRow_v a vector of pointers to the ScanRow s  considered.
 */
void fillHistory(ASDM* ds_p, bool ignoreTime, vector<ScanRow *>& selectedScanRow_v) {
  LOGENTER("fillHistory");
  // Issues :
  // - use executeBlockId for observationId ...to be discussed with Francois.
  // - objectId : not taken into account (it's a string while the MS expects an int).
  try {
    const HistoryTable& historyT = ds_p->getHistory();
    HistoryRow* r = 0;
    int nHistory = historyT.size();
    infostream.str("");
    infostream << "The dataset has " << nHistory << " history(s)...";
    rowsInAScanbyTimeFunctor<HistoryRow> selector(selectedScanRow_v);

    const vector<HistoryRow *>& v = selector(historyT.get(), ignoreTime);;
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected exec blocks / scans ... ";

    info(infostream.str()); 

    for (int i = 0; i < nHistory; i++) {
      r = v.at(i);
      double time =  ((double) r->getTime().get()) / ArrayTime::unitsInASecond ;
      string message     = r->getMessage();
      string priority    = r->getPriority();
      string origin      = r->getOrigin();
      string application = r->getApplication();
      string cliCommand  = r->getCliCommand();
      string appParams   = r->getAppParms();
 
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addHistory(time,
				 r->getExecBlockId().getTagValue(),   
				 message,
				 priority,
				 origin,
				 -1,
				 application,
				 cliCommand,
				 appParams);
      }
    }
    if (nHistory) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->history().nrow() << " history(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillHistory");
}

/** This function fills the MS Observation table.
 * @parameter ds_p a pointer to the ASDM dataset.
 * @parameter selected_eb_scan_m a map which lists the scans to consider for each exec block considered.
 */
void fillObservation(ASDM* ds_p, map<int, set<int> >& selected_eb_scan_m) {
  LOGENTER("fillObservation");
  try {
    const ExecBlockTable& execBlockT = ds_p->getExecBlock(); 
    ExecBlockRow* r = 0;
    int nExecBlock = execBlockT.size();
    infostream.str("");
    infostream << "The dataset has " << nExecBlock << " execBlock(s) ...";

    vector<ExecBlockRow *> temp_v = execBlockT.get();
    vector<ExecBlockRow *> v;
    for (vector<ExecBlockRow *>::iterator iter_v = temp_v.begin(); iter_v != temp_v.end(); iter_v++)
      if ( selected_eb_scan_m.find((*iter_v)->getExecBlockId().getTagValue()) != selected_eb_scan_m.end() )
	v.push_back(*iter_v);
    
    vector<string> schedule; schedule.resize(2);

    infostream << v.size() << " of them in the selected exec blocks / scans ... ";
    info(infostream.str());

    for (unsigned int i = 0; i < v.size(); i++) {
      r = v.at(i);
      
      string telescopeName  = r->getTelescopeName();
      double startTime      = r->getStartTime().getMJD()*86400;
      double endTime        = r->getEndTime().getMJD()*86400;
      string observerName   = r->getObserverName();

      vector<string> observingLog = r->getObservingLog();

      string scheduleType("ALMA");
      schedule[0] = "SchedulingBlock " + ds_p->getSBSummary().getRowByKey(r->getSBSummaryId())->getSbSummaryUID().getEntityId().toString();
      schedule[1] = "ExecBlock " + r->getExecBlockUID().getEntityId().toString();
      string project(r->getProjectUID().getEntityId().toString());
      double releaseDate = r->isReleaseDateExists() ? r->getReleaseDate().getMJD():0.0;

      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addObservation(telescopeName,
				     startTime,
				     endTime,
				     observerName,
				     observingLog,
				     scheduleType,
				     schedule,
				     project,
				     releaseDate
				     );
      }
      if (i==0) { // assume same telescope for all execBlocks 
        if (telescopeName.find("EVLA")!=string::npos) {
          isEVLA=true;
        }
        else if (telescopeName.find("OSF")!=string::npos || 
                 telescopeName.find("AOS")!=string::npos || 
                 telescopeName.find("ATF")!=string::npos) {
          isEVLA=false;
        }     
        string telname = (isEVLA ? "EVLA" : "ALMA");
        infostream.str("");
        infostream << "Telescope Name:" <<telescopeName << ", process as "<<telname<<" data." ; 
        info(infostream.str());
      }
    }
    if (nExecBlock) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->observation().nrow() << " observation(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillObservation");
}

/** This function fills the MS Pointing table.
 * @parameter ds_p a pointer to the ASDM dataset.
 * @parameter withPointingCorrection bool. If true let's add encoder-pointingDirection to the value calculated for direction.
 * @parameter ignoreTime a bool to indicate if a time range must be considered in the selection of the rows to process.
 * @parameter selectedScanRow_v a vector of pointers to the ScanRow s  considered.
 */
void fillPointing(ASDM* ds_p, bool withPointingCorrection, bool ignoreTime, vector<ScanRow *>& selectedScanRow_v) {
  LOGENTER("fillPointing");
  try {
    const PointingTable& pointingT = ds_p->getPointing();
    infostream.str("");
    infostream << "The dataset has " << pointingT.size() << " pointing(s)...";
    rowsInAScanbyTimeIntervalFunctor<PointingRow> selector(selectedScanRow_v);
    
    const vector<PointingRow *>& v = selector(pointingT.get(), ignoreTime);
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected exec blocks / scans ... ";
    
    info(infostream.str());
    int nPointing = v.size();
    
    if (nPointing > 0) {
      
      // Check some assertions.
      //
      // All rows of ASDM-Pointing must have their attribute usePolynomials equal to false
      // and their numTerm attribute equal to 1. Use the opportunity of this check
      // to compute the number of rows to be created in the MS-Pointing by summing
      // all the numSample attributes values.
      //
      int numMSPointingRows = 0;
      for (unsigned int i = 0; i < v.size(); i++) {
	if (v[i]->getUsePolynomials()) {
	  errstream.str("");
	  errstream << "Found usePolynomials equal to true at row #" << i <<". Can't go further.";
	  error(errstream.str());
	}
	
	numMSPointingRows += v[i]->getNumSample();
      }
      
      //
      // Ok now we have verified the assertions and we know the number of rows
      // to be created into the MS-Pointing, we can proceed.
      
      PointingRow* r = 0;
      
      vector<int>	antenna_id_(numMSPointingRows, 0);
      vector<double>	time_(numMSPointingRows, 0.0);
      vector<double>	interval_(numMSPointingRows, 0.0);
      vector<double>	direction_(2 * numMSPointingRows, 0.0);
      vector<double>	target_(2 * numMSPointingRows, 0.0);
      vector<double>	pointing_offset_(2 * numMSPointingRows, 0.0);
      vector<double>	encoder_(2 * numMSPointingRows, 0.0);
      vector<bool>	tracking_(numMSPointingRows, false);
      
      //
      // Let's check if the optional attribute overTheTop is present somewhere in the table.
      //
      unsigned int numOverTheTop = count_if(v.begin(), v.end(), overTheTopExists);
      bool overTheTopExists4All = v.size() == numOverTheTop;
      
      vector<bool> v_overTheTop_ ;
      
      vector<s_overTheTop> v_s_overTheTop_;
      
      if (overTheTopExists4All) 
	v_overTheTop_.resize(numMSPointingRows);
      else if (numOverTheTop > 0) 
	v_overTheTop_.resize(numOverTheTop);
      
      int iMSPointingRow = 0;
      for (int i = 0; i < nPointing; i++) {     // Each row in the ASDM-Pointing ...
	r = v.at(i);
	
	// Let's prepare some values.
	int antennaId = r->getAntennaId().getTagValue();
	
	double time = 0.0, interval = 0.0;
	if (!r->isSampledTimeIntervalExists()) { // If no sampledTimeInterval then
	  // then compute the first value of MS TIME and INTERVAL.
	  interval   = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond / r->getNumSample();
	  // if (isEVLA) {
	  //   time = ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond;
	  // }
	  // else {
	  time = ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval / 2.0;
	  //}
	}
	
	//
	// The size of each vector below 
	// should be checked against numSample !!!
	//
	int numSample = r->getNumSample();
	const vector<vector<Angle> > encoder = r->getEncoder();
	checkVectorSize<vector<Angle> >("encoder", encoder, "numSample", (unsigned int) numSample, "Pointing", (unsigned int)i);
	
	const vector<vector<Angle> > pointingDirection = r->getPointingDirection();
	checkVectorSize<vector<Angle> >("pointingDirection", pointingDirection, "numSample", (unsigned int) numSample, "Pointing", (unsigned int) i);
	
	const vector<vector<Angle> > target = r->getTarget();
	checkVectorSize<vector<Angle> >("target", target, "numSample", (unsigned int) numSample, "Pointing", (unsigned int) i);
	
	const vector<vector<Angle> > offset = r->getOffset();
	checkVectorSize<vector<Angle> >("offset", offset, "numSample", (unsigned int) numSample, "Pointing", (unsigned int) i);
	
	bool   pointingTracking = r->getPointingTracking();
	
	//
	// Prepare some data structures and values required to compute the
	// (MS) direction.
	vector<double> cartesian1(3, 0.0);
	vector<double> cartesian2(3, 0.0);
	vector<double> spherical1(2, 0.0);
	vector<double> spherical2(2, 0.0);
	vector<vector<double> > matrix3x3;
	for (unsigned int ii = 0; ii < 3; ii++) {
	  matrix3x3.push_back(cartesian1); // cartesian1 is used here just as a way to get a 1D vector of size 3.
	}
	double PSI = M_PI_2;
	double THETA;
	double PHI;
	
	vector<ArrayTimeInterval> timeInterval ;
	if (r->isSampledTimeIntervalExists()) timeInterval = r->getSampledTimeInterval();
	
	// Use 'fill' from algorithm for the cases where values remain constant.
	// ANTENNA_ID
	fill(antenna_id_.begin()+iMSPointingRow, antenna_id_.begin()+iMSPointingRow+numSample, antennaId);
	
	// TRACKING 
	fill(tracking_.begin()+iMSPointingRow, tracking_.begin()+iMSPointingRow+numSample, pointingTracking);
	
	// OVER_THE_TOP 
	if (overTheTopExists4All)
	  // it's present everywhere
	  fill(v_overTheTop_.begin()+iMSPointingRow, v_overTheTop_.begin()+iMSPointingRow+numSample,
	       r->getOverTheTop());
	else if (r->isOverTheTopExists()) {
	  // it's present only in some rows.
	  s_overTheTop saux ;
	  saux.start = iMSPointingRow; saux.len = numSample; saux.value = r->getOverTheTop();
	  v_s_overTheTop_.push_back(saux);
	}
	
	// Use an explicit loop for the other values.
	for (int j = 0 ; j < numSample; j++) { // ... must be expanded in numSample MS-Pointing rows.
	  
	  // TIME and INTERVAL
	  if (r->isSampledTimeIntervalExists()) { //if sampledTimeInterval is present use its values.	           
	    // Here the size of timeInterval will have to be checked against numSample !!
	    interval_[iMSPointingRow] = ((double) timeInterval.at(j).getDuration().get()) / ArrayTime::unitsInASecond ;
	    time_[iMSPointingRow] = ((double) timeInterval.at(j).getStart().get()) / ArrayTime::unitsInASecond
	      + interval_[iMSPointingRow]/2;	  
	  }
	  else {                                     // otherwise compute TIMEs and INTERVALs from the first values.
	    interval_[iMSPointingRow]            = interval;
	    time_[iMSPointingRow]                = time + j*interval;
	  }
	  
	  // DIRECTION
	  THETA = target.at(j).at(1).get();
	  PHI   = -M_PI_2 - target.at(j).at(0).get();
	  spherical1[0] = offset.at(j).at(0).get();
	  spherical1[1] = offset.at(j).at(1).get();
	  rect(spherical1, cartesian1);
	  eulmat(PSI, THETA, PHI, matrix3x3);
	  matvec(matrix3x3, cartesian1, cartesian2);
	  spher(cartesian2, spherical2);
	  direction_[2*iMSPointingRow]  = spherical2[0];
	  direction_[2*iMSPointingRow+1]= spherical2[1];
	  if (withPointingCorrection) { // Cf CSV-2878 and ICT-1532
	    direction_[2*iMSPointingRow]   += encoder.at(j).at(0).get() - pointingDirection.at(j).at(0).get();
	    direction_[2*iMSPointingRow+1] += encoder.at(j).at(1).get() - pointingDirection.at(j).at(1).get() ;
	  }

	  // TARGET
	  target_[2*iMSPointingRow]     = target.at(j).at(0).get();
	  target_[2*iMSPointingRow+1]   = target.at(j).at(1).get();
	  
	  // POINTING_OFFSET
	  pointing_offset_[2*iMSPointingRow]   = offset.at(j).at(0).get();
	  pointing_offset_[2*iMSPointingRow+1] = offset.at(j).at(1).get();
	  
	  // ENCODER
	  encoder_[2*iMSPointingRow]           = encoder.at(j).at(0).get();
	  encoder_[2*iMSPointingRow+1]         = encoder.at(j).at(1).get();
	  
	  
	  // increment the row number in MS Pointing.
	  iMSPointingRow++;	
	}
      }
      
      
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addPointingSlice(numMSPointingRows,
				       antenna_id_,
				       time_,
				       interval_,
				       direction_,
				       target_,
				       pointing_offset_,
				       encoder_,
				       tracking_,
				       overTheTopExists4All,
				       v_overTheTop_,
				       v_s_overTheTop_);
      }
      
      if (nPointing) {
	infostream.str("");
	infostream << "converted in " << msFillers.begin()->second->ms()->pointing().nrow() << " pointing(s) in the measurement set." ;
	info(infostream.str()); 
      }
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillPointing");
}

/** This function fills the MS Polarization table.
 * @parameter ds_p a pointer to the ASDM dataset.
 * @parameter polarizationIdx2Idx_v. On exit a vector of int containing a mapping from ASDM polarizations
 * indexes to MS ones
 */
void fillPolarization(ASDM* ds_p, vector<int>& polarizationIdx2Idx_v) {
  LOGENTER("fillPolarization");

  Stokes::StokesTypes linearCorr[] = { Stokes::XX, Stokes::XY, Stokes::YX, Stokes::YY };
  Stokes::StokesTypes circularCorr[] = { Stokes::RR, Stokes::RL, Stokes::LR, Stokes::LL };
  int corrProduct1[] = { 0, 0 };
  int corrProduct2[] = { 0, 0, 1, 1};
  int corrProduct4[] = { 0, 0, 0, 1, 1, 0, 1, 1 };
			 
  int pIdx;
  polarizationIdx2Idx_v.clear();
  try {
    PolarizationTable& polT = ds_p->getPolarization();  
    PolarizationRow* r = 0;
    int nPolarization = polT.size();
    infostream.str("");
    infostream << "The dataset has " << nPolarization << " polarization(s)..."; 
    info(infostream.str());

    for (int i = 0; i < nPolarization; i++) {
      if ((r=polT.getRowByKey(Tag(i, TagType::Polarization))) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the Polarization table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
      
      int numCorr = r->getNumCorr();
      if (numCorr < 1 || numCorr > 4) {
	ostringstream oss ;
	oss << "a polarization row cannot be processed due to  'numCorr = " << numCorr << "'.";
	throw ASDM2MSException(oss.str());
      }

      //Stokes::StokesTypes * corrType;
      //vector<Stokes::StokesTypes> corrType;
      vector<int> corrType;
      StokesMapper stokesMapper;
      if (numCorr != 3) {
	corrType = StokesMapper::toVectorI(r->getCorrType());
      }
      else {
	numCorr  = 4;
	StokesParameterMod::StokesParameter sp = r->getCorrType()[0];
	if ((sp == StokesParameterMod::RR) ||
	    (sp == StokesParameterMod::LL) ||
	    (sp == StokesParameterMod::RL) ||
	    (sp == StokesParameterMod::LR)) {
	  corrType.resize(4);
	  copy (circularCorr, circularCorr+4, corrType.begin());
	}
	else if ((sp == StokesParameterMod::XX) ||
		 (sp == StokesParameterMod::XY) ||
		 (sp == StokesParameterMod::YX) ||
		 (sp == StokesParameterMod::YY)) {
	  corrType.resize(4);
	  copy (linearCorr, linearCorr+4, corrType.begin());
	}
	else {
	  errstream.str("");
	  errstream << " I don't know what to do with the given Stokes parameters for autocorrelation data" << endl;
	  error(errstream.str());
	}
	  
      }
      
      
      /*int* corrProduct = 0;*/
      vector<int> corrProduct; 
      switch (numCorr) {
      case 1: corrProduct.resize(2); copy(corrProduct1, corrProduct1+2, corrProduct.begin()); break;
      case 2: corrProduct.resize(4); copy(corrProduct2, corrProduct2+4, corrProduct.begin()); break;
      case 4: corrProduct.resize(8); copy(corrProduct4, corrProduct4+8, corrProduct.begin()); break;
      }


      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	pIdx = iter->second->addUniquePolarization(numCorr,
						   corrType,
						   corrProduct
						   );
      }
      polarizationIdx2Idx_v.push_back(pIdx);
    }
    if (nPolarization) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->polarization().nrow() << " polarization(s)." ;
      info(infostream.str());
    }
  }
  catch (ASDM2MSException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillPolarization");
}

/** 
 * This function fills the MS Source table.
 *
 * @parameter ds_p a pointer to the ASDM dataset.
 * @param ignoreTime a boolean value to indicate if the selected scans are taken into account or if all the table is going to be processed.
 * @param selectedScanRow_v a vector of pointers on ScanRow used to determine which rows of SysPower are going to be processed.
 */
void fillSource(ASDM* ds_p, bool ignoreTime, vector<ScanRow *>selectedScanRow_v) {
  LOGENTER("fillSource");
  const SourceTable& sourceT = ds_p->getSource();
  try {
    SourceRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << sourceT.size() << " sources(s)...";
    rowsInAScanbyTimeIntervalFunctor<SourceRow> selector(selectedScanRow_v);
    
    const vector<SourceRow *>& v = selector(sourceT.get(), ignoreTime);
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected scans ... ";

    info(infostream.str());
    int nSource = v.size();

    for (int i = 0; i < nSource; i++) {
      r = v.at(i);
      //
      // Check some assertions. 
      // For each row of the Source table, if any of the optional attributes which is an array and depend on numLines for the size of one
      // of its dimensions then the (optional) attribute numLines must be present and the dimensions depending on on numLines must be 
      // consistent with the value of numLines.
      //
      int numLines = r->isNumLinesExists() ? r->getNumLines() : 0;
      
      if (r->isTransitionExists()) {
	if (!r->isNumLinesExists()) {
	  errstream.str("");
	  errstream << "Source row#" << i << ". The attribute 'transition' exists but the attribute 'numLines' which serves to define its shape is missing. Can't go further.";
	  error(errstream.str());
	}

	int transitionSize = r->getTransition().size();
	if (numLines != transitionSize) {
	  errstream.str("");
	  errstream << "The value of 'numLines' (" << numLines << ") is not compatible with the found size of 'transition' (" << transitionSize << "). Can't go further.";
	  error(errstream.str());
	}
      }

      if (r->isRestFrequencyExists()) {
	if (!r->isNumLinesExists()) {
	  errstream.str("");
	  errstream << "Source row#" << i << ". The attribute 'restFrequency' exists but the attribute 'numLines' which serves to define its shape is missing. Cant' go further.";
	  error(errstream.str());
	}
	
	int restFrequencySize = r->getRestFrequency().size();
	if (numLines != restFrequencySize) {
	  errstream.str("");
	  errstream << "The value of 'numLines' (" << numLines << ") is not compatible with the found size of 'restFrequency' (" << restFrequencySize << "). Can't go further.";
	  error(errstream.str());
	}
      }

      if (r->isSysVelExists()) {
	if (!r->isNumLinesExists()) {
	  errstream.str("");
	  errstream << "Source row#" << i << ". The attribute 'sysVel' exists but the attribute 'numLines' which serves to define its shape is missing. Cant' go further.";
	  error(errstream.str());
	}
	
	int sysVelSize = r->getSysVel().size();
	if (numLines != sysVelSize) {
	  errstream.str("");
	  errstream << "The value of 'numLines' (" << numLines << ") is not compatible with the found size of 'sysVel' (" << sysVelSize << "). Can't go further.";
	  error(errstream.str());
	}
      }          

      int sourceId = r->getSourceId();
      // For now we just adapt the types of the time related informations and compute a mid-time.
      //
      double interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  r->getTimeInterval().getStartInMJD()*86400 ;
      // }
      // else {
      time =  r->getTimeInterval().getStartInMJD()*86400 + interval / 2.0 ;
      //}

      int spectralWindowId = swIdx2Idx[r->getSpectralWindowId().getTagValue()];

      string sourceName = r->getSourceName();

      int calibrationGroup = r->isCalibrationGroupExists() ? r->getCalibrationGroup() : 0;
      
      string code = r->getCode();

      vector<double> direction = DConverter::toVectorD(r->getDirection());
 
      DirectionReferenceCodeMod::DirectionReferenceCode dirRefCode = DirectionReferenceCodeMod::J2000;
      if(r->isDirectionCodeExists()){
	dirRefCode = r->getDirectionCode();
	//cout << "found directionCode for source " << sourceName << ": ";
      }
      //else{
      //  cout << "No directionCode in input table. Assuming ";
      //}
      string directionCode = CDirectionReferenceCode::name(dirRefCode);
      //cout << directionCode << endl;

      vector<double> position ;
      if (r->isPositionExists()){
	position = DConverter::toVectorD<Length>(r->getPosition());
      } 
				
      vector<double> properMotion = DConverter::toVectorD(r->getProperMotion());
  
      vector<string> transition;
      if (r->isTransitionExists()) {
	transition = r->getTransition();
      }

      vector<double> restFrequency;
      if (r->isRestFrequencyExists()) {
	restFrequency = DConverter::toVectorD<Frequency>(r->getRestFrequency());
      }

      vector<double> sysVel;
      if (r->isSysVelExists()) {
	sysVel = DConverter::toVectorD<Speed>(r->getSysVel());
      }
   
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addSource(sourceId,
				time,
				interval,
				spectralWindowId,
				numLines,
				sourceName,
				calibrationGroup,
				code,
				direction,
				directionCode,
				position,
				properMotion,
				transition,
				restFrequency,
				sysVel);
      }
    }
    
    if (nSource) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->source().nrow() <<" source(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  }
  catch (IllegalAccessException& e) {
    errstream.str("");
    error(errstream.str());
  }
  LOGEXIT("fillSource");
}

/** 
 * This function fills the MS Spectral Window table.
 * given :
 * @parameter ds_p a pointer to the ASDM dataset.
 */
void fillSpectralWindow(ASDM* ds_p) {
  LOGENTER("fillSpectralWindow");

  try {
    SpectralWindowTable& spwT = ds_p->getSpectralWindow();      
    vector<Tag> reorderedSwIds = reorderSwIds(*ds_p); // The vector of Spectral Window Tags in the order they will be inserted in the MS.
    //for (vector<Tag>::size_type i = 0; i != reorderedSwIds.size() ; i++) cerr<<" reorderedSwIds["<<i<<"]="<<reorderedSwIds[i].getTagValue()<<endl;
 
    for (vector<Tag>::size_type i = 0; i != reorderedSwIds.size() ; i++) swIdx2Idx[reorderedSwIds[i].getTagValue()] = i;

    SpectralWindowRow* r = 0;
    int nSpectralWindow = spwT.size();
    
    infostream.str("");
    infostream << "The dataset has " << nSpectralWindow << " spectral window(s)..."; 
    info(infostream.str());
    
    for (vector<Tag>::size_type i = 0; i < reorderedSwIds.size(); i++) {
      if ((r = spwT.getRowByKey(reorderedSwIds[i])) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the SpectralWindow table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
            
      /* Processing the chanFreqXXX
       *
       * Either (chanFreqStart, chanFreqStep) or (chanFreqArray) must be present
       * with a priority to chanFreqArray if both are present. If chanFreqArray
       * is present it *must* have a size equal to numChan otherwiser EXIT.
       */
      bool chanStartAndStep = r->isChanFreqStartExists() && r->isChanFreqStepExists();
      bool chanArray = r->isChanFreqArrayExists();
      if (!chanStartAndStep && !chanArray) {
	errstream.str("");
	errstream << "Did not find (chanFreqStart, chanFreqStep) nor chanFreqArray. Can't go further.";
	error(errstream.str());
      }
      
      //double* chanFreq1D = (double *) 0;
      //DConverter chanFreqConverter;
      vector<double> chanFreq1D;
      vector<Frequency> chanFreqArray;
      Frequency chanFreqStart, chanFreqStep;
      if (chanArray) { // Frequency channels are specified by an array.
	chanFreqArray = r->getChanFreqArray();
	if (chanFreqArray.size() != (unsigned int)r->getNumChan()) {
	  errstream.str("");
	  errstream << "Size of chanFreqArray ('"
		    << chanFreqArray.size()
		    << "') is not equal to numChan ('"
		    << r->getNumChan()
		    << "'). Can't go further.";
	  error(errstream.str());
	}
      }
      else { // Frequency channels are specified by a (start, step) pair.
	chanFreqStart = r->getChanFreqStart();
	chanFreqStep  = r->getChanFreqStep();
	for (int i = 0; i < r->getNumChan(); i++)
	  chanFreqArray.push_back(chanFreqStart + i * chanFreqStep);
      }
      //chanFreq1D = chanFreqConverter.to1DArray(chanFreqArray);
      chanFreq1D = DConverter::toVectorD<Frequency>(chanFreqArray);
      
      /* Processing the chanWidthXXX
       *
       * Either chanWidth or chanWidthArray must be present
       * with a priority to chanWidthArray if both are present. If chanWidthArray
       * is present it *must* have a size equal to numChan otherwiser EXIT.
       */
      if (!r->isChanWidthExists() && !r->isChanWidthArrayExists()) {
	errstream.str("");
	errstream << "Did not find chanWidth nor chanWidthArray. Can't go further.";
	error(errstream.str());
      }
      
      //double* chanWidth1D = (double *) 0;
      //DConverter chanWidthConverter;
      vector<double> chanWidth1D;
      vector<Frequency> chanWidthArray;
      if (r->isChanWidthArrayExists()) { // Frequency channels widths are specified by an array.
	chanWidthArray = r->getChanWidthArray();
	if (chanWidthArray.size() != (unsigned int) r->getNumChan()) {
	  errstream.str("");
	  errstream << "Size of chanWidthArray ('"
		    << chanWidthArray.size()
		    << "') is not equal to numChan ('"
		    << r->getNumChan()
		    << "'). Can't go further.";
	  error(errstream.str());
	}
      }
      else { // Frequency channels widths are specified by a constant value.
	chanWidthArray.resize(r->getNumChan());
	chanWidthArray.assign(chanWidthArray.size(), r->getChanWidth());
      }
      //chanWidth1D = chanWidthConverter.to1DArray(chanWidthArray);
      chanWidth1D = DConverter::toVectorD<Frequency>(chanWidthArray);

      // To answer JIRA ticket CAS-3265 / Dirk Petry
      if (chanFreq1D[chanFreq1D.size() - 1] < chanFreq1D[0] )
	std::transform(chanWidth1D.begin(), chanWidth1D.end(), chanWidth1D.begin(), negateFunctor<double>());
      
      /* Processing the effectiveBwXXX
       *
       * Either effectiveBw or effectiveBwArray must be present
       * with a priority to effectiveBwArray if both are present. If effectiveBwArray
       * is present it *must* have a size equal to numChan otherwiser EXIT.
       */
      if (!r->isEffectiveBwExists() && !r->isEffectiveBwArrayExists()) {
	errstream.str("");
	errstream << "Did not find effectiveBw nor effectiveBwArray. Can't go further.";
	error(errstream.str());
      }
      
      //double* effectiveBw1D = (double *) 0;
      vector<double> effectiveBw1D;
      //DConverter effectiveBwConverter;
      vector<Frequency> effectiveBwArray;
      if (r->isEffectiveBwArrayExists()) { // Effective BWs are specified by an array.
	effectiveBwArray = r->getEffectiveBwArray();
	if (effectiveBwArray.size() != (unsigned int) r->getNumChan()) {
	  errstream.str("");
	  errstream << "Size of effectiveBwArray ('"
		    << effectiveBwArray.size()
		    << "') is not equal to numChan ('"
		    << r->getNumChan()
		    << "'). Can't go further." ;
	  error(errstream.str());
	}
      }
      else { // Effective BWs are specified by a constant value.
	effectiveBwArray.resize(r->getNumChan());
	effectiveBwArray.assign(effectiveBwArray.size(), r->getEffectiveBw());
      }
      //effectiveBw1D = effectiveBwConverter.to1DArray(effectiveBwArray);
      effectiveBw1D = DConverter::toVectorD<Frequency>(effectiveBwArray);
      
      
      /* Processing the resolutionXXX
       *
       * Either resolution or resolutionArray must be present
       * with a priority to resolutionArray if both are present. If resolutionArray
       * is present it *must* have a size equal to numChan otherwiser EXIT.
       */
      if (!r->isResolutionExists() && !r->isResolutionArrayExists()) {
	errstream.str("");
	errstream << "Did not find resolution nor resolutionArray. Can't go further";
	error(errstream.str());
      }
      
      //double* resolution1D = (double *) 0;
      vector<double> resolution1D;
      //DConverter resolutionConverter;
      vector<Frequency> resolutionArray;
      if (r->isResolutionArrayExists()) { // Resolutions are specified by an array.
	resolutionArray = r->getResolutionArray();
	if (resolutionArray.size() != (unsigned int) r->getNumChan()) {
	  errstream.str("");
	  errstream << "Size of resolutionArray ('"
		    << resolutionArray.size()
		    << "') is not equal to numChan ('"
		    << r->getNumChan()
		    << "'). Can't go further.";
	  error(errstream.str());
	}
      }
      else { // Resolutions are specified by a constant value.
	resolutionArray.resize(r->getNumChan());
	resolutionArray.assign(resolutionArray.size(), r->getResolution());
      }
      //resolution1D = resolutionConverter.to1DArray(resolutionArray);
      resolution1D = DConverter::toVectorD<Frequency>(resolutionArray);

      /*
       * associated spectral windows and and natures.
       */
      
      unsigned int numAssocValues = 0;
      if (r->isNumAssocValuesExists()) numAssocValues = r->getNumAssocValues();
      
      // Test the simultaneous presence or absence of assocNature and assoSpectralWindowId
      if ((r->isAssocNatureExists() && !r->isAssocSpectralWindowIdExists()) ||
	  (!r->isAssocNatureExists() && r->isAssocSpectralWindowIdExists())) {
	errstream.str("");
	errstream << "Only one of the attributes assocSpectralWindowId and assocNature is present. Can't go further."
		  << endl;
	error(errstream.str());
      }
      
      vector<int> assocSpectralWindowId_;
      vector<string> assocNature_ ;

      if (r->isAssocSpectralWindowIdExists()) { // it exists then the assocNature exists also, given the test
	                                        // which is done before.
	vector<Tag> assocSpectralWindowId = r->getAssocSpectralWindowId();
	if (numAssocValues != assocSpectralWindowId.size()) {
	  infostream.str("");
	  infostream << "The size of assocSpectralWindowId ('"
		     << assocSpectralWindowId.size()
		     << "') is not equal to the value announced in numAssocValues ('"
		     << numAssocValues
		     << "'). Ignoring the difference and sending the full vector assocSpectralWindowId to the filler" 
	    ;
	  info(infostream.str());
	}
	numAssocValues = assocSpectralWindowId.size();
	assocSpectralWindowId_ = IConverter::toVectorI(assocSpectralWindowId);

	// Take into account the re ordering of the spectral window indices.
	for (unsigned int iAssocSw = 0; iAssocSw < numAssocValues; iAssocSw++)
	  assocSpectralWindowId_[iAssocSw] =  swIdx2Idx[assocSpectralWindowId_[iAssocSw]];

	vector<SpectralResolutionType> assocNature = r->getAssocNature();

	if (assocNature.size() != assocSpectralWindowId_.size()) {
	  infostream.str("");
	  infostream << "The size of assocNature ('"
		     << assocNature.size() 
		     << "') is not equal to the size of assocSpectralWindowId ('"
		     << assocSpectralWindowId.size()
		     << "'). Ignoring the difference and sending the full assocNature vector to the filler.";
	  info(infostream.str());
	}
	assocNature_ = SConverter::toVectorS<SpectralResolutionType, CSpectralResolutionType>(r->getAssocNature());
      }
      
      int numChan           = r->getNumChan();
      string name           = r->isNameExists()?r->getName():"";
      double refFreq        = r->getRefFreq().get();
      int measFreqRef       = r->isMeasFreqRefExists()?FrequencyReferenceMapper::value(r->getMeasFreqRef()):MFrequency::TOPO;
      double totalBandwidth = r->getTotBandwidth().get();
      int netSideband       = r->getNetSideband();
      int bbcNo             = r->getBasebandName();
      int ifConvChain       = 0;
      int freqGroup         = r->isFreqGroupExists()?r->getFreqGroup():0;
      string freqGroupName  = r->isFreqGroupNameExists()?r->getFreqGroupName().c_str():"";

      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addSpectralWindow(numChan,
					name,
					refFreq,
					chanFreq1D,  
					chanWidth1D, 
					measFreqRef,
					effectiveBw1D, 
					resolution1D,
					totalBandwidth,
					netSideband,
					bbcNo,
					ifConvChain,
					freqGroup,
					freqGroupName,
					numAssocValues,
					assocSpectralWindowId_,
					assocNature_);      
      }      
    }
    if (nSpectralWindow) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->spectralWindow().nrow() << " spectral window(s) in the measurement set(s).";
      info(infostream.str());
    }
  }
  catch (IllegalAccessException& e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillSpectralWindow");
}  

/**
 * This function fills the MS State table.
 * given :
 * @parameter r_p a pointer on the row of the ASDM Main table being processed.
 *
 */

void fillState(MainRow* r_p) {
  LOGENTER("fillState");

  ASDM&			ds	   = r_p -> getTable() . getContainer();
  ScanRow*		scanR_p	   = ds.getScan().getRowByKey(r_p -> getExecBlockId(),	r_p -> getScanNumber());
  vector<ScanIntent>	scanIntent = scanR_p -> getScanIntent();
  SubscanRow*		sscanR_p   = ds.getSubscan().getRowByKey(r_p -> getExecBlockId(),
								 r_p -> getScanNumber(),
								 r_p -> getSubscanNumber());
  if (sscanR_p == 0) {
    errstream.str("");
    errstream << "Could not find a row in the Subscan table for the following key value (execBlockId=" << r_p->getExecBlockId().toString()
	      <<", scanNumber="<< r_p->getScanNumber()
	      <<", subscanNum=" << r_p->getSubscanNumber() << "). Aborting. "
	      << endl;
    throw ASDM2MSException(errstream.str());
  }	  

  SubscanIntent subscanIntent = sscanR_p->getSubscanIntent();
  string obs_mode;
  if (scanIntent.size() > 0) {
    obs_mode = CScanIntent::name(scanIntent.at(0))+"#"+CSubscanIntent::name(subscanIntent);
   
    for (unsigned int iScanIntent = 1; iScanIntent < scanIntent.size(); iScanIntent++) {
      obs_mode += ",";
      obs_mode +=  CScanIntent::name(scanIntent.at(iScanIntent))+"#"+CSubscanIntent::name(subscanIntent);
    }
  }

  const vector<StateRow *>& sRs =  ds.getState().get() ;
  for (unsigned int iState = 0; iState < sRs.size(); iState++) {							     	    
    bool pushed = false;
    
    for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	 iter != msFillers.end();
	 ++iter) {
      int retId = iter->second->addUniqueState(sRs[iState]->getSig(),
					       sRs[iState]->getRef(),
					       0.0, 
					       0.0, 
					       r_p->getSubscanNumber(),
					       obs_mode, 
					       false);
      if (!pushed) {
	stateIdx2Idx[r_p] = retId;
	pushed = true;
      }
    }	    
  }
  LOGEXIT("fillState");
}

/*
** fill the SysCal table.
*
* @parameter ds_p a pointer to the ASDM dataset.
* @param ignoreTime a boolean value to indicate if the selected scans are taken into account or if all the table is going to be processed.
* @param selectedScanRow_v a vector of pointers on ScanRow used to determine which rows of SysPower are going to be processed.
*/
void fillSysCal(ASDM* ds_p, bool ignoreTime, vector<ScanRow*> selectedScanRow_v) {
  LOGENTER("fillSysCal");
  
  const SysCalTable& sysCalT = ds_p->getSysCal();
  try {
    SysCalRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << sysCalT.size() << " sysCal(s)...";
    rowsInAScanbyTimeIntervalFunctor<SysCalRow> selector(selectedScanRow_v);
    
    const vector<SysCalRow *>& v = selector(sysCalT.get(), ignoreTime);
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected scans ... ";
    
    info(infostream.str());
    int nSysCal = v.size();
    
    for (int i = 0; i < nSysCal; i++) {
      r = v.at(i);
      double interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double time;
      // if (isEVLA) {
      //   time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond ;
      // }
      // else {
      time =  ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval / 2.0 ;
      //}
      
      pair<bool, bool> tcal_flag_pair;
      tcal_flag_pair.first   = r->isTcalFlagExists();
      tcal_flag_pair.second  = r->isTcalFlagExists() ? r->getTcalFlag() : false;
      
      pair<bool, vector<float> > tcal_spectrum_pair;
      tcal_spectrum_pair.first  =  r->isTcalSpectrumExists() ;
      if (tcal_spectrum_pair.first)
	tcal_spectrum_pair.second = FConverter::toVectorF<Temperature>(r->getTcalSpectrum(), true);
      
      pair<bool, bool> trx_flag_pair;
      trx_flag_pair.first   = r->isTrxFlagExists();
      trx_flag_pair.second  = r->isTrxFlagExists() ? r->getTrxFlag() : false;
      
      pair<bool, vector<float> > trx_spectrum_pair;
      trx_spectrum_pair.first  =  r->isTrxSpectrumExists() ;
      if (trx_spectrum_pair.first)
	trx_spectrum_pair.second = FConverter::toVectorF<Temperature>(r->getTrxSpectrum(), true);
      
      pair<bool, bool> tsky_flag_pair;
      tsky_flag_pair.first   = r->isTskyFlagExists();
      tsky_flag_pair.second  = r->isTskyFlagExists() ? r->getTskyFlag() : false;
      
      pair<bool, vector<float> > tsky_spectrum_pair;
      tsky_spectrum_pair.first  =  r->isTskySpectrumExists() ;
      if (tsky_spectrum_pair.first)
	tsky_spectrum_pair.second = FConverter::toVectorF<Temperature>(r->getTskySpectrum(), true);

      pair<bool, bool> tsys_flag_pair;
      tsys_flag_pair.first   = r->isTsysFlagExists();
      tsys_flag_pair.second  = r->isTsysFlagExists() ? r->getTsysFlag() : false;

      pair<bool, vector<float> > tsys_spectrum_pair;
      tsys_spectrum_pair.first  =  r->isTsysSpectrumExists() ;
      if (tsys_spectrum_pair.first)
	tsys_spectrum_pair.second = FConverter::toVectorF<Temperature>(r->getTsysSpectrum(), true);

      pair<bool, bool> tant_flag_pair;
      tant_flag_pair.first   = r->isTantFlagExists();
      tant_flag_pair.second  = r->isTantFlagExists() ? r->getTantFlag() : false;

      pair<bool, vector<float> > tant_spectrum_pair;
      tant_spectrum_pair.first  =  r->isTantSpectrumExists() ;
      if (tant_spectrum_pair.first)
	tant_spectrum_pair.second = FConverter::toVectorF(r->getTantSpectrum(), true);
      
      pair<bool, bool> tant_tsys_flag_pair;
      tant_tsys_flag_pair.first   = r->isTantTsysFlagExists();
      tant_tsys_flag_pair.second  = r->isTantTsysFlagExists() ? r->getTantTsysFlag() : false;

      pair<bool, vector<float> > tant_tsys_spectrum_pair;
      tant_tsys_spectrum_pair.first  =  r->isTantTsysSpectrumExists() ;
      if (tant_tsys_spectrum_pair.first)
	tant_tsys_spectrum_pair.second = FConverter::toVectorF(r->getTantTsysSpectrum(), true);

      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addSysCal((int) r->getAntennaId().getTagValue(),
				(int) r->getFeedId(),
				(int) swIdx2Idx[r->getSpectralWindowId().getTagValue()],
				time,
				interval,
				r->getNumReceptor(),
				r->getNumChan(),
				tcal_spectrum_pair,
				tcal_flag_pair,
				trx_spectrum_pair,
				trx_flag_pair,
				tsky_spectrum_pair,
				tsky_flag_pair,
				tsys_spectrum_pair,
				tsys_flag_pair,
				tant_spectrum_pair,
				tant_flag_pair,
				tant_tsys_spectrum_pair,
				tant_tsys_flag_pair);				
      }
    }
    if (nSysCal) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->sysCal().nrow() <<" sysCal(s) in the measurement set(s)." ;
      info(infostream.str());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }  
  LOGEXIT("fillSysCal");
}


unsigned int numberOfRealValuesPerAutoCorrelation( unsigned int numCorr ) {
  LOGENTER("numberOfRealValuesPerAutoCorrelation");
  unsigned int result;
  switch (numCorr) {
  case 1 :
    result = 1;  // 1 real value (XX or YY)
    break;
    
  case 2 :
    result = 2;  // 2 real values (XX and YY)
    break;

  case 3 : 
    result =  4;  // 4 real values == 1 real value + 1 complex value (== 2 real values) + 1 real value (XX, XY, YY)
    break;

  case 4 : 
    result =  6;  // 6 real values == 1 real value + 1 complex value (== 2 real values) + 1 complex value (== 2 real values) + 1 real value (XX, XY, YX, YY)
    break;

  default :
    cout << "Invalid number of correlations '" << numCorr << "'." << endl;
    exit(3);
  }
  LOGEXIT("numberOfRealValuesPerAutoCorrelation");
  return result;
}

template <class T> 
vector<vector<pair<float, float > > > complexifyCrossData(unsigned long int numData, T* data, unsigned int numCorr) {
  LOGENTER("complexifyCrossData");
  vector<vector<pair<float, float > > > result;
  for (unsigned int j = 0; j < numData; j=j+2*numCorr) {
    vector<pair<float, float > >  p_v;
    for (unsigned int i = 0; i < 2*numCorr; i=i+2) 
      p_v.push_back(pair<float, float>(1.0*data[i+j], 1.0*data[i+j+1]));
    result.push_back(p_v);
  }
  LOGEXIT("complexifyCrossData");
  return result;
}

string complexifyCrossData_s (float* data, uint32_t numData, uint32_t numCorr, uint32_t numChar) {
  LOGENTER("complexifyCrossData_s");
  ostringstream oss;
  oss << "(";

  for (uint32_t j = 0; j < numData; j = j + 2 * numCorr) {
    if ( j == 1 ) oss << "," ;
    oss << "(";
    for ( uint32_t i = 0; i < 2 * numCorr; i = i + 2 ) {
      if ( i == 1 ) oss << "," ;
      oss << "(" << data[i+j] << ", " << data[i+j+1] << ")";
    }
    oss << ")";
  }

  oss << ")";
  LOGEXIT("complexifyCrossData_s");
  return oss.str().substr(0, numChar);
}

void complexifyCrossData_f(float* data, uint32_t numData, uint32_t /* numCorr */, vector<float>& result) {
  LOGENTER("complexifyCrossData_f");
  result.resize(numData);
  copy(data, data+numData, result.begin());
  LOGEXIT("complexifyCrossData_f");
}

template <class T>
void listCrossData(int64_t /* midpoint */, int64_t /* interval */, uint32_t numCrossData, T * data_p, vector<AnyValueMap<string> >& bb_spw_v, AnyValueMap<string>& msGroup_avm, uint32_t dataDescIndex, vector<ofstream *>& binaryDataPath_p_v, bool mute) {
  LOGENTER("listCrossData");
  
  LOG("ici_1");
  float			oneOverScaleFactor  = 1.0 / bb_spw_v[dataDescIndex].getValue<float>("scaleFactor");
  vector<float>	cdSlice;
  std::transform ( data_p,
		   data_p+numCrossData,
		   std::back_inserter(cdSlice),
		   [oneOverScaleFactor](int i){return i*oneOverScaleFactor;});
  
  vector<Tag>	antennaId_v	    = msGroup_avm.getValue<vector<Tag> >("antennas");
  unsigned int	blIndex		    = 0;
  unsigned int	numCrossPolProducts = bb_spw_v[dataDescIndex].getValue<vector<StokesParameter> >("crossPolProducts").size();
  // unsigned int  numDD               = msGroup_avm.getValue<vector<Tag> >("ddId").size();

  unsigned int		base		    = 0;
  unsigned int		sliceStart	    = 0;
  unsigned int		sliceEnd	    = 0;
  //cout << "Writing to dataDescIndex = " << dataDescIndex << endl;
  for (unsigned int iAntenna2 = 0; iAntenna2 < antennaId_v.size(); iAntenna2++) {
    for (unsigned int iAntenna1 = 0; iAntenna1 < iAntenna2; iAntenna1++) {
      base				    = 2 * blIndex * msGroup_avm.getValue<uint32_t>("blOffset");	// 2 because we consider complex values.
      LOG("ici_2");
      sliceStart			    = base + 2 * msGroup_avm.getValue<vector<uint32_t> >("cdBBOffset")[dataDescIndex];// 2 because we consider complex values.
      LOG("ici_3");
      sliceEnd				    = sliceStart + 2 * msGroup_avm.getValue<vector<uint32_t> >("numChan")[dataDescIndex] * numCrossPolProducts;// 2 because we consider complex values.
      LOG("ici_4");

      LOG("ici_5");
      
      vector<float> result;
      complexifyCrossData_f(cdSlice.data() + sliceStart, sliceEnd - sliceStart, numCrossPolProducts, result);
      
      // if (! mute) *(binaryDataPath_p_v[dataDescIndex]) << "crossdata - time = " << ArrayTime(midpoint-interval/2).toFITS()
      // 						       << ", scaleFactor = " << 1./oneOverScaleFactor
      // 						       << ", ANTENNA1 = " << antennaId_v[iAntenna1]
      // 						       << ", ANTENNA2 = " << antennaId_v[iAntenna2]
      // 						       << ", base = " << base
      // 						       << ", sliceStart = " << sliceStart
      // 						       << ", sliceEnd = " << sliceEnd
      // 						       << ", bl index = " << blIndex
      // 						       << ", dataDescIndex = " << dataDescIndex << "/" << numDD
      // 						       << ", size of cdSlice = " << cdSlice.size()
      if (! mute) binaryDataPath_p_v[dataDescIndex]->write((char *) &result[0], result.size()*sizeof(char));
      blIndex++;
    }
  }
  LOGEXIT("listCrossData");
}

string complexifyAutoData_s (const AUTODATATYPE* data, uint32_t numAutoData, uint32_t numCorr, uint32_t numChar) {
  LOGENTER("complexifyCrossData_s");
  ostringstream oss;
  oss << "(";
  switch (numCorr) {
  case 1:
  case 2:
    for_each (data, data+numAutoData, [&oss](int i){oss << "(" << i << ")";});
    break;

  case 3 :
    for (uint32_t iAutoData = 0; iAutoData < numAutoData; iAutoData = iAutoData + 4)
      oss << "((" << data[iAutoData] << ", 0.0)"
	  << "(" << data[iAutoData+1] << "," << data[iAutoData+2] << ")"
	  << "(" << data[iAutoData+1] << "," << -data[iAutoData+2] << ")"
	  << "(" << data[iAutoData+3] << ", 0.0))";
    break;

  case 4 :
    for (uint32_t iAutoData = 0; iAutoData < numAutoData; iAutoData = iAutoData + 6)
      oss << "((" << data[iAutoData] << ", 0.0)"
	  << "(" << data[iAutoData+1] << "," << data[iAutoData+2] << ")"
	  << "(" << data[iAutoData+3] << "," << data[iAutoData+4] << ")"
	  << "(" << data[iAutoData+5] << ", 0.0))";   
    
  default :
    cout << "Invalid value for numCorr '" << numCorr << "'." << endl;
    exit(2);
  }

  oss << ")";
  LOGEXIT("complexifyCrossData_s");
  return oss.str().substr(0, numChar);
}

void complexifyAutoData_f (const AUTODATATYPE* data, uint32_t numAutoData, uint32_t numCorr, vector<float>& result) {
  LOGENTER("complexifyCrossData_f");
  
  uint32_t k = 0;
  switch (numCorr) {
  case 1:
  case 2:
    result.resize(numAutoData);
    copy ( data, data+numAutoData, result.begin());
    break;

  case 3 :
    result.resize(2 * numAutoData);
    for (uint32_t iAutoData = 0; iAutoData < numAutoData; iAutoData = iAutoData + 4) {
      result[k] = data[iAutoData]; k++;
      result[k] = 0.0; k++;
      result[k] = data[iAutoData+1]; k++;
      result[k] = data[iAutoData+2]; k++;
      result[k] = data[iAutoData+1]; k++;
      result[k] = -data[iAutoData+1]; k++;
      result[k] = data[iAutoData+3]; k++;
      result[k] = 0.0; k++;      
    }
    break;
    
  case 4 :
    result.resize((numAutoData*4)/3);
    for (uint32_t iAutoData = 0; iAutoData < numAutoData; iAutoData = iAutoData + 6) {
      result[k] = data[iAutoData]; k++;
      result[k] = 0.0; k++;
      result[k] = data[iAutoData+1]; k++;
      result[k] = data[iAutoData+2]; k++;
      result[k] = data[iAutoData+3]; k++;
      result[k] = data[iAutoData+4]; k++;
      result[k] = data[iAutoData+5]; k++;
      result[k] = 0.0; k++;
    }
  default :
    cout << "Invalid value for numCorr '" << numCorr << "'." << endl;
    exit(2);
  }

  LOGEXIT("complexifyCrossData_f");
}

void listAutoData(int64_t /* midpoint */, int64_t /* interval */, uint32_t /* numAutoData */, const AUTODATATYPE * data_p, vector<AnyValueMap<string> >& bb_spw_v, AnyValueMap<string>& msGroup_avm, uint32_t dataDescIndex, vector<ofstream *>& binaryDataPath_p_v, bool mute) {
  LOGENTER("listAutoData");

  vector<Tag>	antennaId_v	    = msGroup_avm.getValue<vector<Tag> >("antennas");
  unsigned int	numSdPolProducts    = bb_spw_v[dataDescIndex].getValue<vector<StokesParameter> >("sdPolProducts").size();
  // uint32_t      numDD               = msGroup_avm.getValue<vector<Tag> >("ddId").size();

  unsigned int sliceStart = 0;
  unsigned int sliceEnd   = 0; 

  for (uint32_t iAntenna = 0; iAntenna < antennaId_v.size(); iAntenna++) {
    sliceStart	= iAntenna * msGroup_avm.getValue<uint32_t>("antOffset") + msGroup_avm.getValue<vector<uint32_t> >("adBBOffset")[dataDescIndex];
    sliceEnd    = sliceStart + msGroup_avm.getValue<vector<uint32_t> >("numChan")[dataDescIndex] * numberOfRealValuesPerAutoCorrelation(numSdPolProducts);
    // string complexAutoData_s = complexifyAutoData_s(data_p + sliceStart, sliceEnd - sliceStart, numSdPolProducts, 100);
    vector<float> result;
    complexifyAutoData_f(data_p + sliceStart, sliceEnd - sliceStart, numSdPolProducts, result);
    // if (! mute) *(binaryDataPath_p_v[dataDescIndex]) << "autodata - time = " << ArrayTime(midpoint-interval/2).toFITS()
    // 						     << ", ANTENNA1 = " << antennaId_v[iAntenna]
    // 						     << ", ANTENNA2 = " << antennaId_v[iAntenna]
    // 						     << ", sliceStart = " << sliceStart
    // 						     << ", sliceEnd = " << sliceEnd
    // 						     << ", dataDescIndex = " << dataDescIndex << "/" << numDD
    // 						     << ", size of autoData = " << numAutoData
    if (! mute) binaryDataPath_p_v[dataDescIndex]->write((char *) &result[0], result.size() * sizeof(char));
  }
  LOGEXIT("listAutoData");
}

string complexifyRadiometricData_s (const AUTODATATYPE* data, uint32_t numRadiometricData, uint32_t numCorr, uint32_t numChar) {
  LOGENTER("complexifyRadiometricData_s");
  ostringstream oss;
  oss << "(";
  for (uint32_t j = 0 ; j < numRadiometricData; j=j+numCorr) {
    oss << "(";
    for (uint32_t i = 0; i < numCorr; i++) 
      oss << "(" << data[i+j] << ", 0)";
    oss << ")";
  }
  oss << ")";
  LOGEXIT("complexifyRadiometricData_s");
  return oss.str().substr(0, numChar);
}

void complexifyRadiometricData_f(const AUTODATATYPE* data, uint32_t numRadiometricData, uint32_t numCorr, vector<float>& result) {
  LOGENTER("complexifyRadiometricData_f");
  result.resize(2 * numRadiometricData);
  uint32_t k = 0;
  for (uint32_t j = 0 ; j < numRadiometricData; j=j+numCorr)
    for (uint32_t i = 0; i < numCorr; i++) { 
      result[k] = data[i+j]; k++;
      result[k] = 0.0; k++;
    }
  LOGEXIT("complexifyRadiometricData_f");
}

void listRadiometricData(int64_t /* midpoint */, int64_t /* interval */, uint32_t numTime, uint32_t /* numRadiometricData */, const AUTODATATYPE * data_p, vector<AnyValueMap<string> >& bb_spw_v, AnyValueMap<string>& msGroup_avm, uint32_t dataDescIndex, vector<ofstream *>& binaryDataPath_p_v, bool mute) {
  LOGENTER("listRadiometricData");

  vector<Tag>	antennaId_v	 = msGroup_avm.getValue<vector<Tag> >("antennas");
  unsigned int	numSdPolProducts = bb_spw_v[dataDescIndex].getValue<vector<StokesParameter> >("sdPolProducts").size();
  // uint32_t      numDD            = msGroup_avm.getValue<vector<Tag> >("ddId").size();
  
  // float deltaTime = interval / ((float) numTime);
  // float startTime = midpoint - (interval - deltaTime) / 2.0;
  
  uint32_t	antOffset	    = msGroup_avm.getValue<uint32_t>("antOffset");
  uint32_t      rmdOffset           = msGroup_avm.getValue<vector<uint32_t> >("rmdBBOffset")[dataDescIndex];
  uint32_t      numChan             = msGroup_avm.getValue<vector<uint32_t> >("numChan")[dataDescIndex];
  uint32_t	allAntennaBlockSize = antennaId_v.size() * antOffset;
  uint32_t	sliceStart	    = 0;
  uint32_t	sliceEnd	    = 0;
  uint32_t	base		    = 0;

  for (uint32_t timeCount = 0; timeCount < numTime; timeCount++) {
    base	 = timeCount *			allAntennaBlockSize;
    for (uint32_t iAntenna = 0; iAntenna < antennaId_v.size(); iAntenna++) {
      sliceStart = base + iAntenna * antOffset + rmdOffset;
      sliceEnd	 = sliceStart + numChan * numSdPolProducts;
      //string complexRadiometricData_s = complexifyRadiometricData_s(data_p + sliceStart, sliceEnd - sliceStart, numSdPolProducts, 100); 
      vector<float> result;
      complexifyRadiometricData_f(data_p + sliceStart, sliceEnd - sliceStart, numSdPolProducts, result); 
      // if (!mute) *(binaryDataPath_p_v[dataDescIndex]) << "radiometricdata - time = " << (ArrayTime(startTime + timeCount * deltaTime)).toFITS()
      // 						      << ", ANTENNA1 = " << antennaId_v[iAntenna]
      // 						      << ", ANTENNA2 = " << antennaId_v[iAntenna]
      // 						      << ", sliceStart = " << sliceStart
      // 						      << ", sliceEnd = " << sliceEnd
      // 						      << ", dataDescIndex = " << dataDescIndex << "/" << numDD
      // 						      << ", size of radiometricData = " << numRadiometricData
      if (!mute) binaryDataPath_p_v[dataDescIndex]->write((char *) &result[0], result.size()*sizeof(char));
    }
  }    
  LOGEXIT("listRadiometricData");
}


void fillMainPar( const string& dsPath, ASDM * ds_p, bool doparallel, bool mute ) {
  LOGENTER("fillMainPar");

  vector<string>	processorType_v	     = CProcessorType::names();
  vector<string>	spectralResolution_v = CSpectralResolutionType::names();

  string processorType("CORRELATOR");
  string spectralResolution("FULL_RESOLUTION");

  map<Tag, AnyValueMap<string> > msGroup_m = msGroup(ds_p);
    
  std::regex crossDataAxesRE("BAL( +BAB)?( +APC)?( +SPP)?( +POL)?");
  std::regex autoDataAxesRE("ANT( +BAB)?( +SPP)?( +POL)?");
  std::regex radiometerAxesRE("TIM ANT( +SPP)?");

  /*
   * Now process the bdf files grouped by configDescriptionId.
   *
   */
  vector<Tag> keys_v;
  std::transform( msGroup_m.begin(),
		  msGroup_m.end(),
		  std::back_inserter(keys_v),
		  std::bind (&std::map<Tag, AnyValueMap<string> >::value_type::first, std::placeholders::_1)
		  );

  // cout << "Size of keys = " << keys_v.size() << endl;
  // copy ( keys_v.begin(), keys_v.end(), std::ostream_iterator<Tag>(std::cout, "\n") );

  //#pragma omp parallel for if(doparallel)
  for (uint32_t iKey = 0; iKey < keys_v.size(); iKey++) {
    Tag configDescriptionId = keys_v[iKey];

    cout << "Beginning of processing for " + configDescriptionId.toString() + "\n";

    AnyValueMap<string>& msGroup_avm = msGroup_m[configDescriptionId];

    ostringstream oss;
    oss << CProcessorType::toString(msGroup_avm.getValue<ProcessorTypeMod::ProcessorType>("processorType"))
	<< "_"
	<< CSpectralResolutionType::toString(msGroup_avm.getValue<SpectralResolutionTypeMod::SpectralResolutionType>("spectralResolution"));
    vector<Tag> dataDescriptionId_v = msGroup_avm.getValue<vector<Tag> >("ddId");
    vector<ofstream *> binaryDataPath_p_v(dataDescriptionId_v.size());
    uint32_t iBinaryDataPath = 0;
    for (Tag dataDescriptionId: dataDescriptionId_v ) {
      binaryDataPath_p_v[iBinaryDataPath] = new ofstream ((oss.str()+"_"+dataDescriptionId.toString()).c_str(), ios_base::binary);
      iBinaryDataPath++;
    }

    ProcessorTypeMod::ProcessorType processorType = msGroup_avm.getValue<ProcessorTypeMod::ProcessorType>("processorType");
    for( string bdfName: msGroup_avm.getValue<vector<string> >("bdfNames") ) {
      try {
	if ( processorType == ProcessorTypeMod::CORRELATOR) {
	  SDMDataObjectStreamReader sdosr;
	  string fullPath;
	  fullPath = dsPath + "/ASDMBinary/" + bdfName;
	  cout << fullPath << endl;
	  sdosr.open(fullPath);
	  /*
	   * Listify the baseband/spectralwindows elements found in the global XML header.
	   */

	  vector<uint32_t> cdBBOffset_v(1, 0);
	  vector<uint32_t> adBBOffset_v(1, 0);
	  vector<uint32_t> numChan_v = msGroup_avm.getValue<vector<uint32_t> >("numChan");

	  vector<AnyValueMap<string> > bb_spw_v;
	  unsigned int ibb_spw = 1;
	  const vector<SDMDataObject::Baseband>& baseband_v = sdosr.dataStruct().basebands();
	  for ( SDMDataObject::Baseband  bb: baseband_v ) {
	    const vector<SDMDataObject::SpectralWindow>& spectralWindow_v = bb.spectralWindows();

	    for ( SDMDataObject::SpectralWindow spectralWindow: spectralWindow_v ) {
	      AnyValueMap<string> attributes_avm;
	      attributes_avm.setValue("sdPolProducts", spectralWindow.sdPolProducts());
	      attributes_avm.setValue("crossPolProducts", spectralWindow.crossPolProducts());
	      attributes_avm.setValue("scaleFactor", spectralWindow.scaleFactor());
	      attributes_avm.setValue("numSpectralPoint", spectralWindow.numSpectralPoint());
	      attributes_avm.setValue("numBin", spectralWindow.numBin());
	         
	      bb_spw_v.push_back(attributes_avm);
	      cdBBOffset_v.push_back(cdBBOffset_v[ibb_spw - 1] +  numChan_v[ibb_spw - 1] * spectralWindow.crossPolProducts().size());
	      adBBOffset_v.push_back(adBBOffset_v[ibb_spw - 1] +  numChan_v[ibb_spw - 1] * numberOfRealValuesPerAutoCorrelation(spectralWindow.sdPolProducts().size()));
	      ibb_spw++;
	    } 
	  }

	  msGroup_avm.setValue("cdBBOffset", cdBBOffset_v);
	  msGroup_avm.setValue("adBBOffset", adBBOffset_v);
	  
	  msGroup_avm.setValue("blOffset", cdBBOffset_v.back());
	  msGroup_avm.setValue("antOffset", adBBOffset_v.back());

	  /*
	   * Then let's traverse the BDF.
	   */
	  
	  while (true) {
	    const vector<SDMDataSubset>& subset_v     = sdosr.nextSubsets(20);
	    uint32_t nSubsets = subset_v.size();
	    cout << "nSubsets = " << nSubsets << endl;

	    if (subset_v.size() == 0) break;

	    vector<int64_t> midpoint_v;
	    vector<int64_t> interval_v;
	    vector<uint32_t> numAutoData_v;
	    vector<const AUTODATATYPE *> autoData_p_v;
	    vector<uint32_t> numCrossData_v;
	    vector<PrimitiveDataTypeMod::PrimitiveDataType> crossDataType_v;
	    vector<void *> crossData_p_v;

	    for (uint32_t subsetIndex = 0; subsetIndex < subset_v.size(); subsetIndex++) {
	      const SDMDataSubset& subset = subset_v[subsetIndex];
	      midpoint_v.push_back(subset.time());
	      interval_v.push_back(subset.interval());

	      autoData_p_v.push_back(0);
	      numAutoData_v.push_back(subset.autoData(autoData_p_v.back()));

	      crossDataType_v.push_back(subset.crossDataType());
	      
	      const SHORTCROSSDATATYPE *	shortData_p;
	      const INTCROSSDATATYPE *		intData_p;
	      const FLOATCROSSDATATYPE *	floatData_p;
	      
	      switch (crossDataType_v.back()) {
	      case PrimitiveDataTypeMod::INT16_TYPE:
		numCrossData_v.push_back(subset.crossData(shortData_p));
		crossData_p_v.push_back((void *) shortData_p);
		break;
		
	      case PrimitiveDataTypeMod::INT32_TYPE:
		numCrossData_v.push_back(subset.crossData(intData_p));
		crossData_p_v.push_back((void *) intData_p);
		break;
		
	      case PrimitiveDataTypeMod::FLOAT32_TYPE:
		numCrossData_v.push_back(subset.crossData(floatData_p));
		crossData_p_v.push_back((void *) floatData_p);
		break;
		
	      default:
		cout << "Impossible data type for cross correlation data ('" << crossDataType_v.back() << "')." << endl; 
	      }	   	      	      
	    }
 
	    /*
	      const SDMDataSubset&	subset	      = sdosr.getSubset();
	      int64_t      		midpoint      = subset.time();
	      int64_t     		interval      = subset.interval();
	      unsigned                  numAutoData   = 0;
	      PrimitiveDataType		crossDataType = subset.crossDataType() ;
	      unsigned long int		numCrossData  = 0;
	    
	      //
	      // Read ONCE
	      //
	    
	      // Autocorrelation data
	      const AUTODATATYPE * autoData_p;
	    
	      numAutoData = subset.autoData(autoData_p);
	    
	      // Cross correlation data

	      const SHORTCROSSDATATYPE *	shortData_p;
	      const INTCROSSDATATYPE *		intData_p;
	      const FLOATCROSSDATATYPE *	floatData_p;

	      switch (crossDataType) {
	      case PrimitiveDataTypeMod::INT16_TYPE:
	      numCrossData	      = subset.crossData(shortData_p);
	      break;
		
	      case PrimitiveDataTypeMod::INT32_TYPE:
	      numCrossData = subset.crossData(intData_p);
	      break;
		
	      case PrimitiveDataTypeMod::FLOAT32_TYPE:
	      numCrossData = subset.crossData(floatData_p);
	      break;
		
	      default:
	      cout << "Impossible data type for cross correlation data ('" << crossDataType << "')." << endl; 
	      }	    
	      
	    */

	    /**
	     * Write MANY
	     */

	    const SHORTCROSSDATATYPE *	shortData_p;
	    const INTCROSSDATATYPE *	intData_p;
	    const FLOATCROSSDATATYPE *	floatData_p;
	    const AUTODATATYPE *        autoData_p;

#pragma omp parallel for schedule(dynamic) if(doparallel)
	    for (uint32_t dataDescIndex = 0; dataDescIndex < dataDescriptionId_v.size(); dataDescIndex++) {
	      
	      //	      if (omp_get_thread_num() == 1) 
	      //cout << "Number of threads=" + TO_STRING(omp_get_num_threads()) + "\n";

	      for (uint32_t subsetIndex = 0; subsetIndex < subset_v.size(); subsetIndex++) {
		int64_t midpoint = midpoint_v[subsetIndex];
		int64_t interval = interval_v[subsetIndex];
		uint32_t numAutoData = numAutoData_v[subsetIndex];
		PrimitiveDataTypeMod::PrimitiveDataType crossDataType = crossDataType_v[subsetIndex];
		autoData_p = autoData_p_v[subsetIndex];
		listAutoData(midpoint, interval, numAutoData, autoData_p, bb_spw_v, msGroup_avm, dataDescIndex, binaryDataPath_p_v, mute);

		/* Cross correlation data */
		uint32_t numCrossData = numCrossData_v[subsetIndex];
		switch (crossDataType) {
		case PrimitiveDataTypeMod::INT16_TYPE:
		  shortData_p = ( SHORTCROSSDATATYPE *) crossData_p_v[subsetIndex];
		  listCrossData(midpoint, interval, numCrossData, shortData_p, bb_spw_v, msGroup_avm,  dataDescIndex, binaryDataPath_p_v, mute);
		  break;
		  
		case PrimitiveDataTypeMod::INT32_TYPE:
		  intData_p = ( INTCROSSDATATYPE *) crossData_p_v[subsetIndex];
		  listCrossData(midpoint, interval, numCrossData, intData_p, bb_spw_v, msGroup_avm, dataDescIndex, binaryDataPath_p_v, mute);
		  break;
		  
		case PrimitiveDataTypeMod::FLOAT32_TYPE:
		  shortData_p = ( SHORTCROSSDATATYPE *) crossData_p_v[subsetIndex];
		  listCrossData(midpoint, interval, numCrossData, floatData_p, bb_spw_v, msGroup_avm, dataDescIndex, binaryDataPath_p_v, mute);
		  break;
		  
		default:
		  cout << "Impossible data type for cross correlation data ('" << crossDataType << "')." << endl; 
		}	    
	      }
	    }
	  }
	  sdosr.close();     
	}
	else if (processorType == ProcessorTypeMod::RADIOMETER) {
	  SDMDataObjectReader sdor;
	  string fullPath;
	  fullPath = dsPath + "/ASDMBinary/" + bdfName;
	  const SDMDataObject& sdo = sdor.read(fullPath);

	  vector<uint32_t> rmdBBOffset_v(1, 0);
	  vector<uint32_t> numChan_v = msGroup_avm.getValue<vector<uint32_t> >("numChan");
	  
	  vector<AnyValueMap<string> > bb_spw_v;
	  unsigned int ibb_spw = 1;
	  const vector<SDMDataObject::Baseband>& baseband_v = sdo.dataStruct().basebands();
	  for ( SDMDataObject::Baseband  bb: baseband_v ) {
	    const vector<SDMDataObject::SpectralWindow>& spectralWindow_v = bb.spectralWindows();

	    for ( SDMDataObject::SpectralWindow spectralWindow: spectralWindow_v ) {
	      AnyValueMap<string> attributes_avm;
	      attributes_avm.setValue("sdPolProducts", spectralWindow.sdPolProducts());
	      attributes_avm.setValue("numSpectralPoint", spectralWindow.numSpectralPoint());
	      attributes_avm.setValue("numBin", spectralWindow.numBin());
	         
	      bb_spw_v.push_back(attributes_avm);
	      rmdBBOffset_v.push_back(rmdBBOffset_v[ibb_spw - 1] +  numChan_v[ibb_spw - 1] * spectralWindow.sdPolProducts().size());
	      ibb_spw++;
	    } 
	  }
	  msGroup_avm.setValue("rmdBBOffset", rmdBBOffset_v);
	  msGroup_avm.setValue("antOffset", rmdBBOffset_v.back());
	  
	  const SDMDataSubset& subset = sdo.tpDataSubset();
	  const AUTODATATYPE* radiometricData_p;
	  uint32_t numRadiometricData = subset.autoData(radiometricData_p);
	  int64_t  midpoint = subset.time();
	  int64_t  interval = subset.interval();
	  uint32_t numTime  = sdo.numTime();
#pragma omp parallel for schedule(dynamic) if(doparallel)
	  for (uint32_t dataDescIndex = 0; dataDescIndex < dataDescriptionId_v.size(); dataDescIndex++) {
	    listRadiometricData(midpoint, interval, numTime, numRadiometricData, radiometricData_p, bb_spw_v, msGroup_avm, dataDescIndex, binaryDataPath_p_v, mute);
	  }
	  sdor.done();
	}
      }
      catch (SDMDataObjectReaderException e) {
	cout << e.getMessage() << endl;
      }
      catch (SDMDataObjectStreamReaderException e) {
	cout << e.getMessage() << endl;
      }
    }
    for (unsigned int iBinaryDataPath = 0; iBinaryDataPath < binaryDataPath_p_v.size(); iBinaryDataPath++)
      binaryDataPath_p_v[iBinaryDataPath]->close();
    cout << "End of processing for " + configDescriptionId.toString() + "\n";
  }

  LOGEXIT("fillMainPar");
}

/**
 * This function fills the MS Main table from an ASDM Main table which refers to correlator data.
 *
 * given:
 * @parameter rowNum an integer expected to contain the number of the row being processed.
 * @parameter r_p a pointer to the MainRow being processed.
 * @parameter sdmBinData a reference to the SDMBinData containing a lot of information about the binary data being processed. Useful to know the requested ordering of data.
 * @parameter uvwCoords a reference to the UVW calculator.
 * @parameter complexData a bool which says if the DATA is going to be filled (true) or if it will be the FLOAT_DATA (false).
 * @parameter mute if the value of this parameter is false then nothing is written in the MS .
 *
 * !!!!! One must be carefull to the fact that fillState must have been called before fillMain. During the execution of fillState , the global vector<int> msStateID
 * is filled and will be used by fillMain.
 */ 
void fillMain(int		/* rowNum */,
	      MainRow*		r_p,
	      SDMBinData&	sdmBinData,
	      const VMSData*	vmsData_p,
	      UvwCoords&	uvwCoords,
	      bool		complexData,
	      bool              mute) {
  
  if (debug) cout << "fillMain : entering" << endl;

  // ASDM & ds = r_p -> getTable() . getContainer();

  // Then populate the Main table.
  ComplexDataFilter filter; // To process the case numCorr == 3
  
  if (vmsData_p->v_antennaId1.size() == 0) {
    infostream.str("");
    infostream << "No MS data produced for the current row." << endl;
    info(infostream.str());
    return;
  }

  vector<vector<unsigned int> > filteredShape = vmsData_p->vv_dataShape;
  for (unsigned int ipart = 0; ipart < vmsData_p->vv_dataShape.size(); ipart++) {
    if (filteredShape.at(ipart).at(0) == 3)
      filteredShape.at(ipart).at(0) = 4;
  }
	  
  vector<int> filteredDD;
  for (unsigned int idd = 0; idd < vmsData_p->v_dataDescId.size(); idd++){
    filteredDD.push_back(dataDescriptionIdx2Idx.at(vmsData_p->v_dataDescId.at(idd)));
  }


  vector<float *> uncorrectedData;
  vector<float *> correctedData;

  /* compute the UVW */
  vector<double> uvw(3*vmsData_p->v_time.size());
	  
  vector<casacore::Vector<casacore::Double> > vv_uvw;
#if DDPRIORITY
  uvwCoords.uvw_bl(r_p, sdmBinData.timeSequence(), e_query_cm, 
		   sdmbin::SDMBinData::dataOrder(),
		   vv_uvw);
#else
  uvwCoords.uvw_bl(r, vmsData_p->v_timeCentroid, e_query_cm, 
		   sdmbin::SDMBinData::dataOrder(),
		   vv_uvw);
#endif
  int k = 0;
  for (unsigned int iUvw = 0; iUvw < vv_uvw.size(); iUvw++) {
    uvw[k++] = vv_uvw[iUvw](0); 
    uvw[k++] = vv_uvw[iUvw](1);
    uvw[k++] = vv_uvw[iUvw](2);
  } 

  // Here we make the assumption that the State is the same for all the antennas and let's use the first State found in the vector stateId contained in the ASDM Main Row
  // int asdmStateIdx = r_p->getStateId().at(0).getTagValue();  
  vector<int> msStateId(vmsData_p->v_m_data.size(), stateIdx2Idx[r_p]);

  ComplexDataFilter cdf;
  map<AtmPhaseCorrectionMod::AtmPhaseCorrection, float*>::const_iterator iter;

  vector<double>	correctedTime;
  vector<int>		correctedAntennaId1;
  vector<int>		correctedAntennaId2;
  vector<int>		correctedFeedId1;
  vector<int>		correctedFeedId2;
  vector<int>		correctedFieldId;
  vector<int>           correctedFilteredDD;
  vector<double>	correctedInterval;
  vector<double>	correctedExposure;
  vector<double>	correctedTimeCentroid;
  vector<int>		correctedMsStateId(msStateId);
  vector<double>	correctedUvw ;
  vector<unsigned int>	correctedFlag;

  // Do we have to fill an MS with uncorrected data + radiometric data (radiometric data are considered as uncorrected data)  ?
  if (msFillers.find(AP_UNCORRECTED) != msFillers.end()) {
    for (unsigned int iData = 0; iData < vmsData_p->v_m_data.size(); iData++) {
      if ((iter=vmsData_p->v_m_data.at(iData).find(AtmPhaseCorrectionMod::AP_UNCORRECTED)) != vmsData_p->v_m_data.at(iData).end()){
	uncorrectedData.push_back(cdf.to4Pol(vmsData_p->vv_dataShape.at(iData).at(0),
					     vmsData_p->vv_dataShape.at(iData).at(1),
					     iter->second));
      }
    }
  }
	    
  // Do we have to fill an MS with corrected data + radiometric data ( again considered as uncorrected data) ?
  if (msFillers.find(AP_CORRECTED) != msFillers.end() ) {

    // Are we with radiometric data ? Then we assume that the data are labelled AP_UNCORRECTED.
    if (sdmBinData.processorType(r_p) == RADIOMETER) {
      for (unsigned int iData = 0; iData < vmsData_p->v_m_data.size(); iData++) {
	if ((iter=vmsData_p->v_m_data.at(iData).find(AtmPhaseCorrectionMod::AP_UNCORRECTED)) != vmsData_p->v_m_data.at(iData).end()){
	  correctedTime.push_back(vmsData_p->v_time.at(iData));
	  correctedAntennaId1.push_back(vmsData_p->v_antennaId1.at(iData));
	  correctedAntennaId2.push_back(vmsData_p->v_antennaId2.at(iData));
	  correctedFeedId1.push_back(vmsData_p->v_feedId1.at(iData));
	  correctedFeedId2.push_back(vmsData_p->v_feedId2.at(iData));
	  correctedFilteredDD.push_back(filteredDD.at(iData));
	  correctedFieldId.push_back(vmsData_p->v_fieldId.at(iData));
	  correctedInterval.push_back(vmsData_p->v_interval.at(iData));
	  correctedExposure.push_back(vmsData_p->v_exposure.at(iData));
	  correctedTimeCentroid.push_back(vmsData_p->v_timeCentroid.at(iData));
	  correctedUvw.push_back(vv_uvw.at(iData)(0));
	  correctedUvw.push_back(vv_uvw.at(iData)(1));
	  correctedUvw.push_back(vv_uvw.at(iData)(2));
	  correctedData.push_back(cdf.to4Pol(vmsData_p->vv_dataShape.at(iData).at(0), // Force radiometric data to go
					     vmsData_p->vv_dataShape.at(iData).at(1), // into correctedData.
					     iter->second));
	  correctedFlag.push_back(vmsData_p->v_flag.at(iData));
	}
      }
    }
    
    else {
      // we are in front of non radiometric data , but data supposedly with AP_CORRECTED data.
      for (unsigned int iData = 0; iData < vmsData_p->v_m_data.size() ; iData++) {
	if  (vmsData_p->v_antennaId1.at(iData) == vmsData_p->v_antennaId2.at(iData) ) {
	  /*
	  ** do not forget to prepend the autodata copied from the uncorrected data, because the lower layers of the software do not put the (uncorrected) autodata in the
	  ** corrected data.
	  */
	  correctedTime.push_back(vmsData_p->v_time.at(iData));
	  correctedAntennaId1.push_back(vmsData_p->v_antennaId1.at(iData));
	  correctedAntennaId2.push_back(vmsData_p->v_antennaId2.at(iData));
	  correctedFeedId1.push_back(vmsData_p->v_feedId1.at(iData));
	  correctedFeedId2.push_back(vmsData_p->v_feedId2.at(iData));
	  correctedFilteredDD.push_back(filteredDD.at(iData));
	  correctedFieldId.push_back(vmsData_p->v_fieldId.at(iData));
	  correctedInterval.push_back(vmsData_p->v_interval.at(iData));
	  correctedExposure.push_back(vmsData_p->v_exposure.at(iData));
	  correctedTimeCentroid.push_back(vmsData_p->v_timeCentroid.at(iData));
	  correctedUvw.push_back(vv_uvw.at(iData)(0));
	  correctedUvw.push_back(vv_uvw.at(iData)(1));
	  correctedUvw.push_back(vv_uvw.at(iData)(2));
	  correctedData.push_back(uncorrectedData.at(iData));    // <-------- Here we re-use the autodata already present in the uncorrected data.
	  correctedFlag.push_back(vmsData_p->v_flag.at(iData));	    
	}
	else {
	  /*
	  ** And now finally the correlation corrected data.
	  */
	  correctedTime.push_back(vmsData_p->v_time.at(iData));
	  correctedAntennaId1.push_back(vmsData_p->v_antennaId1.at(iData));
	  correctedAntennaId2.push_back(vmsData_p->v_antennaId2.at(iData));
	  correctedFeedId1.push_back(vmsData_p->v_feedId1.at(iData));
	  correctedFeedId2.push_back(vmsData_p->v_feedId2.at(iData));
	  correctedFilteredDD.push_back(filteredDD.at(iData));
	  correctedFieldId.push_back(vmsData_p->v_fieldId.at(iData));
	  correctedInterval.push_back(vmsData_p->v_interval.at(iData));
	  correctedExposure.push_back(vmsData_p->v_exposure.at(iData));
	  correctedTimeCentroid.push_back(vmsData_p->v_timeCentroid.at(iData));
	  correctedUvw.push_back(vv_uvw.at(iData)(0));
	  correctedUvw.push_back(vv_uvw.at(iData)(1));
	  correctedUvw.push_back(vv_uvw.at(iData)(2));
	  iter=vmsData_p->v_m_data.at(iData).find(AtmPhaseCorrectionMod::AP_CORRECTED);
	  float* theData = cdf.to4Pol(vmsData_p->vv_dataShape.at(iData).at(0),
				      vmsData_p->vv_dataShape.at(iData).at(1),
				      iter->second);
	  correctedData.push_back(theData);
	  correctedFlag.push_back(vmsData_p->v_flag.at(iData));
	}
      }
    }
  }
	  
  if (uncorrectedData.size() > 0 && (msFillers.find(AP_UNCORRECTED) != msFillers.end())) {
    if (! mute) {
      // when found, missing weight and sigma arguments, added here to set default values to check with compiler
      vector<double> weight_v(vmsData_p->v_time.size(),1.0);
      vector<double> sigma_v(vmsData_p->v_time.size(),1.0);
      msFillers[AP_UNCORRECTED]->addData(complexData,
					 (vector<double>&) vmsData_p->v_time, // this is already time midpoint
					 (vector<int>&) vmsData_p->v_antennaId1,
					 (vector<int>&) vmsData_p->v_antennaId2,
					 (vector<int>&) vmsData_p->v_feedId1,
					 (vector<int>&) vmsData_p->v_feedId2,
					 filteredDD,
					 vmsData_p->processorId,
					 (vector<int>&)vmsData_p->v_fieldId,
					 (vector<double>&) vmsData_p->v_interval,
					 (vector<double>&) vmsData_p->v_exposure,
					 (vector<double>&) vmsData_p->v_timeCentroid,
					 (int) r_p->getScanNumber(), 
					 0,                                               // Array Id
					 (int) r_p->getExecBlockId().getTagValue(), // Observation Id
					 (vector<int>&)msStateId,
					 uvw,
					 filteredShape, // vmsData_p->vv_dataShape after filtering the case numCorr == 3
					 uncorrectedData,
					 (vector<unsigned int>&)vmsData_p->v_flag,
					 weight_v,
					 sigma_v);
    }
  }

  if (correctedData.size() > 0 && (msFillers.find(AP_CORRECTED) != msFillers.end())) {
    if (! mute) {
      // missing correctedWeight and correctedSigma, added default values here to 
      // check compiler. Not correct
      vector<double> correctedWeight(correctedTime.size(),1.0);
      vector<double> correctedSigma(correctedTime.size(),1.0);
      msFillers[AP_CORRECTED]->addData(complexData,
				       correctedTime, // this is already time midpoint
				       correctedAntennaId1, 
				       correctedAntennaId2,
				       correctedFeedId1,
				       correctedFeedId2,
				       correctedFilteredDD,
				       vmsData_p->processorId,
				       correctedFieldId,
				       correctedInterval,
				       correctedExposure,
				       correctedTimeCentroid,
				       (int) r_p->getScanNumber(), 
				       0,                                               // Array Id
				       (int) r_p->getExecBlockId().getTagValue(), // Observation Id
				       correctedMsStateId,
				       correctedUvw,
				       filteredShape, // vmsData_p->vv_dataShape after filtering the case numCorr == 3
				       correctedData,
				       correctedFlag,
				       correctedWeight,
				       correctedSigma);
    }
  }
  if (debug) cout << "fillMain : exiting" << endl;
}


void testFunc(string& tstr) {
  cerr<<tstr<<endl;
}

/**
 * compute the UVW (put in a method to keep sepearate from fillMain for
 parallel case) returns a Matrix, mat_uvw for data passing thread-safe
 way
*/
void calcUVW(MainRow* r_p, 
             SDMBinData& sdmBinData, 
             const VMSData* vmsData_p, 
             UvwCoords& uvwCoords,
             casacore::Matrix< casacore::Double >& mat_uvw) {

  vector< casacore::Vector<casacore::Double> > vv_uvw;
  mat_uvw.resize(3,vmsData_p->v_time.size());

#if DDPRIORITY
  uvwCoords.uvw_bl(r_p, sdmBinData.timeSequence(), e_query_cm,
                   sdmbin::SDMBinData::dataOrder(),
                   vv_uvw);
#else
  uvwCoords.uvw_bl(r, vmsData_p->v_timeCentroid, e_query_cm,
                   sdmbin::SDMBinData::dataOrder(),
                   vv_uvw);
#endif

  //put in a Matrix
  for (unsigned int iUvw = 0; iUvw < vv_uvw.size(); iUvw++) {
    mat_uvw(iUvw, 0) = vv_uvw[iUvw](0);
    mat_uvw(iUvw, 1) = vv_uvw[iUvw](1);
    mat_uvw(iUvw, 2) = vv_uvw[iUvw](2);
  }
} 


void fillProcessor(ASDM* ds_p) {
  LOGENTER("fillProcessor");
  try {
    ProcessorTable& processorT = ds_p->getProcessor();
    ProcessorRow* r = 0;
    int nProcessor = processorT.size();
    
    infostream.str("");
    infostream << "The dataset has " << nProcessor << " processor(s)...";
    info(infostream.str());
    
    for (int i = 0; i < nProcessor; i++) {
      if ((r=processorT.getRowByKey(Tag(i, TagType::Processor))) == 0) {
	errstream.str("");
	(errstream << "Problem while reading the Processor table, the row with key = Tag(" << i << ") does not exist.Aborting." << endl);
	error(errstream.str());
      }
      
      string processorType    = CProcessorType::name(r->getProcessorType());
      string processorSubType = CProcessorSubType::name(r->getProcessorSubType());
      
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addProcessor(processorType,
				   processorSubType,
				   -1,    // Since there is no typeId in the ASDM.
				   r->getModeId().getTagValue());
      }  
    }
    if (nProcessor) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->processor().nrow() << " processor(s) in the measurement set." ;
      info(infostream.str());
    } 
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillProcessor");
}

void fillSysPower_aux (const vector<SysPowerRow *>& sysPowers, map<AtmPhaseCorrection, ASDM2MSFiller*>& msFillers_m) {
  if (debug) cout << "fillSysPower_aux : entering" << endl;
  vector<int>		antennaId;
  vector<int>		spectralWindowId;
  vector<int>		feedId;
  vector<double>	time;
  vector<double>	interval;
  vector<int>		numReceptor;
  vector<float>		switchedPowerDifference;
  vector<float>		switchedPowerSum;
  vector<float>		requantizerGain;
  
  if (debug) cout << "fillSysPower_aux : resizing the arrays (" << sysPowers.size() << ") to populate the columns of the MS SYSPOWER table." << endl;

  antennaId.resize(sysPowers.size());
  spectralWindowId.resize(sysPowers.size());
  feedId.resize(sysPowers.size());
  time.resize(sysPowers.size());
  interval.resize(sysPowers.size());
  
  /*
   * Prepare the mandatory attributes.
   */
  if (debug) cout << "fillSysPower_aux : filling the arrays to populate the columns of the MS SYSPOWER table." << endl;

  std::transform(sysPowers.begin(), sysPowers.end(), antennaId.begin(), sysPowerAntennaId);
  std::transform(sysPowers.begin(), sysPowers.end(), spectralWindowId.begin(), sysPowerSpectralWindowId);
  std::transform(sysPowers.begin(), sysPowers.end(), feedId.begin(), sysPowerFeedId);
  std::transform(sysPowers.begin(), sysPowers.end(), time.begin(), sysPowerMidTimeInSeconds);
  std::transform(sysPowers.begin(), sysPowers.end(), interval.begin(), sysPowerIntervalInSeconds);
  
  /*
   * Prepare the optional attributes.
   */
  if (debug) cout << "fillSysPower_aux : working on the optional attributes." << endl;

  unsigned int numReceptor0 = sysPowers[0]->getNumReceptor();
  if (debug) cout << "fillSysPower_aux : numReceptor = " << numReceptor0 << endl;
 
  bool switchedPowerDifferenceExists0 = sysPowers[0]->isSwitchedPowerDifferenceExists();
  if (switchedPowerDifferenceExists0) {
    switchedPowerDifference.resize(numReceptor0 * sysPowers.size());
    for_each(sysPowers.begin(), sysPowers.end(), sysPowerSwitchedPowerDifference(switchedPowerDifference.begin()));
  }

  bool switchedPowerSumExists0 = sysPowers[0]->isSwitchedPowerSumExists();
  if (switchedPowerSumExists0) {
    switchedPowerSum.resize(numReceptor0 * sysPowers.size());
    for_each(sysPowers.begin(), sysPowers.end(), sysPowerSwitchedPowerSum(switchedPowerSum.begin()));
  }
  
  bool requantizerGainExists0 = sysPowers[0]->isRequantizerGainExists();
  if (requantizerGainExists0) {
    requantizerGain.resize(numReceptor0 * sysPowers.size());
    for_each(sysPowers.begin(), sysPowers.end(), sysPowerRequantizerGain(requantizerGain.begin()));  
  }

  if (debug) cout << "fillSysPower_aux : about to append a slice to the MS SYSPOWER table." << endl;
 
  for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator msIter = msFillers_m.begin();
       msIter != msFillers_m.end();
       ++msIter) {
    msIter->second->addSysPowerSlice(antennaId.size(),
				     antennaId,
				     spectralWindowId,
				     feedId,
				     time,
				     interval,
				     (unsigned int) numReceptor0,
				     switchedPowerDifference,
				     switchedPowerSum,
				     requantizerGain);
  }
  infostream << "Appended " << sysPowers.size() << " rows to the MS SYSPOWER table." << endl;
  if (debug) cout << "fillSysPower_aux : exiting" << endl;
}

/**
 * This function fills the MS SysPower table from an ASDM SysPower table.
 *
 * @param ds the ASDM dataset the ASDM SysPower table belongs to.
 * @param ignoreTime a boolean value to indicate if the selected scans are taken into account or if all the table is going to be processed.
 * @param selectedScanRow_v a vector of pointers on ScanRow used to determine which rows of SysPower are going to be processed.
 * @param msFillers_m a map of ASDM2MSFillers depending on AtmosphericPhaseCorrection.
 *
 */
void fillSysPower(const string asdmDirectory, ASDM* ds_p, bool ignoreTime, const vector<ScanRow *>& selectedScanRow_v, map<AtmPhaseCorrection, ASDM2MSFiller*>& msFillers_m) {
  LOGENTER("fillSysPower");

  const SysPowerTable& sysPowerT = ds_p->getSysPower();

  infostream.str("");
  infostream << "The dataset has " << sysPowerT.size() << " syspower(s).";
  info(infostream.str()); 
  
  if (sysPowerT.size() > 0 ) {
    try {
      // Prepare a row filter based on the time intervals of the selected scans.
      rowsInAScanbyTimeIntervalFunctor<SysPowerRow> selector(selectedScanRow_v);

      //
      // We can assume that there is an SysPower table , but we don't know yet if it's stored in a binary or an XML file.
      // 
      if (file_exists(uniqSlashes(asdmDirectory + "/SysPower.bin"))) {

	if (debug) cout << "fillSysPower : working with SysPower.bin by successive slices." << endl;

	TableStreamReader<SysPowerTable, SysPowerRow> tsrSysPower;
	tsrSysPower.open(asdmDirectory);

	// We can process the SysPower table by slice when it's stored in a binary file so let's do it.
	while (tsrSysPower.hasRows()) {
	  cout << "About to read rows until 5000000 bytes have been consumed" << endl;
	  const vector<SysPowerRow*>&	sysPowerRows = tsrSysPower.untilNBytes(50000000);
	  infostream.str("");
	  infostream << "(considering the next " << sysPowerRows.size() << " rows of the SysPower table. ";

	  if (debug) cout << "fillSysPower : determining which rows are in the selected scans." << endl;
	  const vector<SysPowerRow *>& sysPowers = selector(sysPowerRows, ignoreTime);
  
	  if (!ignoreTime) 
	    infostream << sysPowers.size() << " of them are in the selected exec blocks / scans";
  
	  infostream << ")";	     
	  info(infostream.str());
	  
	  infostream.str("");
	  errstream.str("");
  
	  if (sysPowers.size() > 0)
	    fillSysPower_aux(sysPowerRows, msFillers_m);
	}
	tsrSysPower.close();
      }
      
      else if (file_exists(uniqSlashes(asdmDirectory + "/SysPower.xml"))) {

	if (debug) cout << "fillSysPower : working with SysPower.xml read with a TableSAXReader" << endl;

	//
	// Instantiate a TableSAXReader functor with T==SysPowerTable, R==SysPowerRow 
	// and RFilter==rowsInAScanbyTimeIntervalFunctor<SysPowerRow>.
	//
	TableSAXReader<SysPowerTable, SysPowerRow, rowsInAScanbyTimeIntervalFunctor<SysPowerRow> >
	  tableSAXReader(verbose,
			 selector,
			 &fillSysPower_aux,
			 msFillers_m);
	
	// Execute the functor
	tableSAXReader(asdmDirectory, ignoreTime);

      }
      else 
	throw ConversionException ("fillSysPower: no file found for SysPower", "SysPower");

      unsigned int numMSSysPowers =  (const_cast<casacore::MeasurementSet*>(msFillers_m.begin()->second->ms()))->rwKeywordSet().asTable("SYSPOWER").nrow();
      if (numMSSysPowers > 0) {
	infostream.str("");
	infostream << "converted in " << numMSSysPowers << " syspower(s) in the measurement set.";
	info(infostream.str());
      }
      
    } // end of filling SysPower by slice.
    
    catch (ConversionException e) {
      errstream.str("");
      errstream << e.getMessage();
      error(errstream.str());
    }
    catch ( std::exception & e) {
      errstream.str("");
      errstream << e.what();
      error(errstream.str());      
    }
  }
  LOGEXIT("fillSysPower");
}

/**
 * @param ds the ASDM dataset the ASDM SysPower table belongs to.
 * @param ignoreTime a boolean value to indicate if the selected scans are taken into account or if all the table is going to be processed.
 * @param selectedScanRow_v a vector of pointers on ScanRow used to determine which rows of SysPower are going to be processed.
 */
void fillWeather(ASDM* ds_p, bool ignoreTime, vector<ScanRow *>& selectedScanRow_v) {
  LOGENTER("fillWeather");
  const WeatherTable& weatherT = ds_p->getWeather();

  try {
    WeatherRow* r = 0;
    infostream.str("");
    infostream << "The dataset has " << weatherT.size() << " weather(s)...";
    rowsInAScanbyTimeIntervalFunctor<WeatherRow> selector(selectedScanRow_v);
    
    const vector<WeatherRow *>& v = selector(weatherT.get(), ignoreTime);
    if (!ignoreTime) 
      infostream << v.size() << " of them in the selected scans ... ";

    info(infostream.str());
    int nWeather = v.size();

    infostream.str("");
    infostream << "The dataset has " << nWeather << " weather(s)...";
    info(infostream.str());
    
    pair<bool, float>
      pressureOpt,
      relHumidityOpt,
      temperatureOpt,
      windDirectionOpt,
      windSpeedOpt,
      dewPointOpt;
    
#define OPT_ATTR_PAIR( rowPtr, AttributeName ) rowPtr -> is ## AttributeName ## Exists() ? make_pair ( true, rowPtr -> get ## AttributeName ().get()) : make_pair( false, 0.)        
    
    for (int i = 0; i < nWeather; i++) {
      r			 = v.at(i);      
      double	interval = ((double) r->getTimeInterval().getDuration().get()) / ArrayTime::unitsInASecond ;
      double	time	 = ((double) r->getTimeInterval().getStart().get()) / ArrayTime::unitsInASecond + interval / 2.0;
      
      pressureOpt			   = OPT_ATTR_PAIR(r, Pressure);
      pressureOpt.second		  /= 100. ;	// We consider that ASDM stores Pascals & MS expects hectoPascals
      relHumidityOpt			   = OPT_ATTR_PAIR(r, RelHumidity);
      temperatureOpt			   = OPT_ATTR_PAIR(r, Temperature);
      windDirectionOpt			   = OPT_ATTR_PAIR(r, WindDirection);
      windSpeedOpt			   = OPT_ATTR_PAIR(r, WindSpeed);
      dewPointOpt			   = OPT_ATTR_PAIR(r, DewPoint);
      int		wxStationId        = r->getStationId().getTagValue();
      vector<double>	wxStationPosition  = DConverter::toVectorD(r->getStationUsingStationId()->getPosition());
    
      for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
	   iter != msFillers.end();
	   ++iter) {
	iter->second->addWeather(-1,
				 time,
				 interval,
				 pressureOpt,
				 relHumidityOpt,
				 temperatureOpt,
				 windDirectionOpt,
				 windSpeedOpt,
				 dewPointOpt,
				 wxStationId,
				 wxStationPosition);
      }
    }

    if (nWeather) {
      infostream.str("");
      infostream << "converted in " << msFillers.begin()->second->ms()->weather().nrow() <<" weather(s) in the measurement set." ;
      info(infostream.str());
    }
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch ( std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());      
  }
  LOGEXIT("fillWeather");
}

//
// Create and fill the "toplevel" Measurement Set.
void fillTopMS(ASDM*					ds_p,
	       const string&                            dsPath,
	       map<AtmPhaseCorrection, string>&	        msNames_m,
	       bool					complexData,
	       bool					withCompression,
	       const string&				telescopeName,
	       int					maxNumCorr,
	       int					maxNumChan,
	       bool                                     ignoreTime,
	       vector<ScanRow *>&                       selectedScanRow_v,
	       map<int, set<int> >&                     selected_eb_scan_m,
	       bool	                                processCalDevice,
	       bool                                     processEphemeris,
	       bool                                     processPointing,
	       bool                                     withPointingCorrection,
	       bool                                     processSysPower,
	       bool                                     asis,
	       const vector<string>&                    asisTablenames_v) {
  LOGENTER("fillTopMS");
  // Create the measurement set(s). 
  try {
    for (map<AtmPhaseCorrection, string>::const_iterator iter = msNames_m.begin(); iter != msNames_m.end(); ++iter) {
      info("About to create a filler for the measurement set '" + msNames_m[iter->first] + "'");
      msFillers[iter->first] = new ASDM2MSFiller(msNames_m[iter->first],
						 0.0,
						 false,
						 complexData,
						 withCompression,
						 telescopeName,
						 maxNumCorr,
						 maxNumChan,
						 false,
						 lazy);
    }
  }
  catch(AipsError & e) {
    errstream.str("");
    errstream << e.getMesg();
    error(errstream.str());
  }
  catch (std::exception & e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());
  }
  
  /* ASDM2MSFiller* msFiller_p = msFillers.begin()->second; */
 
  //
  // Firstly convert the basic tables.
  //
  // For purpose of simplicity we assume that in all ASDM basic tables having a Tag identifier
  // these Tag values are forming a sequence of integer 0 -> table.size() - 1. If that's not the case, the
  // program aborts.
  //

  //
  // Process the Antenna table.
  fillAntenna(ds_p);

  //
  // Process the SpectralWindow table.
  fillSpectralWindow(ds_p);

  //
  // Process the Polarization table
  vector<int> polarizationIdx2Idx_v;
  fillPolarization(ds_p, polarizationIdx2Idx_v);
  
  //
  // Process the DataDescription table.
  fillDataDescription(ds_p, polarizationIdx2Idx_v);

  // 
  // Process the ExecBlock table in order to build the MS Observation table.
  fillObservation(ds_p, selected_eb_scan_m);
  
  //
  // Process the Feed table
  fillFeed(ds_p, ignoreTime, selectedScanRow_v);

  //
  // Process the Ephemeris table.
  // Create and fill the MS ephemeris table(s) with a time interpolation time step set to 86400000000 nanoseconds ( 1/1000 day).
  if (processEphemeris) 
    fillEphemeris(ds_p, 86400000000LL);

  //
  // Process the Field table.
  // Now it respects the degree of the polynomials but it ignores the ephemerisId.
  // The ephemerisId will be processed during the call to fillEphemeris.
  fillField(ds_p, processEphemeris);
   
  //
  // Process the FlagCmd table.
  fillFlagCmd(ds_p, ignoreTime, selectedScanRow_v);

  //
  // Process the History table.
  fillHistory(ds_p, ignoreTime, selectedScanRow_v);


  //
  // Process the Pointing table.
  // Issues :
  // - pointingModelId , phaseTracking, sourceOffset and overTheTop not taken into account.
  if (processPointing) 
    fillPointing(ds_p, withPointingCorrection, ignoreTime, selectedScanRow_v);

  //
  // Process the processor table.
  fillProcessor(ds_p);
 
  //
  // Process the Source table.
  fillSource(ds_p, ignoreTime, selectedScanRow_v);

  //
  // Process the SysCal table.
  fillSysCal(ds_p, ignoreTime, selectedScanRow_v);
  
  //
  // Process the CalDevice table.
  if (processCalDevice && ds_p->getCalDevice().size() > 0) 
    fillCalDevice(ds_p, ignoreTime, selectedScanRow_v);

  //
  // Process the SysPower table.
  if ( processSysPower )
    fillSysPower(dsPath, ds_p, ignoreTime, selectedScanRow_v, msFillers);

  //
  // Process the weather table
  fillWeather(ds_p, ignoreTime, selectedScanRow_v);


  //
  // And then finally process the state and the main table.
  //fillMainPar( dsPath, ds_p, doparallel, mute );
  
  // Do we also want to store the verbatim copies of some tables of the ASDM dataset ?
  if (asis) {
    ASDMVerbatimFiller avf(const_cast<casacore::MS*>(msFillers.begin()->second->ms()), Name2Table::find(asisTablenames_v, verbose));
    avf.fill(*ds_p);
  }

  for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
       iter != msFillers.end();
       ++iter)
    iter->second->end(0.0);

  for (map<AtmPhaseCorrection, ASDM2MSFiller*>::iterator iter = msFillers.begin();
       iter != msFillers.end();
       ++iter)
    delete iter->second;

  LOGEXIT("fillTopMS");
}

// ------ data partition function
/*
void partitionMS(vector<int> SwIds, 
                 //vector< map<AtmPhaseCorrection,ASDM2MSFiller* > >&  msFillers_vec,
                 //map<AtmPhaseCorrection,string>  msFillers,
                 map<AtmPhaseCorrection,string> msNames,
                 bool complexData,
                 bool withCompression,
                 string telName,
                 int maxNumCorr,
                 int maxNumChan) {
  LOGENTER("partitionMS");
  for (int i=0; i<SwIds.size(); i++) {
    ostringstream oss;
    oss<< SwIds.at(i);
    string msname_suffix = ".SpW"+oss;
    //cerr<<"msname_prefix="<<msname_suffix<<endl;
    for (map<AtmPhaseCorrection, string>::iterator iter = msNames.begin(); iter != msNames.end(); ++iter) {
      msFillers[iter->first] = new ASDM2MSFiller(msNames[iter->first]+msname_suffix,
						 0.0,
						 false,
						 complexData,
						 withCompression,
						 telName,
						 maxNumCorr,
						 maxNumChan,
						 false,
						 lazy
						 );
      info("About to create a filler for the measurement set '" + msNames[iter->first] + msname_suffix + "'");
    }
    //vector<std::pair<const AtmPhaseCorrection,ASDM2MSFiller*> > msFillers_vec(msFillers.begin(),msFillers.end());
    msFillers_v.push_back(msFillers);
    //store ms names locally
  }
  //cerr<<"msFillers_v.size="<<msFillers_v.size()<<endl;
  LOGEXIT("partitionMS");
}
*/

class MSMainRowsInSubscanChecker {
public:
  MSMainRowsInSubscanChecker();
  virtual ~MSMainRowsInSubscanChecker();
  void check(const VMSData* vmsData_p, MainRow* mainRow_p, unsigned int mainRowIndex, const string& BDFName);
  const vector<string>& report() const;
  void reset();

private:
  vector<string> report_v;
};

MSMainRowsInSubscanChecker::MSMainRowsInSubscanChecker() {;}
MSMainRowsInSubscanChecker::~MSMainRowsInSubscanChecker() {;}
void MSMainRowsInSubscanChecker::reset() {
  LOGENTER("MSMainRowsInSubscanChecker::reset");
  report_v.clear();
  LOGEXIT("MSMainRowsInSubscanChecker::reset");
}

void MSMainRowsInSubscanChecker::check( const VMSData* vmsData_p,
					MainRow* mainRow_p,
					unsigned int mainRowIndex,
					const string& BDFName ) {
  LOGENTER("MSMainRowsInSubscanChecker::check");
  SubscanTable & subscanTable = mainRow_p->getTable().getContainer().getSubscan();

  SubscanRow* subscanRow_p = subscanTable.getRowByKey(mainRow_p->getExecBlockId(),
						      mainRow_p->getScanNumber(),
						      mainRow_p->getSubscanNumber());
  if (subscanRow_p == NULL) {
    infostream.str("");
    infostream << "Could not find a row in the subscan table with the key 'execBlockId = "<< mainRow_p->getExecBlockId()
	       << ", scanNumber = " << mainRow_p->getScanNumber()
	       << ", subscanNumber = " << mainRow_p->getSubscanNumber()
	       << "'. I can't check if the BDF contents is in the subscan's time range.";
    info(infostream.str());
    LOGEXIT("MSMainRowsInSubscanChecker::check");
    return;
  }

  // We make the assumption that the content pointed by vmsData_p is ordered by time.
  double subscanStartTime = subscanRow_p->getStartTime().getMJD()*86400.0;
  double subscanEndTime   = subscanRow_p->getEndTime().getMJD()*86400.0;

  //
  // Now detect one of two abnormal situations : the 1st data time is anterior to the subscan start time or the last data time
  // is posterior to the subscan end time. 
  if ( (vmsData_p->v_time[0] < subscanStartTime) || (subscanEndTime < vmsData_p->v_time[vmsData_p->v_time.size() - 1])) {
    ostringstream oss;
    oss << "Main row #" << mainRowIndex
	<< " - The BDF '" << BDFName << "' contained data not in the time range of scan=" << mainRow_p->getScanNumber()
	<< ", subscan=" << mainRow_p->getSubscanNumber() << ".";
    string s = oss.str();
    if (!(report_v.size() > 0 && s == report_v.back()))
      report_v.push_back(s);
  }
  LOGEXIT("MSMainRowsInSubscanChecker::check");
}

const vector<string>& MSMainRowsInSubscanChecker::report() const {
  return report_v;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * The main function.
 */
int main(int argc, char *argv[]) {

  string dsPath;
  string msNamePrefix;
  string msNameExtension;

  appName = string(argv[0]);
  ofstream ofs;

  // LogSinkInterface& lsif = LogSink::globalSink();

  uint64_t bdfSliceSizeInMb = 0; // The default size of the BDF slice hold in memory.

  bool mute = false;

  bool doparallel = false;

  vector<string> asisTablenames_v;  

  //   Process command line options and parameters.
  po::variables_map vm;

  try {


    // Declare the supported options.

    po::options_description generic("Converts an ASDM dataset into a CASA measurement set.\n"
				    "Usage : " + appName +" [options] asdm-directory [ms-directory-prefix]\n\n"
				    "Command parameters: \n"
				    " asdm-directory : the pathname to the ASDM dataset to be converted \n"
				    " ms-directory-prefix : the prefix of the pathname(s) of the measurement set(s ) to be created,\n"
				    " this prefis is completed by a suffix to form the name(s) of the resulting measurement set(s), \n"
				    " this suffix depends on the selected options (see options compression and wvr-correction) \n"
				    ".\n\n"
				    "Allowed options:");
    generic.add_options()
      ("help", "produces help message.")
      ("icm",  po::value<string>()->default_value("all"), "specifies the correlation mode to be considered on input. A quoted string containing a sequence of 'ao' 'co' 'ac' 'all' separated by whitespaces is expected")
      ("isrt", po::value<string>()->default_value("all"), "specifies the spectral resolution type to be considered on input. A quoted string containing a sequence of 'fr' 'ca' 'bw' 'all' separated by whitespaces is expected")
      ("its",  po::value<string>()->default_value("all"), "specifies the time sampling (INTEGRATION and/or SUBINTEGRATION)  to be considered on input. A quoted string containing a sequence of 'i' 'si' 'all' separated by whitespaces is expected")  
      ("ocm",  po::value<string>()->default_value("ca"),  "output data for correlation mode AUTO_ONLY (ao) or CROSS_ONLY (co) or CROSS_AND_AUTO (ca)")
      ("compression,c", "produces compressed columns in the resulting measurement set (not set by default). When this option is selected the string '-compressed' is inserted in the pathname of the resulting measurement set.")
      ("asis", po::value<string>(), "creates verbatim copies of the ASDM tables in the output measurement set. The value given to this option must be a quoted string containing a list of table names separated by space characters; the wildcard character '*' is allowed in table names.")
      ("wvr-corrected-data", po::value<string>()->default_value("no"),  "specifies wich values are considered in the ASDM binary data to fill the DATA column in the MAIN table of the MS. Expected values for this option are 'no' for the uncorrected data (this is the default), 'yes' for the corrected data and 'both' for corrected and uncorrected data. In the latter case, two measurement sets are created, one containing the uncorrected data and the other one, whose name is suffixed by '-wvr-corrected', containing the corrected data.")
      ("scans,s", po::value<string>(), "processes only the scans specified in the option's value. This value is a semicolon separated list of scan specifications. A scan specification consists in an exec bock index followed by the character ':' followed by a comma separated list of scan indexes or scan index ranges. A scan index is relative to the exec block it belongs to. Scan indexes are 1-based while exec blocks's are 0-based. \"0:1\" or \"2:2~6\" or \"0:1,1:2~6,8;2:,3:24~30\" \"1,2\" are valid values for the option. \"3:\" alone will be interpreted as 'all the scans of the exec block#3'. An scan index or a scan index range not preceded by an exec block index will be interpreted as 'all the scans with such indexes in all the exec blocks'.  By default all the scans are considered.")
      ("logfile,l", po::value<string>(), "specifies the log filename. If the option is not used then the logged informations are written to the standard error stream.")
      ("verbose,v", "logs numerous informations as the filler is working.")
      ("revision,r", "logs information about the revision of this application.")
      ("dry-run,m", "does not fill the MS MAIN table.")
      ("ignore-time,t", "all the rows of the tables Feed, History, Pointing, Source, SysCal, CalDevice, SysPower and Weather are processed independently of the time range of the selected exec block / scan.")
      ("no-caldevice", "The CalDevice table will be ignored.")
      ("no-ephemeris", "The ephemeris table will be ignored.")
      ("no-syspower", "the SysPower table will be  ignored.")
      ("no-pointing", "The Pointing table will be ignored.")
      ("check-row-uniqueness", "The row uniqueness constraint will be checked in the tables where it's defined")
      ("bdf-slice-size", po::value<uint64_t>(&bdfSliceSizeInMb)->default_value(500),  "The maximum amount of memory expressed as an integer in units of megabytes (1024*1024) allocated for BDF data. The default is 500 (megabytes)") 
      ("parallel", "run with multithreading mode.")
      ("with-pointing-correction", "add (ASDM::Pointing::encoder - ASDM::Pointing::pointingDirection) to the value to be written in MS::Pointing::direction - (related with JIRA tickets CSV-2878 and ICT-1532))");

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("asdm-directory", po::value< string >(), "asdm directory")
      ;
    hidden.add_options()
      ("ms-directory-prefix", po::value< string >(), "ms directory prefix")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    
    po::positional_options_description p;
    p.add("asdm-directory", 1);
    p.add("ms-directory-prefix", 1);
    
    // Parse the command line and retrieve the options and parameters.
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

    po::notify(vm);
    
    // Where do the log messages should go ?
    if (vm.count("logfile")) {
      //LogSinkInterface *theSink;
#if 0     
      // Replaced the change at the LogSink level ...
      ofs.open(vm["logfile"].as<string>().c_str(), ios_base::app);
      LogSinkInterface *theSink = new casacore::StreamLogSink(&ofs);
      LogSink::globalSink(theSink);
#else
      // ... with a change at the cerr (stderr) level since by default global logs are going to cerr (stderr).
      freopen(vm["logfile"].as<string>().c_str(), "a", stderr);
#endif
    }

    // Help ? displays help's content and don't go further.

    if (vm.count("help")) {
      errstream.str("");
      errstream << generic << "\n" ;
      error(errstream.str());
    }

    // Verbose or quiet ?
    verbose = vm.count("verbose") > 0;
   
    // Revision ? displays revision's info and don't go further if there is no dataset to process otherwise proceed....
    string revision = "$Id: asdm2MS.cpp,v 1.84 2011/10/25 14:56:48 mcaillat Exp $\n";
    if (vm.count("revision")) {
      if (!vm.count("asdm-directory")) {
	errstream.str("");
	errstream << revision ;
	error(errstream.str());
      }
      else {
	infostream.str("");
	infostream << revision ;
	info(infostream.str());
      }
    }

    // Selection of correlation mode of data to be considered on input.
    istringstream iss;
    string token;

    string icm_opt = vm["icm"].as< string >();
    iss.clear();
    iss.str(icm_opt);
    
    while (iss >> token) {
      if (token.compare("co") == 0) 
	es_cm.fromString("CROSS_ONLY", false);
      else if (token.compare("ao") == 0)
	es_cm.fromString("AUTO_ONLY", false);
      else if (token.compare("ac") == 0)
	es_cm.fromString("CROSS_AND_AUTO", false);
      else if (token.compare("all") == 0)
	es_cm.fromString("CROSS_ONLY AUTO_ONLY CROSS_AND_AUTO", false);
      else {
	errstream.str("");
	errstream << "Token '" << token << "' invalid for --icm option." << endl;
	errstream << generic << endl;
	error(errstream.str());
      }
    }

    // Selection of spectral resolution type of data to be considered.

    string isrt_opt = vm["isrt"].as< string >();
    iss.clear();
    iss.str(isrt_opt);

    while (iss >> token) {
      if (token.compare("fr") == 0)
	es_srt.fromString("FULL_RESOLUTION", false);
      else if (token.compare("ca") == 0)
	es_srt.fromString("CHANNEL_AVERAGE", false);
      else if (token.compare("bw") == 0)
	es_srt.fromString("BASEBAND_WIDE", false);
      else if (token.compare("all") == 0)
	es_srt.fromString("FULL_RESOLUTION CHANNEL_AVERAGE BASEBAND_WIDE", false);
      else { 
	errstream.str("");
	errstream << "Token '" << token << "' invalid for --isrt option." << endl;
	errstream << generic;
	error(errstream.str());
      }
    }


    // Selection of the time sampling of data to be considered (integration and/or subintegration)

    string its_opt = vm["its"].as < string >();
    iss.clear();
    iss.str(its_opt);

    while ( iss >> token ) {
      if (token.compare("i") == 0)
	es_ts.fromString("INTEGRATION",false);
      else if (token.compare("si") == 0)
	es_ts.fromString("SUBINTEGRATION", false);
      else if (token.compare("all") == 0)
	es_ts.fromString("INTEGRATION SUBINTEGRATION", false);
      else {
	errstream.str("");
	errstream << "Token '" << token << "' invalid for its option." << endl;
	errstream << generic ;
	error(errstream.str());
      }
    }

    // Selection of the correlation mode of data to be produced in the measurement set.

    string ocm_opt = vm["ocm"].as< string >();
    if ( ocm_opt.compare("co") == 0 )
      e_query_cm = CROSS_ONLY;
    else if ( ocm_opt.compare("ao") == 0 )
      e_query_cm = AUTO_ONLY;
    else if ( ocm_opt.compare("ca") == 0 )
      e_query_cm = CROSS_AND_AUTO;
    else {
      errstream.str("");
      errstream << "Token '" << ocm_opt << "' invalid for ocm option." << endl;
      errstream << generic ;
      error(errstream.str());
    }

    if (vm.count("asdm-directory")) {
      string dummy = vm["asdm-directory"].as< string >();
      dsPath = lrtrim(dummy) ;
      if (dsPath.back()=='/') dsPath.erase(dsPath.size()-1);
    }
    else {
      errstream.str("");
      errstream << generic ;
      error(errstream.str());
    }

    if (vm.count("ms-directory-prefix")) {
      string dummyMSName = vm["ms-directory-prefix"].as< string >();
      dummyMSName = lrtrim(dummyMSName);
      if (dummyMSName.back()=='/') dummyMSName.erase(dummyMSName.size()-1);
      Path msPath(dummyMSName);
      string msDirectory = msPath.dirName();
      if (msDirectory.size() == 0) msDirectory = ".";
      // extract the prefix and extension. Prefix is everything before any final "."
      // etension is everything after any final ".", including that final dot.
      string msPathBasename = msPath.baseName();
      size_t rdot = msPathBasename.find_last_of('.');
      if (rdot != std::string::npos) {
	msNameExtension = msPathBasename.substr(rdot,msPathBasename.size()-rdot);
	msNamePrefix = msPathBasename.substr(0,rdot);
      } else {
	msNameExtension = "";
	msNamePrefix = msPathBasename;
      }
    }
    else {
      msNamePrefix = dsPath;
      msNameExtension = ".ms";
    }
    
    // Does the user want compressed columns in the resulting MS ?
    if ((withCompression = (vm.count("compression") != 0))) {
      infostream.str("");
      infostream << "Compressed columns in the resulting MS(s) : Yes" ;
      info(infostream.str());
    }
    else {
      infostream.str("");
      infostream << "Compressed columns in the resulting MS(s) : No" ;
      info(infostream.str());
    }

    // WVR uncorrected and|or corrected data required ?
    string wvr_corrected_data = vm["wvr-corrected-data"].as<string>();
    if (wvr_corrected_data.compare("no") == 0)
      es_query_apc.fromString("AP_UNCORRECTED");
    else if (wvr_corrected_data.compare("yes") == 0)
      es_query_apc.fromString("AP_CORRECTED");
    else if (wvr_corrected_data.compare("both") == 0)
      es_query_apc.fromString("AP_CORRECTED AP_UNCORRECTED");
    else {
      errstream.str("");
      errstream << "Token '" << wvr_corrected_data << "' invalid for wvr-corrected-data." << endl;
      errstream << generic;
      error(errstream.str());
    }
    
    // Do we want an MS Main table to be filled or not ?
    mute = vm.count("dry-run") != 0;
    if (mute) {
      infostream.str("");
      infostream << "option dry-run is used, the MS Main table will not be filled" << endl;
      info(infostream.str());
    }

    // What is the amount of memory allocated to the BDF slices.
    infostream.str("");
    infostream << "the BDF slice size is set to " << bdfSliceSizeInMb << " megabytes." << endl;
    info(infostream.str());

    // Do we process in parallel ?
    doparallel = vm.count("parallel") != 0;
    if (doparallel) {
      infostream.str("");
      infostream << "run in multithreading mode" << endl;
      info(infostream.str());
    }

    lazy = false;

    // Do we want verbatim copies of some tables ?
    if (vm.count("asis")) {
      istringstream iss;
      iss.str(vm["asis"].as<string>());
      string word;
      while (iss>>word)
	asisTablenames_v.push_back(word);
    } 
  }
  catch (std::exception& e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());
  }
  // this just a dummy number for now (innthread may
  // come from user input in future..)
  // Also setting environment variable, OMP_NUM_THREADS=1
  // one can excute multiwrite part in a single thread.
  int innthread = 4;
  if (doparallel && innthread > 1) {
    doparallel = true;
  }

  //if(doparallel) cerr<<"DO PARALLEL...."<<endl;
  //
  // Try to open an ASDM dataset whose name has been passed as a parameter on the command line
  //
  if ( (dsPath.size() > 0) && dsPath.at(dsPath.size()-1) == '/' ) dsPath.erase(dsPath.size()-1);


  double cpu_time_parse_xml  = 0.0;
  double real_time_parse_xml = 0.0;
  int mode;
  mode = 0; myTimer(&cpu_time_parse_xml, &real_time_parse_xml, &mode);

  ASDM* ds_p = new ASDM();
  bool  checkRowUniqueness = vm.count("check-row-uniqueness") != 0;
  infostream.str("");

  if (checkRowUniqueness) 
    infostream << "Row uniqueness constraint will be applied." << endl;
  else
    infostream << "Row uniqueness constraint will be ignored." << endl;

  info(infostream.str());

  try {
    infostream.str("");
    infostream << "Input ASDM dataset : " << dsPath << endl;
    info(infostream.str());
    
    ds_p->setFromFile(dsPath, ASDMParseOptions().loadTablesOnDemand(true).checkRowUniqueness(checkRowUniqueness));
  }
  catch (ConversionException e) {
    errstream.str("");
    errstream << e.getMessage();
    error(errstream.str());
  }
  catch (std::exception e) {
    errstream.str("");
    errstream << e.what();
    error(errstream.str());
  }
  catch (...) {
    errstream.str("");
    errstream << "Uncaught exception !" << endl;
    error(errstream.str());
  }
  
  mode = 1; myTimer(&cpu_time_parse_xml, &real_time_parse_xml, &mode);
  infostream.str("");
  infostream << "Time spent parsing the ASDM medata : " << cpu_time_parse_xml << " s.";
  info(infostream.str());

  //
  // What are the apc literals present in the binary data.
  //
  es_apc = apcLiterals(*ds_p);
 
  //
  // Determine what kind of data complex (DATA column) or float (FLOAT_DATA) will be
  // stored in the measurement set by using the method isDataComplex on the first row of 
  // the ASDM Dataset. 
  // This method called on all remaining rows of the ASDM should return the same result
  // otherwise the filling process will stop.
  //
  SDMBinData sdmBinData(ds_p, dsPath);

  // Define the SDM Main table subset of rows to be accepted
  sdmBinData.select( es_cm, es_srt, es_ts);   

  // Define the subset of data to be extracted from the BLOBs
  sdmBinData.selectDataSubset(e_query_cm, es_query_apc);

  //
  // Selection of the scans to consider.
  //
  vector<ScanRow *>	scanRow_v	       = ds_p->getScan().get();
  map<int, set<int> > all_eb_scan_m;
  for (vector<ScanRow *>::size_type i = 0; i < scanRow_v.size(); i++)
    all_eb_scan_m[scanRow_v[i]->getExecBlockId().getTagValue()].insert(scanRow_v[i]->getScanNumber());

  vector<ScanRow *>	selectedScanRow_v;
  map<int, set<int> >   selected_eb_scan_m;

  string scansOptionInfo;
  if (vm.count("scans")) {
    string scansOptionValue = vm["scans"].as< string >();

    map<int, set<int> > eb_scan_m;
    int status = scansParser(scansOptionValue, eb_scan_m);

    if (status == 0) {
      errstream.str("");
      errstream << "'" << scansOptionValue << "' is an invalid scans selection." << endl;
      error(errstream.str());
    }

    vector<ScanRow *> scanRow_v = ds_p->getScan().get();
    map<int, set<int> >::iterator iter_m = eb_scan_m.find(-1);

    if (iter_m != eb_scan_m.end())
      for (map<int, set<int> >::iterator iterr_m = all_eb_scan_m.begin(); iterr_m != all_eb_scan_m.end(); iterr_m++)
	if ((iter_m->second).empty())
	  selected_eb_scan_m[iterr_m->first] = iterr_m->second;
	else
	  selected_eb_scan_m[iterr_m->first] = SetAndSet<int>(iter_m->second, iterr_m->second);

    for (map<int, set<int> >::iterator iterr_m = all_eb_scan_m.begin(); iterr_m != all_eb_scan_m.end(); iterr_m++) {
      if ((iter_m=eb_scan_m.find(iterr_m->first)) != eb_scan_m.end()) {
	if ((iter_m->second).empty()) {
	  selected_eb_scan_m[iterr_m->first].insert((iterr_m->second).begin(), (iterr_m->second).end());
	} else {
	  set<int> s = SetAndSet<int>(iter_m->second, iterr_m->second);
	  selected_eb_scan_m[iterr_m->first].insert(s.begin(), s.end());
	}
      }
    }

    ostringstream	oss;
    oss << "The following scans will be processed : " << endl;
    for (map<int, set<int> >::const_iterator iter_m = selected_eb_scan_m.begin(); iter_m != selected_eb_scan_m.end(); iter_m++) {
      oss << "eb#" << iter_m->first << " -> " << displaySet<int>(iter_m->second) << endl;

      Tag		execBlockTag  = Tag(iter_m->first, TagType::ExecBlock);
      for (set<int>::const_iterator iter_s = iter_m->second.begin();
	   iter_s		     != iter_m->second.end();
	   iter_s++)
	selectedScanRow_v.push_back(ds_p->getScan().getRowByKey(execBlockTag, *iter_s));

    }

    scansOptionInfo = oss.str();
  }
  else {
    selectedScanRow_v = ds_p->getScan().get();
    selected_eb_scan_m = all_eb_scan_m;
    scansOptionInfo = "All scans of all exec blocks will be processed \n";
  }

  bool	ignoreTime		= vm.count("ignore-time") != 0;
  bool	processCalDevice	= vm.count("no-caldevice") == 0;
  bool  processEphemeris	= vm.count("no-ephemeris") == 0;
  bool  processPointing		= vm.count("no-pointing") == 0;
  bool  withPointingCorrection	= vm.count("with-pointing-correction") != 0;
  bool	processSysPower		= vm.count("no-syspower") == 0;
  bool  asis                    = vm.count("asis");
  //
  // Report the selection's parameters.
  //
  infostream.str("");
  infostream << "Correlation modes requested : " << e_query_cm.str() << endl;
  infostream << "Spectral resolution types requested : " << es_srt.str() << endl;
  infostream << "Time sampling requested : " << es_ts.str() << endl;
  infostream << "WVR uncorrected and|or corrected data requested : " << es_query_apc.str() << endl;
  if (selectedScanRow_v.size() == 0) { 
    errstream.str("");
    errstream << "No scan number corresponding to your request. Can't go further.";
    error(errstream.str());
  }

  infostream << scansOptionInfo;

  if (ignoreTime)
    infostream << "All rows of the tables depending on time intervals will be processed independently of the selected exec block / scan.";
  info(infostream.str());

  infostream.str("");
  if (!processSysPower)   infostream << "The SysPower table will not be processed." << endl;
  if (!processCalDevice)  infostream << "The CalDevice table will not be processed." << endl;
  if (!processPointing)   infostream << "The Pointing table will not be processed." << endl;
  if (processPointing && withPointingCorrection ) infostream << "The correction (encoder - pointingDirection) will be applied" << endl;
  if (!processEphemeris)  infostream << "The Ephemeris table will not be processed." << endl;

  info(infostream.str());
  //
  // Shall we have Complex or Float data ?
  //
  bool complexData = true;
  try {
    complexData =  sdmBinData.isComplexData();
  }
  catch (Error & e) {
    errstream.str("");
    errstream << e.getErrorMessage();
    error(errstream.str());
  }

  //
  // Prepare a map AtmPhaseCorrection -> name of measurement set.
  // Three cases are possible :
  // only AP_CORRECTED -> MS name is suffixed with "-wvr-corrected",
  // only AP_UNCORRECTED -> MS name has no particular suffix,
  // AP_CORRECTED and AP_UNCORRECTED -> 2 MSs whith names defined by the two above rules.
  //
  map<AtmPhaseCorrection, string> msNames_m;
  if (hasCorrectedData(es_apc) && es_query_apc[AP_CORRECTED]) {
    msNames_m[AP_CORRECTED] = msNamePrefix + "-wvr-corrected";
  }
  if (hasUncorrectedData(es_apc) && es_query_apc[AP_UNCORRECTED]) {
    msNames_m[AP_UNCORRECTED] = msNamePrefix;
  }
      
  if (msNames_m.size() == 0) {
    //
    // no MS can be produced due to the selection parameters values.
    // 
    infostream.str("");
    infostream << "No measurement set can be produced with  your selection criteria on '" << dsPath << "'" << endl;
    info(infostream.str());
    delete ds_p;
    exit(1);
  }
  else {
    //
    // OK, we are going to produce at least one MS.
    // If the '--compression' option has been used then append a suffix ".compressed" to
    // the MS name.
    // And eventually always suffix with ".ms".
    //
    for (map<AtmPhaseCorrection, string>::iterator iter=msNames_m.begin(); iter != msNames_m.end(); ++iter) {
      if (withCompression)
	iter->second = iter->second + ".compressed";
      iter->second +=  msNameExtension;
    }
  }

  infostream.str("");
  infostream << "The resulting measurement set will contain a '" << ((complexData) ? "DATA" : "FLOAT_DATA") << "' column" << endl; 
  info(infostream.str());

#if DDPRIORITY
  // The data returned by getDataCols will be ordered in priority by Data Descriptions.
  sdmBinData.setPriorityDataDescription();
#endif

  //get numCorr, numChan, telescope name for setupMS
  
  ExecBlockTable& temp_execBlockT = ds_p->getExecBlock();
  //take first row of the table (assuming telescope name is all the same)
  ExecBlockRow* temp_ebtrow = temp_execBlockT.get()[0];
  string telescopeName  = temp_ebtrow->getTelescopeName();
  //cout<<"telName="<<telName<<endl;


  //Determine max numPol
  int maxNumCorr=1;
  const vector<PolarizationRow *>& polRs_v = ds_p->getPolarization().get();

  for(PolarizationRow* polR_p: polRs_v) {
    maxNumCorr=max(maxNumCorr, polR_p->getNumCorr());
  }

  //Determine max NumChan
  int maxNumChan=1;
  const vector<SpectralWindowRow *>& spwRs_v = ds_p->getSpectralWindow().get();

  for(SpectralWindowRow* spwR_p: spwRs_v) {
    maxNumChan = max(maxNumChan, spwR_p->getNumChan());
  }

  // // Create the measurement set(s). 
  // try {
  //    for (map<AtmPhaseCorrection, string>::iterator iter = msNames_m.begin(); iter != msNames_m.end(); ++iter) {
  //     info("About to create a filler for the measurement set '" + msNames_m[iter->first] + "'");
  //     msFillers[iter->first] = new ASDM2MSFiller(msNames_m[iter->first],
  // 						 0.0,
  // 						 false,
  // 						 complexData,
  // 						 withCompression,
  // 						 telescopeName,
  // 						 maxNumCorr,
  // 						 maxNumChan,
  // 						 false,
  // 						 lazy);
  //   }
  // }
  // catch(AipsError & e) {
  //   errstream.str("");
  //   errstream << e.getMesg();
  //   error(errstream.str());
  // }
  // catch (std::exception & e) {
  //   errstream.str("");
  //   errstream << e.what();
  //   error(errstream.str());
  // }
  
  // msFiller = msFillers.begin()->second;
 
  fillTopMS(ds_p,
	    dsPath,
	    msNames_m,
	    complexData,
	    withCompression,
	    telescopeName,
	    maxNumCorr,
	    maxNumChan,
	    ignoreTime,
	    selectedScanRow_v,
	    selected_eb_scan_m,
	    processCalDevice,
	    processEphemeris,
	    processPointing,
	    withPointingCorrection,
	    processSysPower,
	    asis,
	    asisTablenames_v);

  // Create as many subdirectories of the toplevel MS as there will be subMSes.
  // 
  // Partition the ASDM.
  map<Tag, AnyValueMap<string> > msGroup_avm = msGroup(ds_p);

  // Collect the configuration description ids in a vector.
  vector<Tag> keys_v;
  std::transform( msGroup_avm.begin(),
		  msGroup_avm.end(),
		  std::back_inserter(keys_v),
		  std::bind (&std::map<Tag, AnyValueMap<string> >::value_type::first, std::placeholders::_1)
		  );

  // Iterate over the vector elements.
  vector<Path> subMSpath_v;
  for(Tag configDescriptionId: keys_v) {
    ostringstream oss;
    oss << msNames_m.begin()->second << "/";
    oss << CProcessorType::toString(msGroup_avm[configDescriptionId].getValue<ProcessorTypeMod::ProcessorType>("processorType"))
	<< "_"
	<< CSpectralResolutionType::toString(msGroup_avm[configDescriptionId].getValue<SpectralResolutionTypeMod::SpectralResolutionType>("spectralResolution"));
    vector<Tag> dataDescriptionId_v = msGroup_avm[configDescriptionId].getValue<vector<Tag> >("ddId");
    for (Tag dataDescriptionId: dataDescriptionId_v ) {
      subMSpath_v.push_back(Path(oss.str()+"_"+dataDescriptionId.toString()));
    }
  }
  
  for(Path path: subMSpath_v) {
    File subMSfile(path);
    if (subMSfile.isDirectory()) {
      Directory subMSdir(path);
      int nRemoved = subMSdir.nEntries()+1; // count itself
      subMSdir.removeRecursive();
      cout << "removing directory " << path.originalName() << " (" << nRemoved << " files deleted.)" << endl;
    } else {
      cout << "directory " << path.originalName() << " does not exist." << endl;
    }
    
    Directory subMSdir(path);
    subMSdir.create();
    cout << "successfully created directory " << path.originalName() << endl;
  }

  //
  // Let's create symbolic links to the toplevel MS subtables in each subMS directory.
  //
  vector<Path> subTablePath_v;
  Directory msDir(msNames_m.begin()->second);
  DirectoryIterator msDirIter(msDir);
  while (!msDirIter.pastEnd()) {
    if (msDirIter.file().isDirectory()) {
      Path tabPath(msDirIter.file().path());
      tabPath.append("table.dat");
      File tabDat(tabPath);
      if (tabDat.exists()) {
	subTablePath_v.push_back(msDirIter.file().path());
	cout << subTablePath_v.back().originalName() << endl;
      }
      msDirIter++;
    }
  }

  for(Path subMSPath: subMSpath_v) { 
    for(Path subTablePath: subTablePath_v) {
      Path newSubTabPath(subMSPath);
      newSubTabPath.append(subTablePath.baseName());
      SymLink(newSubTabPath).create(subTablePath);
    }
  }

  //
  // Let's copy the toplevel MS main table into each subMS directory.
  //
  
  //
  // Determine the MS main table files.
  DirectoryIterator dirIter(Path(msNames_m.begin()->second));
  while (!dirIter.pastEnd()) {
    string filename = dirIter.name();
    string tabPrefix = "table";
    if (filename.compare(0, tabPrefix.size(), tabPrefix)) {
      for(Path subMSPath: subMSpath_v) {
	cout << "Trying to copy " << dirIter.file().path().originalName() << " to " << subMSPath.originalName() << endl;
	RegularFile origFile(dirIter.file());
	Path newCopyPath(subMSPath);
	newCopyPath.append(filename);
	origFile.copy(newCopyPath);
      }
    }
    dirIter++;
  }
}
