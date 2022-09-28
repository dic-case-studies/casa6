
/***
 * Framework independent implementation file for utils...
 *
 * Implement the utils component here.
 *
 * // TODO: WRITE YOUR DESCRIPTION HERE!
 *
 * @author
 * @version
 ***/

#include <iostream>
#include <fstream>
#include <stdcasa/record.h>
#include <stdcasa/version.h>
#include <utils_cmpt.h>
#include <tools/utils/stdBaseInterface.h>
#include <climits>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/OS/DOos.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableUtil.h>
#include <casacore/casa/System/Aipsrc.h>
#include <casacore/casa/OS/HostInfo.h>
#ifndef NO_CRASH_REPORTER
#include <stdcasa/StdCasa/CrashReporter.h>
#endif
#include <stdlib.h>
#include <signal.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <casacore/casa/Quanta/UnitMap.h>
#include <casatools/Config/State.h>
#include <fftw3.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <casacore/scimath/Mathematics/FFTW.h>
#include <asdmstman/AsdmStMan.h>
#include <casacore/derivedmscal/DerivedMC/Register.h>
#include <toolversion.h>

using namespace std;
using namespace casacore;
using namespace casa;

using namespace casacore;
namespace casac {

utils::utils()
{
  myConstraints = 0;
  itsLog = new casacore::LogIO;
}

utils::~utils()
{
  if(myConstraints)
     delete myConstraints;
  delete itsLog;
}

std::string
utils::getrc(const std::string& rcvar)
{
  String rstat1;
  if(!rcvar.length()){
	  rstat1 = Aipsrc::aipsRoot();
  } else {
	  if(!Aipsrc::find(rstat1, rcvar))
		  rstat1 = "Unknown value";
  }
  string rstat(rstat1.c_str());
  return rstat;
}

bool
utils::removetable(const std::vector<std::string> &tablenames)
{
  bool rstat(true);
  try {
     *itsLog << LogOrigin("utils", "removetable");
     for(vector<std::string>::const_iterator iter = tablenames.begin();
		     iter != tablenames.end(); iter++){
       String fileName(*iter);
       if (fileName.empty()) {
          *itsLog << LogIO::WARN << "Empty filename" << LogIO::POST;
          rstat = false;
       }
       File f(fileName);
       if (! f.exists()) {
           *itsLog << LogIO::WARN << fileName << " does not exist." << LogIO::POST;
          rstat = false;
       }

// Now try and blow it away.  If it's open, tabledelete won't delete it.
       String message;
       if(rstat && Table::isReadable(fileName)){
          if (TableUtil::canDeleteTable(message, fileName, true)) {
             TableUtil::deleteTable(fileName, true);
          } else {
             *itsLog << LogIO::WARN << "Cannot delete file " << fileName
             << " because " << message << LogIO::POST;
          }
       } else {
           *itsLog << LogIO::WARN << "Cannot delete file " << fileName
           << " because it's not a table." << LogIO::POST;
       }
  }
    } catch (AipsError x) {
       *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
       << LogIO::POST;
       RETHROW(x);
  }
  return rstat;
}

::casac::record *utils::tableinfo(const std::string &tablename) {
    Vector<Int> info = casacore::DOos::lockInfo(tablename);
    ::casac::record *result = new record( );

    switch( info[0] ) {
    case 3:
	result->insert("lockstatus", "write");
	break;
    case 2:
	result->insert("lockstatus", "read");
	break;
    case 1:
	result->insert("lockstatus", "open");
	break;
    case 0:
	result->insert("lockstatus", "not in use");
	break;
    default:
	result->insert("lockstatus", "unknown");
    }

    result->insert("lockpid", (long) info[1]);
    result->insert("lockperm", info[2] ? true : false);
    return result;
}

std::vector<std::string> utils::lockedtables( ) {
    Vector<String> locks = Table::getLockedTables( );
    std::vector<std::string> result;
    for (unsigned int x = 0; x < locks.nelements(); ++x ) {
	result.push_back(locks[x]);
    }
    return result;
}


typedef long NUMCAST;
::casac::record *utils::hostinfo( ) {
    ::casac::record *result = new record( );

    ::casac::record *swap = new record( );
    swap->insert( "total", (NUMCAST) HostInfo::swapTotal( ) );
    swap->insert( "used", (NUMCAST) HostInfo::swapUsed( ) );
    swap->insert( "free", (NUMCAST) HostInfo::swapFree( ) );
    result->insert( "swap", swap );

    ::casac::record *memory = new record( );
    memory->insert( "total", (NUMCAST) HostInfo::memoryTotal( ) );
    memory->insert( "available", (NUMCAST) HostInfo::memoryTotal(true) );
    memory->insert( "used", (NUMCAST) HostInfo::memoryUsed( ) );
    memory->insert( "free", (NUMCAST) HostInfo::memoryFree( ) );
    result->insert( "memory", memory );

    ::casac::record *cpus = new record( );
    cpus->insert( "total", (NUMCAST) HostInfo::numCPUs( ) );
    cpus->insert( "available", (NUMCAST) HostInfo::numCPUs(true) );
    result->insert( "cpus", cpus );

    result->insert( "endian", HostInfo::bigEndian( ) ? "big" : "little" );
    result->insert( "hostname", HostInfo::hostName( ) );
    result->insert( "pid", (NUMCAST) HostInfo::processID( ) );

    result->insert( "seconds", (NUMCAST) HostInfo::secondsFrom1970( ) );

    return result;
}

std::string
utils::c_exception ()
{
  String lastMessage, lastStackTrace;
  AipsError::getLastInfo (lastMessage, lastStackTrace);

  String result = lastMessage + "\n" + lastStackTrace;

  return result;
}

void
utils::c_exception_clear ()
{
  AipsError::clearLastInfo ();
}

void bogusHandler (int, siginfo_t *, void *)
{
    // Do nothing
}

string
utils::_crash_reporter_initialize (const string & crashDirectory,
                                   const string & crashPosterApplication,
                                   const string & crashPostingUrl,
				   const string & logFile)
{
#ifndef NO_CRASH_REPORTER
    // *NOTE*: Not intended for casual use!

    string status = casa::CrashReporter::initialize(crashDirectory, crashPosterApplication,
                                                    crashPostingUrl, logFile);

    return status;
#else
    return "no-op";
#endif
}

bool
utils::_trigger_segfault (long faultType)
{
    // *NOTE*: Not intended for casual use!

    switch (faultType) {

    case 0:{
	bool * p;
        long zero = 0;
	p = (bool *) zero;
	return * p;
	break;
    }

    default:
    case 1:{
	throw exception();
	break;
    }

    }

    return false;
}

// ------------------------------------------------------------
// -------------------- initialize CASAtools ------------------

static std::vector<std::string> default_data_path;
static std::string python_path;

long utils::maxint( ) { return INT_MAX; }
long utils::minint( ) { return INT_MIN; }
long utils::maxlong( ) { return LONG_MAX; }
long utils::minlong( ) { return LONG_MIN; }

double utils::tryit(const ::casac::record &input) {
    auto found = input.find("value");
    fprintf(stderr,"\t\t<<%s>>\n",(*found).second.typeString( ).c_str( ));
    return (*found).second.toDouble( );
}

// CASA 6
bool utils::initialize( const std::string &pypath, 
                        const std::string &distro_data,
                        const std::vector<std::string> &default_path,
                        bool nogui,
                        bool agg,
                        bool pipeline) {
    static bool initialized = false;
    if ( initialized ) return false;
    default_data_path = default_path;
    python_path = pypath;
    casatools::get_state( ).setDataPath(default_data_path);
    casatools::get_state( ).setDistroDataPath(distro_data);
    casatools::get_state( ).setPythonPath(python_path);
    casatools::get_state( ).setNoGui(nogui);
    casatools::get_state( ).setAgg(agg);
    casatools::get_state( ).setPipeline(pipeline);
    // configure quanta/measures customizations...
    UnitMap::putUser( "pix", UnitVal(1.0), "pixel units" );

    casa::AsdmStMan::registerClass( );
    register_derivedmscal();

    // --- --- --- configure fftw --- CAS-13342 --- --- --- --- --- --- --- --- 
    casacore::FFTW init_casacore_fftw;
#ifdef _OPENMP
    int numCPUs = omp_get_max_threads();
#else
    int numCPUs = HostInfo::numCPUs();
#endif
    int nthreads = 1;
    if (numCPUs > 1) {
        nthreads = numCPUs;
    }
    fftwf_plan_with_nthreads(nthreads);
    fftw_plan_with_nthreads(nthreads);
    initialized = true;
    return true;
}

// ------------------------------------------------------------
// -------------------- handling rundata path -----------------
std::string utils::rundata( ) {
    return casatools::get_state( ).measuresDir( );
}

void utils::setrundata( const std::string &data ) {
    casatools::get_state( ).setDistroDataPath(data);
}

// ------------------------------------------------------------
// -------------------- handling data path --------------------
std::vector<std::string> utils::defaultpath( ) {
    return default_data_path;
}

bool utils::setpath(const std::vector<std::string> &dirs) {
    casatools::get_state( ).setDataPath(dirs);
    return casatools::get_state( ).dataPath( ).size( ) == dirs.size( );
}

std::vector<std::string> utils::getpath( ) {
    std::vector<std::string> result;
    const std::list<std::string> &path = casatools::get_state( ).dataPath( );
    std::copy( path.begin( ), path.end( ), std::back_inserter(result) );
    return result;
}

std::string utils::getpython( ) {
    return casatools::get_state( ).pythonPath( );
}

void utils::clearpath( ) {
    casatools::get_state( ).clearDataPath( );
}

std::string utils::resolve(const std::string &subdir) {
    return casatools::get_state( ).resolve(subdir);
}
// ------------------------------------------------------------

// ------------------------------------------------------------
// -------------- handling service registry -------------------
::casac::record *utils::registry( ) {
    casac::record *regrec = new casac::record;
    regrec->insert("uri",casatools::get_state( ).registryURI( ));
    return regrec;
}

::casac::record *utils::services( ) {
    std::list<casatools::ServiceId> servs = casatools::get_state( ).services( );
    casac::record *regrec = new casac::record;
    unsigned int count = 1;
    for ( std::list<casatools::ServiceId>::const_iterator it=servs.begin( ); it != servs.end( ); ++it ) {
        casac::record *sub = new casac::record;
        sub->insert("id",it->id( ));
        sub->insert("uri",it->uri( ));
        sub->insert("types",std::vector<std::string>(it->types( ).begin( ),it->types( ).end( )));
        sub->insert("priority",(NUMCAST) it->priority( ));
        regrec->insert(std::to_string(count++),sub);
    }
    return regrec;
}

bool utils::remove_service( const std::string& uri ) {
    std::list<casatools::ServiceId> servs = casatools::get_state( ).services( );
    for ( std::list<casatools::ServiceId>::const_iterator it=servs.begin( ); it != servs.end( ); ++it ) {
        if ( uri == it->uri( ) ) {
            *itsLog << LogOrigin("utils","remove_service") <<
                "removing service " <<
                it->uri( ) << "/" << it->id( ) << std::endl;
            auto id = it->id( );
            return casatools::get_state( ).removeService(id);
        }
    }
    return false;
}

void utils::shutdown( ) {
    casatools::get_state( ).shutdown( );
    // this will result in the deletion of casacore state object
    casacore::AppStateSource::initialize(0);
}

// ------------------------------------------------------------

std::vector<long>
utils::version( ) {
    std::vector<long> result = {
        VersionInfo::major( ),
        VersionInfo::minor( ),
        VersionInfo::patch( ),
        VersionInfo::feature( )
    };
    return result;
}

std::string
utils::version_desc( ) { return VersionInfo::desc( ); }

std::string
utils::version_variant( ) { return VersionInfo::variant( ); }


std::string
utils::version_info( ) { return VersionInfo::info( ); }

std::string
utils::version_string( ) { return VersionInfo::str( ); }

bool utils::compare_version(const  string& comparitor,  const std::vector<long>& vec) {
    return VersionInfo::compare(comparitor,vector<int>(vec.begin(),vec.end()));
}

std::vector<long>
utils::toolversion( ) {
    std::vector<long> result = {
        ToolVersionInfo::major( ),
        ToolVersionInfo::minor( ),
        ToolVersionInfo::patch( ),
        ToolVersionInfo::feature( ),
    };
    return result;
}

std::string
utils::toolversion_string( ) {
    return ToolVersionInfo::version( );
}

// ------------------------------------------------------------
// -------------------- Other configuration params ------------

bool utils::getnogui( ) {
    return casatools::get_state( ).noGui( );
}
bool utils::getagg( ) {
    return casatools::get_state( ).agg( );
}
bool utils::getpipeline( ) {
    return casatools::get_state( ).pipeline( );
}

} // casac namespace
