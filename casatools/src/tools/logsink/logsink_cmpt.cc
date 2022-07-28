
/***
 * Framework independent implementation file for logsink...
 *
 * Implement the logsink component here.
 *
 * // TODO: WRITE YOUR DESCRIPTION HERE!
 *
 * @author
 * @version
 ***/

#include <iostream>
#include <fstream>
#include <logsink_cmpt.h>
#include <casacore/casa/Logging/LogSink.h>
#include <casacore/casa/Logging/LogFilter.h>
#include <parallel/Logging/LogFilterParallel.h>
#include <TSLogSink.h>
#include <casacore/casa/Logging/NullLogSink.h>
#include <casacore/casa/Logging/StreamLogSink.h>
#include <stdcasa/version.h>
#include <unistd.h>
#include <sys/param.h>
#include <CasapyWatcher.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/BasicMath/Random.h>
#include <casacore/casa/System/Aipsrc.h>
#include <casacore/casa/OS/HostInfo.h>
#include <chrono>
#include <ctime>
#include <mutex>
#include <sys/file.h>

#include <casatools/Config/State.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>

using namespace std;
using namespace casacore;
using namespace casa;

using namespace casacore;
namespace casac {

	/*
static ofstream *logfile = new ofstream("logfile");
static LogSinkInterface *thelogsink
  = new StreamLogSink(logfile);
  */
static string theLogName;
static string telemetryLog;
static bool telemetryEnabled;
string telemetryLoggerPid;

//logsink::logsink():thelogsink(0)
logsink::logsink(const std::string &filename, bool telemetrytoggle, const std::string &inTelemetryfilename)
{
  if( ! theLogName.size( ) ){
     char *buff = NULL;
     if ( filename.at(0) == '/' )
         theLogName = filename;
     else {
         char *mybuff = getcwd(buff, MAXPATHLEN);
         theLogName = string(mybuff) + "/" + filename;
     }
     telemetryEnabled = telemetrytoggle;
  }

  if (telemetryEnabled) {
     setstatslogfile(inTelemetryfilename);
  }

  // jagonzal: Set task and processor name
  taskname = "casa";
  processor_name = "casa";

  //cout << "thelogsink=" << thelogsink << endl;
  thelogsink = new casa::TSLogSink(theLogName);
  setlogfile(theLogName);
  itsorigin = new LogOrigin("casa");
  logLevel =  LogMessage::NORMAL;
  thelogsink->postLocally(LogMessage("", *itsorigin, LogMessage::NORMAL));
  globalsink = false;
  logname = theLogName ;
  filterMsgList.clear();
}

std::string logsink::version(){
  if(!thelogsink){
       thelogsink = &LogSink().globalSink();
  }
  ostringstream os1, os2;
  os1 << "CASA Version ";
  casa::VersionInfo::report(os1);
  std::string mymess = os1.str();

  thelogsink->postLocally(LogMessage(mymess, *itsorigin, LogMessage::NORMAL));
  thelogsink->postLocally(LogMessage(os2.str(), *itsorigin, LogMessage::NORMAL));
  return mymess;
}

logsink::~logsink()
{
	delete itsorigin;
	if(!globalsink)
	   delete thelogsink;
	thelogsink=0;
	itsorigin=0;
}

bool logsink::origin(const std::string &fromwhere)
{
    bool rstat(true);
    if(!thelogsink){
       thelogsink = &LogSink().globalSink();
    }
    delete itsorigin;
    itsorigin = new LogOrigin(processor_name);
    //String* taskname = new String(fromwhere);
    taskname = fromwhere;
    itsorigin->taskName(taskname);
    thelogsink->setTaskName(taskname);
    LogSink().globalSink().setTaskName(taskname);
    return rstat;
}

// jagonzal: Allow to set the processor origin (normally "casa"
// but in the MPI case we use the hostname and rank involved)
bool logsink::processorOrigin(const std::string &fromwhere)
{
    bool rstat = true;

    processor_name = fromwhere;

    return rstat;
}

bool logsink::filter(const std::string &level)
{
	bool rstat = true;

	logLevel = getLogLevel(level);

	installLogFilter();

	return rstat;
}

void logsink::filterMsg(const std::vector<std::string>& msgList)
{
	// By default SWIG create a vector of one empty string
	if ((msgList.size() >= 1) and (msgList[0].size() > 0))
	{
		for(std::vector<std::string>::const_iterator  it = msgList.begin(); it != msgList.end(); ++it)
		{
			if (it->size() > 0) filterMsgList.push_back(*it);
		}

		installLogFilter();
	}

	return;
}

void logsink::clearFilterMsgList()
{
	filterMsgList.clear();
	installLogFilter();

	return;
}

void logsink::installLogFilter()
{
    if (!thelogsink)
    {
       thelogsink = &LogSink().globalSink();
    }

	// Use a LogFilterParallel to filter out the specified messages
	if (filterMsgList.size() >= 1)
	{
		LogFilterParallel filter(logLevel);
		for(std::vector<std::string>::const_iterator  it = filterMsgList.begin(); it != filterMsgList.end(); ++it)
		{
			if (it->size() > 0) filter.filterOut(it->c_str());
		}

		thelogsink->filter(filter);
	}
	// Normal operation filtering only on priority
	else
	{
		LogFilter filter(logLevel);
		thelogsink->filter(filter);
	}

	// Also set for any watchers.
	CasapyWatcher::logChanged_(logLevel);
}

LogMessage::Priority logsink::getLogLevel(const std::string &level)
{
	LogMessage::Priority priority = LogMessage::NORMAL;

	if (level == "DEBUG" || level == "DEBUGGING")
		priority = LogMessage::DEBUGGING;
	else if (level == "DEBUG1")
		priority = LogMessage::DEBUG1;
	else if (level == "DEBUG2")
		priority = LogMessage::DEBUG2;
	else if (level == "NORMAL")
		priority = LogMessage::NORMAL;
	else if (level == "NORMAL1")
		priority = LogMessage::NORMAL1;
	else if (level == "NORMAL2")
		priority = LogMessage::NORMAL2;
	else if (level == "NORMAL3")
		priority = LogMessage::NORMAL3;
	else if (level == "NORMAL4")
		priority = LogMessage::NORMAL4;
	else if (level == "NORMAL5")
		priority = LogMessage::NORMAL5;
	else if (level == "INFO" || level == "NORMAL")
		priority = LogMessage::NORMAL;
	else if (level == "INFO1")
		priority = LogMessage::NORMAL1;
	else if (level == "INFO2")
		priority = LogMessage::NORMAL2;
	else if (level == "INFO3")
		priority = LogMessage::NORMAL3;
	else if (level == "INFO4")
		priority = LogMessage::NORMAL4;
	else if (level == "INFO5")
		priority = LogMessage::NORMAL5;
	else if (level == "WARN")
		priority = LogMessage::WARN;
	else if (level == "ERROR" || level == "SEVERE")
		priority = LogMessage::SEVERE;

	return priority;
}

bool logsink::post(const std::string& message,
		   const std::string& priority,
		   const std::string& origin)
{
	return postLocally(message, priority, origin);
}

bool
logsink::postLocally(const std::string& message,
                     const std::string& priority,
                     const std::string& origin)
{
    if(!thelogsink) thelogsink = &LogSink().globalSink();

    if(!itsorigin) itsorigin = new LogOrigin(processor_name);

    // Set origin
    itsorigin->className(origin);
    thelogsink->setTaskName(origin);

    // Get priority
    LogMessage::Priority messagePriority = getLogLevel(priority);

    // Recursive call to LogSinkInterface::postLocally
    return thelogsink->postLocally(LogMessage(message,*itsorigin,messagePriority));
}

std::mutex _stat_mutex;
bool logsink::poststat(const std::string& message,
                       const std::string& /*origin*/)
{
    int fd = open(telemetryLog.c_str(), O_RDWR | O_CREAT, 0666);
    int telemetryFileLocked = flock(fd, LOCK_EX | LOCK_NB);
    if (!telemetryFileLocked) {
        if (telemetryEnabled ) {
            try {
                std::lock_guard<std::mutex> lock(_stat_mutex);
                ofstream tlog;
                tlog.open (telemetryLog, ios_base::app);
                auto timestamp = std::chrono::system_clock::now();
                std::time_t convertedtime = std::chrono::system_clock::to_time_t(timestamp);
                char formattedtime[20];
                strftime(formattedtime, 20, "%Y-%m-%d %H:%M:%S", localtime(&convertedtime));
                std::string telemetryOrigin =" ";
                if(itsorigin) {
                    std::string origOrigin = itsorigin->toString();
                    std::remove_copy(origOrigin.begin(), origOrigin.end(), std::back_inserter(telemetryOrigin), ':');
                }
                tlog << formattedtime << " :: " << telemetryLoggerPid<< " ::"<< telemetryOrigin << " :: " << message << "\n";
                tlog.close();
            } catch (const std::exception& e) { // caught by reference to base
                std::cout << "Writing stats failed: '" << e.what() << "'\n";
            }
        }
        //else {
        //    cout << "Telemetry is disabled. Won't write stats.";
        //}
    }
    else {
        cout << "Telemetry file is locked by another process. Won't write stats.";
    }
    flock(fd, LOCK_UN);
    close(fd);
    return true;
}

std::string
logsink::localId()
{
  return thelogsink->localId();
}

std::string
logsink::id()
{
  return thelogsink->id();
}

bool logsink::setglobal(const bool isglobal)
{
   bool rstat(true);

   if(isglobal){
      if (globalsink==isglobal && thelogsink)
         return rstat;
      LogSink().globalSink(thelogsink);
      globalsink = isglobal;
   } //else {
   //   LogSinkInterface *dummy = new StreamLogSink(LogMessage::NORMAL, &cerr);
   //   LogSink().globalSink(dummy);
   //}
   return rstat;
}
//
bool logsink::setlogfile(const std::string& filename)
{
   bool rstat(true);

   string tmpname = filename ;
   if(!filename.size()){
      String logfileKey="logfile.no.default";
      String logname2;
      if(!Aipsrc::find(logname2, logfileKey)){
         tmpname = "casa.log";
      } else {
         //tmpname = logname2;
         ACG g(7326458, 98);
         tmpname = "/tmp/"+String::toString(g.asuInt());
         //tmpname = "null";
      }
   }
   if(tmpname != "null") {
      casacore::File filein( tmpname ) ;
      logname = filein.path().absoluteName() ;
      if ( ! thelogsink )
          thelogsink = &LogSink().globalSink();
      static_cast<TSLogSink*>(thelogsink)->setLogSink(logname);
   }
   else
      thelogsink = new NullLogSink();

   //
   // Also set for any watchers.
   CasapyWatcher::logChanged_(logname);
   casatools::get_state( ).setLogPath(logname);

   return rstat;
}

bool logsink::setstatslogfile(const std::string& filename)
{
  //cout << "Setting stats file to " << filename << "\n";
  char *buff = NULL;
  if ( filename.at(0) == '/' )
     telemetryLog = filename;
  else {
     char *mybuff = getcwd(buff, MAXPATHLEN);
     telemetryLog = string(mybuff) + "/" + filename;
  }
  telemetryLoggerPid = std::to_string(getpid());
  return true;
}

string logsink::getstatslogfile()
{
   return telemetryLog ;
}


bool
logsink::showconsole(const bool onconsole)
{
   bool rstat(true);
   if(thelogsink)
      thelogsink->cerrToo(onconsole);
   return rstat;
}
std::string
logsink::logfile()
{
  //return theLogName;
  return logname ;
}

long
logsink::ompNumThreadsTest()
{
	Int num_thread(-1);

	#ifdef _OPENMP
		num_thread = 0;
		#pragma omp parallel reduction( + : num_thread )
		{
			num_thread = num_thread + 1;
		}
	#endif

	return num_thread;
}

long
logsink::ompGetNumThreads()
{
	Int num_threads(-1);

	#ifdef _OPENMP
		#pragma omp parallel
		{
			#pragma omp master
			{
				num_threads = omp_get_num_threads();
			}
		}
	#endif

	return num_threads;
}

bool
logsink::ompSetNumThreads(long numThreads)
{
	Int res(true);

	#ifdef _OPENMP
		omp_set_num_threads(numThreads);
	#endif

	return res;
}

long
logsink::setMemoryTotal(long memory)
{
	return HostInfo::setMemoryTotal(memory);
}

long
logsink::setMemoryFraction(long memfrac)
{
	return HostInfo::setMemoryFraction(memfrac);
}

long
logsink::setNumCPUs(long numCPUs)
{
	return HostInfo::setNumCPUs(numCPUs);
}

long
logsink::getMemoryTotal(bool use_aipsrc)
{
	return HostInfo::memoryTotal(use_aipsrc);
}

long
logsink::getNumCPUs(bool use_aipsrc)
{
	return HostInfo::numCPUs(use_aipsrc);
}



} // casac namespace
