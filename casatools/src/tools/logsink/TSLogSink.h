//# TSLogSink.h: Save distributed log messages
//# Copyright (C) 2005
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#ifndef CASA_TSLOGSINK_H
#define CASA_TSLOGSINK_H

//# Includes
#include <casa/aips.h>
#include <casa/Logging/LogSink.h>
#include <casa/Containers/Block.h>
#include <casa/BasicSL/String.h>

#include <fstream>

namespace casacore{

	class LogSinkInterface;
	class LogFilterInterface;
}

namespace casa { //# NAMESPACE CASA - BEGIN

//# Forward Declarations
#ifndef AIPS_LOG4CPLUS
#endif

// <summary>
// Save distributed log messages.
// </summary>

// <use visibility=export>

// <reviewed reviewer="UNKNOWN" date="before2004/08/25" tests="tLogging.cc" demos="dLogging.cc">
// </reviewed>

// <prerequisite>
//   <li> <linkto class=casacore::LogSink>LogSink</linkto>
// </prerequisite>
//
// <synopsis>
// This class posts messages which pass the filter to 
// a distributed logger.
// </synopsis>
//
// <example>
// See <linkto file="Logging.h">Logging.h</linkto>.
// </example>
//
// <motivation>
// For a distributed systems log messages.
// </motivation>
//
//# <todo asof="2005/08/15">
//# </todo>

class TSLogSink : public casacore::LogSinkInterface
{
public:
  // Create an empty sink without a filter.
    TSLogSink(const std::string &name);

  // Create an empty sink with the given filter.
  // <group>
  explicit TSLogSink (casacore::LogMessage::Priority filter);
  explicit TSLogSink (const casacore::LogFilterInterface& filter);
  // </group>

  // Copy constructor (copy semantics).
  TSLogSink (const TSLogSink& other);

  // Assignment (copy semantics).
  TSLogSink& operator= (const TSLogSink& other);
  
  virtual ~TSLogSink();

  // If the message passes the filter, write it
  virtual casacore::Bool postLocally (const casacore::LogMessage& message);

  // Returns the id for this class...
  static casacore::String localId( );
  // Returns the id of the casacore::LogSink in use...
  casacore::String id( ) const;
  const casacore::LogFilterInterface& filter() const;
  casacore::LogSinkInterface& filter(const casacore::LogFilterInterface &newfilter);
#ifndef AIPSLOG4CPLUS
  void cerrToo(casacore::Bool);
  void setLogSink(casacore::String logname = "");
#endif

private:
  // Avoid duplicating code in copy ctor and assignment operator
  void copy_other (const TSLogSink& other);
#ifndef AIPS_LOG4CPLUS
  casacore::LogSinkInterface *logsink;
  casacore::Bool send2cerr;
  std::ofstream *logfile;
#endif

};



} //# NAMESPACE CASA - END

#endif
