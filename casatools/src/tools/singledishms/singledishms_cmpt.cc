/***
 * Tool bindings of the test single dish module that work on an MS
 *
 * @author kana
 * @version
 ***/
#include <singledishms_cmpt.h>
#include <string>
#include <iostream>

#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogOrigin.h>
#include <casacore/casa/Exceptions/Error.h>

#include <singledish/SingleDish/SingleDishMS.h>

#include <casacore/casa/namespace.h> // using casa namespace

using namespace std;

#define _ORIGIN LogOrigin("singledishms", __func__, WHERE)

using namespace casacore;
using namespace casa;

namespace casac {

singledishms::singledishms()
{
  itsSd = 0;
  itsLog = new LogIO();
}

singledishms::~singledishms()
{
  if (itsSd != 0) delete itsSd;
  delete itsLog;
}

bool
singledishms::open(string const& ms_name)
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    // In case already open, close it!
    close();
    // create instanse
    itsSd = new SingleDishMS(ms_name);
    if (itsSd != 0) rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::close()
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    if(itsSd != 0) delete itsSd;
    itsSd = 0;
    rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::done()
{
   return close();
}

void
singledishms::assert_valid_ms()
{
  if (itsSd == 0)
    throw(AipsError("No MeasurementSet has been assigned, please run open."));
}

Record
singledishms::get_time_averaging_record(const bool& timeaverage,
					const string& timebin,
					const string& timespan)
{
      Record average_param;
      average_param.define("timeaverage", timeaverage);
      if (timeaverage) {
	String average_string;
	average_string = toCasaString(timebin);
	if (average_string != "") {
	  average_param.define("timebin", average_string);
	}
	average_string = toCasaString(timespan);
	if (average_string != "") {
	  average_param.define("timespan", average_string);
	}
      }
      return average_param;
}

string
singledishms::name()
{
  *itsLog << _ORIGIN;
  try {
    assert_valid_ms();
    return itsSd->name();
  } catch (AipsError x) {
    *itsLog << LogIO::SEVERE << "MS is not yet assigned." << LogIO::POST;
    RETHROW(x);
  }
  return "";
}

bool
singledishms::subtract_baseline(string const& datacolumn,
                                string const& outfile,
                                string const& bloutput,
                                bool const dosubtract,
                                ::casac::variant const& spw,
                                bool const updateweight,
                                string const& sigmavalue,
                                string const& blfunc,
                                long const order,
                                float const clip_threshold_sigma,
                                long const num_fitting_max,
                                bool const linefinding,
                                float const threshold,
                                long const avg_limit,
                                long const minwidth,
                                vector<long> const& edge)
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    assert_valid_ms();
    itsSd->subtractBaseline(datacolumn,
                            outfile,
                            bloutput,
                            dosubtract,
                            toCasaString(spw),
                            updateweight,
                            sigmavalue,
                            blfunc,
                            order,
                            clip_threshold_sigma,
                            num_fitting_max,
                            linefinding,
                            threshold,
                            avg_limit,
                            minwidth,
                            vector<int>(edge.begin(),edge.end()));
    rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::subtract_baseline_cspline(string const& datacolumn,
                                        string const& outfile,
                                        string const& bloutput,
                                        bool const dosubtract,
                                        ::casac::variant const& spw,
                                        bool const updateweight,
                                        string const& sigmavalue,
                                        long const npiece,
                                        float const clip_threshold_sigma,
                                        long const num_fitting_max,
                                        bool const linefinding,
                                        float const threshold,
                                        long const avg_limit,
                                        long const minwidth,
                                        vector<long> const& edge)
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    assert_valid_ms();
    itsSd->subtractBaselineCspline(datacolumn,
                                   outfile,
                                   bloutput,
                                   dosubtract,
                                   toCasaString(spw),
                                   updateweight,
                                   sigmavalue,
                                   npiece,
                                   clip_threshold_sigma,
                                   num_fitting_max,
                                   linefinding,
                                   threshold,
                                   avg_limit,
                                   minwidth,
                                   vector<int>(edge.begin(),edge.end()));
    rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::subtract_baseline_sinusoid(string const& datacolumn,
                                         string const& outfile,
                                         string const& bloutput,
                                         bool const dosubtract,
                                         ::casac::variant const& spw,
                                         bool const updateweight,
                                         string const& sigmavalue,
                                         string const& addwn,
                                         string const& rejwn,
                                         bool const applyfft,
                                         string const& fftmethod,
                                         ::casac::variant const& fftthresh,
                                         float const clip_threshold_sigma,
                                         long const num_fitting_max,
                                         bool const linefinding,
                                         float const threshold,
                                         long const avg_limit,
                                         long const minwidth,
                                         vector<long> const& edge)
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    assert_valid_ms();
    itsSd->subtractBaselineSinusoid(datacolumn,
                                    outfile,
                                    bloutput,
                                    dosubtract,
                                    toCasaString(spw),
                                    updateweight,
                                    sigmavalue,
                                    addwn,
                                    rejwn,
                                    applyfft,
                                    fftmethod,
                                    toCasaString(fftthresh),
                                    clip_threshold_sigma,
                                    num_fitting_max,
                                    linefinding,
                                    threshold,
                                    avg_limit,
                                    minwidth,
                                    vector<int>(edge.begin(),edge.end()));
    rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported1: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::subtract_baseline_variable(string const& datacolumn,
                                         string const& outfile,
                                         string const& bloutput,
                                         bool const dosubtract,
                                         ::casac::variant const& spw,
                                         bool const updateweight,
                                         string const& sigmavalue,
                                         string const& blparam,
					 bool const verbose)
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    assert_valid_ms();
    itsSd->subtractBaselineVariable(datacolumn,
                                    outfile,
                                    bloutput,
                                    dosubtract,
                                    toCasaString(spw),
                                    updateweight,
                                    sigmavalue,
                                    blparam,
				    verbose);
    rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::apply_baseline_table(string const& bltable,
				   string const& datacolumn,
				   ::casac::variant const& spw,
                                   bool const updateweight,
                                   string const& sigmavalue,
				   string const& outfile)
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    assert_valid_ms();
    itsSd->applyBaselineTable(datacolumn,
			      bltable,
			      toCasaString(spw),
			      updateweight,
			      sigmavalue,
			      outfile);
    rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::fit_line(string const& datacolumn,
		       ::casac::variant const& spw,
		       ::casac::variant const& pol,
                       string const& timebin,
                       string const& timespan,
                       string const& polaverage,
		       string const& fitfunc,
		       string const& nfit,
		       bool const linefinding,
		       float const threshold,
		       long const avg_limit,
		       long const minwidth,
		       vector<long> const& edge,
		       string const& tempfile,
		       string const& tempoutfile)
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    assert_valid_ms();
    if (timebin != "") {
      Record average_param = get_time_averaging_record(true, timebin,
						       timespan);
      itsSd->setAverage(average_param);
    }

    if (polaverage != "") {
      Record average_param;
      average_param.define("polaverage", true);
      average_param.define("polaveragemode", polaverage);
      itsSd->setPolAverage(average_param, True);
    }

    itsSd->fitLine(datacolumn, toCasaString(spw), toCasaString(pol),
		   fitfunc, nfit, linefinding, threshold, avg_limit,
		   minwidth, vector<int>(edge.begin(),edge.end()), tempfile, tempoutfile);
    rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::set_selection(::casac::variant const& spw,
		    ::casac::variant const& field,
		    ::casac::variant const& antenna,
		    ::casac::variant const& timerange,
		    ::casac::variant const& scan,
		    ::casac::variant const& observation,
		    ::casac::variant const& polarization,
		    ::casac::variant const& beam,
		    ::casac::variant const& intent,
        ::casac::variant const& feed,
		    string const& taql,
			bool const reindex)
{
  bool rstat(false);
  *itsLog << _ORIGIN;
  try {
    assert_valid_ms();

    // make selection string to record.
    // Doing this here to make future extention easier.
    Record selection;
    String selection_string;
    // spw
    selection_string = toCasaString(spw);
    if (selection_string != "") {
      selection.define("spw", selection_string);
    }
    // field
    selection_string = toCasaString(field);
    if (selection_string != "")
      selection.define("field", selection_string);
    // antenna
    selection_string = toCasaString(antenna);
    if (selection_string != "")
      selection.define("antenna", selection_string);
    // time
    selection_string = toCasaString(timerange);
    if (selection_string != "")
      selection.define("timerange", selection_string);
    // scan
    selection_string = toCasaString(scan);
    if (selection_string != "")
      selection.define("scan", selection_string);
    // observation
    selection_string = toCasaString(observation);
    if (selection_string != "")
      selection.define("observation", selection_string);
    // polarization
    selection_string = toCasaString(polarization);
    if (selection_string != "")
      selection.define("correlation", selection_string);
    // beam
    selection_string = toCasaString(beam);
    if (selection_string != "")
      *itsLog << LogIO::WARN << "Beam selection is not yet supported. Ignoring beam selection" << LogIO::POST;
      //selection.define("beam", selection_string);
    // intent
    selection_string = toCasaString(intent);
    if (selection_string != "")
      selection.define("intent", selection_string);
    // feed
    selection_string = toCasaString(feed);
    if (selection_string != "")
      selection.define("feed", selection_string);
    // taql
    selection_string = toCasaString(taql);
    if (selection_string != "")
      selection.define("taql", selection_string);

    selection.define("reindex", Bool(reindex));

    itsSd->setSelection(selection);
    rstat = true;
  } catch  (AipsError x) {
    *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
	    << LogIO::POST;
    RETHROW(x);
  }
  return rstat;
}

bool
singledishms::smooth(string const &type, float const width,
        string const &datacolumn, string const &outfile)
{
    bool rstat(false);
    *itsLog << _ORIGIN;
    try {
      assert_valid_ms();
      itsSd->smooth(type, width, datacolumn, outfile);
      rstat = true;
    } catch  (AipsError x) {
      *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
          << LogIO::POST;
      RETHROW(x);
    }
    return rstat;
}

bool
singledishms::atmcor(const ::casac::record& config,
        string const &datacolumn, string const &outfile)
{
    bool rstat(false);
    *itsLog << _ORIGIN;
    try {
      assert_valid_ms();
      std::unique_ptr<Record> atmCorConfig(toRecord(config));
      itsSd->atmcor(*atmCorConfig, datacolumn, outfile);
      rstat = true;
    } catch  (AipsError x) {
      *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
          << LogIO::POST;
      RETHROW(x);
    }
    return rstat;
}

bool
singledishms::importasap(string const &infile, string const &outfile, bool const parallel)
{
    bool rstat(false);
    *itsLog << _ORIGIN;
    try {
      rstat = SingleDishMS::importAsap(infile, outfile, parallel);
    } catch  (AipsError x) {
      *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
          << LogIO::POST;
      RETHROW(x);
    }
    return rstat;
}

bool
singledishms::importnro(string const &infile, string const &outfile, bool const parallel)
{
    bool rstat(false);
    *itsLog << _ORIGIN;
    try {
      rstat = SingleDishMS::importNRO(infile, outfile, parallel);
    } catch  (AipsError x) {
      *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
          << LogIO::POST;
      RETHROW(x);
    }
    return rstat;
}

} // end of casac namespace
