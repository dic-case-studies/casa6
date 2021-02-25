//# PlotMSSelection.cc: MS Selection parameters.
//# Copyright (C) 2009
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
//# $Id: $
#include <plotms/PlotMS/PlotMSSelection.h>

#include <ms/MSSel/MSSelectionTools.h>
#include <synthesis/CalTables/NewCalTable.h>
#include <synthesis/CalTables/CTInterface.h>
#include <synthesis/CalTables/CTSelection.h>
#include <QDebug>

using namespace casacore;
namespace casa {

/////////////////////////////////
// PLOTMSSELECTION DEFINITIONS //
/////////////////////////////////

// Static //

String PlotMSSelection::defaultValue(Field /*f*/) { return ""; }


// Non-Static //

PlotMSSelection::PlotMSSelection() {
    initDefaults(); }

PlotMSSelection::PlotMSSelection(const PlotMSSelection& copy) {
    operator=(copy); }

PlotMSSelection::~PlotMSSelection() { }


void PlotMSSelection::fromRecord(const RecordInterface& record) {
    const vector<String>& fields = fieldStrings();
    for(unsigned int i = 0; i < fields.size(); i++)
        if(record.isDefined(fields[i])&&record.dataType(fields[i]) == TpString)
            setValue(field(fields[i]), record.asString(fields[i]));
    if (record.isDefined("forceNew")&&record.dataType("forceNew")==TpInt)
      setForceNew(record.asInt("forceNew"));
}



Record PlotMSSelection::toRecord() const {
    Record record(Record::Variable);
    
    const vector<Field>& f = fields();
    for(unsigned int i = 0; i < f.size(); i++)
        record.define(field(f[i]), getValue(f[i]));
    record.define("forceNew",forceNew());
    
    return record;
}

bool PlotMSSelection::isEmpty() const {
    bool emptySelection = true;
    const vector<Field>& f = fields();
    for(unsigned int i = 0; i < f.size(); i++){
        String fieldValue = getValue( f[i] );
        if ( fieldValue.length() > 0 ){
            emptySelection = false;
            break;
        }
    }
    return emptySelection;
}

String PlotMSSelection::toStringShort() const {
    stringstream ss;
    ss.precision(6);
    if ( !isEmpty() ){
        ss << "  Selection: ";
        const vector<Field>& f = fields();
        bool valueWritten = false;
        for(unsigned int i = 0; i < f.size(); i++){
            String fieldValue = getValue( f[i] );
            if ( fieldValue.length() > 0 ){
                if ( valueWritten ){
                    ss << ", ";
                } else {
                    valueWritten = true;
                }
                ss << PlotMSSelection::field(f[i]) << ": "<<fieldValue;
            }
        }
    }
    return ss.str();
}

void PlotMSSelection::apply(MeasurementSet& ms, MeasurementSet& selMS,
    Vector<Vector<Slice> >& chansel,
    Vector<Vector<Slice> >& corrsel){
    // Set the selected MeasurementSet to be the same initially
    // as the input MeasurementSet
    selMS = ms;
    if (!isEmpty()) {
        MSSelection mss;
        String spwstr = spw(); // for errormsg
        try {
            mssSetData2(ms, selMS, chansel,corrsel, "", 
               timerange(), antenna(), field(), spwstr,
               uvrange(), msselect(), corr(), scan(), array(),
               intent(), observation(), feed(), 1, &mss );
        } catch(AipsError& x) {
            String errormsg = x.getMesg();
            if (errormsg.startsWith("Spw Expression: No match found") && (spwstr[0] != '"') && (spwstr.find('-') != std::string::npos)) {
                errormsg += "\nTIP: For a name match (particularly names with a hyphen), add double quotes around the name in the spw string."; 
            }
            throw(AipsError(errormsg));
        }

        selectedAnt1.resize();
        selectedAnt2.resize();
        if ( antenna().length() > 0 ){
            selectedAnt1 = mss.getAntenna1List();
            selectedAnt2 = mss.getAntenna2List();
        }
    }
}

Vector<int> PlotMSSelection::getSelectedAntennas1() {
    return selectedAnt1;
}

Vector<int> PlotMSSelection::getSelectedAntennas2() {
    return selectedAnt2;
}

Matrix<int> PlotMSSelection::getSelectedChannels() {
    return selectedChan;
}

void PlotMSSelection::apply(NewCalTable& ct, NewCalTable& selCT) {
  // Trap unsupported selections
  if (uvrange().length()>0)
    throw(AipsError("Selection by uvrange not supported for NewCalTable"));
  if (array().length()>0)
    throw(AipsError("Selection by array not supported for NewCalTable"));
  if (intent().length()>0)
    throw(AipsError("Selection by intent not supported for NewCalTable"));
  if (feed().length()>0)
    throw(AipsError("Selection by feed not supported for NewCalTable"));

  if (isEmpty()) {
    // Set the selected NewCalTable to be the same as the input NewCalTable
    selCT = ct;
  } else {
    // set up CTSelection with expressions
    CTSelection cts;
    cts.setTimeExpr(timerange());
    cts.setAntennaExpr(antenna());
    cts.setFieldExpr(field());
    cts.setSpwExpr(spw());
    cts.setTaQLExpr(msselect());
    // poln selection handled in loadCalAxis (uses getParSlice)
    cts.setScanExpr(scan());
    cts.setObservationExpr(observation());
    // do selection
    CTInterface cti(ct);
    TableExprNode ten = cts.toTableExprNode(&cti);

    try {
      selCT = ct(ten);

      // Store antenna selection
      selectedAnt1.resize();
      selectedAnt2.resize();
      if (antenna().length() > 0) {
        selectedAnt1 = cts.getAntenna1List();
        selectedAnt2 = cts.getAntenna2List();
      }

      // Store channel selection
      selectedChan.resize();
      if (spw().length() > 0) {
          selectedChan = cts.getChanList();
      }
    } catch(AipsError& x) {
      throw(AipsError("Error selecting on caltable:\n" + x.getMesg()));
    }

    selectedAnt1.resize();
    selectedAnt2.resize();
    if (antenna().length() > 0) {
      selectedAnt1 = cts.getAntenna1List();
      selectedAnt2 = cts.getAntenna2List();
    }
  }
}

void PlotMSSelection::apply(CalTable& ct, CalTable& selectedCT,
    Vector<Vector<Slice> >& chansel,
    Vector<Vector<Slice> >& corrsel) {
  // Set the selected CalTable to be the same initially 
  // as the input CalTable
  selectedCT = ct;
  selectedAnt1.resize(); // set in getAntTaql if selected
  selectedAnt2.resize(); // not applicable
  selectedChan.resize();

  if (!isEmpty()) {
    // do taql selection and check result
    String calselect(getTaql(ct));
    selectedCT = ct.select(calselect);
    if (selectedCT.nRowMain()==0) 
      throw(AipsError("NullSelection: The selected table has zero rows."));

    // do channel selection to set chansel vector
    String spwExpr(spw());
    if (!spwExpr.empty() && spwExpr.contains(":")) {
      MSSelection mss;
      MeasurementSet ms(getMSName(ct)), selMS;
      mssSetData2(ms, selMS, chansel,corrsel, "", "", "", "", spwExpr,
        "", "", "", "", "", "", "", "", 1, &mss );
    }
  }
}

String PlotMSSelection::getTaql(CalTable& ct) {
  // convert selection expressions to taql expressions
  String taqlExpr("");
  // Use MSSelection with MS to get selected lists
  MSSelection mss;
  MeasurementSet ms(getMSName(ct));
  // add taql for each selected Field
  const vector<Field>& f = fields();
  for(unsigned int i = 0; i < f.size(); i++) {
    switch (f[i]) {
      case FIELD: {
        if (!field().empty()) {
          mss.setFieldExpr(field());
          mss.toTableExprNode(&ms);
          Vector<Int> fieldlist = mss.getFieldList();
          // add taql
          if (!taqlExpr.empty()) taqlExpr += " && ";
          String fieldlistStr = mss.indexExprStr(fieldlist);
          taqlExpr += "FIELD_ID IN [" + fieldlistStr + "]";
        }
        break;
      }
      case SPW: {
        if (!spw().empty()) {
          mss.setSpwExpr(spw()); // CAL_DESC table, not main
          mss.toTableExprNode(&ms);
          Vector<Int> spwlist = mss.getSpwList();
          Vector<Int> calDesclist;
          // Find cal desc id for each spw
          for (uInt ispw=0; ispw<spwlist.size(); ++ispw) {
            for (Int cdrow=0; cdrow<ct.nRowDesc(); ++cdrow) {
              Record rec = ct.getRowDesc(cdrow);
              Vector<Int> spwsPerDesc = rec.asArrayInt("SPECTRAL_WINDOW_ID");
              for (uInt ispwd=0; ispwd<spwsPerDesc.size(); ++ispwd) {
                // Is spw in this row?
                if (spwsPerDesc(ispwd)==spwlist(ispw)) {
                  uInt idsize = calDesclist.size();
                  calDesclist.resize(idsize+1, True);
                  calDesclist(idsize) = cdrow;
                  break;
                }
              }
            }
          }

          if (calDesclist.empty())
            throw(AipsError("Spw Expression: no CAL_DESC_ID match found."));
          // add taql
          if (!taqlExpr.empty()) taqlExpr += " && ";
          String calDesclistStr = mss.indexExprStr(calDesclist);
          taqlExpr += "CAL_DESC_ID IN [" + calDesclistStr + "]";
        }
        break;
      }
      case TIMERANGE: {
        if (!timerange().empty()) {
          mss.setTimeExpr(timerange());
          mss.toTableExprNode(&ms);
          // timelist is [start, end, dT (exposure)] in columns
          Matrix<Double> timelist = mss.getTimeList();
          for (uInt icol=0; icol<timelist.ncolumn(); ++icol) {
            Vector<Double> timesel = timelist.column(icol);
            // add taql
            if (!taqlExpr.empty()) taqlExpr += " && ";
            taqlExpr += "TIME BETWEEN " + String::toString(timesel(0)) + " AND " + 
              String::toString(timesel(1));
          }
        }
        break;
      }
      case UVRANGE: {
        if (!uvrange().empty())
          throw(AipsError("Selection by uvrange not supported for CalTable"));
        break;
      }
      case ANTENNA: {
        if (!antenna().empty()) {
          // add taql
          if (!taqlExpr.empty()) taqlExpr += " && ";
          taqlExpr += getAntTaql(mss, ms, antenna());
        }
        break;
      }
      case SCAN: {
        if (!scan().empty()) {
          mss.setScanExpr(scan());
          mss.toTableExprNode(&ms);
          Vector<Int> scanlist = mss.getScanList();
          // add taql
          if (!taqlExpr.empty()) taqlExpr += " && ";
          String scanListStr = mss.indexExprStr(scanlist);
          taqlExpr += "SCAN_NUMBER IN [" + scanListStr + "]";
        }
        break;
      }
      case CORR:  // done in getParSlice() to slice cubes
        break;
      case ARRAY: {
        if (!array().empty())
          throw(AipsError("Selection by array not supported for CalTable"));
        break;
      }
      case OBSERVATION: {
        if (!observation().empty()) {
          mss.setObservationExpr(observation());
          mss.toTableExprNode(&ms);
          Vector<Int> obslist = mss.getObservationList();
          // add taql
          if (!taqlExpr.empty()) taqlExpr += " && ";
          String obsListStr = mss.indexExprStr(obslist);
          taqlExpr += "OBSERVATION_ID IN [" + obsListStr + "]";
        }
        break;
      }
      case INTENT: {
        if (!intent().empty())
          throw(AipsError("Selection by intent not supported for CalTable"));
        break;
      }
      case FEED: {
        if (!feed().empty())
          throw(AipsError("Selection by feed not supported for CalTable"));
        break;
      }
      case MSSELECT: {
        if (!msselect().empty()) {
          // just add taql
          if (!taqlExpr.empty()) taqlExpr += " && ";
          taqlExpr += msselect();
        }
        break;
      }
    } // switch
  } // for
  return taqlExpr;
}

String PlotMSSelection::getMSName(CalTable& ct) {
  // Get ms filename from caltable path + msname
  // Throw exception if no MS (cannot use classes which calculate solutions)
  Path calpath(ct.tableName());
  String caldir(calpath.dirName());
  if (!caldir.empty()) caldir += "/";
  String msname(ct.getRowDesc(0).asString("MS_NAME"));
  String fullmsname(caldir+msname);
  if (msname.empty() || !Table::isReadable(fullmsname))
    throw(AipsError("Associated MS is not available, cannot select or plot solutions."));
  return fullmsname;
}

String PlotMSSelection::getAntTaql(MSSelection& mss, MeasurementSet& ms, String antExpr) {
  // convert negative values to NOT if necessary
  String taql;
  String beginAnt("ANTENNA1 IN ["), beginNotAnt("ANTENNA1 NOT IN [");
  String selAnt(""), selNotAnt("");
  Vector<Int> selAntIds, selNotAntIds;

  Vector<String> selections;
  selections = split(antExpr, ';', selections);

  for (uInt expr=0; expr<selections.size(); ++expr) {
    String subexpr = selections(expr);
    bool notselected = subexpr.contains("!"); // needed for zero
    mss.setAntennaExpr(subexpr);
    mss.toTableExprNode(&ms);
    Vector<Int> ant1list = mss.getAntenna1List();

    if (notselected) {
      ant1list *= -1;

      // Add to Vector
      selNotAntIds = casacore::set_union(selNotAntIds, ant1list);
      // Add to String
      if (!selNotAnt.empty()) selNotAnt += ",";
      selNotAnt += mss.indexExprStr(ant1list);
    } else  {
      // Add to Vector
      selAntIds = casacore::set_union(selAntIds, ant1list);
      // Add to String
      if (!selAnt.empty()) selAnt += ",";
      selAnt += mss.indexExprStr(ant1list);
    }
  }

  // Set selected ant1
  selectedAnt1.resize();
  selectedAnt2.resize();
  selectedAnt1 = selAntIds;
  // Append de-selected ant1
  size_t sel_size(selAntIds.size()), notsel_size(selNotAntIds.size());
  selectedAnt1.resize(sel_size + notsel_size, true);
  for (size_t i = 0; i < notsel_size; ++i) {
    selectedAnt1(sel_size + i) = selNotAntIds(i);
  } 

  // Set and return taqlExpr
  taql = (selAnt.empty() ? "" : beginAnt + selAnt + "]");
  if (!taql.empty() && !selNotAnt.empty()) taql += " && ";
  taql += (selNotAnt.empty() ? "" : beginNotAnt + selNotAnt + "]");
  return taql;
}

const String& PlotMSSelection::getValue(Field f) const {
    return const_cast<map<Field,String>&>(itsValues_)[f]; }
void PlotMSSelection::setValue(Field f, const String& value) {
    itsValues_[f] = value; }


bool PlotMSSelection::operator==(const PlotMSSelection& other) const {    

  // Check forceNew counter first
  //  not equal (first reset forceNew so that it isn't sticky)
  if (forceNew()!=other.forceNew())  {
    return false;
  }

  return fieldsEqual(other);

}

bool PlotMSSelection::fieldsEqual(const PlotMSSelection& other) const {    
  vector<Field> f = fields();
  for(unsigned int i = 0; i < f.size(); i++)
    if(getValue(f[i]) != other.getValue(f[i])) return false;
  
  return true;
}

PlotMSSelection& PlotMSSelection::operator=(const PlotMSSelection& copy) {
    itsValues_ = copy.itsValues_;    
    forceNew_ = copy.forceNew_;
    return *this;
}


void PlotMSSelection::initDefaults() {
    vector<Field> f = fields();
    for(unsigned int i = 0; i < f.size(); i++)
        itsValues_[f[i]] = defaultValue(f[i]);

    // forceNew_ is a counter, start it at zero
    forceNew_=0;
}

}
