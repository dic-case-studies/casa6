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

#include <ms/MeasurementSets/MeasurementSet.h>
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
		 ss << " Sel: ";
		 const vector<Field>& f = fields();
		 bool valueWritten = false;
		 for(unsigned int i = 0; i < f.size(); i++){
			 String fieldValue = getValue( f[i] );
			 if ( fieldValue.length() > 0 ){
				 if ( valueWritten ){
					 ss << ", ";
				 }
				 else {
					 valueWritten = true;
				 }
				 ss << f[i] << ": "<<fieldValue;
			 }
		 }
	 }
	 return ss.str();
}

void PlotMSSelection::apply(MeasurementSet& ms, MeasurementSet& selMS,
			    Vector<Vector<Slice> >& chansel,
			    Vector<Vector<Slice> >& corrsel){
    // Set the selected MeasurementSet to be the same initially as the input
    // MeasurementSet
    selMS = ms;
    MSSelection mss;
    String spwstr = spw(); // for errormsg
    try {
        mssSetData2(ms, selMS, chansel,corrsel, "", 
           timerange(), antenna(), field(), spwstr,
           uvrange(), msselect(), corr(), scan(), array(),
           intent(), observation(), feed(), 1, &mss );
    } catch(AipsError x) {
        String errormsg = x.getMesg();
        if (errormsg.startsWith("Spw Expression: No match found") && (spwstr[0] != '"') && (spwstr.find('-') != std::string::npos)) {
            errormsg += "\nTIP: For a name match (particularly names with a hyphen), add double quotes around the name in the spw string."; 
        }
        throw(AipsError(errormsg));
    }

    selAnts1.resize(0);
    selAnts2.resize(0);
    if ( antenna().length() > 0 ){
    	selAnts1 = mss.getAntenna1List();
    	selAnts2 = mss.getAntenna2List();
    }
}

Vector<int> PlotMSSelection::getSelectedAntennas1(){
	return selAnts1;
}

Vector<int> PlotMSSelection::getSelectedAntennas2(){

	return selAnts2;
}

void PlotMSSelection::apply(NewCalTable& ct, NewCalTable& selCT,
  		            Vector<Vector<Slice> >& /*chansel*/,
  		            Vector<Vector<Slice> >& /*corrsel*/) {
  // Trap unsupported selections

  if (uvrange().length()>0)
    throw(AipsError("Selection by uvrange not supported for NewCalTable"));
  if (array().length()>0)
    throw(AipsError("Selection by array not supported for NewCalTable"));
  if (feed().length()>0)
    throw(AipsError("Selection by feed not supported for NewCalTable"));

  // Set the selected NewCalTable to be the same initially as the input
  // NewCalTable
  selCT = ct;
  // set up CTSelection with expressions
  CTSelection cts;
  cts.setTimeExpr(timerange());
  cts.setAntennaExpr(antenna());
  cts.setFieldExpr(field());
  cts.setSpwExpr(spw());
  cts.setTaQLExpr(msselect());
  // poln selection handled in loadCalAxis (uses getParSlice)
  cts.setScanExpr(scan());
  cts.setStateExpr(intent());
  cts.setObservationExpr(observation());
  // do selection
  CTInterface cti(ct);
  TableExprNode ten = cts.toTableExprNode(&cti);
  try {
    getSelectedTable(selCT, ct, ten, "");
  } catch(AipsError x) {
      throw(AipsError("Error selecting on caltable:\n"+
        x.getMesg()));
  }
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
