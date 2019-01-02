//# PlotMSFlagging.h: Flagging parameters.
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
#ifndef PLOTMSFLAGGING_H_
#define PLOTMSFLAGGING_H_

#include <plotms/PlotMS/PlotMSConstants.h>
#include <plotms/PlotMS/PlotMSSelection.h>

namespace casa {

// Specifies flagging parameters (including possibly flag extension) for an MS.
class PlotMSFlagging {    
public:
    // Static //
    
    // Enum and methods to define the different fields for an casacore::MS flagging.  All
    // fields have a bool flag for on/off, and some of them also have a double
    // value (see fieldHasValue()) or a PlotMSSelection value (see
    // fieldHasSelectionValue()).  Most fields are off by default (see
    // fieldDefault()), with a default double value of 0 or empty/default
    // PlotMSSelection value if applicable.  Some fields are in mutually
    // exclusive groups (see fieldMutuallyExclusiveGroup()).
    // **If these are changed, also update: convenience methods below,
    // xmlcasa/implement/plotms/plotms*, xmlcasa/tasks/plotms.xml,
    // xmlcasa/scripts/task_plotms.py.**
    // <group>
    PMS_ENUM1(Field, fields, fieldStrings, field,
              EXTEND, CORR, CORR_ALL,
              CORR_POLN_DEP, CHANNEL, SPW, ANTENNA,
              ANTENNA_ANTENNA, ANTENNA_BASELINES, TIME,
              SCANS, FIELD, SEL_SELECTED, SEL_ALTERNATE)
    PMS_ENUM2(Field, fields, fieldStrings, field,
              "extend", "correlation", "correlation:all",
              "correlation:poln-dep", "channel", "spw", "antenna",
              "antenna:antenna-based", "antenna:all-baselines", "time",
              "scans", "field", "selection:selected", "selection:alternate")
    // </group>
              
    // Returns whether the given field has a double value associated with it or
    // not.
    static bool fieldHasValue(Field f);
    
    // Returns whether the given field has a PlotMSSelection value associated
    // with it or not.
    static bool fieldHasSelectionValue(Field f);

    // Returns whether the given field is on or off by default.  This is
    // important for mutually exclusive groups, which need to know which field
    // in that group is on to begin with.
    static bool fieldDefault(Field f);
    
    // Returns the list of fields, NOT including the given, with which the
    // given field is mutually exclusive.  In a mutually exclusive group, only
    // and exactly one field can be turned on at a given time.  NOTE: this is
    // different from mutually exclusive groups in PlotMSAveraging, because
    // here it is NOT allowed to have all fields in a group turned off.
    static const std::vector<Field>& fieldMutuallyExclusiveGroup(Field f);

    // Returns true if the given field is in a mutually exclusive group, false
    // otherwise.  See fieldMutuallyExclusiveGroup().
    static bool fieldIsInMutuallyExclusiveGroup(Field f) {
        return fieldMutuallyExclusiveGroup(f).size() > 0; }
    
    
    // Non-Static //
    
    // Default constructor.
    PlotMSFlagging();
    
    // Destructor.
    ~PlotMSFlagging();
        
    // Converts this object to/from a record.  Each field will have a key that
    // is its enum name, with a bool value for its flag value.  Fields that
    // also have double values will have an additional key that is its enum
    // name + "Value" for double values, or its enum name + "SelectionValue"
    // for PlotMSSelection values.  PlotMSSelection values are stored in Record
    // form (see PlotMSSelection::fromRecord()).  The casacore::MS objects (the casacore::MS, the
    // selected casacore::MS, and the vis set) are NOT included in the record.
    // <group>
    void fromRecord(const casacore::RecordInterface& record);
    casacore::Record toRecord(bool useStrings = false) const;
    // </group>
    
    // Gets/Sets the on/off flag for the given field.
    // <group>
    bool getFlag(Field f) const;
    void getFlag(Field f, bool& flag) const { flag = getFlag(f); }
    void setFlag(Field f, bool on);
    // </group>
    
    // Gets/Sets the double value for the given field, if applicable.
    // <group>
    double getValue(Field f) const;
    void getValue(Field f, double& value) const { value = getValue(f); }
    void setValue(Field f, double value);
    // </group>
    
    // Gets/Sets the value for the given field as a String.  Only applicable
    // for special cases (correlation, antenna).
    // Correlation: "", "all", or "poln-dep".
    // Antenna: "", "all", or antenna-based value.
    // <group>
    casacore::String getValueStr(Field f) const;
    void getValue(Field f, casacore::String& value) const { value = getValueStr(f); }
    void setValue(Field f, const casacore::String& value);
    // </group>
    
    // Gets/Sets the selection value for the given field, if applicable.
    // <group>
    PlotMSSelection getSelectionValue(Field f) const;
    void getSelectionValue(Field f, PlotMSSelection& value) const {
        value = getSelectionValue(f); }
    void setSelectionValue(Field f, const PlotMSSelection& value);
    // </group>
    
    
    // Convenience methods for returning the standard field values.
    // <group>
    bool extend() const { return getFlag(EXTEND); }
    bool corr() const { return (extend() && getFlag(CORR)); }
    bool corrAll() const { return (corr() && getFlag(CORR_ALL)); }
    bool corrPolnDep() const { return (corr() && getFlag(CORR_POLN_DEP)); }
    casacore::String corrStr() const { return getValueStr(CORR); }
    bool channel() const { return (extend() && getFlag(CHANNEL)); }
    bool spw() const { return (extend() && getFlag(SPW)); }
    bool antenna() const { return (extend() && getFlag(ANTENNA)); }
    bool antennaAntennaBased() const { return (antenna() && getFlag(ANTENNA_ANTENNA)); }
    double antennaAntennaBasedValue() const{ return getValue(ANTENNA_ANTENNA);}
    bool antennaBaselinesBased() const { return (antenna() && getFlag(ANTENNA_BASELINES)); }
    casacore::String antennaStr() const { return getValueStr(ANTENNA); }
    bool time() const { return (extend() && getFlag(TIME)); }
    bool scans() const { return (time() && getFlag(SCANS)); }
    bool field() const { return (time() && getFlag(FIELD)); }
    bool selectionSelected() const { return (extend() && getFlag(SEL_SELECTED)); }
    bool selectionAlternate() const { return (extend() && getFlag(SEL_ALTERNATE)); }
    PlotMSSelection selectionAlternateSelection() const {
        return getSelectionValue(SEL_ALTERNATE); }
    // </group>
    
    // Convenience methods for setting the standard field values.
    // <group>
    void setExtend(bool flag) { setFlag(EXTEND, flag); }
    void setCorr(bool flag) { setFlag(CORR, flag); }
    void setCorrAll(bool flag) { setFlag(CORR_ALL, flag); }
    void setCorrPolnDep(bool flag) { setFlag(CORR_POLN_DEP, flag); }
    void setCorr(const casacore::String& value) { setValue(CORR, value); }
    void setChannel(bool flag) { setFlag(CHANNEL, flag); }
    void setSpw(bool flag) { setFlag(SPW, flag); }
    void setAntenna(bool flag) { setFlag(ANTENNA, flag); }
    void setAntennaAntennaBased(bool flag) { setFlag(ANTENNA_ANTENNA, flag); }
    void setAntennaAntennaBasedValue(double value) {
        setValue(ANTENNA_ANTENNA, value); }
    void setAntennaBaselinesBased(bool flag){setFlag(ANTENNA_BASELINES, flag);}
    void setAntenna(const casacore::String& value) { setValue(ANTENNA, value); }
    void setTime(bool flag) { setFlag(TIME, flag); }
    void setScans(bool flag) { setFlag(SCANS, flag); }
    void setField(bool flag) { setFlag(FIELD, flag); }
    void setSelectionSelected(bool flag) { setFlag(SEL_SELECTED, flag); }
    void setSelectionAlternate(bool flag) { setFlag(SEL_ALTERNATE, flag); }
    void setSelectionAlternateSelection(const PlotMSSelection& value) {
        setSelectionValue(SEL_ALTERNATE, value); }
    // </group>
    
    
    // Equality operators.
    // <group>
    bool operator==(const PlotMSFlagging& other) const;
    bool operator!=(const PlotMSFlagging& other) const {
        return !(operator==(other)); }
    // </group>
    
private:
    
    // Flagging field flags.
    std::map<Field, bool> itsFlags_;
    
    // Flagging field double values.
    std::map<Field, double> itsValues_;
    
    // Flagging field selection values.
    std::map<Field, PlotMSSelection> itsSelectionValues_;

    
    // Sets the default values.
    void setDefaults();
    
    
    // casacore::String constant for what to append to the enum name in the record to
    // get the key for the value.
    // <group>
    static const casacore::String RKEY_VALUE;
    static const casacore::String RKEY_SELVALUE;
    // </group>
};

}

#endif /* PLOTMSFLAGGING_H_ */
