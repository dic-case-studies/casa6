//# PlotMSLabelFormat.h: Class for generating labels based on a format.
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
#ifndef PLOTMSLABELFORMAT_H_
#define PLOTMSLABELFORMAT_H_

#include <plotms/PlotMS/PlotMSConstants.h>

namespace casa {

//# Forward declarations.
class PlotMSPlotParameters;


// Class for generating labels based upon axes, units, etc.  A format is a
// casacore::String that consists of tags and non-tags.  Tags are special substrings
// differentiated by being surrounded by a separator; tags are replaced with
// given values when the label is built.  For example: the %%axis%% tag is
// replaced with the name of the Axis that is passed in when generating the
// String.  Non-tags are copied exactly into the label.  An example of a format
// is: "%%xaxis%% vs. %%yaxis%%", when passed in the axes TIME and AMP, would
// return the label "time vs. amp".  Currently, labels can either be single
// axis (like for an axis label) or double axes (like for a plot title).
class PlotMSLabelFormat {
public:
    // Static //
    
    // Separator that goes before and after tags.
    static const casacore::String& TAGSEPARATOR();
    
    // Tag for axes, both the single case and the double case.  This tag is
    // replaced with the axis name (see PMS::axis()) during generation.
    // <group>
    static const casacore::String& TAG_AXIS();
    static const casacore::String& TAG_XAXIS();
    static const casacore::String& TAG_YAXIS();
    // </group>
    
    // Tag for axes units, both the single case and the double case.  This tag
    // is replaced with the axis unit name (see PMS::axisUnit()) during
    // generation.
    // <group>
    static const casacore::String& TAG_UNIT();
    static const casacore::String& TAG_XUNIT();
    static const casacore::String& TAG_YUNIT();
    // </group>
    
    // Tags for if/endif for axes units.  Parts of the format that are
    // surrounded with an IF at the beginning and ENDIF at the end are only
    // copied to the label if the given unit is NOT PMS::UNONE.  For example,
    // "TAG_IF_UNIT(sample text)TAG_ENDIF_UNIT" would copy over "(sample text)"
    // only if the given unit was NOT PMS::UNONE.
    // <group>
    static const casacore::String& TAG_IF_UNIT();
    static const casacore::String& TAG_IF_XUNIT();
    static const casacore::String& TAG_IF_YUNIT();
    static const casacore::String& TAG_ENDIF_UNIT();
    static const casacore::String& TAG_ENDIF_XUNIT();
    static const casacore::String& TAG_ENDIF_YUNIT();
    // </group>
    
    // Tags for axis reference values, both the single case and the double
    // case.  This tag is replaced with the reference value for the given axis,
    // if there is one, during generation.  Note: if the associated axis is a
    // date, the value will be in date format.
    // <group>
    static const casacore::String& TAG_REFVALUE();
    static const casacore::String& TAG_XREFVALUE();
    static const casacore::String& TAG_YREFVALUE();
    // </group>
    
    // Tags for if/endif for axis reference values.  Parts of the format that
    // are surrounded with an IF at the beginning and ENDIF at the end are only
    // copied to the label if the given axis has a reference value set.  For
    // example, "TAG_IF_REFVALUE(sample text)TAG_ENDIF_REFVALUE" would copy
    // over "(sample text)" only if the given axis had a set reference value.
    // <group>
    static const casacore::String& TAG_IF_REFVALUE();
    static const casacore::String& TAG_IF_XREFVALUE();
    static const casacore::String& TAG_IF_YREFVALUE();
    static const casacore::String& TAG_ENDIF_REFVALUE();
    static const casacore::String& TAG_ENDIF_XREFVALUE();
    static const casacore::String& TAG_ENDIF_YREFVALUE();
    // </group>
    
    // Convenience method to surround the given tag with the separator.
    static casacore::String TAG(const casacore::String& tag);
    
    
    // Non-Static //
    
    // Constructor which takes an optional starting format.
    PlotMSLabelFormat(const casacore::String& format = "");
    
    // Copy constructor.  See operator=().
    PlotMSLabelFormat(const PlotMSLabelFormat& copy);
    
    // Destructor.
    ~PlotMSLabelFormat();
    
    
    // Format.
    casacore::String format;
    
    
    // Generates a label, using the given single axis and reference value.  If
    // any double axes tags are in the format, the given axis will be used for
    // them.
    casacore::String getLabel(PMS::Axis axis, bool refValueSet = false,
            double refValue = 0, PMS::DataColumn data = PMS::DATA,
            bool polnRatio = false) const;
    
    casacore::String getLabel(std::vector<PMS::Axis> axes, std::vector<bool> refValueSets,
    		std::vector<double> refValues, std::vector<PMS::DataColumn> datas,
            bool polnRatio = false ) const;

    // Generates a label, using the given double axes and reference values.  If
    // any single axes tags are in the format, the x axis will be used for it.
    casacore::String getLabel(PMS::Axis xAxis, std::vector<PMS::Axis> yAxes,
            bool xRefValueSet = false, double xRefValue = 0,
            std::vector<bool> yRefValueSets=std::vector<bool>(1,false), 
            std::vector<double> yRefValues =std::vector<double>(1,0.0),
            PMS::DataColumn xData=PMS::DATA,
            std::vector<PMS::DataColumn> yData = std::vector<PMS::DataColumn>(1,PMS::DATA),
            bool polnRatio = false ) const;
    
    // Equality operators.
    // <group>
    bool operator==(const PlotMSLabelFormat& other) const;
    bool operator!=(const PlotMSLabelFormat& other) const {
        return !(operator==(other)); }
    // </group>
    
    // Copy operator.
    PlotMSLabelFormat& operator=(const PlotMSLabelFormat& copy);
    
private:
    // Generates a label using the given format, single axis, and double axes.
    static casacore::String getLabel(const casacore::String& format, PMS::Axis axis,
            PMS::Axis xAxis, std::vector<PMS::Axis> yAxes, bool refValueSet,
            double refValue, bool xRefValueSet, double xRefValue,
            std::vector<bool> yRefValueSets, std::vector<double> yRefValues,
            PMS::DataColumn xData, const std::vector<PMS::DataColumn>& ydatas,
            bool polnRatio);
    
    //Adds the data label such as DATA, MODEL, CORRECTED, etc to the tag
    static void addDataToTag( casacore::String& tag, PMS::Axis axis, PMS::DataColumn column );
    // Adds "POLN Ratio" to tag for cal table ratio plots
    static void addPolnRatioToTag( casacore::String& tag, PMS::Axis axis );

    // Helper method for getLabel() which gets the next token in the format.
    // Returns true if a token was returned; false if the end of the format was
    // reached.  The given format will be automatically shortened as tokens are
    // taken out.  The next token is put in the token parameter, and the
    // tokenWasTag parameter is set to true if the token was a tag, false
    // otherwise.
    static bool nextToken(casacore::String& format, casacore::String& token, bool& tokenWasTag);
    
    // casacore::Format for when reference values are dates.
    static const casacore::String REFERENCE_DATE_FORMAT;
};

}

#endif /* PLOTMSLABELFORMAT_H_ */
