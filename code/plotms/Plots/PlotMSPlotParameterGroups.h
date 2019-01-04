//# PlotMSPlotParameterGroups.h: Implementations of plot subparameter groups.
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
#ifndef PLOTMSPLOTPARAMETERGROUPS_H_
#define PLOTMSPLOTPARAMETERGROUPS_H_

#include <plotms/Plots/PlotMSPlotParameters.h>

#include <plotms/PlotMS/PlotMSAveraging.h>
#include <plotms/PlotMS/PlotMSExportParam.h>
#include <plotms/PlotMS/PlotMSIterParam.h>
#include <plotms/PlotMS/PlotMSSelection.h>
#include <plotms/PlotMS/PlotMSTransformations.h>
#include <plotms/PlotMS/PlotMSCalibration.h>
#include <plotms/PlotMS/PlotMSLabelFormat.h>
#include <plotms/PlotMS/PlotMSPageHeaderParam.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/TableRecord.h>

namespace casa {

// Container class to hold constants for groups.
class PMS_PP {

public:
	// Update flag for redrawing.
	// <group>
	static const casacore::String UPDATE_REDRAW_NAME;
	static const int UPDATE_REDRAW;
	// </group>

	// Update flag for casacore::MS data group.
	// <group>
	static const casacore::String UPDATE_MSDATA_NAME;
	static const int UPDATE_MSDATA;
	// </group>

	// Update flag for cache group.
	// <group>
	static const casacore::String UPDATE_CACHE_NAME;
	static const int UPDATE_CACHE;
	// </group>

	// Update flag for axes group.
	// <group>
	static const casacore::String UPDATE_AXES_NAME;
	static const int UPDATE_AXES;
	// </group>

	// Update flag for canvas group.
	// <group>
	static const casacore::String UPDATE_CANVAS_NAME;
	static const int UPDATE_CANVAS;
	// </group>

	// Update flag for display group.
	// <group>
	static const casacore::String UPDATE_DISPLAY_NAME;
	static const int UPDATE_DISPLAY;
	// </group>

	// Update flag for page header group.
	// <group>
	static const casacore::String UPDATE_PAGEHEADER_NAME;
	static const int UPDATE_PAGEHEADER;
	// </group>

	// Update flag for iteration group.
	// <group>
	static const casacore::String UPDATE_ITERATION_NAME;
	static const int UPDATE_ITERATION;
	// </group>

	// Update flag for log group.
	// <group>
	static const casacore::String UPDATE_LOG_NAME;
	static const int UPDATE_LOG;
	// </group>
	//
	// Update flag for plotms_options group.
	// <group>
	static const casacore::String UPDATE_PLOTMS_OPTIONS_NAME;
	static const int UPDATE_PLOTMS_OPTIONS;
	// </group>

private:
	// Disable constructor.
	PMS_PP() {
	}
};




// Subclass of PlotMSPlotParameters::Group to handle subparameters for casacore::MS data.
// Currently includes:
// * filename
// * selection
// * averaging
//
class PMS_PP_MSData : public PlotMSPlotParameters::Group {

public:
	/* Constructor which takes a factory */
	PMS_PP_MSData (PlotFactoryPtr factory);

	/* Copy constructor.  See operator=(). */
	PMS_PP_MSData (const PMS_PP_MSData & copy);
	~PMS_PP_MSData();


	/* Overrides PlotMSPlotParameters::Group::operator=(). */
	Group *clone() const {
		return new PMS_PP_MSData (*this);
	}

	/* Implements PlotMSPlotParameters::Group::name(). */
	const casacore::String & name() const {
		static casacore::String groupName = PMS_PP::UPDATE_MSDATA_NAME;
		return groupName;
	}

	/* Implements PlotMSPlotParameters::Group::toRecord(). */
	casacore::Record toRecord() const;

	/* Implements PlotMSPlotParameters::Group::fromRecord(). */
	void fromRecord (const casacore::Record & record);

	/* Implements PlotMSPlotParameters::Group::requiresRedrawOnChanged(). */
	bool requiresRedrawOnChange() const {
		return true;
	}

	/* Overrides the real assignment operator, operator=(). */
	PMS_PP_MSData& operator=(const PMS_PP_MSData& other);

	/* Overrides PlotMSPlotParameters::Group::operator=(). */
	PMS_PP_MSData & operator= (const Group & other);

	/* Overrides PlotMSPlotParameters::Group::operator==(). */
	bool operator== (const Group & other) const;


	bool isSet() const {
		return !itsFilename_.empty();
	}

	const casacore::String & filename() const {
		return itsFilename_;
	}
	void setFilename (const casacore::String & value) {
		if (itsFilename_ != value) {
			itsFilename_ = value;
			updated();
		}
	}

    // based on PlotMSCacheBase::Type enum {MS, CAL}
    const casacore::Int & type() const {
        return itsType_;
    }
	void setType (const casacore::Int & value) {
		if (itsType_ != value) {
			itsType_ = value;
			updated();
		}
	}

	const PlotMSSelection & selection() const {
		return itsSelection_;
	}
	void setSelection (const PlotMSSelection & value) {
		if (itsSelection_ != value) {
			itsSelection_ = value;
			updated();
		}
	}

	const PlotMSAveraging & averaging() const {
		return itsAveraging_;
	}
	void setAveraging (const PlotMSAveraging & value) {
		if (itsAveraging_ != value) {
			itsAveraging_ = value;
			updated();
		}
	}

	const PlotMSTransformations & transformations() const {
		return itsTransformations_;
	}
	void setTransformations (const PlotMSTransformations & value) {
		if (itsTransformations_ != value) {
			itsTransformations_ = value;
			updated();
		}
	}

	const PlotMSCalibration & calibration() const {
		return itsCalibration_;
	}
	void setCalibration (const PlotMSCalibration & value) {
		if (itsCalibration_ != value) {
			itsCalibration_ = value;
			updated();
		}
	}


private:
	//Does the work of the operator=()s.
	PMS_PP_MSData& assign(const PMS_PP_MSData* other);

	/* Parameters' values */
	casacore::String itsFilename_;
	casacore::Int itsType_;
	PlotMSSelection itsSelection_;
	PlotMSAveraging itsAveraging_;
	PlotMSTransformations itsTransformations_;
	PlotMSCalibration itsCalibration_;

	/* Key strings for casacore::Record */
	static const casacore::String REC_FILENAME;
	static const casacore::String REC_TYPE;
	static const casacore::String REC_SELECTION;
	static const casacore::String REC_AVERAGING;
	static const casacore::String REC_TRANSFORMATIONS;
	static const casacore::String REC_CALIBRATION;

	void setDefaults();
};






// Subclass of PlotMSPlotParameters::Group to handle cache parameters.
// Currently includes:
// * x and y axes
// * x and y data columns
// Parameters are vector-based, on a per-plot basis.
//
class PMS_PP_Cache : public PlotMSPlotParameters::Group {

public:
	/* Constructor which takes a factory */
	PMS_PP_Cache (PlotFactoryPtr factory);

	/* Copy constructor.  See operator=(). */
	PMS_PP_Cache (const PMS_PP_Cache & copy);

	~PMS_PP_Cache();

	/* Implements PlotMSPlotParameters::Group::clone(). */
	Group *clone() const {
		return new PMS_PP_Cache (*this);
	}

	/* Implements PlotMSPlotParameters::Group::name(). */
	const casacore::String & name() const {
		static casacore::String groupName = PMS_PP::UPDATE_CACHE_NAME;
		return groupName;
	}

	/* Implements PlotMSPlotParameters::Group::toRecord(). */
	casacore::Record toRecord() const;

	/* Implements PlotMSPlotParameters::Group::fromRecord(). */
	void fromRecord (const casacore::Record & record);

	/* Implements PlotMSPlotParameters::Group::requiresRedrawOnChanged(). */
	bool requiresRedrawOnChange() const {
		return true;
	}

	/* Overrides the real assignment operator=().*/
	PMS_PP_Cache& operator=(const PMS_PP_Cache& other);

	/* Overrides PlotMSPlotParameters::Group::operator=(). */
	PMS_PP_Cache & operator= (const Group & other);

	/* Overrides PlotMSPlotParameters::Group::operator==(). */
	bool operator== (const Group & other) const;

	// Gets how many axes and data columns there are.
	// <group>
	unsigned int numXAxes() const;
	unsigned int numYAxes() const;
	// </group>


	// Setting the data and data columns for the x- and y-axes
	// <group>
	void setXAxis (const PMS::Axis & axis, const PMS::DataColumn & data,
			unsigned int index = 0);
	void setYAxis (const PMS::Axis & axis, const PMS::DataColumn & data,
			unsigned int index = 0);
	void setAxes (const PMS::Axis & xAxis, const PMS::Axis & yAxis,
			const PMS::DataColumn & xData,
			const PMS::DataColumn & yData, unsigned int index = 0);
	//<group>

	const std::vector<PMS::Axis> &xAxes() const {
		return itsXAxes_;
	}
	void setXAxes (const std::vector<PMS::Axis> &value) {
		if (itsXAxes_ != value) {
			itsXAxes_ = value;
			updated();
		}
	}
	PMS::Axis xAxis (unsigned int index = 0) const {
		if (index >= itsXAxes_.size())
			const_cast< std::vector<PMS::Axis>& >(itsXAxes_).resize (index + 1);
		return itsXAxes_[index];
	}
	void setXAxis (const PMS::Axis & value, unsigned int index = 0) {
		if (index >= itsXAxes_.size())
			itsXAxes_.resize (index + 1);
		if (itsXAxes_[index] != value) {
			itsXAxes_[index] = value;
			updated();
		}
	}


	const std::vector<PMS::Axis>&yAxes() const {
		return itsYAxes_;
	}
	void setYAxes (const std::vector<PMS::Axis> &value) {
		if (itsYAxes_ != value) {
			itsYAxes_ = value;
			updated();
		}
	}
	PMS::Axis yAxis (unsigned int index = 0) const {
		if (index >= itsYAxes_.size())
			const_cast< std::vector<PMS::Axis> &>(itsYAxes_).resize (index + 1);
		return itsYAxes_[index];
	}
	void setYAxis (const PMS::Axis & value, unsigned int index = 0) {
		if (index >= itsYAxes_.size())
			itsYAxes_.resize (index + 1);
		if (itsYAxes_[index] != value)   {
			itsYAxes_[index] = value;
			updated();
		}
	}



	const std::vector<PMS::DataColumn> &xDataColumns() const {
		return itsXData_;
	}
	void setXDataColumns (const vector < PMS::DataColumn > &value) {
		if (itsXData_ != value)   {
			itsXData_ = value;
			updated();
		}
	}
	PMS::DataColumn xDataColumn (unsigned int index = 0) const {
		if (index >= itsXData_.size())
			const_cast < vector < PMS::DataColumn >
		&>(itsXData_).resize (index + 1);
		return itsXData_[index];
	}
	void setXDataColumn (const PMS::DataColumn & value, unsigned int index =
			0) {
		if (index >= itsXData_.size())
			itsXData_.resize (index + 1);
		if (itsXData_[index] != value)   {
			itsXData_[index] = value;
			updated();
		}
	}


	const vector < PMS::DataColumn > &yDataColumns() const {
		return itsYData_;
	}
	void setYDataColumns (const vector < PMS::DataColumn > &value) {
		if (itsYData_ != value) {
			itsYData_ = value;
			updated();
		}
	}
	PMS::DataColumn yDataColumn (unsigned int index = 0) const {
		if (index >= itsYData_.size())
			const_cast < vector < PMS::DataColumn >
		&>(itsYData_).resize (index + 1);
		return itsYData_[index];
	}
	void setYDataColumn (const PMS::DataColumn & value, unsigned int index =
			0) {
		if (index >= itsYData_.size())
			itsYData_.resize (index + 1);
		if (itsYData_[index] != value)   {
			itsYData_[index] = value;
			updated();
		}
	}

	bool showAtm() const {
		return itsShowAtm_;
	}
	void setShowAtm (const bool & value) {
	    if (itsShowAtm_!= value) {
		    itsShowAtm_ = value;
		    updated();
	    }
    }

    bool showTsky() const {
		return itsShowTsky_;
	}
	void setShowTsky (const bool & value) {
	    if (itsShowTsky_!= value) {
		    itsShowTsky_ = value;
		    updated();
	    }
    }

	void resize( int count );

private:
	//Does the work for the operator=()s.
	PMS_PP_Cache& assign(const PMS_PP_Cache* o);

	/* Parameters' values */
	std::vector<PMS::Axis> itsXAxes_;
	std::vector<PMS::Axis> itsYAxes_;
	std::vector<PMS::DataColumn> itsXData_;
	std::vector<PMS::DataColumn> itsYData_;
	bool itsShowAtm_;
	bool itsShowTsky_;

	/* Key strings for casacore::Record */
	static const casacore::String REC_XAXES;
	static const casacore::String REC_YAXES;
	static const casacore::String REC_XDATACOLS;
	static const casacore::String REC_YDATACOLS;
	static const casacore::String REC_SHOWATM;
	static const casacore::String REC_SHOWTSKY;

	void setDefaults();
};


// Subclass of PlotMSPlotParameters::Group to handle axes parameters.
// Currently includes:
// * canvas attach axes
// * axes ranges, if any
// Parameters are vector-based, on a per-plot basis.
//
class PMS_PP_Axes : public PlotMSPlotParameters::Group {

public:
	/* Constructor which takes a factory */
	PMS_PP_Axes (PlotFactoryPtr factory);

	/* Copy constructor.  See operator=(). */
	PMS_PP_Axes (const PMS_PP_Axes & copy);

	~PMS_PP_Axes();


	/* Implements PlotMSPlotParameters::Group::clone(). */
	Group *clone() const {
		return new PMS_PP_Axes (*this);
	}

	/* Implements PlotMSPlotParameters::Group::name(). */
	const casacore::String & name() const {
		static casacore::String groupName = PMS_PP::UPDATE_AXES_NAME;
		return groupName;
	}

	/* Implements PlotMSPlotParameters::Group::toRecord(). */
	casacore::Record toRecord() const;

	/* Implements PlotMSPlotParameters::Group::fromRecord(). */
	void fromRecord (const casacore::Record & record);

	/* Implements PlotMSPlotParameters::Group::requiresRedrawOnChanged(). */
	bool requiresRedrawOnChange() const {
		return true;
	}

	/* Overrides the real assignment operator, operator= */
	PMS_PP_Axes& operator=(const PMS_PP_Axes& other);

	/* Overrides PlotMSPlotParameters::Group::operator=(). */
	PMS_PP_Axes & operator= (const Group & other);

	/* Overrides PlotMSPlotParameters::Group::operator==(). */
	bool operator== (const Group & other) const;


	// Gets how many axes there are.
	// <group>
	unsigned int numXAxes() const;
	unsigned int numYAxes() const;
	// </group>


	// Sets single versions of the parameters for the given index.
	// <group>
	void setAxes (const PlotAxis & xAxis, const PlotAxis & yAxis,
			unsigned int index = 0);
	void setXRange (const bool & set, const prange_t & range,
			unsigned int index = 0) {
		setRanges (set, yRangeSet (index), range, yRange (index), index);
	}
	void setYRange (const bool & set, const prange_t & range,
			unsigned int index = 0) {
		setRanges (xRangeSet (index), set, xRange (index), range, index);
	}
	void setRanges (const bool & xSet, const bool & ySet,
			const prange_t & xRange, const prange_t & yRange,
			unsigned int index = 0);
	// </group>


	const vector < PlotAxis > &xAxes() const {
		return itsXAxes_;
	}
	void setXAxes (const vector < PlotAxis > &value) {
		if (itsXAxes_ != value) {
			itsXAxes_ = value;
			updated();
		}
	}


	PlotAxis xAxis (unsigned int index = 0) const {
		if (index >= itsXAxes_.size())
			const_cast < vector < PlotAxis > &>(itsXAxes_).resize (index + 1);
		return itsXAxes_[index];
	}
	void setXAxis(const PlotAxis & value, unsigned int index = 0) {
		if (index >= itsXAxes_.size()){
			itsXAxes_.resize (index + 1);

		}
		if (itsXAxes_[index] != value) {

			itsXAxes_[index] = value;
			updated();
		}
	}


	const std::vector<PlotAxis> &yAxes() const {
		return itsYAxes_;
	}
	void setYAxes (const vector < PlotAxis > &value) {
		if (itsYAxes_ != value) {
			itsYAxes_ = value;
			updated();
		}
	}


	PlotAxis yAxis (unsigned int index = 0) const {
		if (index >= itsYAxes_.size())
			const_cast < vector < PlotAxis > &>(itsYAxes_).resize (index + 1);
		return itsYAxes_[index];
	}
	void setYAxis (const PlotAxis & value, unsigned int index = 0) {
		if (index >= itsYAxes_.size())
			itsYAxes_.resize (index + 1);
		if (itsYAxes_[index] != value) {
			itsYAxes_[index] = value;
			updated();
		}
	}

	void setYAxis( casacore::String& value, unsigned int index = 0 ){
		PlotAxis axisLocation = Y_LEFT;
		if ( value == "right"){
			axisLocation = Y_RIGHT;
		}
		setYAxis( axisLocation, index );
	}



	const std::vector<bool> &xRangesSet() const {
		return itsXRangesSet_;
	}
	void setXRanges (const vector < bool > &value) {
		if (itsXRangesSet_ != value) {
			itsXRangesSet_ = value;
			updated();
		}
	}


	bool xRangeSet (unsigned int index = 0) const {
		if (index >= itsXRangesSet_.size())
			const_cast < vector < bool > &>(itsXRangesSet_).resize (index + 1);
		return itsXRangesSet_[index];
	}
	void setXRange (const bool & value, unsigned int index = 0) {
		if (index >= itsXRangesSet_.size())
			itsXRangesSet_.resize (index + 1);
		if (itsXRangesSet_[index] != value) {
			itsXRangesSet_[index] = value;
			updated();
		}
	}


	const vector < bool > &yRangesSet() const {
		return itsYRangesSet_;
	}
	void setYRanges (const vector < bool > &value) {
		if (itsYRangesSet_ != value) {
			itsYRangesSet_ = value;
			updated();
		}
	}


	bool yRangeSet (unsigned int index = 0) const {
		if (index >= itsYRangesSet_.size())
			const_cast < vector < bool > &>(itsYRangesSet_).resize (index + 1);
		return itsYRangesSet_[index];
	}
	void setYRange (const bool & value, unsigned int index = 0) {
		if (index >= itsYRangesSet_.size())
			itsYRangesSet_.resize (index + 1);
		if (itsYRangesSet_[index] != value) {
			itsYRangesSet_[index] = value;
			updated();
		}
	}



	const vector < prange_t > &xRanges() const {
		return itsXRanges_;
	}
	void setXRanges (const vector < prange_t > &value) {
		if (itsXRanges_ != value) {
			itsXRanges_ = value;
			updated();
		}
	}


	const prange_t & xRange (unsigned int index = 0) const {
		return itsXRanges_[index];
	}
	void setXRange (const prange_t & value, unsigned int index = 0) {
		if (itsXRanges_[index] != value) {
			itsXRanges_[index] = value;
			updated();
		}
	}


	const vector < prange_t > &yRanges() const {
		return itsYRanges_;
	}
	void setYRanges (const vector < prange_t > &value) {
		if (itsYRanges_ != value) {
			itsYRanges_ = value;
			updated();
		}
	}


	const prange_t & yRange (unsigned int index = 0) const {
		return itsYRanges_[index];
	}
	void setYRange (const prange_t & value, unsigned int index = 0) {
		if (itsYRanges_[index] != value) {
			itsYRanges_[index] = value;
			updated();
		}
	}

	//Change the size of the vectors.
	void resize( int count, bool copyValues=False );

private:

	//Does the work for operator=()s.
	PMS_PP_Axes& assign(const PMS_PP_Axes* o);

	/* Parameters' values */
	std::vector<PlotAxis> itsXAxes_;
	std::vector<PlotAxis> itsYAxes_;
	std::vector<bool> itsXRangesSet_;
	std::vector<bool> itsYRangesSet_;
	std::vector<prange_t> itsXRanges_;
	std::vector<prange_t> itsYRanges_;

	/* Key strings for casacore::Record */
	static const casacore::String REC_XAXES;
	static const casacore::String REC_YAXES;
	static const casacore::String REC_XRANGESSET;
	static const casacore::String REC_YRANGESSET;
	static const casacore::String REC_XRANGES;
	static const casacore::String REC_YRANGES;

	void setDefaults();
};



// Subclass of PlotMSPlotParameters::Group to handle canvas parameters.
// Currently includes:
// * axes label formats
// * whether to show the canvas axes or not
// * whether to show the legend or not, and its position
// * canvas title label format
// * whether to show grid lines, and their properties
// Parameters are vector-based, on a per-canvas basis.
//
class PMS_PP_Canvas : public PlotMSPlotParameters::Group {

public:
	/* Constructor which takes a factory */
	PMS_PP_Canvas (PlotFactoryPtr factory);

	/* Copy constructor.  See operator=(). */
	PMS_PP_Canvas (const PMS_PP_Canvas & copy);

	~PMS_PP_Canvas();

	/* Implements PlotMSPlotParameters::Group::clone(). */
	Group *clone() const {
		return new PMS_PP_Canvas (*this);
	}

	/* Implements PlotMSPlotParameters::Group::name(). */
	const casacore::String & name() const {
		static casacore::String groupName = PMS_PP::UPDATE_CANVAS_NAME;
		return groupName;
	}

	/* Implements PlotMSPlotParameters::Group::toRecord(). */
	casacore::Record toRecord() const;

	/* Implements PlotMSPlotParameters::Group::fromRecord(). */
	void fromRecord (const casacore::Record & record);

	/* Implements PlotMSPlotParameters::Group::requiresRedrawOnChanged(). */
	bool requiresRedrawOnChange() const {
		return true;
	}

	/* Overrides PlotMSPlotParameters::Group::operator=(). */
	PMS_PP_Canvas& operator= (const Group & other);

	/* Overrides the actual operator=(). */
	PMS_PP_Canvas& operator=(const PMS_PP_Canvas& other );

	/* Overrides PlotMSPlotParameters::Group::operator==(). */
	bool operator== (const Group & other) const;


	// Gets how many canvases there are.
	unsigned int numCanvases() const;

	// Sets single versions of the parameters for the given index.
	// <group>
	void setLabelFormats (const PlotMSLabelFormat & xFormat,
			const PlotMSLabelFormat & yFormat,
			unsigned int index = 0);
	void showAxes (const bool & xShow, const bool & yShow,
			unsigned int index = 0);
	void showLegend (const bool & show,
			const PlotCanvas::LegendPosition & pos,
			unsigned int index = 0);

	void showLegend( const bool& show, const casacore::String& pos, unsigned int index=0);

	void showGridMajor (const bool & show, const PlotLinePtr & line,
			unsigned int index = 0) {
		showGrid (show, gridMinorShown (index), line, gridMinorLine (index),
				index);
	}
	void showGridMinor (const bool & show, const PlotLinePtr & line,
			unsigned int index = 0) {
		showGrid (gridMajorShown (index), show, gridMajorLine (index), line,
				index);
	}
	void showGrid (const bool & showMajor, const bool & showMinor,
			const PlotLinePtr & majorLine,
			const PlotLinePtr & minorLine, unsigned int index = 0);
	// </group>



	const vector < PlotMSLabelFormat > &xLabelFormats() const {
		return itsXLabels_;
	}
	void setXLabelFormats (const vector < PlotMSLabelFormat > &value) {
		if (itsXLabels_ != value) {
			itsXLabels_ = value;
			updated();
		}
	}
	const PlotMSLabelFormat & xLabelFormat (unsigned int index = 0) const {
		return itsXLabels_[index];
	}
	void setXLabelFormat (const PlotMSLabelFormat & value,
			unsigned int index = 0) {
		if (itsXLabels_[index] != value) {
			itsXLabels_[index] = value;
			updated();
		}
	}


	const vector < PlotMSLabelFormat > &yLabelFormats() const {
		return itsYLabels_;
	}
	void setYLabelFormats (const vector < PlotMSLabelFormat > &value) {
		if (itsYLabels_ != value) {
			itsYLabels_ = value;
			updated();
		}
	}
	const PlotMSLabelFormat & yLabelFormat (unsigned int index = 0) const {
		return itsYLabels_[index];
	}
	void setYLabelFormat (const PlotMSLabelFormat & value,
			unsigned int index = 0) {
		if (itsYLabels_[index] != value) {
			itsYLabels_[index] = value;
			updated();
		}
	}

	const vector < bool > &xFontsSet() const {
		return itsXFontsSet_;
	}
	void setXFontsSet (const vector < bool > &value) {
		if (itsXFontsSet_ != value) {
			itsXFontsSet_ = value;
			updated();
		}
	}
	bool xFontSet (unsigned int index = 0) const {
		return itsXFontsSet_[index];
	}
	void setXFontSet (const bool & value, unsigned int index = 0) {
		if (itsXFontsSet_[index] != value) {
			itsXFontsSet_[index] = value;
			updated();
		}
	}
    const vector < bool > &yFontsSet() const {
		return itsYFontsSet_;
	}
	void setYFontsSet (const vector < bool > &value) {
		if (itsYFontsSet_ != value) {
			itsYFontsSet_ = value;
			updated();
		}
	}
	bool yFontSet (unsigned int index = 0) const {
		return itsYFontsSet_[index];
	}
	void setYFontSet (const bool & value, unsigned int index = 0) {
		if (itsYFontsSet_[index] != value) {
			itsYFontsSet_[index] = value;
			updated();
		}
	}

    const vector < casacore::Int > &xAxisFonts() const {
		return itsXAxisFonts_;
	}
	void setXAxisFonts (const vector < casacore::Int > &value) {
		if (itsXAxisFonts_ != value) {
			itsXAxisFonts_ = value;
			updated();
		}
	}
	const casacore::Int & xAxisFont (unsigned int index = 0) const {
		return itsXAxisFonts_[index];
	}
	void setXAxisFont (const casacore::Int value, unsigned int index = 0) {
		if (itsXAxisFonts_[index] != value) {
			itsXAxisFonts_[index] = value;
			updated();
		}
	}

    const vector < casacore::Int > &yAxisFonts() const {
		return itsYAxisFonts_;
	}
	void setYAxisFonts (const vector < casacore::Int > &value) {
		if (itsYAxisFonts_ != value) {
			itsYAxisFonts_ = value;
			updated();
		}
	}
	const casacore::Int & yAxisFont (unsigned int index = 0) const {
		return itsYAxisFonts_[index];
	}
	void setYAxisFont (const casacore::Int value, unsigned int index = 0) {
		if (itsYAxisFonts_[index] != value) {
			itsYAxisFonts_[index] = value;
			updated();
		}
	}


	const vector < bool > &xAxesShown() const {
		return itsXAxesShown_;
	}
	void showXAxes (const vector < bool > &value) {
		if (itsXAxesShown_ != value) {
			itsXAxesShown_ = value;
			updated();
		}
	}
	bool xAxisShown (unsigned int index = 0) const {
		if (index >= itsXAxesShown_.size())
			const_cast < vector < bool > &>(itsXAxesShown_).resize (index + 1);
		return itsXAxesShown_[index];
	}
	void showXAxis (const bool & value, unsigned int index = 0) {
		if (index >= itsXAxesShown_.size())
			itsXAxesShown_.resize (index + 1);
		if (itsXAxesShown_[index] != value) {
			itsXAxesShown_[index] = value;
			updated();
		}
	}


	const vector < bool > &yAxesShown() const {
		return itsYAxesShown_;
	}
	void showYAxes (const vector < bool > &value) {
		if (itsYAxesShown_ != value) {
			itsYAxesShown_ = value;
			updated();
		}
	}
	bool yAxisShown (unsigned int index = 0) const {
		if (index >= itsYAxesShown_.size())
			const_cast < vector < bool > &>(itsYAxesShown_).resize (index + 1);
		return itsYAxesShown_[index];
	}
	void showYAxis (const bool & value, unsigned int index = 0) {
		if (index >= itsYAxesShown_.size())
			itsYAxesShown_.resize (index + 1);
		if (itsYAxesShown_[index] != value) {
			itsYAxesShown_[index] = value;
			updated();
		}
	}


	const vector < bool > &legendsShown() const {
		return itsLegendsShown_;
	}
	void showLegends (const vector < bool > &value) {
		if (itsLegendsShown_ != value) {
			itsLegendsShown_ = value;
			updated();
		}
	}
	bool legendShown (unsigned int index = 0) const {
		if (index >= itsLegendsShown_.size())
			const_cast < vector < bool > &>(itsLegendsShown_).resize (index + 1);
		return itsLegendsShown_[index];
	}
	void showLegend (const bool & value, unsigned int index = 0) {
		if (index >= itsLegendsShown_.size())
			itsLegendsShown_.resize (index + 1);
		if (itsLegendsShown_[index] != value) {
			itsLegendsShown_[index] = value;
			updated();
		}
	}


	const vector < PlotCanvas::LegendPosition > &legendPositions() const {
		return itsLegendsPos_;
	}
	void setLegendPositions (const vector < PlotCanvas::LegendPosition > &value) {
		if (itsLegendsPos_ != value) {
			itsLegendsPos_ = value;
			updated();
		}
	}
	PlotCanvas::LegendPosition legendPosition (unsigned int index = 0) const {
		if (index >= itsLegendsPos_.size())
			const_cast < vector < PlotCanvas::LegendPosition >
		&>(itsLegendsPos_).resize (index + 1);
		return itsLegendsPos_[index];
	}
	void setLegendPosition (const PlotCanvas::LegendPosition & value,
			unsigned int index = 0) {
		if (index >= itsLegendsPos_.size())
			itsLegendsPos_.resize (index + 1);
		if (itsLegendsPos_[index] != value) {
			itsLegendsPos_[index] = value;
			updated();
		}
	}


	const vector < PlotMSLabelFormat > &titleFormats() const {
		return itsTitles_;
	}
	void setTitleFormats (const vector < PlotMSLabelFormat > &value) {
		if (itsTitles_ != value) {
			itsTitles_ = value;
			updated();
		}
	}
	const PlotMSLabelFormat & titleFormat (unsigned int index = 0) const {
		return itsTitles_[index];
	}
	void setTitleFormat (const PlotMSLabelFormat & value, unsigned int index =
			0) {
		if (itsTitles_[index] != value) {
			itsTitles_[index] = value;
			updated();
		}
	}
	const vector < bool > &titleFontsSet() const {
		return itsTitleFontsSet_;
	}
	void setTitleFontsSet (const vector < bool > &value) {
		if (itsTitleFontsSet_ != value) {
			itsTitleFontsSet_ = value;
			updated();
		}
	}
	bool titleFontSet (unsigned int index = 0) const {
		return itsTitleFontsSet_[index];
	}
	void setTitleFontSet (const bool & value, unsigned int index = 0) {
		if (itsTitleFontsSet_[index] != value) {
			itsTitleFontsSet_[index] = value;
			updated();
		}
	}
    const vector < casacore::Int > &titleFonts() const {
		return itsTitleFonts_;
	}
	void setTitleFonts (const vector < casacore::Int > &value) {
		if (itsTitleFonts_ != value) {
			itsTitleFonts_ = value;
			updated();
		}
	}
	const casacore::Int & titleFont (unsigned int index = 0) const {
		return itsTitleFonts_[index];
	}
	void setTitleFont (const casacore::Int value, unsigned int index =
			0) {
		if (itsTitleFonts_[index] != value) {
			itsTitleFonts_[index] = value;
			updated();
		}
	}

	const vector < bool > &gridMajorsShown() const {
		return itsGridMajsShown_;
	}
	void showGridMajors (const vector < bool > &value) {
		if (itsGridMajsShown_ != value) {
			itsGridMajsShown_ = value;
			updated();
		}
	}
	bool gridMajorShown (unsigned int index = 0) const {
		if (index >= itsGridMajsShown_.size())
			const_cast < vector < bool > &>(itsGridMajsShown_).resize (index + 1);
		return itsGridMajsShown_[index];
	}
	void showGridMajor (const bool & value, unsigned int index = 0) {
		if (index >= itsGridMajsShown_.size())
			itsGridMajsShown_.resize (index + 1);
		if (itsGridMajsShown_[index] != value) {
			itsGridMajsShown_[index] = value;
			updated();
		}
	}


	const vector < bool > &gridMinorsShown() const {
		return itsGridMinsShown_;
	}
	void showGridMinors (const vector < bool > &value) {
		if (itsGridMinsShown_ != value) {
			itsGridMinsShown_ = value;
			updated();
		}
	}
	bool gridMinorShown (unsigned int index = 0) const {
		if (index >= itsGridMinsShown_.size())
			const_cast < vector < bool > &>(itsGridMinsShown_).resize (index + 1);
		return itsGridMinsShown_[index];
	}
	void showGridMinor (const bool & value, unsigned int index = 0) {
		if (index >= itsGridMinsShown_.size())
			itsGridMinsShown_.resize (index + 1);
		if (itsGridMinsShown_[index] != value) {
			itsGridMinsShown_[index] = value;
			updated();
		}
	}


	const vector < PlotLinePtr > &gridMajorLines() const {
		return itsGridMajLines_;
	}
	void setGridMajorLines (const vector < PlotLinePtr > &value) {
		if (itsGridMajLines_ != value) {
			itsGridMajLines_ = value;
			updated();
		}
	}
	PlotLinePtr gridMajorLine (unsigned int index = 0) const {
		if (index >= itsGridMajLines_.size())
			const_cast < vector < PlotLinePtr >
		&>(itsGridMajLines_).resize (index + 1);
		return itsGridMajLines_[index];
	}
	void setGridMajorLine (const PlotLinePtr & value, unsigned int index = 0) {
		if (index >= itsGridMajLines_.size())
			itsGridMajLines_.resize (index + 1);
		if (itsGridMajLines_[index] != value) {
			itsGridMajLines_[index] = value;
			updated();
		}
	}


	const vector < PlotLinePtr > &gridMinorLines() const {
		return itsGridMinLines_;
	}
	void setGridMinorLines (const vector < PlotLinePtr > &value) {
		if (itsGridMinLines_ != value) {
			itsGridMinLines_ = value;
			updated();
		}
	}
	PlotLinePtr gridMinorLine (unsigned int index = 0) const {
		if (index >= itsGridMinLines_.size())
			const_cast < vector < PlotLinePtr >
		&>(itsGridMinLines_).resize (index + 1);
		return itsGridMinLines_[index];
	}
	void setGridMinorLine (const PlotLinePtr & value, unsigned int index = 0) {
		if (index >= itsGridMinLines_.size())
			itsGridMinLines_.resize (index + 1);
		if (itsGridMinLines_[index] != value) {
			itsGridMinLines_[index] = value;
			updated();
		}
	}


private:

	//Does the work for the operator=()s.
	PMS_PP_Canvas& assign(const PMS_PP_Canvas* o );

	/* Parameters' values */
	std::vector<PlotMSLabelFormat> itsXLabels_;
	std::vector<bool> itsXFontsSet_;
	std::vector<casacore::Int> itsXAxisFonts_;
	std::vector<PlotMSLabelFormat> itsYLabels_;
	std::vector<bool> itsYFontsSet_;
	std::vector<casacore::Int> itsYAxisFonts_;
	std::vector<bool> itsXAxesShown_;
	std::vector<bool> itsYAxesShown_;
	std::vector<bool> itsLegendsShown_;
	std::vector<PlotCanvas::LegendPosition > itsLegendsPos_;
	std::vector<PlotMSLabelFormat> itsTitles_;
	std::vector<bool> itsTitleFontsSet_;
	std::vector<casacore::Int> itsTitleFonts_;
	std::vector<bool> itsGridMajsShown_;
	std::vector<bool> itsGridMinsShown_;
	std::vector<PlotLinePtr> itsGridMajLines_;
	std::vector<PlotLinePtr> itsGridMinLines_;

	/* Key strings for casacore::Record */
	static const casacore::String REC_XLABELS;
	static const casacore::String REC_XFONTSSET;
	static const casacore::String REC_XAXISFONTS;
	static const casacore::String REC_YLABELS;
	static const casacore::String REC_YFONTSSET;
	static const casacore::String REC_YAXISFONTS;
	static const casacore::String REC_SHOWXAXES;
	static const casacore::String REC_SHOWYAXES;
	static const casacore::String REC_SHOWLEGENDS;
	static const casacore::String REC_LEGENDSPOS;
	static const casacore::String REC_TITLES;
	static const casacore::String REC_TITLEFONTSSET;
	static const casacore::String REC_TITLEFONTS;
	static const casacore::String REC_SHOWGRIDMAJS;
	static const casacore::String REC_SHOWGRIDMINS;
	static const casacore::String REC_GRIDMAJLINES;
	static const casacore::String REC_GRIDMINLINES;


	void setDefaults();
};





// Subclass of PlotMSPlotParameters::Group to handle display parameters.
// Currently includes:
// * flagged and unflagged symbols
// * plot title format
// * colorize flag and axis
// Parameters are vector-based, on a per-plot basis.
//
class PMS_PP_Display : public PlotMSPlotParameters::Group {

public:
	/* Constructor which takes a factory */
	PMS_PP_Display (PlotFactoryPtr factory);

	/* Copy constructor.  See operator=(). */
	PMS_PP_Display (const PMS_PP_Display & copy);

	~PMS_PP_Display();


	/* Implements PlotMSPlotParameters::Group::clone(). */
	Group *clone() const {
		return new PMS_PP_Display (*this);
	}

	/* Implements PlotMSPlotParameters::Group::name(). */
	const casacore::String & name() const {
		static casacore::String groupName = PMS_PP::UPDATE_DISPLAY_NAME;
		return groupName;
	}

	/* Implements PlotMSPlotParameters::Group::toRecord(). */
	casacore::Record toRecord() const;

	/* Implements PlotMSPlotParameters::Group::fromRecord(). */
	void fromRecord (const casacore::Record & record);

	/* Implements PlotMSPlotParameters::Group::requiresRedrawOnChanged(). */
	bool requiresRedrawOnChange() const {
		return true;
	}

	/* Overrides PlotMSPlotParameters::Group::operator=(). */
	PMS_PP_Display & operator= (const Group & other);

	/*  Overrides the actual assignment operator=() */
	PMS_PP_Display& operator=(const PMS_PP_Display& other);

	/* Overrides PlotMSPlotParameters::Group::operator==(). */
	bool operator== (const Group & other) const;


	void setColorize (const bool & colorize, const PMS::Axis & axis,
			unsigned int index = 0);


	void resizeVectors (unsigned int newSize);


	const vector < PlotSymbolPtr > &unflaggedSymbols() const {
		return itsUnflaggedSymbols_;
	}
	void setUnflaggedSymbols (const vector < PlotSymbolPtr > &value) {
		if (itsUnflaggedSymbols_ != value) {
			itsUnflaggedSymbols_ = value;
			updated();
		}
	}
	PlotSymbolPtr unflaggedSymbol (unsigned int index = 0) const {
		if (index >= itsUnflaggedSymbols_.size()){
			int newSize = index+1;
			std::vector<PlotSymbolPtr> & unflaggedSymbols = const_cast < vector <PlotSymbolPtr > &>(itsUnflaggedSymbols_);
			unflaggedSymbols.resize ( newSize);
			for ( int j = 0; j < newSize; j++ ){
				if ( unflaggedSymbols[j].null() ){
					unflaggedSymbols[j]=PMS::DEFAULT_UNFLAGGED_SYMBOL(factory());
				}
			}
		}
		return itsUnflaggedSymbols_[index];
	}
	void setUnflaggedSymbol (const PlotSymbolPtr & value, unsigned int index =0);

	const vector < PlotSymbolPtr > &flaggedSymbols() const {
		return itsFlaggedSymbols_;
	}
	void setFlaggedSymbols (const vector < PlotSymbolPtr > &value) {
		if (itsFlaggedSymbols_ != value) {
			itsFlaggedSymbols_ = value;
			updated();
		}
	}
	PlotSymbolPtr flaggedSymbol (unsigned int index = 0) const {
		if (index >= itsFlaggedSymbols_.size()){
			int newSize = index + 1;
			std::vector<PlotSymbolPtr> & flaggedSymbols = const_cast < vector <PlotSymbolPtr > &>(itsFlaggedSymbols_);
			flaggedSymbols.resize( newSize );
			for ( int j = 0; j < newSize; j++ ){
				if ( flaggedSymbols[j].null()){
					flaggedSymbols[j] = PMS::DEFAULT_FLAGGED_SYMBOL(factory());
				}
			}
		}
		return itsFlaggedSymbols_[index];
	}
	void setFlaggedSymbol (const PlotSymbolPtr & value, unsigned int index =0);

	const vector < PlotMSLabelFormat > &titleFormats() const {
		return itsTitleFormats_;
	}
	void setTitleFormats (const vector < PlotMSLabelFormat > &value) {
		if (itsTitleFormats_ != value) {
			itsTitleFormats_ = value;
			updated();
		}
	}
	const PlotMSLabelFormat & titleFormat (unsigned int index = 0) const {
		return itsTitleFormats_[index];
	}
	void setTitleFormat (const PlotMSLabelFormat & value, unsigned int index =
			0) {
		if (itsTitleFormats_[index] != value) {
			itsTitleFormats_[index] = value;
			updated();
		}
	}


	const vector < bool > &colorizeFlags() const {
		return itsColorizeFlags_;
	}
	void setColorize (const vector < bool > &value) {
		if (itsColorizeFlags_ != value) {
			itsColorizeFlags_ = value;
			updated();
		}
	}
	bool colorizeFlag (unsigned int index = 0) const {
		if (index >= itsColorizeFlags_.size())
			const_cast < vector < bool > &>(itsColorizeFlags_).resize (index + 1);
		return itsColorizeFlags_[index];
	}
	void setColorize (const bool & value, unsigned int index = 0);

	const vector < PMS::Axis > &colorizeAxes() const {
		return itsColorizeAxes_;
	}
	void setColorize (const vector < PMS::Axis > &value) {
		if (itsColorizeAxes_ != value) {
			itsColorizeAxes_ = value;
			updated();
		}
	}
	PMS::Axis colorizeAxis (unsigned int index = 0) const {
		if (index >= itsColorizeAxes_.size())
			const_cast < vector < PMS::Axis >
		&>(itsColorizeAxes_).resize (index + 1);
		return itsColorizeAxes_[index];
	}
	void setColorize (const PMS::Axis & value, unsigned int index = 0) {
		if (index >= itsColorizeAxes_.size())
			itsColorizeAxes_.resize (index + 1);
		if (itsColorizeAxes_[index] != value) {
			itsColorizeAxes_[index] = value;
			updated();
		}
	}




private:

	/* Does the work for both versions of operator=() */
	PMS_PP_Display& assign( const PMS_PP_Display* o );


	/* Parameters' values */
	std::vector<PlotSymbolPtr> itsUnflaggedSymbols_;
	std::vector<PlotSymbolPtr> itsFlaggedSymbols_;
	std::vector<PlotMSLabelFormat> itsTitleFormats_;
	std::vector<bool> itsColorizeFlags_;
	std::vector<PMS::Axis> itsColorizeAxes_;


	/* Key strings for casacore::Record */
	static const casacore::String REC_UNFLAGGEDS;
	static const casacore::String REC_FLAGGEDS;
	static const casacore::String REC_TITLES;
	static const casacore::String REC_COLFLAGS;
	static const casacore::String REC_COLAXES;

	void setDefaults();
};

// Subclass of PlotMSPlotParameters::Group to handle page header parameters.
// Currently includes:
// * page header items

// Parameters are vector-based, on a per-plot basis.
//
class PMS_PP_PageHeader : public PlotMSPlotParameters::Group {

public:
	/* Constructor which takes a factory */
	PMS_PP_PageHeader (PlotFactoryPtr factory);

	/* Copy constructor.  See operator=(). */
	PMS_PP_PageHeader (const PMS_PP_PageHeader & copy);

	~PMS_PP_PageHeader();

	/* Implements PlotMSPlotParameters::Group::clone(). */
	Group *clone() const {
		return new PMS_PP_PageHeader (*this);
	}

	/* Implements PlotMSPlotParameters::Group::name(). */
	const casacore::String & name() const {
		static casacore::String groupName = PMS_PP::UPDATE_PAGEHEADER_NAME;
		return groupName;
	}

	/* Implements PlotMSPlotParameters::Group::toRecord(). */
	casacore::Record toRecord() const;

	/* Implements PlotMSPlotParameters::Group::fromRecord(). */
	void fromRecord (const casacore::Record & record);

	/* Overrides PlotMSPlotParameters::Group::operator=(). */
	PMS_PP_PageHeader & operator= (const Group & other);

	/*  Overrides the actual assignment operator=() */
	PMS_PP_PageHeader& operator=(const PMS_PP_PageHeader& other);

	/* Overrides PlotMSPlotParameters::Group::operator==(). */
	bool operator== (const Group & other) const;

	/* Implements PlotMSPlotParameters::Group::requiresRedrawOnChanged(). */
	bool requiresRedrawOnChange() const {
		return false;
	}

	void setPageHeaderItems (const PageHeaderItems & value) {
		if (itsPageHeaderItems_ != value) {
			itsPageHeaderItems_ = value;
			updated();
		}
	}

	const PageHeaderItems& pageHeaderItems() const { return itsPageHeaderItems_;  }

private:
	/* Parameters' values */
	PageHeaderItems itsPageHeaderItems_;

	/* Key strings for casacore::Record */
	static const casacore::String REC_ITEMS;

	//Does the work for the operator=()s.
	PMS_PP_PageHeader& assign( const PMS_PP_PageHeader* o );
	void setDefaults();

};


// Subclass of PlotMSPlotParameters::Group to handle iteration parameters.
// Currently includes:
// * which axis to use
// * rows, columns to fit onto each page
// Parameters apply to whole set of plots
//
class PMS_PP_Iteration : public PlotMSPlotParameters::Group {



public:

	/* Constructor which takes a factory */
	PMS_PP_Iteration (PlotFactoryPtr factory);

	/* Copy constructor.  See operator=(). */
	PMS_PP_Iteration (const PMS_PP_Iteration & copy);

	~PMS_PP_Iteration();


	/* Implements PlotMSPlotParameters::Group::clone(). */
	Group *clone() const {
		return new PMS_PP_Iteration (*this);
	}

	/* Implements PlotMSPlotParameters::Group::name(). */
	const casacore::String & name() const {
		static casacore::String groupName = PMS_PP::UPDATE_ITERATION_NAME;
		return groupName;
	}

	/* Implements PlotMSPlotParameters::Group::toRecord(). */
	casacore::Record toRecord() const;

	/* Implements PlotMSPlotParameters::Group::fromRecord(). */
	void fromRecord (const casacore::Record & record);

	/* Implements PlotMSPlotParameters::Group::requiresRedrawOnChanged(). */
	bool requiresRedrawOnChange() const {
		return true;
	}

	/* Overrides the real operator=().  */
	PMS_PP_Iteration& operator=(const PMS_PP_Iteration& other);

	/* Overrides PlotMSPlotParameters::Group::operator=(). */
	PMS_PP_Iteration & operator= (const Group & other);

	/* Overrides PlotMSPlotParameters::Group::operator==(). */
	bool operator== (const Group & other) const;

	//Returns whether or not we are iterating on an axis.
	bool isIteration() const;

	const PlotMSIterParam& iterParam() const {
		return itsIterParam_;
	}
	void setIterParam(PlotMSIterParam iterparam) {
		if (itsIterParam_ != iterparam) {
			itsIterParam_=iterparam;
			updated();
		}
	}

	PMS::Axis  iterationAxis() const {
		return itsIterParam_.iterAxis();
	}

	void setIterationAxis (const PMS::Axis & value) {
		if (iterationAxis()!=value) {
			itsIterParam_.setIterAxis(value);
			updated();
		}
	}


	int getGridRow() const  {
		return itsIterParam_.getGridRow();
	}
	void setGridRow(const int &value)  {
		if (getGridRow() != value) {
			itsIterParam_.setGridRow(value);
			updated();
		}
	}

	int getGridCol() const  {
		return itsIterParam_.getGridCol();
	}
	void setGridCol(const int &value)  {
		if (getGridCol() != value) {
			itsIterParam_.setGridCol(value);
			updated();
		}
	}

	casacore::Bool isCommonAxisX() const {
		return itsIterParam_.isCommonAxisX();
	}
	void setCommonAxisX( bool commonAxis ){
		if ( isCommonAxisX() != commonAxis ){
			bool validValue = false;
			if ( isGlobalScaleX()){
				validValue = true;
			}
			else if ( !commonAxis ){
				validValue = true;
			}

			if ( validValue ){
				itsIterParam_.setCommonAxisX( commonAxis );
				updated();
			}
		}
	}
	casacore::Bool isCommonAxisY() const {
		return itsIterParam_.isCommonAxisY();
	}
	void setCommonAxisY( bool commonAxis ){
		if ( isCommonAxisY() != commonAxis ){
			bool validValue = false;
			if ( isGlobalScaleY() ){
				validValue = true;
			}
			else if ( !commonAxis ){
				validValue = true;
			}
			if ( validValue ){
				itsIterParam_.setCommonAxisY( commonAxis );
				updated();
			}
		}
	}
	casacore::Bool isGlobalScaleX() const {
			return itsIterParam_.isGlobalAxisX();
		}
		void setGlobalScaleX( bool globalAxis ){
			if ( isGlobalScaleX() != globalAxis ){
				itsIterParam_.setGlobalScaleX( globalAxis );
				updated();
			}
		}
		casacore::Bool isGlobalScaleY() const {
			return itsIterParam_.isGlobalAxisY();
		}
		void setGlobalScaleY( bool globalAxis ){
			if ( isGlobalScaleY() != globalAxis ){
				itsIterParam_.setGlobalScaleY( globalAxis );
				updated();
			}
		}
private:

	//Does the work for the operator=()s.
	PMS_PP_Iteration& assign(const PMS_PP_Iteration* o);

	/* Parameters' values */
	PlotMSIterParam itsIterParam_;
	void setDefaults();
};




}

#endif /* PLOTMSPLOTPARAMETERGROUPS_H_ */
