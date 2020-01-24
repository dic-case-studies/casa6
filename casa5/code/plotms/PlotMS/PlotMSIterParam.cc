//# PlotMSIterateParams.cc: Iteration paremters container implementation
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
#include <plotms/PlotMS/PlotMSIterParam.h>
#include <QDebug>

using namespace casacore;
namespace casa {

///////////////////////////////////////
// PLOTMSITERATIONPARM DEFINITIONS //
///////////////////////////////////////
const String PlotMSIterParam::ITER_AXIS = "iterAxis";
const String PlotMSIterParam::GLOBAL_SCALE_X = "globalScaleX";
const String PlotMSIterParam::GLOBAL_SCALE_Y = "globalScaleY";
const String PlotMSIterParam::COMMON_AXIS_X = "commonAxisX";
const String PlotMSIterParam::COMMON_AXIS_Y = "commonAxisY";
const String PlotMSIterParam::ROW_INDEX = "rowIndex";
const String PlotMSIterParam::COL_INDEX = "colIndex";

// Non-Static //

PlotMSIterParam::PlotMSIterParam()
{ setDefaults(); }
PlotMSIterParam::~PlotMSIterParam() { }


void PlotMSIterParam::fromRecord(const RecordInterface& record) {    

	// Set defaults first (in case any missing bits in record)
	setDefaults();

	// Set from the Record
	if (record.isDefined(ITER_AXIS)){
		setIterAxis(record.asString(ITER_AXIS));
	}

	if (record.isDefined(ROW_INDEX)){
		setGridRow(record.asInt(ROW_INDEX));
	}

	if (record.isDefined(COL_INDEX)){
		setGridCol(record.asInt(COL_INDEX));
	}

	if (record.isDefined(GLOBAL_SCALE_X)){
		setGlobalScaleX(record.asBool(GLOBAL_SCALE_X));
	}

	if (record.isDefined(GLOBAL_SCALE_Y)){
		setGlobalScaleY(record.asBool(GLOBAL_SCALE_Y));
	}

	if (record.isDefined(COMMON_AXIS_X)){
		setCommonAxisX(record.asBool(COMMON_AXIS_X));
	}

	if (record.isDefined(COMMON_AXIS_Y)){
		setCommonAxisY(record.asBool(COMMON_AXIS_Y));
	}

}

Record PlotMSIterParam::toRecord() const {

	// Fill a record
	Record rec(Record::Variable);
	rec.define(ITER_AXIS,iterAxisStr());
	rec.define(COMMON_AXIS_X, commonAxisX_);
	rec.define(COMMON_AXIS_Y, commonAxisY_);
	rec.define(GLOBAL_SCALE_X, globalScaleX_);
	rec.define(GLOBAL_SCALE_Y, globalScaleY_);
	rec.define(ROW_INDEX,gridRow);
	rec.define(COL_INDEX,gridCol);

	// Return it
	return rec;

}



bool PlotMSIterParam::operator==(const PlotMSIterParam& other) const {

	return (iterAxis_ == other.iterAxis_ &&
			gridRow == other.gridRow &&
			gridCol == other.gridCol &&
			commonAxisX_ == other.commonAxisX_ &&
			commonAxisY_ == other.commonAxisY_ &&
			globalScaleX_ == other.globalScaleX_ &&
			globalScaleY_ == other.globalScaleY_);
}


void PlotMSIterParam::setDefaults() {    

	setIterAxis(PMS::NONE);  // No iteration
	setGridRow(0);
	setGridCol(0);
	setGlobalScaleX(false);
	setGlobalScaleY(false);
	setCommonAxisX(false);
	setCommonAxisY(false);
}

bool PlotMSIterParam::isIteration() const {
	bool iterationSet = false;
	if ( iterAxis_ != PMS::NONE ){
		iterationSet = true;
	}
	return iterationSet;
}

String PlotMSIterParam::summary() const {

	stringstream ss;
	ss << boolalpha;

	ss << "Iteration parameters:" << endl;
	ss << " Iteration Axis = " << iterAxisStr() << endl;

	ss << " Grid Row        = " << gridRow << endl;
	ss << " Grid Col        = " << gridCol << endl;
	ss << " Global Scale X  = " << globalScaleX_ << endl;
	ss << " Global Scale Y  = " << globalScaleY_ << endl;
	ss << " Common Axis X   = " << commonAxisX_ << endl;
	ss << " Common Axis Y   = " << commonAxisY_ << endl;
	return ss.str();
}


}
