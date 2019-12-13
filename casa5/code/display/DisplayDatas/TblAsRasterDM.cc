//# TblAsRasterDM.cc:  Display Method for raster displays of data from tables
//# Copyright (C) 2000,2001
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
//# $Id$

#include <casa/aips.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/ColumnDesc.h>
#include <tables/TaQL/TableParse.h>
#include <casa/Arrays/Vector.h>
#include <casa/Containers/Record.h>
#include <display/Display/Attribute.h>
#include <casa/Utilities/DataType.h>
#include <tables/Tables/ArrayColumn.h>
#include <casa/Arrays/ArrayMath.h>
#include <display/DisplayDatas/TblAsRasterDD.h>
#include <display/DisplayDatas/TblAsRasterDM.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

// constructor
	TblAsRasterDM::TblAsRasterDM(WorldCanvas *worldCanvas,
	                             AttributeBuffer *wchAttributes,
	                             AttributeBuffer *ddAttributes,
	                             CachingDisplayData *dd) :
		CachingDisplayMethod(worldCanvas, wchAttributes, ddAttributes, dd) {
	}

// destructor
	TblAsRasterDM::~TblAsRasterDM() {
		cleanup();
	}

// cleanup function
	void TblAsRasterDM::cleanup() {
	}

	Bool TblAsRasterDM::drawIntoList(Display::RefreshReason /*reason*/,
	                                 WorldCanvasHolder &wcHolder) {

		// which world canvas do we draw on?
		WorldCanvas *wc = wcHolder.worldCanvas();

		// can we do a dynamic cast to the correct data type?
		TblAsRasterDD *parent = dynamic_cast<TblAsRasterDD *>
		                        (parentDisplayData());
		if (!parent) {
			throw(AipsError("invalid parent of TblAsRasterDM"));
		}

		// get the table column data type
		TableDesc tdesc(parent->table()->tableDesc());
		DataType type =
		    tdesc.columnDesc(parent->itsXColumnName->value()).trueDataType();

		// array to contain data to be plotted
		Array<float> data;
		Table *theTable = parent->table();

		if (type == TpArrayDouble) {
			Array<double> typedata;
			// read the column into an array
			ArrayColumn<double>
			dataCol(*theTable,parent->itsXColumnName->value());
			dataCol.getColumn(typedata,true);
			// now convert array to type float
			data.resize(typedata.shape());
			convertArray(data,typedata);
		}
		if (type == TpArrayFloat) {
			ArrayColumn<float>
			dataCol(*theTable,parent->itsXColumnName->value());
			dataCol.getColumn(data,true);
		}
		if (type == TpArrayUShort) {
			Array<ushort> typedata;
			ArrayColumn<ushort>
			dataCol(*theTable,parent->itsXColumnName->value());
			dataCol.getColumn(typedata,true);
			data.resize(typedata.shape());
			convertArray(data,typedata);
		}
		if (type == TpArrayInt) {
			Array<int> typedata;
			ArrayColumn<int>
			dataCol(*theTable,parent->itsXColumnName->value());
			dataCol.getColumn(typedata,true);
			data.resize(typedata.shape());
			convertArray(data,typedata);
		}
		if (type == TpArrayUInt) {
			Array<uInt> typedata;
			ArrayColumn<uInt>
			dataCol(*theTable,parent->itsXColumnName->value());
			dataCol.getColumn(typedata,true);
			data.resize(typedata.shape());
			convertArray(data,typedata);
		}
		// put the data into the matrix
		const Matrix<float> theData = data;

		// define several things for drawing
		Bool usePixelEdges = true;

		// now plot the data
		wc->drawImage(parent->itsLinblc, parent->itsLintrc, theData,
		              usePixelEdges);
		return true;
	}

// (required) default constructor
	TblAsRasterDM::TblAsRasterDM() {
	}

// (required) copy constructor
	TblAsRasterDM::TblAsRasterDM(const TblAsRasterDM & other) :
		CachingDisplayMethod(other) {
	}

// (required) copy assignment
	void TblAsRasterDM::operator=(const TblAsRasterDM &) {
	}


} //# NAMESPACE CASA - END

