//# TblAsRasterDD.h: Display casacore::Data for raster displays of data from a table
//# Copyright (C) 2000,2001,2002
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

#ifndef TRIALDISPLAY_TBLASRASTERDD_H
#define TRIALDISPLAY_TBLASRASTERDD_H

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <display/Display/DParameterRange.h>
#include <display/Display/DParameterChoice.h>
#include <display/Display/DParameterString.h>
#include <coordinates/Coordinates/LinearCoordinate.h>
#include <display/DisplayDatas/ActiveCaching2dDD.h>

namespace casacore{

	class Table;
	class Regex;
}

namespace casa { //# NAMESPACE CASA - BEGIN

	class TblAsRasterDM;

//# Forward Declarations

// <summary>
//Class for displaying data within a table as a raster image.
// </summary>

// <use visibility=local>   or   <use visibility=export>

// <reviewed reviewer="" date="" tests="" demos="">
// </reviewed>

// <prerequisite>
//   <li> ActiveCaching2dDD
//   <li> CachingDisplayData
//   <li> Table
// </prerequisite>
//
// <etymology>
// "TblAsRasterDD" is a implementation of a <linkto class=ActiveCaching2dDD>
// ActiveCaching2dDD </linkto> which provides for the display of data held
// within a table to be displayed as a raster image in an environment where
// individual depictions of the data are automatically cached.
// </etymology>
//
// <synopsis>
// This class adds to the interface defined in <linkto
// class=DisplayData>DisplayData </linkto>.  It adds the capability to
// display vector/array data from a <linkto class=casacore::Table>Table</linkto>
// column as a raster image.  It is assumed that the Y axis is defined
// to be either the row number of the table column being displayed or
// the scalar value from the same row number in a different table
// column (e.g. plotting intensity as a function of frequency against
// row number or time determined from a different column of the table).
// The X axis is assumed to be a one dimensional array or vector of
// of data for each row in the column being displayed.  It is assumed
// that the length of this array/vector does not change throughout the
// column of the table.
// </synopsis>
//
// <example>
// A TblAsRasterDD object could be construed and used as follows:
// <srcblock>
//    TblAsRasterDD *tardd1 = 0;
//    tardd1 = new TblAsRasterDD("tablename");
//    Colormap cmap1("Hot Metal 2");
//    tardd1->setColormap(&cmap1, 1.0);
//    wcHolder->addDisplayData((DisplayData *)tardd1);
// </srcblock>
// </example>
//
// <motivation>
// To allow the display of data from a table as a raster image.
// </motivation>
//
// <templating arg=T>
// </templating>
//
// <thrown>
// </thrown>
//
// <todo asof="2000/10/30">
//   <li> make sure complex data are handled correctly
//   <li> make sure table column units can be determined properly
//   <li> make sure a scalar table column can be used for y axis
//   <li> extend to n-dimensional arrays in table column
//   <li> when constructed with *table we need to keep table from being deleted
//   <li> handle movie axis once ActiveCachingNDim exists
//   <li> implement showValue()
// </todo>

	class TblAsRasterDD : public ActiveCaching2dDD {

	public:

		// constructors
		// given an already constructed table
		TblAsRasterDD(casacore::Table *table);

		// given a string which gives the full pathname and filename of a table
		// on disk
		TblAsRasterDD(const casacore::String tablename);

		// Destructor
		virtual ~TblAsRasterDD();

		// format the table value at the give world position
		virtual casacore::String showValue(const casacore::Vector<casacore::Double> &world);

		// get the data unit
		virtual const casacore::Unit dataUnit(const casacore::String column);
		virtual const casacore::Unit dataUnit();

		// install the default options for this DisplayData
		virtual void setDefaultOptions();

		// Apply options stored in <src>rec</src> to the DisplayData.  A
		// return value of <src>true</src> means a refresh is needed.
		// <src>recOut</src> contains any fields which were implicitly
		// changed as a result of the call to this function.
		virtual casacore::Bool setOptions(casacore::Record &rec, casacore::Record &recOut);

		// Retrieve the current and default options and parameter types.
		virtual casacore::Record getOptions( bool scrub=false ) const;

		// Return the type of this DisplayData.
		virtual Display::DisplayDataType classType() {
			return Display::Raster;
		}

		// Create a new TblAsRasterDM for drawing on the given
		// WorldCanvas when the AttributeBuffers are suitably matched to the
		// current state of this DisplayData and of the WorldCanvas/Holder.
		// The tag is a unique number used to identify the age of the newly
		// constructed CachingDisplayMethod.
		virtual CachingDisplayMethod *newDisplayMethod(WorldCanvas *worldCanvas,
		        AttributeBuffer *wchAttributes,
		        AttributeBuffer *ddAttributes,
		        CachingDisplayData *dd);

		// Return the current options of this DisplayData as an
		// AttributeBuffer.
		virtual AttributeBuffer optionsAsAttributes();

		//provide read-only access to the table
		casacore::Table *table();

		// Clean up (ie. delete any existing cached display list).
		virtual void cleanup();

	protected:

		// (Required) default constructor.
		TblAsRasterDD();

		// (Required) copy constructor.
		TblAsRasterDD(const TblAsRasterDD &other);

		// (Required) copy assignment.
		void operator=(const TblAsRasterDD &other);

// Get the value of the named keyword, or the first keyword matching
		// <src>regex</src>, and return it in <src>value</src>.  The return
		// value is <src>true</src> for success, and <src>false</src> for
		// failure, which is the result if the wrong type <src>T</src> is
		// requested.
		// <group>
		template <class T> casacore::Bool getTableKeyword(T &value,
		                                        const casacore::String keyword) const;
		template <class T> casacore::Bool getTableKeyword(T &value, const casacore::Regex &regex) const;
		// </group>

		// Get the value of the named keyword, or the first keyword matching
		// <src>regex</src> for the named column, and return it in
		// <src>value</src>. The return value is <src>true</src> for
		// success, and <src>false</src> for failure, which is the result if           // the wrong type <src>T</src> is requested, or if the keyword
		// doesn't exist.
		// <group>
		template <class T> casacore::Bool getColumnKeyword(T &value, const casacore::String column,
		        const casacore::String keyword) const;
		template <class T> casacore::Bool getColumnKeyword(T &value, const casacore::String column,
		        const casacore::Regex &regex) const;
		// </group>
	private:

		friend class TblAsRasterDM;

		// The table to be displayed
		casacore::Table *itsTable;

		// The result from a table query
		casacore::Table *itsQueryTable;

		// store all the table column names
		casacore::Vector<casacore::String> itsColumnNames;

		// what columns are we displaying and do we have a movie axis available
		DParameterChoice *itsXColumnName;
		DParameterChoice *itsYColumnName;
		DParameterChoice *itsMColumnName;
		DParameterChoice *itsMColumnSet;

		// options - what is the query string and is it unset?
		casacore::String itsOptQueryString;
		casacore::Bool itsOptQueryStringUnset;

		// set the default options for this display data
		void installDefaultOptions();

		// Arrange the query table (called after changing an option).
		casacore::Bool arrangeQueryTable();

		// holder for the current coordinate system
		DisplayCoordinateSystem itsCoord;
		casacore::Vector<casacore::Double> itsLinblc, itsLintrc;

		// update/set the coordinate system
		void getCoordinateSystem();
		void setCoordinateSystem();

		// get all of the table columnNames
		void getTableColumnNames();

		// get the table column world coordinate range
		casacore::Vector<double> columnStatistics(const casacore::String& columnName);

		// get all of the table columnNames with a certain data type
		casacore::Vector<casacore::String> getColumnNamesOfType(const casacore::Bool isarray);

		// Construct and destruct the parameter set.
		// <group>
		void constructParameters();
		void destructParameters();
		// </group>

	};



} //# NAMESPACE CASA - END

#ifndef AIPS_NO_TEMPLATE_SRC
#include <display/DisplayDatas/TblAsRasterDDTemplates.tcc>
#endif //# AIPS_NO_TEMPLATE_SRC
#endif


