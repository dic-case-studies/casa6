//# TBArray.h: Holds a potentially multi-dimensional array.
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
//# $Id: $
#ifndef TBARRAY_H_
#define TBARRAY_H_

#include <casa/BasicSL/String.h>

#include <vector>

namespace casa {

//# Forward Declarations
class TBTable;

// <summary>
// Holds a potentially multi-dimensional array.
// <summary>
//
// <synopsis>
// A TBArray holds a array object as used by the table browser.  The data
// is represented by casacore::String values, while the array structure is represented by
// vectors of void*s.  For all but the "last" dimension, the void*s point to
// other vector<void*>s, and on the last dimension the void*s point to Strings.
// NOTE: this class is mostly obsolete now that the browser uses the TBData
// structure.
// </synopsis>

class TBArray {
public:
    // Constructor for a data array (in other words, an array found in the data
    // of a table).  Takes as input the row and column of the array in table,
    // the type of table, and the casacore::String holding the values of the array which
    // need to be parsed.
    TBArray(int row, int col, casacore::String type, casacore::String array);

    // Constructor for a table keyword array (in other words, an array found in
    // the table keywords).  Takes as input the index of the keyword in the
    // table keywords list, the type of the array, and the casacore::String holding the
    // values of the array which need to be parsed.
    TBArray(int keywordIndex, casacore::String type, casacore::String array);

    // Constructor for a field keyword array (in other words, an array found in
    // the field keywords).  Takes as input the name of the field, the index of
    // the keyword in the field keywords list, the type of the array, and the
    // casacore::String holding the values of the array which need to be parsed.
    TBArray(casacore::String col, int index, casacore::String type, casacore::String array);
    
    // Constructor for a non-specific array.  Takes as input the type of the
    // array and the casacore::String holding the values of the array which need to be
    // parsed.
    TBArray(casacore::String type, casacore::String array);

    ~TBArray();

    
    // Returns true if the array is valid, false otherwise.
    bool isValid();
    
    // Returns the dimensions of the array in list format.
    std::vector<int> getDimensions();

    // Returns the dimensionality of the array.  For example, a 4x4 array
    // would return 2 while a 4x4x4 array would return 3.
    unsigned int dim();

    // If the dimensions are thought of as a list (e.g., 2x4x2), then this
    // method returns the ith dimension in the list.
    int dimensionAt(unsigned int i);

    // Returns true if the array is one-dimensional, false otherwise.
    bool isOneDimensional();

    // Returns the type of the array.
    casacore::String getType();

    // Returns the data representation.  In all but the last dimension, the
    // void*s point to vector<void*>s; in the last dimension the void*s point
    // to Strings.
    std::vector<void*>* getData();

    // Returns the row of the table where the array is located.  This is only
    // valid for data arrays.
    int getRow();

    // Returns the column of the table where the array is located.  This is
    // only valid for data arrays.
    int getCol();
    
    
    // Returns true if the given array is in the same location as this array,
    // false otherwise.
    // For data arrays: true if the row and col values are equal;
    // for table keyword arrays: true if they refer to the same keyword index;
    // for field keyword arrays: true if the fields are the same and the
    // keyword indices are the same.
    bool sameLocationAs(TBArray* array);

    // Returns the name of this array, assuming that it belongs to the given
    // table. For data arrays: "[table name][[row],[col]]";
    // for table keyword arrays: "[table name] [keyword name]";
    // for field keyword arrays: "[table name] [field name, keyword name]".
    casacore::String getName(TBTable* table);

    // Returns the data at the given coordinates, or blank if the coordinates
    // are invalid.
    casacore::String dataAt(std::vector<int> d);

    // Sets the data at the given coordinates to the given value.  This call
    // does NOT write through to the underlying table; it only updates the data
    // representation.
    void setDataAt(std::vector<int> d, casacore::String newVal);

    // Returns true if the given coordinates are valid for this array, false
    // otherwise.
    bool dimensionIsValid(std::vector<int> d);

    // Returns a "flattened" casacore::String representation of this array.  Each cell
    // is appended to the casacore::String separated by a space.
    casacore::String toFlattenedString();

    // Returns true if this array contains the given value, false otherwise.
    bool contains(casacore::String value);

    // Returns true if this array contains any value that is between the two
    // given values, false otherwise.
    bool containsBetween(casacore::String value1, casacore::String value2);

    // Returns true if this array contains any value that is less than the
    // given value, false otherwise.
    bool containsLessThan(casacore::String value);
    
    // Returns true if this array contains any value that is greater than the
    // given value, false otherwise.
    bool containsGreaterThan(casacore::String value);

private:
    // Holds the dimensions of this array.
    std::vector<int> dimensions;

    // casacore::Data representation.
    std::vector<void*> data;

    // Indicates whether the array is valid or not.
    bool valid;

    // The type of the array.
    casacore::String type;

    // The row of the array for data arrays, or the keyword index for other
    // arrays.
    int row;

    // The column of the array for data arrays, invalid for other arrays.
    int col;

    // Indicates whether this array is one-dimensional or not.
    bool oneDim;

    // Indicates whether this is a data array or not.
    bool isData;

    // Indicates whether this is a field keyword array or not.
    bool isColKeyword;

    // Holds the field name for a field keyword array, empty otherwise.
    casacore::String field;

    
    // Parses the given casacore::String into the array.
    void parseArray(casacore::String* table);

    // Helper for parseArray().  Parses a single row into the given vector.
    casacore::String parseRow(casacore::String& str, std::vector<void*>* r, std::vector<int> d, int x);

    // Helper for parseArray().  Parses a table with dimension > 1.
    void parseMultidimensionalTable(casacore::String str);

    // Helper for parseArray().  Creates placeholder objects (such as empty
    // Strings and vectors) into the given row.
    void insertPlaceholders(std::vector<void*>* r, std::vector<int> d, int x);

    // Helper method for toFlattenedString();
    casacore::String toFlattenedString(std::vector<void*>* row, int d);

    // Helper method for contains().
    bool contains(std::vector<void*>* data, int n, casacore::String v);

    // Helper method for containsBetween().
    bool containsBetween(std::vector<void*>* data, int n, double v1, double v2);

    // Helper method for containsLessThan().
    bool containsLessThan(std::vector<void*>* data, int n, double v);

    // Helper method for containsGreaterThan().
    bool containsGreaterThan(std::vector<void*>* data, int n, double v);
    
    // Deletes the data in the given row.
    void deleteData(std::vector<void*>* data, int n);
};

}

#endif /* TBARRAY_H_ */
