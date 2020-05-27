//# TBConstants.h: Constants, defaults, and common functions for the browser.
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
#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <sstream>
#include <vector>
#include <map>

#include <QFrame>

#include <casa/BasicSL/String.h>
#include <casa/Utilities/DataType.h>
#include <casa/Containers/Record.h>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>

/*
#ifdef AIPS_HAS_QWT
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#endif
*/
#include <graphics/GenericPlotter/PlotFactory.h>

using namespace xercesc;

namespace casa {

//# Forward Declarations
class TBTableTabs;

 namespace tb {
    // Enum listing the possible types of comparisons that can be made in, for
    // example, a filter rule.
    enum Comparator {
        EQUALS, CONTAINS, BETWEEN, CONTAINSBT, LESSTHAN, CONTAINSLT, GREATERTHAN,
        CONTAINSGT
    };

    // Enum listing the format for boolean values: "true/false", "t/f", "1/0", etc.
    enum BooleanFormat {
        DEFAULT, TRUEFALSE, TF, B10
    };

    // Enum to indicate the different driver types.
    enum Driver {
        DIRECT//, XML
    };

    // Enum to indicate the different parsers available for XML drivers.
    enum Parser {
        HOME, XERCES_DOM, XERCES_SAX
    };
 }


// <summary>
// Convenience class for a casacore::String/bool tuple.
// <summary>
//
// <synopsis>
// Result is nothing more than a casacore::String/bool pair.  Semantically, it can be
// used to return the result of an operation: the casacore::String is the value and the
// bool indicates whether that value is valid or the operation was successful.
// </synopsis>

class Result {
public:
    // Construct a Result with the given values.
    Result(casacore::String r= "", bool v = true): result(r), valid(v) { }

    ~Result() { }

    // Result string.
    casacore::String result;
    
    // Result bool.
    bool valid;
};


// <summary>
// Results of a row locate on at least one table.
// </summary>
//
// <synopsis>
// A TBLocatedRows basically consists of a map from TBTableTabs to vectors of
// ints that represent row numbers.
// </synopsis>

class TBLocatedRows {
public:
    // Default Constructor.
    TBLocatedRows();
    
    ~TBLocatedRows();
    
    // Returns a list of all the tables in this TBLocatedRows.
    std::vector<TBTableTabs*> tables();
    
    // Associates the given list of row numbers with the given TBTableTabs.
    void put(TBTableTabs* tt, std::vector<int>* r);
    
    
    // Results.
    std::map<TBTableTabs*, std::vector<int>*> results;
};


// <summary>
// Constants, defaults, and commonly-used functions for the table browser.
// <summary>
//
// <synopsis>
// TBConstants contains definitions that may be common or useful to multiple
// classes.
// </synopsis>

class TBConstants {
public:
    /* Debugging Constants/Defaults */
    
    // The current debug level.
    static int debugThreshold;
    
    // Debug levels.
    // <group>
    static const int DEBUG_OFF;
    static const int DEBUG_LOW;
    static const int DEBUG_MED;
    static const int DEBUG_HIGH;
    // </group>

    // Print the given message to standard out with the given indentation level,
    // IF the message's level is at or below the current debugThreshold.
    static void dprint(int level, casacore::String message, int indentLevel = 0);

    
    /* Browsing Constants/Defaults */

    // The default number of rows to load from the table at a time.
    static const unsigned int DEFAULT_SELECT_NUM;

    // The default number of rows to load at a time while exporting.
    static const unsigned int DEFAULT_EXPORT_NUM;

    // The maximum number of rows to load from the table at a time.
    static const unsigned int MAX_SELECT_NUM;

    // The maximum number of actions to keep in action lists.
    static const unsigned int MAX_ACTION_BUFFER;

    // The default row interval for plotting.
    static const unsigned int DEFAULT_ROW_INTERVAL;

    // Kinds of queries that can be sent to the table through the TBXMLDriver.
    // <group>
    static const casacore::String QUERY_QUERY;
    static const casacore::String QUERY_ARRAY;
    static const casacore::String QUERY_UPDATE;
    static const casacore::String QUERY_FULL;
    // </group>

    // Returns the display value of the given Comparator.
    static casacore::String compToString(tb::Comparator c);

    // Returns the Comparator corresponding to the given display value.
    static tb::Comparator stringToComp(casacore::String str);

    // Returns the display value of the given BooleanFormat.
    static casacore::String bFormToString(tb::BooleanFormat bf);

    // Returns the BooleanFormat corresponding to the given display value.
    static tb::BooleanFormat stringToBForm(casacore::String str);

    // Constants used for an array slicer.  See TBSlicer.
    // <group>
    static const int SLICER_ROW_AXIS;
    static const int SLICER_COL_AXIS;
    // </group>
    
    // Number of rows in an embedded record widget to show by default.
    static const unsigned int DEFAULT_RECORD_VISIBLE_ROWS;

    // Returns the path of the CASA top-level directory.  This value is
    // retrieved from the environment variables the first time it is called,
    // and then returns the stored result on subsequent calls.
    static casacore::String aipspath();

    // Returns the absolute path of the ~/.casa directory.  This value is
    // retrieved from the environment variables the first time it is called,
    // and then returns the stored result on subsequent calls.
    static casacore::String dotCasapyDir();
    
    // Converts a vector<int> to an IPosition.
    static casacore::IPosition ipos(std::vector<int>& d);
    
    // Converts an casacore::IPosition to a vector<int>.
    static std::vector<int> ipos(casacore::IPosition& d);
    
    // Increments the given dimension, using the given shape as a maximum.
    static bool increment(casacore::IPosition& shape, casacore::IPosition& d);
    
    // Increments the given dimension, using the given shape as a maximum.
    static bool increment(std::vector<int>& shape, std::vector<int>& d);
    
    // Inserts the given widget into the given placeholder frame.
    static void insert(QFrame* frame, QWidget* widget);
    
    // Inserts the given layout into the given placeholder frame.
    static void insert(QFrame* frame, QLayout* layout);
    
    
    /* casacore::Data Types and Related Methods */

    // casacore::Table data types.
    // <group>
    static const casacore::String TYPE_STRING;
    static const casacore::String TYPE_DOUBLE;
    static const casacore::String TYPE_FLOAT;
    static const casacore::String TYPE_INT;
    static const casacore::String TYPE_UINT;
    static const casacore::String TYPE_BOOL;
    static const casacore::String TYPE_CHAR;
    static const casacore::String TYPE_UCHAR;
    static const casacore::String TYPE_SHORT;
    static const casacore::String TYPE_COMPLEX;
    static const casacore::String TYPE_DCOMPLEX;
    static const casacore::String TYPE_TABLE;
    static const casacore::String TYPE_RECORD;
    static const casacore::String TYPE_DATE;

    static const casacore::String TYPE_ARRAY_STRING;
    static const casacore::String TYPE_ARRAY_DOUBLE;
    static const casacore::String TYPE_ARRAY_FLOAT;
    static const casacore::String TYPE_ARRAY_INT;
    static const casacore::String TYPE_ARRAY_UINT;
    static const casacore::String TYPE_ARRAY_BOOL;
    static const casacore::String TYPE_ARRAY_CHAR;
    static const casacore::String TYPE_ARRAY_UCHAR;
    static const casacore::String TYPE_ARRAY_SHORT;
    static const casacore::String TYPE_ARRAY_COMPLEX;
    static const casacore::String TYPE_ARRAY_DCOMPLEX;
    // </group>

    // Unsupported types
    // static const casacore::String TYPE_USHORT;
    // static const casacore::String TYPE_ARRAY_USHORT;
    // static const casacore::String TYPE_OTHER;
    
    // Return a vector containing all the valid data types.
    static std::vector<casacore::String>* allTypes();

    // Return a vector containing all the valid array data types.
    static std::vector<casacore::String>* arrayTypes();

    // Return a vector containing all the valid non-array data types.
    static std::vector<casacore::String>* nonArrayTypes();

    // The comments for casacore::Double fields that indicate that they should be
    // interpreted as a date.
    static const casacore::String COMMENT_DATE;
    static const casacore::String COMMENT_TIMP;
    static const casacore::String COMMENT_TIMP2;

    // Convert a CASA DataType enum to its casacore::String equivalent.
    static casacore::String typeName(casacore::DataType dt);

    // Returns a human-readable name for the given type.
    static casacore::String typeName(casacore::String type);

    // Converts the given type into its corresponding VOTable type.  See
    // http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC11
    static casacore::String VOType(casacore::String& type);

    // Returns true if the given type is a table type, false otherwise.
    static bool typeIsTable(casacore::String& type);

    // Returns true if the given type is an array type, false otherwise.
    static bool typeIsArray(casacore::String& type);

    // Returns true if the given type is plottable, false otherwise.
    static bool typeIsPlottable(casacore::String& type);

    // Returns true if the given type is numberable, false otherwise.
    static bool typeIsNumberable(casacore::String& type);

    // Returns true if the given type can be used as an index (i.e.,
    // integer-like), false otherwise.
    static bool typeIsIndexable(casacore::String& type);

    // Retruns true if the given type is a complex or array of complex types,
    // false otherwise.
    static bool typeIsComplex(casacore::String& type);

    // Returns the value type from the given array type.  For example,
    // arrayType(TYPE_ARRAY_STRING) = TYPE_STRING.
    static casacore::String arrayType(casacore::String& at);

    // Returns true if the given value is valid for the given type, false
    // otherwise.
    static bool valueIsValid(casacore::String& value, casacore::String& type);

    // Formats the given value for the given data type and returns the
    // formatted value.
    static casacore::String formatValue(casacore::String& value, casacore::String& type);

    // Converts the given value of the given type to a double and return it.
    // Returns -1 for invalid doubles.
    static double valueToDouble(casacore::String& value, casacore::String& type);
    
    // Converts the given casacore::String in date format (year-month-day-hour:min:sec)
    // to its double representation.  The double returned is in modified
    // julian seconds.
    static double date(casacore::String d);

    // Converts the given double into a casacore::String in date format
    // (year-month-day-hour:min:sec).  The input double should be in modified
    // julian seconds.
    static casacore::String date(const double d, const int numDecimals = -1);
    
    // Holds the default number of decimals displayed in a number.
    static const int DEFAULT_DECIMALS;
    
    // Holds the default number of decimals displayed in a date.
    static const int DEFAULT_DATE_DECIMALS;
    
    // Returns true if the given casacore::String contains a valid date format, false
    // otherwise.  A valid date format has %y, %m, %d, %n, and %s appear
    // exactly once.  %y = year, %m = month, %d = day of month, %n = minute
    // %s = second; all fields are integers.
    static bool dateFormatIsValid(casacore::String& d);

    // Converts a casacore::String in complex format (a,b) into a pair of doubles.
    static std::pair<double, double> toComplex(casacore::String str);
    
    // casacore::String used in parsing arrays.
    static const casacore::String ARRAY_AXES_LENGTHS;


    /* XML Driver Constants/Defaults */

    // XML token names that are used in table parsing.
    // <group>
    static const casacore::String XML_VOTABLE;
    static const casacore::String XML_RESOURCE;
    static const casacore::String XML_TABLE;
    static const casacore::String XML_TOTAL;
    static const casacore::String XML_FIELD;
    static const casacore::String XML_KEYWORD;
    static const casacore::String XML_COLUMNKW;
    static const casacore::String XML_RWINFO;
    static const casacore::String XML_DATA;
    static const casacore::String XML_TABLEDATA;
    static const casacore::String XML_TR;
    static const casacore::String XML_TD;
    static const casacore::String XML_INSERTROW;
    static const casacore::String XML_REMOVEROW;
    static const casacore::String XML_TRUE;
    static const casacore::String XML_FALSE;
    static const casacore::String XML_ROW;
    static const casacore::String XML_NAME;

    static const casacore::String XML_ID;
    static const casacore::String XML_KEYWORD_NAME;
    static const casacore::String XML_KEYWORD_TYPE;
    static const casacore::String XML_KEYWORD_VAL;
    static const casacore::String XML_COLKW_COL;
    static const casacore::String XML_FIELD_NAME;
    static const casacore::String XML_FIELD_TYPE;
    static const casacore::String XML_FIELD_UNIT;
    static const casacore::String XML_FIELD_UCD;
    static const casacore::String XML_FIELD_REF;
    static const casacore::String XML_FIELD_PRECISION;
    static const casacore::String XML_FIELD_WIDTH;
    // </group>

    // Error text that is thrown/returned on an empty table or row.
    static const casacore::String ERROR_EMPTY;

    // Converts a char* to an XMLCh*.  See XMLString::transcode().
    static const XMLCh* xstr(char* str);

    // Converts a casacore::String to an XMLCh*.  See XMLString::transcode().
    static const XMLCh* xstr(casacore::String str);

    // Converts an XMLCh* to a String.  See XMLString::transcode().
    static casacore::String xstr(XMLCh* str);

    // Converts a const XMLCh* to a String.  See XMLString::transcode().
    static casacore::String xstr(const XMLCh* str);


    /* casacore::String and Number Methods */

    // Converts the given int into its casacore::String representation and returns it.
    static casacore::String itoa(int n);
    
    // Converts the given unsigned int into its casacore::String representation and
    // returns it.
    static casacore::String uitoa(unsigned int n);
    
    // Converts the given short int into its casacore::String representation and
    // returns it.
    static casacore::String sitoa(short int n);

    // Converts the given float into its casacore::String representation and returns it.
    static casacore::String ftoa(float f, int decimals = -1, bool scientific = false);

    // Converts the given double into its casacore::String representation and returns it.
    static casacore::String dtoa(double d, int decimals = -1, bool scientific = false);

    // Converts the given casacore::String into its int representation.  Returns the
    // result of sscanf.
    static int atoi(casacore::String& str, int* i);
    
    // Converts the given casacore::String into its unsigned int representation.  Returns
    // the result of sscanf.
    static int atoui(casacore::String& str, unsigned int* i);
    
    // Converts the given casacore::String into its short int representation.  Returns
    // the result of sscanf.
    static int atosi(casacore::String& str, short int* i);

    // Converts the given casacore::String into its float representation.  Returns the
    // result of sscanf.
    static int atof(casacore::String& str, float* f);

    // Converts the given casacore::String into its double representation.  Returns the
    // result of sscanf.
    static int atod(casacore::String& str, double* d);
    
    // Rounds the given double to the nearest int and returns it.
    static int round(double d);

    // Returns true if the given char is a whitespace, false otherwise.
    static bool isWS(char c);

    // Trims the given casacore::String of its leading and trailing whitespace using
    // casacore::String::erase().
    static void strtrim(casacore::String& str);

    // Converts all upper-case letters to lower-case letters in the given
    // String.
    static void toLower(casacore::String& str);

    // Returns the index of the next instance of a whitespace character in the
    // given casacore::String from the given index, or the string's length if there is none.
    static unsigned int findWS(casacore::String& str, int index = 0);

    // Returns the index of the next instance of the find casacore::String in the given
    // casacore::String starting from the given index, or str's length if there is none.
    static unsigned int findIgnoreCase(casacore::String& str, casacore::String& find,
                                       int index = 0);

    // Returns true if the two Strings are equal (ignoring case), false
    // otherwise.
    static bool equalsIgnoreCase(casacore::String str1, casacore::String str2);

    // Returns the name of the file from the given path (the last part of the
    // path).
    static casacore::String nameFromPath(casacore::String& path);

    // Returns the directory from the given path (all but the last part of the
    // path).
    static casacore::String dirFromPath(casacore::String& path);


    /* TBConnection Constants/Defaults */

    // Help descriptions used by TBConnection.
    // <group>
    static const casacore::String OPEN_TEXT_LOCAL;
    static const casacore::String OPEN_TEXT_REMOTE;
    static const casacore::String OPEN_TEXT_HOST;
    static const casacore::String OPEN_TEXT_PORT;
    static const casacore::String OPEN_TEXT_LOCATION;
    static const casacore::String OPEN_TEXT_DIRECT;
    static const casacore::String OPEN_TEXT_XML;
    static const casacore::String OPEN_TEXT_HOME;
    static const casacore::String OPEN_TEXT_DOM;
    static const casacore::String OPEN_TEXT_SAX;
    static const casacore::String OPEN_TEXT_START;
    static const casacore::String OPEN_TEXT_NUM;
    // </group>


    /* Help Constants/Defaults */

    // Constants for the help system.
    // <group>
    static const casacore::String HELP_DIR;
    static const casacore::String HELP_INDEX;
    static const casacore::String HELP_XML;
    static const casacore::String HELP_XML_HELP;
    static const casacore::String HELP_XML_CATEGORY;
    static const casacore::String HELP_XML_NAME;
    static const casacore::String HELP_XML_GROUP;
    static const casacore::String HELP_XML_ITEM;
    // </group>


    /* Plotting Constants/Defaults */

    // Plotting defaults.
    // <group>
    
    // The default plotter implementation for the browser.
    static const Plotter::Implementation defaultPlotterImplementation;
    
    // Returns a default line to use with new plots, using the given factory
    // (should match the implementation of defaultPlotterImplementation).
    static PlotLinePtr defaultPlotLine(PlotFactoryPtr factory);
    
    // Returns a default plot symbol to use with new plots, using the given
    // factory (should match the implementation of
    // defaultPlotterImplementation).
    static PlotSymbolPtr defaultPlotSymbol(PlotFactoryPtr factory);


    /* View Constants/Defaults */

    // View constants.
    // <group>
    static const casacore::String VIEW_SAVE_LOC;

    static const casacore::String VIEW_DOCUMENT;
    static const casacore::String VIEW_VIEW;
    static const casacore::String VIEW_LASTDIR;
    static const casacore::String VIEW_HISTLIMIT;
    static const casacore::String VIEW_TABLE;
    static const casacore::String VIEW_LOCATION;
    static const casacore::String VIEW_SELECTED;
    static const casacore::String VIEW_TAQL;
    static const casacore::String VIEW_HIDDEN;
    static const casacore::String VIEW_HIDDEN_LENGTH;
    static const casacore::String VIEW_VISIND;
    static const casacore::String VIEW_ROWS;
    static const casacore::String VIEW_ROWS_FROM;
    static const casacore::String VIEW_ROWS_NUM;
    static const casacore::String VIEW_FILTER;
    static const casacore::String VIEW_FILTER_RULE;
    static const casacore::String VIEW_FILTER_RULE_FIELD;
    static const casacore::String VIEW_FILTER_RULE_COMPARATOR;
    static const casacore::String VIEW_FILTER_RULE_VALUE;
    static const casacore::String VIEW_FILTER_RULE_VALUE2;
    static const casacore::String VIEW_FILTER_RULE_NOT;
    static const casacore::String VIEW_FILTER_RULE_ANY;
    static const casacore::String VIEW_FILTER_RULE_TYPE;
    static const casacore::String VIEW_FORMATS;
    static const casacore::String VIEW_FORMAT;
    static const casacore::String VIEW_FORMAT_COL;
    static const casacore::String VIEW_FORMAT_DECIMALS;
    static const casacore::String VIEW_FORMAT_SFORMAT;
    static const casacore::String VIEW_FORMAT_BFORMAT;
    static const casacore::String VIEW_FORMAT_DFORMAT;
    static const casacore::String VIEW_FORMAT_VTHRESHOLD;
    static const casacore::String VIEW_FORMAT_ALLFONT;
    static const casacore::String VIEW_FORMAT_FONT;
    static const casacore::String VIEW_FORMAT_COLOR;
    static const casacore::String VIEW_FORMAT_FAMILY;
    static const casacore::String VIEW_FORMAT_SIZE;
    static const casacore::String VIEW_FORMAT_BOLD;
    static const casacore::String VIEW_FORMAT_ITALICS;
    static const casacore::String VIEW_FORMAT_ULINE;
    static const casacore::String VIEW_FORMAT_STRIKE;
    static const casacore::String VIEW_SORT;
    static const casacore::String VIEW_NAME;
    static const casacore::String VIEW_SORT_ASCENDING;
    // </group>
    
private:
    // casacore::Path to the CASA top-level directory.
    static casacore::String AIPS_PATH;

    // Common string in the TBConnection defaults.
    static const casacore::String OPEN_PAGE;

    // Holds the absolute location of the ~/.casa directory.
    static casacore::String DOT_CASAPY_DIR;
};

}

#endif /*CONSTANTS_H_*/
