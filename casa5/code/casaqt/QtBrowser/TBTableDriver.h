//# TBTableDriver.h: Driver for interacting with the table on disk.
//# Copyright (C) 2007-2008
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
#ifndef TBTABLEDRIVER_H_
#define TBTABLEDRIVER_H_

#include <casaqt/QtBrowser/TBConstants.h>
#include <casaqt/QtBrowser/TBTable.h>

#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>

#include <vector>

namespace casacore{

class TableRecord;
class Table;
}

namespace casa {

//# Forward Declarations
class TableParams;
class TBField;
class TBKeyword;
class TBArray;
class TBParser;
class DriverParams;
class ProgressHelper;
class TBData;
class TBArrayData;

// <summary>
// Driver for interacting with the table on disk.
// </summary>
//
// <synopsis>
// TBTableDriver is an abstract superclass that defines how any implementing
// subclass must behave in order to be used by a TBTable.  Any implementing
// subclass should update the table parameters that the driver has references
// to via the TableParams passed into the constructor.
// </synopsis>

class TBTableDriver {
public:
    // Constructor that takes the table and the table parameters.
    TBTableDriver(TableParams* tp, TBTable* table);

    virtual ~TBTableDriver();

    
    // Sets whether the driver should print debug information or not.
    void setPrintDebug(bool pdb);
    
    
    // canRead() must be implemented by any subclass.
    // Returns whether or not the underlying table can be read (i.e., the read
    // lock is set).  Does not necessary guarantee that the next read operation
    // will succeed, but a sufficiently close check before a read should
    // suffice in most situations.
    virtual bool canRead() = 0;
    
    // canWrite() must be implemented by any subclass.
    // Returns whether or not the underlying table can be written (i.e., the
    // write lock is set).  Does not necessary guarantee that the next write
    // operation will succeed, but a sufficiently close check before a write
    // should suffice in most situations.
    virtual bool canWrite() = 0;
    
    // tryWriteLock() must be implemented by any subclass.
    // Tries to reopen the table using a write lock, and returns whether the
    // operation succeeded or not.  This MUST be done before any editing
    // operations, and the lock should be released afterwards using
    // releaseWriteLock().
    virtual bool tryWriteLock() = 0;
    
    // releaseWriteLock() must be implemented by any subclass.
    // Releases the write lock if needed.
    virtual bool releaseWriteLock() = 0;
    
    // loadRows() must be implemented by any subclass.
    // Loads the given rows into the table.  See TBTable::loadRows();
    virtual Result loadRows(int start = 0,
                            int num = TBConstants::DEFAULT_SELECT_NUM,
                            bool full = false, std::vector<casacore::String>* fields = NULL,
                            bool parsedata = true,
                            ProgressHelper* pp = NULL) = 0;
    
    // loadArray() must be implemented by any subclass.
    // Loads the data into the given array at the given coordinates.
    virtual void loadArray(TBArrayData* d, unsigned int row,
                           unsigned int col) = 0;

    // dimensionsOf() must be implemented by any subclass.
    // Returns the array dimensions of the given field.  See
    // TBTable::dimensionsOf().
    virtual std::vector<int> dimensionsOf(unsigned int col) = 0;
    
    // editData() must be implemented by any subclass.
    // Updates the cell at the given coordinates to have the given value.
    // If the data is an array, the array coordinates are provided in d.
    // Returns the result of the operation.
    virtual Result editData(unsigned int row, unsigned int col, TBData* newVal,
                            std::vector<int>* d = NULL) = 0;
                                
    // totalRowsOf() must be implemented by any subclass.
    // Returns the total rows of the table at the given location.  See
    // TBTable::totalRowsOf().
    virtual int totalRowsOf(casacore::String location) = 0;

    // insertRows() must be implemented by any subclass.
    // Inserts the given number of rows at the end of the table.  See
    // TBTable::insertRows().
    virtual Result insertRows(int n) = 0;

    // deleteRows() must be implemented by any subclass.
    // Deletes the given rows.  See TBTable::deleteRows().
    virtual Result deleteRows(std::vector<int> r) = 0;

protected:
    // The location of the table.
    casacore::String location;
    
    // Reference to the table's insertRow parameter.
    bool& insertRow;
    
    // Reference to the table's removeRow parameter.
    bool& removeRow;
    
    // Reference to the table's data parameter.    
    std::vector<std::vector<TBData*>*>& data;
    
    // Reference to the table's fields parameter.
    std::vector<TBField*>& fields;
    
    // Reference to the table's keywords parameter.
    std::vector<TBKeyword*>& keywords;
    
    // Reference to the table's subtableRows parameter.
    std::vector<int>& subtableRows;
    
    // Reference to the table's totalRows parameter.
    int& totalRows;
    
    // Reference to the table's loadedRows parameter.
    int& loadedRows;
    
    // Reference to the table's writable parameter.
    std::vector<bool>& writable;
    
    // Reference to the table's taql parameter.
    bool& taql;
    
    // Copy of the table's driver parameters.
    DriverParams* dp;
    
    // Indicates whether the driver should print debug information or not.
    bool printdebug;
    
    // Reference to the table.
    TBTable* table;

    
    // Sends a direct query through the XMLDriver and returns the Result.
    // Result query(casacore::String type, casacore::String query);
};

// <summary>
// TBTableDriver implementation that directly accesses the table on disk.
// </summary>
//
// <synopsis>
// TBTableDriverDirect is a table driver that accesses the table on disk via
// the CASA tables code module.  It is therefore faster than the XML driver
// and is the default for the table browser.  See casa::Table.
// </synopsis>

class TBTableDriverDirect : public TBTableDriver {
public:
    // Constructor that takes a table and its parameters.
    TBTableDriverDirect(TableParams* tp, TBTable* table);

    virtual ~TBTableDriverDirect();

    
    // Implements TBTableDriver::canRead().
    bool canRead();
    
    // Implements TBTableDriver::canWrite().
    bool canWrite();
    
    // Implements TBTableDriver::tryWriteLock().
    bool tryWriteLock();

    // Implements TBTableDriver::releaseWriteLock().
    bool releaseWriteLock();
    
    // Implements TBTableDriver::loadRows().
    Result loadRows(int start, int num, bool full, std::vector<casacore::String>* fields,
                    bool parsedata, ProgressHelper* progressPanel);
    
    // Implements TBTableDriver::loadArray().
    void loadArray(TBArrayData* d, unsigned int row, unsigned int col);

    // Implements TBTableDriver::dimensionsOf().
    std::vector<int> dimensionsOf(unsigned int col);
    
    // Implements TBTableDriver::editData().
    Result editData(unsigned int row, unsigned int col, TBData* newVal,
                    std::vector<int>* d = NULL);
                                
    // Implements TBTableDriver::totalRowsOf().
    int totalRowsOf(casacore::String location);

    // Implements TBTableDriver::insertRows().
    Result insertRows(int n);

    // Implements TBTableDriver::deleteRows().
    Result deleteRows(std::vector<int> r);

    // Converts keywords in a casacore::TableRecord to a vector of TBKeywords.
    static std::vector<TBKeyword*>* getKeywords(casacore::RecordInterface& kws);
    
private:
    // Reference to table on disk.
    casacore::Table m_table;
};

// NOTE: the TBTableDriverXML has been disabled.  If it is to be used in the
// future, the problems with the new TBData infrastructure must be addressed.

/*
// <summary>
// TBTableDriver implementation that accesses tables through XML.
// </summary>
//
// <synopsis>
// TBTableDriverXML is a table driver that collects table information into
// a casacore::String of information in XML format and then parses the XML.  The XML
// is generated by the TBXMLDriver (which is what was used with the old
// Java table browser) and then parsed by a subclass of TBParser that is
// provided in the TableParams.  This approach has the advantage of
// abstraction and portability; for example, in the future a server/client
// architecture could be implemented to browse remote tables.  However, it
// also as the disadvantage of being slower than the direct approach,
// especially when dealing with large amounts of data such as plotting.
// See TBParser.
// </synopsis>

class TBTableDriverXML : public TBTableDriver {
public:
    // Constructor that takes table parameters.
    TBTableDriverXML(TableParams* tp, TBTable* table);

    ~TBTableDriverXML();

    
    // Implements TBTableDriver::loadRows().
    Result loadRows(int start, int num, bool full, std::vector<casacore::String>* fields,
                    bool parsedata, ProgressHelper* progressPanel);
    
    void loadArray(TBArrayData* d, unsigned int row, unsigned int col);

    // Implements TBTableDriver::dimensionsOf().
    std::vector<int> dimensionsOf(unsigned int col);
    
    Result editData(unsigned int row, unsigned int col, TBData* newVal,
                        std::vector<int>* d = NULL);
                                
    // Implements TBTableDriver::totalRowsOf().
    int totalRowsOf(casacore::String location);

    // Overrides TBTableDriver::setPrintDebug().
    void setPrintDebug(bool pdb);

    // Implements TBTableDriver::insertRows().
    Result insertRows(int n);

    // Implements TBTableDriver::deleteRows().
    Result deleteRows(std::vector<int> r);
    
private:
    // XML parser.
    TBParser* parser;
};
*/

}

#endif /* TBTABLEDRIVER_H_ */
