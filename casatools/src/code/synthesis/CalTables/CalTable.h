//# CalTable.h: Calibration table access and creation
//# Copyright (C) 1996,1997,1998,1999,2001,2002,2003
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
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#ifndef CALIBRATION_CALTABLE_H
#define CALIBRATION_CALTABLE_H

#include <casacore/casa/aips.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Utilities/Sort.h>
#include <synthesis/CalTables/CalTableDesc.h>
#include <synthesis/CalTables/CalMainRecord.h>
#include <synthesis/CalTables/CalDescRecord.h>
#include <synthesis/CalTables/CalHistRecord.h>

namespace casa { //# NAMESPACE CASA - BEGIN

// <summary> 
// CalTable: Calibration table access and creation
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="" tests="" demos="">

// <prerequisite>
//   <li> <linkto class="CalTableDesc">CalTableDesc</linkto> module
// </prerequisite>
//
// <etymology>
// From "calibration" and "table".
// </etymology>
//
// <synopsis>
// The CalTable classes provide basic access to calibration tables. This 
// includes the creation of calibration tables, sorting and selection
// capabilities, and basic data access. Specializations for baseline-based,
// time-variable and solvable VisJones types, and sub-types, are provided 
// through inheritance.
// </etymology>
//
// <example>
// <srcblock>
// </srcblock>
// </example>
//
// <motivation>
// This class is used by other calibration table accessors and iterators.
// </motivation>
//
// <todo asof="98/01/01">
// (i) Generic synthesis selection
// </todo>

class CalTable
{
 public:
   // Default null constructor, and destructor
   CalTable();
   virtual ~CalTable();

   // Construct from a specified table name, calibration table descriptor 
   // and table access option. Used for creating new tables.
   CalTable (const casacore::String& tableName, CalTableDesc& ctableDesc,
	     casacore::Table::TableOption access = casacore::Table::New);

   // Construct from a specified table name, and access option. Used
   // for accessing existing tables.
   CalTable (const casacore::String& tableName, casacore::Table::TableOption access = casacore::Table::Old);

   // Construct from an existing table object
   CalTable (const casacore::Table& table);

   // Copy constructor
   CalTable (const CalTable& other);

   // Assignment operator
   CalTable& operator= (const CalTable& other);

   // casacore::Sort (cal_main)
   CalTable sort (const casacore::Block <casacore::String>& columnNames, 
		  casacore::Sort::Order order = casacore::Sort::Ascending,
		  casacore::Sort::Option option = casacore::Sort::HeapSort); 

   // casacore::Sort (cal_main) IN PLACE
   void sort2 (const casacore::Block <casacore::String>& columnNames, 
	      casacore::Sort::Order order = casacore::Sort::Ascending,
	      casacore::Sort::Option option = casacore::Sort::HeapSort); 

   // Apply selection to the calibration table 
   CalTable select (const casacore::String& calSelect);

   // Apply selection to the calibration table  IN PLACE
   void select2 (const casacore::String& calSelect);

   // Return number of rows in cal_main, cal_desc or cal_history
   casacore::Int nRowMain() const;
   casacore::Int nRowDesc() const;
   casacore::Int nRowHistory() const;
   casacore::Int nRowObservation() const;

   // Add rows to cal_main, cal_desc or cal_history
   void addRowMain (casacore::uInt nrrow = 1, casacore::Bool initialize = false) 
     {itsMainTable->addRow(nrrow, initialize);};
   void addRowDesc (casacore::uInt nrrow = 1, casacore::Bool initialize = false)
     {itsDescTable->addRow(nrrow, initialize);};
   void addRowHistory (casacore::uInt nrrow = 1, casacore::Bool initialize = false)
     {itsHistoryTable->addRow(nrrow, initialize);};

   // Get a row from cal_main, cal_desc or cal_history
   casacore::Record getRowMain (const casacore::Int& jrow);
   casacore::Record getRowDesc (const casacore::Int& jrow);
   casacore::Record getRowHistory (const casacore::Int& jrow);
   casacore::Record getRowObservation (const casacore::Int& jrow);

   // Put a row to cal_main, cal_desc or cal_history
   void putRowMain (const casacore::Int& jrow, CalMainRecord& tableRec);
   void putRowDesc (const casacore::Int& jrow, CalDescRecord& tableRec);
   void putRowHistory (const casacore::Int& jrow, CalHistoryRecord& tableRec);

   // Return the maximum antenna number in the table
   virtual casacore::Int maxAntenna();

   // Return the number of unique time slots in the table
   casacore::Int numberTimeSlots (const casacore::Double& fracError);

   // Return the number of rows per CalDescId
   void rowsPerCalDescId(casacore::Vector<casacore::Int>& rowspercdi);

   const casacore::String tableName() { return calMainAsTable().tableName(); }

   const casacore::String type() { return calMainAsTable().tableInfo().subType(); }

 protected:
   // Create a new table
   void createCalTable (const casacore::String& tableName, CalTableDesc& ctableDesc,
			casacore::Table::TableOption access = casacore::Table::New);

   // Open an existing table
   void openCalTable (const casacore::String& tableName, 
		      casacore::Table::TableOption access = casacore::Table::Old);

   // Return cal_main and sub-tables as casacore::Table references
   casacore::Table& calMainAsTable() {return *itsMainTable;};
   casacore::Table& calDescAsTable() {return *itsDescTable;};
   casacore::Table& calHistoryAsTable() {return *itsHistoryTable;};
   casacore::Table& calObservationAsTable() {return *itsObservationTable;};
   const casacore::Table& calMainAsTable() const {return *itsMainTable;};
   const casacore::Table& calDescAsTable() const {return *itsDescTable;};
   const casacore::Table& calHistoryAsTable() const {return *itsHistoryTable;};
   const casacore::Table& calObservationAsTable() const {return *itsObservationTable;};

   // Friend class access from the ROCalMainColumns, ROCalDescColumns
   // and ROCalHistoryColumns class hierarchies
   friend class ROCalMainColumns;
   friend class CalMainColumns;
   friend class ROCalDescColumns;
   friend class CalDescColumns;
   friend class ROCalHistoryColumns;
   friend class CalHistoryColumns;

 private:
   // Pointer to underlying cal_main casacore::Table object
   casacore::Table* itsMainTable;

   // Pointers to the cal_desc and cal_history sub-tables
   casacore::Table* itsDescTable;
   casacore::Table* itsHistoryTable;
   casacore::Table* itsObservationTable;

   // Check if has optional OBSERVATION Table
   bool hasObsTable() const { 
     return (itsMainTable->keywordSet().fieldNumber("OBSERVATION") != -1); }
 };


} //# NAMESPACE CASA - END

#endif
   
