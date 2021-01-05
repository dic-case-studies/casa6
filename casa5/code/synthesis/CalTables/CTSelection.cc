//# CTSelection.cc: Implementation of CTSelection class
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2003
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
//----------------------------------------------------------------------------

#include <synthesis/CalTables/CTSelection.h>
#include <synthesis/CalTables/CTColumns.h>
#include <synthesis/CalTables/CTInterface.h>
#include <casacore/ms/MSSel/MSSelectionTools.h>
#include <casacore/casa/Utilities/GenSort.h>

namespace casa { //# NAMESPACE CASA - BEGIN

  CTSelection::CTSelection() {
    msSelection_p = new casacore::MSSelection();
  };

  CTSelection::CTSelection(
			const NewCalTable& ct,
			const casacore::MSSelection::MSSMode& mode,
			const casacore::String& timeExpr,
			const casacore::String& antennaExpr,
			const casacore::String& fieldExpr,
			const casacore::String& spwExpr,
			const casacore::String& taqlExpr,
			const casacore::String& scanExpr,
			const casacore::String& stateExpr,
			const casacore::String& observationExpr) {
    msSelection_p = new casacore::MSSelection();
    setTimeExpr(timeExpr);
    setAntennaExpr(antennaExpr);
    setFieldExpr(fieldExpr);
    setSpwExpr(spwExpr);
    setTaQLExpr(taqlExpr);
    setScanExpr(scanExpr);
    setStateExpr(stateExpr);
    setObservationExpr(observationExpr);

    if (mode==casacore::MSSelection::PARSE_NOW) {
        CTInterface cti(ct);
        toTableExprNode(&cti);
    }
  }
  
  CTSelection::CTSelection (const CTSelection& other) {
    if (this != &other) {
        msSelection_p = other.msSelection_p;
    }
  }

  CTSelection& CTSelection::operator= (const CTSelection& other) {
    if (this != &other) {
        msSelection_p = other.msSelection_p;
    }
    
    return *this;
  }

  CTSelection::~CTSelection() {
    if (msSelection_p) {
        delete msSelection_p;
        msSelection_p = NULL;
    }
  }

  //---------------------------------------------------------------------------

  void CTSelection::reset(casacore::MSSelectableTable& msLike,
			  const casacore::MSSelection::MSSMode& mode,
			  const casacore::String& timeExpr,
			  const casacore::String& antennaExpr,
			  const casacore::String& fieldExpr,
			  const casacore::String& spwExpr,
			  const casacore::String& taqlExpr,
			  const casacore::String& scanExpr,
			  const casacore::String& stateExpr,
			  const casacore::String& observationExpr) {
    clear();
    setTimeExpr(timeExpr);
    setAntennaExpr(antennaExpr);
    setFieldExpr(fieldExpr);
    setSpwExpr(spwExpr);
    setTaQLExpr(taqlExpr);
    setScanExpr(scanExpr);
    setStateExpr(stateExpr);
    setObservationExpr(observationExpr);
    if (mode==casacore::MSSelection::PARSE_NOW)
        toTableExprNode(&msLike);
  };

  //---------------------------------------------------------------------------

  casacore::TableExprNode CTSelection::toTableExprNode(
        casacore::MSSelectableTable* msLike) {
    // Convert the CT selection to a TableExprNode object, 
    // representing a TaQL selection in C++.
    // Input:
    //    msLike           const MSSelectableTable&  CalTable to bind TaQL
    // Output:
    //    toTableExprNode  TableExprNode             Table expression node
    //
    // Interpret all expressions and produce a consolidated TEN.  
    //

    // Specialize antenna selection if necessary
    casacore::String antExpr =
      msSelection_p->getExpr(casacore::MSSelection::ANTENNA_EXPR);
    if (!antExpr.empty()) {
      casacore::MSSelectableTable::MSSDataType type = msLike->dataType();
      switch (type) {
        case casacore::MSSelectableTable::PURE_ANTENNA_BASED:
          doAntenna1Selection(antExpr, msLike);
          break;
        case casacore::MSSelectableTable::REF_ANTENNA_BASED:
          doRefAntennaSelection(antExpr, msLike);
          break;
        case casacore::MSSelectableTable::BASELINE_BASED:
          // Use MSSelection antenna/baseline selection
          break;
      }
    }

    return msSelection_p->toTableExprNode(msLike);
  }

  //---------------------------------------------------------------------------

  void CTSelection::doAntenna1Selection(
    const casacore::String& antennaExpr, casacore::MSSelectableTable* msLike) {
    // Use antenna/baseline selections for ANTENNA1 only
    casacore::Vector<casacore::String> selections;
    selections = split(antennaExpr, ';', selections);

    // List of antennas to select or exclude
    casacore::String ant1(""), notAnt1("");

    for (auto& selection : selections) {
      // Check and strip "!" from selection
      bool neg = isSelectionExcluded(selection);

      casacore::Vector<casacore::Int> ant1List = getAntennaList(selection, msLike);
      casacore::String ant1ListStr = casacore::MSSelection::indexExprStr(ant1List);

      if (neg) {
        notAnt1 += (notAnt1.empty() ? "" : ",") + ant1ListStr;
      } else {
        ant1 += (ant1.empty() ? "" : ",") + ant1ListStr;
      }
    }

    // Create taql strings for selection/exclusion
    casacore::String antTaqlExpr("");
    if (!ant1.empty()) {
      antTaqlExpr = "ANTENNA1 IN [" + ant1 + "]";
    }
    if (!notAnt1.empty()) {
      casacore::String sep = (antTaqlExpr.empty() ? "" : " && ");
      antTaqlExpr += sep + "ANTENNA1 NOT IN [" + notAnt1 + "]";
    }

    // Reset antenna expr
    msSelection_p->setAntennaExpr("");

    // Append taql expr
    casacore::String taqlExpr = msSelection_p->getExpr(casacore::MSSelection::TAQL_EXPR);
    casacore::String sep = (taqlExpr.empty() ? "" : " && ");
    taqlExpr += sep + antTaqlExpr;
    msSelection_p->setTaQLExpr(taqlExpr);
  }

  bool CTSelection::isSelectionExcluded(casacore::String& selection) {
    // Finds and strips exclude '!' syntax
    // Returns true if found
    bool exclude = (selection[0] == '!');
    if (exclude) {
      selection = selection.after('!');
    }
    return exclude;
  }

  casacore::Vector<casacore::Int> CTSelection::getAntennaList(
    const casacore::String& selection, casacore::MSSelectableTable* msLike,
    bool getAnt1) {
    // Use MSSelection to get Antenna1 (getAnt1=true) or Antenna2 ids
    NewCalTable nct = NewCalTable(*msLike->table());
    CTInterface cti(nct);
    casacore::MSSelection mssel;
    mssel.setAntennaExpr(selection);
    casacore::TableExprNode ten = mssel.toTableExprNode(&cti);
    casacore::Vector<casacore::Int> antlist;

    if (getAnt1) {
      antlist = mssel.getAntenna1List();
    } else {
      antlist = mssel.getAntenna2List();
    }

    return antlist;
  }

  void CTSelection::doRefAntennaSelection(const casacore::String& antennaExpr,
      casacore::MSSelectableTable* msLike) {
    // Convert MSSelection antenna selection for refant-based cal tables to TaQL.
    // Needs special handling when expr ANT1 or ANT2 is ref antenna.
    casacore::String antTaqlExpr("");

    NewCalTable nct = NewCalTable(*msLike->table());
    casacore::rownr_t nAnt = nct.antenna().nrow();
    casacore::Vector<casacore::Int> refAntlist = getRefAntIds(msLike);

    // Separate antenna selections
    casacore::Vector<casacore::String> selections;
    selections = split(antennaExpr, ';', selections);

    for (auto& selection : selections) {
      // antenna 1
      casacore::Vector<casacore::Int> ant1list = getAntennaList(selection, msLike);
      casacore::String ant1select, ant1exclude;
      setAntSelectExclude(ant1list, selection, msLike, ant1select, ant1exclude);

      // antenna 2
      casacore::Vector<casacore::Int> ant2list = getAntennaList(selection, msLike, false);
      if (ant2list.size() == nAnt) {
        ant2list.resize();
        ant2list = refAntlist;
        if (selection.startsWith("!")) {
          ant2list *= -1;
        }
      }
      casacore::String ant2select, ant2exclude;
      setAntSelectExclude(ant2list, selection, msLike, ant2select, ant2exclude);

      // TaQL for antenna1 and 2
      bool select(!ant1select.empty()), exclude(!ant1exclude.empty());
      casacore::String selectTaql, excludeTaql, antTaql;
      if (select) {
        selectTaql = "(ANTENNA1 IN [" + ant1select + "]";
        selectTaql += " && ANTENNA2 IN [" + ant2select + "])";
        antTaql = (exclude ? "(" : "") + selectTaql;
      }
      if (!ant1exclude.empty()) {
        // !(ANT1 & ANT2) == !ANT1 || !ANT2
        excludeTaql += "(ANTENNA1 NOT IN [" + ant1exclude + "]";
        excludeTaql += " || ANTENNA2 NOT IN [" + ant2exclude + "])";
        antTaql += (select ? " && " : "") + excludeTaql + (select ? ")" : "");
      }

      antTaqlExpr += (antTaqlExpr.empty() ? "" : " || ") + antTaql;
    }

    // Reset antenna expr
    msSelection_p->setAntennaExpr("");

    // Append taql expr
    casacore::String taqlExpr = msSelection_p->getExpr(casacore::MSSelection::TAQL_EXPR);
    taqlExpr = (taqlExpr.empty() ? "" : " && ") + antTaqlExpr;
    msSelection_p->setTaQLExpr(taqlExpr);
  }

  casacore::Vector<casacore::Int> CTSelection::getRefAntIds(
      casacore::MSSelectableTable* msLike) {
    // get unique antenna2 ids from antenna columns
    NewCalTable nct = NewCalTable(*msLike->table());
    ROCTMainColumns ctmain(nct);
    casacore::Vector<casacore::Int> ant2 = ctmain.antenna2().getColumn();
    casacore::uInt nval = genSort(ant2, casacore::Sort::Ascending, 
      (casacore::Sort::QuickSort | casacore::Sort::NoDuplicates));
    ant2.resize(nval, casacore::True);
    return ant2;
  }

  void CTSelection::setAntSelectExclude(
    const casacore::Vector<casacore::Int>& antlist,
    const casacore::String& selection, casacore::MSSelectableTable* msLike,
    casacore::String& antselect, casacore::String& antexclude) {
    // Returns string list of selected and excluded antennas
    std::vector<casacore::Int> select, exclude;

    for (auto ant : antlist) {
      if (ant > 0) {
        select.push_back(ant);
      } else if (ant < 0) {
        exclude.push_back(abs(ant));
      } else { // == 0
        if (zeroIsExcluded(selection, msLike)) {
          exclude.push_back(0);
        } else {
          select.push_back(0);
        }
      }
    }
    antselect = casacore::MSSelection::indexExprStr(select);
    antexclude = casacore::MSSelection::indexExprStr(exclude);
  }

  bool CTSelection::zeroIsExcluded(const casacore::String& antennaExpr,
    casacore::MSSelectableTable* msLike) {
    // Check if antenna Id 0 is excluded in antennaExpr
    // Normally designated by -Id
    if (!antennaExpr.contains("!")) {
      return false;
    }

    // Check for !name
    NewCalTable nct = NewCalTable(*msLike->table());
    ROCTColumns ctcol(nct);
    casacore::String ant0name = ctcol.antenna().name().get(0);
    return antennaExpr.contains("0") || antennaExpr.contains(ant0name);
  }

} //# NAMESPACE CASA - END

