//# CTSelection_GTest.cc: Google Test program for selecting on a NewCalTable
//# Copyright (C) 2017
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#

#include <synthesis/CalTables/NewCalTable.h>
#include <synthesis/CalTables/CTColumns.h>
#include <synthesis/CalTables/CTInterface.h>
#include <synthesis/CalTables/CTSelection.h>
#include <ms/MSSel/MSSelectionTools.h>
#include <casa/Arrays/ArrayLogical.h>

#include <gtest/gtest.h>

// <summary>
// Google Test program for CTSelection class.
// Adapted from tCTSelection.cc legacy test.
// </summary>

using namespace casacore;
using namespace casa;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


TEST( CTSelectionTest, CTFieldSpwSelection) {
  // Select fields and spws, then reset selection
  // Make a testing NewCalTable (Table::Memory)
  uInt nFld(20), nAnt(10), nSpw(4), nObs(1), nScan(20),nTime(10);
  Vector<Int> nChan(nSpw, 32);
  Double refTime(4832568000.0); // 2012 Jan 06 @ noon
  Double tint(60.0);
  NewCalTable tnct("tCTSelection1.ct","Complex", nObs, nScan, nTime, nAnt,
    nSpw, nChan, nFld, refTime, tint);
  // some sanity checks on the test NewCalTable
  EXPECT_EQ( tnct.tableType(), Table::Memory );
  EXPECT_EQ( tnct.nrow(), nObs * nScan * nTime * nSpw * nAnt);

  // Make selection strings
  ROCTColumns ctc(tnct);
  // Field selection by name
  Vector<String> fldnames = ctc.field().name().getColumn();
  Vector<Int> fldids(3);
  fldids(0)=3; fldids(1)=10; fldids(2)=17;
  String fieldsel("");
  for (uInt i=0; i<fldids.nelements(); ++i) {
    if (i>0) fieldsel += ",";
    fieldsel += fldnames(fldids(i));
  }
  // Spw selection
  Vector<Int> spwids(2);
  spwids(0)=1; spwids(1)=3;
  String spwsel("");
  for (uInt i=0; i<spwids.nelements(); ++i) {
    if (i>0) spwsel += ",";
    spwsel += String::toString(spwids(i));
  }

  // Create CTSelection and set selections
  CTSelection cts;
  cts.setFieldExpr(fieldsel);
  cts.setSpwExpr(spwsel);
  CTInterface cti(tnct);
  TableExprNode ten = cts.toTableExprNode(&cti);

  // check CTSelection accessors for lists of selected items 
  EXPECT_TRUE( allEQ(cts.getFieldList(), fldids) );
  EXPECT_TRUE( allEQ(cts.getSpwList(), spwids) );

  // Create selected table and check nrow
  NewCalTable selnct(tnct);
  getSelectedTable(selnct,tnct,ten,"");
  ASSERT_EQ( selnct.nrow(), 
    nAnt * nTime * spwids.nelements() * fldids.nelements());

  // get columns from selected table and check elements
  ROCTMainColumns ctmc(selnct);
  Vector<Int> selfieldcol, selspwcol;
  ctmc.fieldId().getColumn(selfieldcol);
  ctmc.spwId().getColumn(selspwcol);

  Vector<Bool> fldok(selfieldcol.nelements(), false);
  for (uInt i=0; i<fldids.nelements(); ++i) {
    fldok |= (selfieldcol==fldids(i));
  }
  Vector<Bool> spwok(selspwcol.nelements(), false);
  for (uInt i=0; i<spwids.nelements(); ++i) {
    spwok |= (selspwcol==spwids(i));
  }

  EXPECT_TRUE( allEQ(fldok, true) );
  EXPECT_TRUE( allEQ(spwok, true) );
}

TEST( CTSelectionTest, CTRefAntSelection ) {
  // Antenna selection for refant-based table.
  // Includes ref antenna (1), to check that not all rows (baselines) were
  // selected. Cannot use virtual table with no antenna subtable for this test.
  String casapath = getenv("CASAPATH");
  String datapath = casapath.substr(0,casapath.find(" "));
  String caldata = "/casatestdata/unittest/gaincal/ngc5921.ref2a.gcal";
  //std::cout << "Test refant-based caltable: " << datapath << caldata << endl; 
  NewCalTable tnct(datapath + caldata);
  uInt nTime(12); 

  // Make antenna selection string: 0, 1, 6, 8, 9
  ROCTColumns ctc(tnct);
  Vector<String> antnames=ctc.antenna().name().getColumn();
  Vector<Int> antids(5);
  antids(0)=0; antids(1)=1; antids(2)=6; antids(3)=8; antids(4)=9;
  String antsel("");
  for (uInt i=0; i<antids.nelements(); ++i) {
    if (i>0) antsel += ",";
    antsel += antnames(antids(i));
  }

  CTSelection cts;
  cts.setAntennaExpr(antsel);
  CTInterface cti(tnct);
  TableExprNode ten = cts.toTableExprNode(&cti);

  // check CTSelection accessors for lists of selected antennas 
  Vector<Int> ant2ids(1, 1); // one refant, ID=1
  EXPECT_TRUE( allEQ(cts.getAntenna1List(), antids) );
  EXPECT_TRUE( allEQ(cts.getAntenna2List(), ant2ids) );

  // Create selected table and check nrow
  NewCalTable selnct(tnct);
  getSelectedTable(selnct, tnct, ten, "");
  ASSERT_EQ( selnct.nrow(), nTime * antids.nelements());

  // get antenna1 column from selected table and check elements
  ROCTMainColumns ctmc(selnct);
  Vector<Int> selantcol;
  ctmc.antenna1().getColumn(selantcol);
  Vector<Bool> antok(selantcol.nelements(), false);
  for (uInt i=0; i<antids.nelements(); ++i) {
    antok |= (selantcol==antids(i));
  }
  EXPECT_TRUE(allEQ(antok, true));

  // Repeat test using CTSelection constructor with selection expressions
  CTSelection* newcts = new CTSelection(tnct, casacore::MSSelection::PARSE_NOW,
    "", antsel, "", "", "", "", "", "");
  // check CTSelection accessors for lists of selected items 
  EXPECT_TRUE( allEQ(newcts->getAntenna1List(), antids) );
  EXPECT_TRUE( allEQ(newcts->getAntenna2List(), ant2ids) );

  // Create selected table and check nrow
  ten = newcts->getTEN();
  getSelectedTable(selnct, tnct, ten, "");
  ASSERT_EQ(selnct.nrow(), nTime * antids.nelements());
}

TEST( CTSelectionTest, CTPureAntSelection ) {
  // Antenna selection for pure ant-based table.
  // Cannot use virtual table with no antenna subtable for this test.
  String casapath = getenv("CASAPATH");
  String datapath = casapath.substr(0,casapath.find(" "));
  String caldata = "/casatestdata/unittest/gaincal/ngc5921.ref1a.gcal";
  //std::cout << "Test antenna-based caltable: " << datapath << caldata << endl;
  NewCalTable tnct(datapath + caldata);
  uInt nTime(7); 

  // Make antenna selection string: 0, 1, 6, 8, 9
  ROCTColumns ctc(tnct);
  Vector<String> antnames=ctc.antenna().name().getColumn();
  Vector<Int> antids(5);
  antids(0)=0; antids(1)=1; antids(2)=6; antids(3)=8; antids(4)=9;
  String antsel("");
  for (uInt i=0; i<antids.nelements(); ++i) {
    if (i>0) antsel += ",";
    antsel += antnames(antids(i));
  }

  CTSelection cts;
  cts.setAntennaExpr(antsel);
  CTInterface cti(tnct);
  TableExprNode ten = cts.toTableExprNode(&cti);

  // check CTSelection accessors for lists of selected antennas 
  Vector<Int> ant2ids; // No antenna2 list
  EXPECT_TRUE( allEQ(cts.getAntenna1List(), antids) );
  EXPECT_TRUE( allEQ(cts.getAntenna2List(), ant2ids) );

  // Create selected table and check nrow
  NewCalTable selnct(tnct);
  getSelectedTable(selnct, tnct, ten, "");
  ASSERT_EQ(selnct.nrow(), nTime * antids.nelements());

  // get columns from selected table and check elements
  ROCTMainColumns ctmc(selnct);
  Vector<Int> selantcol;
  ctmc.antenna1().getColumn(selantcol);
  Vector<Bool> antok(selantcol.nelements(), false);
  for (uInt i=0; i<antids.nelements(); ++i) {
    antok |= (selantcol==antids(i));
  }
  
  EXPECT_TRUE(allEQ(antok, true));

  // Repeat test using CTSelection constructor with selection expressions
  CTSelection* newcts = new CTSelection(tnct, casacore::MSSelection::PARSE_NOW,
    "", antsel, "", "", "", "", "", "");
  // check CTSelection accessors for lists of selected items 
  EXPECT_TRUE( allEQ(newcts->getAntenna1List(), antids) );
  EXPECT_TRUE( allEQ(newcts->getAntenna2List(), ant2ids) );

  // Create selected table and check nrow
  ten = newcts->getTEN();
  getSelectedTable(selnct, tnct, ten, "");
  ASSERT_EQ(selnct.nrow(), nTime * antids.nelements());
}

TEST( CTSelectionTest, CTBaselineAntSelection ) {
  // Antenna selection for baseline-based table.
  // Cannot use virtual table with no antenna subtable for this test.
  String casapath = getenv("CASAPATH");
  String datapath = casapath.substr(0,casapath.find(" "));
  String caldata = "/casatestdata/unittest/plotms/a_mueller.uvcont.tbl";
  //std::cout << "Test antenna-based caltable: " << datapath << caldata << endl; 
  NewCalTable tnct(datapath + caldata);

  // Make antenna selection string: 0, 1, 6, 8, 9
  ROCTColumns ctc(tnct);
  Vector<String> antnames=ctc.antenna().name().getColumn();
  uInt nAnt = antnames.nelements();

  Vector<Int> antids(5);
  antids(0)=0; antids(1)=1; antids(2)=6; antids(3)=8; antids(4)=9;
  String antsel("");
  for (uInt i=0; i<antids.nelements(); ++i) {
    if (i>0) antsel += ",";
    antsel += antnames(antids(i));
  }

  CTSelection cts;
  cts.setAntennaExpr(antsel);
  cts.setSpwExpr("16");
  cts.setFieldExpr("3");
  cts.setScanExpr("7");
  CTInterface cti(tnct);
  TableExprNode ten = cts.toTableExprNode(&cti);

  // check CTSelection accessors for lists of selected antennas 
  EXPECT_TRUE( allEQ(cts.getAntenna1List(), antids) );
  ASSERT_EQ( cts.getAntenna2List().size(), nAnt );
}

