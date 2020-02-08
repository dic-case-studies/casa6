//# tGridFT.cc: Tests the Synthesis model data serving
//# Copyright (C) 2016
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
//# $Id$
#include <measures/Measures/Stokes.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/StokesCoordinate.h>
#include <coordinates/Coordinates/Projection.h>
#include <casa/Arrays/ArrayMath.h>
#include <images/Images/ImageInterface.h>
#include <images/Images/PagedImage.h>
#include <images/Images/TempImage.h>
#include <components/ComponentModels/ComponentList.h>
#include <components/ComponentModels/ComponentShape.h>
#include <components/ComponentModels/Flux.h>
#include <tables/TaQL/ExprNode.h>
#include <measures/Measures/MeasTable.h>
#include <ms/MSSel/MSSelection.h>
#include <synthesis/TransformMachines2/BriggsCubeWeightor.h>
#include <synthesis/TransformMachines2/FTMachine.h>
#include <synthesis/TransformMachines2/GridFT.h>
#include <synthesis/TransformMachines2/SetJyGridFT.h>
#include <synthesis/TransformMachines2/WProjectFT.h>
#include <synthesis/TransformMachines2/SimpleComponentFTMachine.h>
#include <msvis/MSVis/VisImagingWeight.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <casa/namespace.h>
#include <casa/OS/Directory.h>
#include <casa/Utilities/Regex.h>
#include <synthesis/TransformMachines2/test/MakeMS.h>
using namespace casacore;
using namespace casa;
using namespace casacore;
using namespace casa::refim;
using namespace casacore;
using namespace casa::test;

FTMachine* makeFTM(String ftmachine, MPosition loc, MDirection dir){
  if(ftmachine=="GridFT")
    return new refim::GridFT(1000000, 16, "SF", loc, 1.0, false, false);
  if(ftmachine=="SetJyGridFT")
    return new refim::SetJyGridFT(1000000, 16, "SF", loc, dir, 1.0, false, false);
  if(ftmachine=="WProjectFT")
    return new refim::WProjectFT(50, loc, 1000000, 16);
  return nullptr;

}

Int main(/*int argc, char **argv*/){
  try{
    MDirection thedir(Quantity(20.0, "deg"), Quantity(20.0, "deg"));
    String msname("Test.ms");
    MakeMS::makems(msname, thedir);
    MeasurementSet thems(msname, Table::Update);
    thems.markForDelete();
    vi::VisibilityIterator2 vi2(thems,vi::SortColumns(),true);
    vi::VisBuffer2 *vb=vi2.getVisBuffer();
    VisImagingWeight viw("natural");
    vi2.useImagingWeight(viw);
    ComponentList cl;
    SkyComponent otherPoint(ComponentType::POINT);
    otherPoint.flux() = Flux<Double>(6.66e-2, 0.0, 0.0, 0.00000);
    otherPoint.shape().setRefDirection(thedir);
    cl.add(otherPoint);
    Matrix<Double> xform(2,2);
    xform = 0.0;
    xform.diagonal() = 1.0;
    DirectionCoordinate dc(MDirection::J2000, Projection::SIN, Quantity(20.0,"deg"), Quantity(20.0, "deg"),
			   Quantity(0.5, "arcsec"), Quantity(0.5,"arcsec"),
			   xform, 50.0, 50.0, 999.0, 
			   999.0);
    Vector<Int> whichStokes(1, Stokes::I);
    StokesCoordinate stc(whichStokes);
    SpectralCoordinate spc(MFrequency::LSRK, 1.5e9, 1e6, 0.0 , 1.420405752E9);
    CoordinateSystem cs;
    cs.addCoordinate(dc); cs.addCoordinate(stc); cs.addCoordinate(spc);
    TempImage<Complex> im(IPosition(4,100,100,1,100), cs);//, "gulu.image");
    std::vector<String> ftmachines={"GridFT"};
    FTMachine * ftm=nullptr;
    for (uInt k=0; k < ftmachines.size(); ++k){
     cerr << "########### Test imagingweight" <<  << endl; 
      im.set(Complex(0.0));
      MPosition loc;
      MeasTable::Observatory(loc, MSColumns(thems).observation().telescopeName()(0));
      ftm=makeFTM(ftmachines[k], loc, thedir);
      Matrix<Float> weight;
      vi2.originChunks();
      vi2.origin();
      ftm->initializeToSky(im, weight, *vb);
	  String err;
	  Record ftmRec;
	  ftm->toRecord(err, ftmRec);
      refim::SimpleComponentFTMachine cft;
	  BriggsCubeWeightor brwgt("norm", Quantity(0.0, "Jy"), 0.0);
	  BriggsCubeWeightor brwgt2("norm", Quantity(0.0, "Jy"), 0.0);
	  brwgt.init(vi2, im, ftmRec);
	  brwgt2.initImgWeightCol(vi2, im, ftmRec);
	  Matrix<Float> wgt1;
	  Matrix<Float> wgt2;
      for (vi2.originChunks();vi2.moreChunks(); vi2.nextChunk()){
	for (vi2.origin(); vi2.more(); vi2.next()){
		brwgt.weightUniform(wgt1, *vb);
		brwgt2.weightUniform(wgt2, *vb);
		//cerr << "Diff " << max(wgt2-wgt1) << "   max " << max(wgt2) << endl; 
	  //cft.get(*vb, cl);
	  //vb->setVisCube(vb->visCubeModel());
	  //ftm->put(*vb);	
		AlwaysAssertExit(near(0.0, max(wgt2-wgt1), 1.0e-6));
		AlwaysAssertExit(near(0.0, min(wgt2-wgt1), 1.0e-6));
	}
      }
      //ftm->finalizeToSky();
      
	}
     } catch (AipsError x) {
    cout << "Caught exception " << endl;
    cout << x.getMesg() << endl;
    return(1);
  }
  cerr <<"OK" << endl;
  exit(0);
}
