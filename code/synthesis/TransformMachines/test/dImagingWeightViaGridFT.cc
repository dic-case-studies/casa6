//# tVisModelData.cc: Tests the Synthesis model data serving
//# Copyright (C) 2011
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

#include <casa/Arrays/ArrayMath.h>
#include <components/ComponentModels/ComponentList.h>
#include <components/ComponentModels/ComponentShape.h>
#include <components/ComponentModels/Flux.h>
#include <tables/TaQL/ExprNode.h>
#include <measures/Measures/MeasTable.h>
#include <measures/Measures/MRadialVelocity.h>
#include <ms/MSSel/MSSelection.h>
#include <synthesis/MeasurementEquations/Imager.h>
#include <synthesis/TransformMachines/VisModelData.h>
#include <synthesis/TransformMachines/FTMachine.h>
#include <synthesis/TransformMachines/GridFT.h>
#include <msvis/MSVis/VisibilityIterator.h>
#include <msvis/MSVis/VisBuffer.h>
#include <casa/OS/Timer.h>
#include <casa/namespace.h>
#include <casa/OS/Directory.h>
#include <casa/Utilities/Regex.h>
#include <iomanip>

using namespace casa;
using namespace casacore;
using namespace std;

int main(int argc, char **argv){


  if (argc<2) {
    cout <<"Usage: dImagingWeightViaGridFT ms-table-name "<<endl;
    exit(1);
  }
  try{
    
    MeasurementSet myms(argv[1],Table::Old);
    
    String spwsel=String("*");
    MSSelection selector;
    selector.setSpwExpr(spwsel);
    TableExprNode exprNode=selector.toTableExprNode(&myms);
    cerr << "is expr null " << exprNode.isNull() << endl;
    MeasurementSet mssel(myms(exprNode));
    
    Int nx=3000;
    Quantity cellx(30, "arcsec");


    
    Block<Int> sort(0); // create empty block with sortColumns
    //VisibilityIterator vi(mssel,sort);
    VisibilityIterator vi(mssel,sort);
    VisBuffer vb(vi);
    MDirection myDir=vi.msColumns().field().phaseDirMeas(0);

    ///Let's make a complex image for use in the FTMachine
    ///
    String cImageName="TempComplex.image";
    String rImageName="TempReal.image";
    {
      Imager im(mssel, false, false);
      im.defineImage(nx, nx, cellx, cellx,
		     "I", myDir, 0, "mfs", 1, 0, 1, MFrequency(), MRadialVelocity(), Quantity(0.0, "km/s"), Vector<Int>(1,0));
      im.makeimage("psf", rImageName, cImageName);

    }

    
		 
 
   
    // tm.mark();
    /*    VisImagingWeight vWght(vi, "norm", Quantity(0.0, "Jy"), -1.0,  nx, nx,  cellx, cellx,0, 0);
    vi.origin();
    Float sumval=0;
    Float elmax, elmin;
    IPosition posmax, posmin;
    vi.useImagingWeight(vWght);
    for (vi.originChunks();vi.moreChunks(); vi.nextChunk()){
      for (vi.origin(); vi.more(); vi++){
	
	sumval+=sum(vb.imagingWeight());
	//vm.getModelVis(vb);

	minMax(elmin, elmax, posmin, posmax, vb.imagingWeight());
	cerr << "field " << vb.fieldId() << "  spw " << vb.spectralWindow() <<"max " <<  elmax << " " << posmax  << "  min " << elmin << "  " << posmin << endl;
	//cerr << vb.modelVisCube().xyPlane(40) << endl;
	//vb.visCube();
      }
    }
    Block<Matrix<Float> > density;
    vWght.getWeightDensity(density);
    cerr << "sum of grid "<< sum(density[0]) << " shape " << density[0].shape() << endl;
    cerr << "sum of weight "<< sumval << endl;
    */
    VisImagingWeight vWghtNat("natural");
    vi.useImagingWeight(vWghtNat);
    GridFT ft(Long(1000000), Int(200), "BOX",1.0, false, true);
    PagedImage<Complex> theImage("TempComplex.image");
    cerr << "SHAPE " << theImage.shape() << endl;
    ////////////////
    //PagedImage<Float> laloo(theImage.shape(), theImage.coordinates(),"density1.image");
    //laloo.put(density[0].reform((theImage.shape())));
    ///////

    Matrix<Float> dummy;
    // Initialize the ft machine with the backward grid image
    ft.initializeToSky(theImage,dummy,vb);
    Vector<Double> convFunc(2, 1.0);
    ft.modifyConvFunc(convFunc, 0, 1);
     for (vi.originChunks();vi.moreChunks();vi.nextChunk()) {
      for (vi.origin(); vi.more(); vi++) {
	 ft.put(vb, -1, true, FTMachine::PSF);
      }
     }
     Array<Float> griddedWght;
     // get the gridded weights
     ft.getGrid(griddedWght);
     cerr << std::setprecision(12) << "sum of grid "<< sum(griddedWght) << " shape " << griddedWght.shape() << endl;
     /*
     ////////////////
    PagedImage<Float> laloo2(theImage.shape(), theImage.coordinates(),"density2.image");
    laloo2.put(griddedWght.reform((theImage.shape())));
    ///////
     Block<Matrix<Float> > grids(1);
     grids[0].assign(griddedWght.reform(IPosition(2, nx, nx)));
     //Setup new VisImaging weight with weight density 
     VisImagingWeight vWght2(vi, "norm", Quantity(0.0, "Jy"), -1.0,  nx, nx,  cellx, cellx,0, 0);
     //VisImagingWeight vWght2(vi, grids, "abs", Quantity(0.0, "Jy"), -1,  Quantity(1.0, "arcsec"), Quantity(1.0, "arcsec"));
     float sumval2=0.0;
     vWght2.setWeightDensity(grids);
     vi.origin();
     vi.useImagingWeight(vWght2);
     for (vi.originChunks();vi.moreChunks(); vi.nextChunk()){
      for (vi.origin(); vi.more(); vi++){
	sumval2+=sum(vb.imagingWeight());
	elmax=0;
	elmin=0;
	minMax(elmin, elmax, posmin, posmax, vb.imagingWeight());
	cerr << "field " << vb.fieldId() << "  spw " << vb.spectralWindow() <<"max " <<  elmax << " " << posmax  << "  min " << elmin << "  " << posmin << endl;
      }
    }
    cerr << "sum of weight "<< sumval2 << endl;
     */

    } catch (AipsError x) {
    cout << "Caught exception " << endl;
    cout << x.getMesg() << endl;
    return(1);
  }


  cout << "Done" << endl;
  exit(0);

}
