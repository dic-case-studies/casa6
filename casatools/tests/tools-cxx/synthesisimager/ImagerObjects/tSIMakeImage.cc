/*
 * dMPIMakeCube.cc: demonstrator of Cube imaging using MPI
//# Copyright (C) 2019
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by
//# the Free Software Foundation; either version 3 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
//# License for more details.
//#
//# https://www.gnu.org/licenses/
//#
//# Queries concerning CASA should be submitted at
//#        https://help.nrao.edu
//#
//#        Postal address: CASA Project Manager 
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
 *  Created on: Oct 9, 2019
 *      Author: kgolap
 */
#include <synthesis/ImagerObjects/CubeMakeImageAlgorithm.h>
#include <synthesis/TransformMachines2/test/MakeMS.h>
#include <synthesis/ImagerObjects/SynthesisImagerVi2.h>
#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>

extern casa::Applicator casa::applicator;
int main()
{
	using namespace casa;
	using namespace casacore;
	using namespace casa::test;
	
	
	
	String msname("TestCube.ms");
	const Int numchan=10;
	MDirection phasecenter;
	Quantity freqBeg;
	Quantity freqWidth;
	
	///Lets make the ms
    try{   
	{
		MDirection thedir(Quantity(20.0, "deg"), Quantity(20.0, "deg"));
		
		const Int numchan=10;
		MakeMS::makems(msname, thedir, 1.0e9, 1e8, numchan, 20);
		//MakeMS::makems(msname, thedir, 1.5e9, 1e6, numchan, 20);
		MeasurementSet thems(msname, Table::Update);
		MSColumns(thems).data().fillColumn(Matrix<Complex>(4,numchan, Complex(6.66e-2)));
		MSColumns(thems).correctedData().fillColumn(Matrix<Complex>(4,numchan, Complex(6.66e-2)));
		phasecenter=MSFieldColumns(thems.field()).phaseDirMeas(0,0.0);
		freqBeg=MSSpWindowColumns(thems.spectralWindow()).chanFreqQuant()(0)(IPosition(1,0));
		Int ndataChan=MSSpWindowColumns(thems.spectralWindow()).numChan()(0);
		freqWidth=MSSpWindowColumns(thems.spectralWindow()).chanFreqQuant()(0)(IPosition(1,ndataChan-1));
		freqWidth-=freqBeg;
		freqWidth /= Double(numchan);
		
		thems.flush();
	}
	SynthesisParamsSelect datSel;
	{datSel.msname=msname;datSel.spw="0";datSel.freqbeg="";datSel.freqend="";
			datSel.freqframe=MFrequency::LSRK;datSel.field="0";datSel.antenna=""; datSel.timestr="";datSel.scan="";datSel.obs="";datSel.state="";datSel.uvdist="";
			datSel.taql="";datSel.usescratch=false;datSel.readonly=true;datSel.incrmodel=false;
			datSel.datacolumn="corrected";
	}
	SynthesisParamsImage imSel;
	{
		imSel.imageName = "SIMICube";
		imSel.imsize.resize(2); imSel.imsize.set(100);
		imSel.cellsize.resize(2); imSel.cellsize.set( Quantity(10.0,"arcsec") );
		imSel.stokes="I";imSel.phaseCenter=phasecenter;imSel.phaseCenterFieldId=-1;
		imSel.projection=Projection::SIN;imSel.useNCP=False;imSel.startModel=Vector<String>(0);
		imSel.overwrite=False;imSel.pseudoi=False;imSel.nchan=numchan;imSel.mode="cube";
		imSel.start="";imSel.step="";imSel.chanStart=0;imSel.chanStep=1;    imSel.freqStart=freqBeg;imSel.freqStep=freqWidth;imSel.velStart=Quantity(0,"");
        imSel.velStep=Quantity(0,"");imSel.veltype=String("radio");
        imSel.restFreq=Vector<Quantity>(1,Quantity(1.420, "GHz"));
        imSel.refFreq = Quantity(0,"Hz");imSel.frame = ""; imSel.freqFrame=MFrequency::LSRK;imSel.sysvel="";imSel.sysvelframe=""; imSel.sysvelvalue=Quantity(0.0,"m/s");imSel.nTaylorTerms=1;imSel.deconvolver="hogbom";
	}	
	SynthesisParamsGrid gridSel;		
	{
		gridSel.ftmachine="gridft"; 
		gridSel.gridder = "gridft";
		gridSel.imageName = imSel.imageName;
	 
	}
	{
	// Lets build the full imager and images
	SynthesisImagerVi2 imgr;
	imgr.selectData(datSel);
	//We do not use chanchunking in this model
	gridSel.chanchunks = 1;
	imgr.defineImage(imSel,gridSel);
	imgr.weight("natural");
        imgr.makeImage("observed", "dirty.image", "dirty_complex.image");
        imgr.makeImage("psf", "psf.real", "psf.cmplx");
        }

        PagedImage<Float> dirty( "dirty.image");
        ////////////////////////////////////
        CoordinateSystem csys=dirty.coordinates();
        IPosition shp=dirty.shape();
        PagedImage<Bool> boolim(shp, csys, "boolean.image");
        boolim.set(False);
        boolim.putAt(True, shp/2);
        PagedImage<Int> intim(shp, csys, "integer.image");
        boolim.set(22);
        boolim.putAt(-23, shp/2);

        ////////////////////////////////////
        LatticeExprNode LENMaxRes = max(dirty);
        AlwaysAssertExit(near(6.66e-2, LENMaxRes.getFloat(), 1.0e-5));
        PagedImage<Float> psf( "psf.real");
        LatticeExprNode LENMaxPsf = max(psf);
        AlwaysAssertExit(near(1.0, LENMaxPsf.getFloat(), 1.0e-5));
  
   }catch( AipsError e ){
    cout << "Exception ocurred." << endl;
    cout << e.getMesg() << endl;
    exit(-1);
  }   
   
  
        cout <<"OK"<< endl;
        exit(0); 
}
