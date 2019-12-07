/*
 * dMPICubeImager.cc: demonstrator of Cube imaging using MPI
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
#include <synthesis/ImagerObjects/CubeMajorCycleAlgorithm.h>
#include <synthesis/TransformMachines2/test/MakeMS.h>
#include <synthesis/ImagerObjects/SynthesisImagerVi2.h>
#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>

extern casa::Applicator casa::applicator;
int main(int argc, char **argv)
{
	using namespace casa;
	using namespace casacore;
	using namespace casa::test;
 	for (int lala=0; lala <2; ++lala){
	CubeMajorCycleAlgorithm *cmc=new CubeMajorCycleAlgorithm();
	//casa::applicator.defineAlgorithm(cmc);
	casa::applicator.init(argc, argv);
	Int numprocs = applicator.numProcs(); 
	cerr << "Number of procs: " << numprocs << endl;
	
	String msname("TestCube.ms");
	const Int numchan=10;
	MDirection phasecenter;
	Quantity freqBeg;
	Quantity freqWidth;
	cerr << "#####CONTROLLER ?" <<  applicator.isController() <<  " worker? " <<  applicator.isWorker() <<  endl;
	///Lets make the ms
	if(casa::applicator.isController())
 {
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
	Record selparsRec = datSel.toRecord();
	Record vecSelParsRec, vecImParsRec, vecGridParsRec;
	vecSelParsRec.defineRecord(String::toString(0), selparsRec);
	SynthesisParamsImage imSel;
	{
		imSel.imageName = "mpiCube_procs_"+String::toString(numprocs);
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
	Record imparsRec = imSel.toRecord();
	imparsRec.define("polrep", 0);
	vecImParsRec.defineRecord(String::toString(0), imparsRec);
	SynthesisParamsGrid gridSel;		
	{
		gridSel.ftmachine="gridft"; 
		gridSel.gridder = "gridft";
		gridSel.imageName = imSel.imageName;
	 
	}
	
	// Lets build the full imager and images
	SynthesisImagerVi2 imgr;
	imgr.selectData(datSel);
	//We do not use chanchunking in this model
	gridSel.chanchunks = 1;
	imgr.defineImage(imSel,gridSel);
	imgr.weight("natural");
	CountedPtr<SIImageStore> imstor = imgr.imageStore(0);
	// checking that psf,  residual and sumwt is allDone
	cerr << "shapes " <<  imstor->residual()->shape() <<  " " <<  imstor->psf()->shape() <<  " " <<  imstor->sumwt()->shape() <<  endl;
	//imstor->psf()->set(0.0);
	imstor->psf()->unlock();
	//imstor->residual()->set(0.0);
	//imstor->residual()->flush();
	imstor->residual()->unlock();
	imstor->sumwt()->unlock();
	Record gridparsRec = gridSel.toRecord();
	vecGridParsRec.defineRecord(String::toString(0), gridparsRec);
	Record controlRecord;   
	Record weightRecord=SynthesisUtilMethods::fillWeightRecord();
	
	controlRecord.define("dividebyweight", True);
	Vector<String>weightnames(1,"");
	controlRecord.define("weightnames", weightnames);
	Int rank(0);
	Bool assigned; //(casa::casa::applicator.nextAvailProcess(pwrite, rank));
	Bool allDone(false);
	Bool dopsf=False;
	Vector<Int> chanRange(2);
	//if(casa::applicator.isController()) {
        for (Int k=0; k < numchan; ++k) {
            assigned=casa::applicator.nextAvailProcess(*cmc, rank);
            cerr << "assigned "<< assigned << endl;
            while(!assigned) {
			    cerr << "SErial ? " << casa::applicator.isSerial() << endl;
                rank = casa::applicator.nextProcessDone(*cmc, allDone);
                cerr << "while rank " << rank << endl;
                Bool status;
                casa::applicator.get(status);
                if(status)
                    cerr << k << " rank " << rank << " successful " << endl;
                else
                    cerr << k << " rank " << rank << " failed " << endl;
                assigned = casa::applicator.nextAvailProcess(*cmc, rank);

            }

            ///send process info
            // put data sel params #1
			applicator.put(vecSelParsRec);
			// put image sel params #2
			applicator.put(vecImParsRec);
			// put gridders params #3
			applicator.put(vecGridParsRec);
			// put which channel to process #4
			chanRange.set(k);
			applicator.put(chanRange);
			// psf or residual CubeMajorCycleAlgorithm #5
			if(lala==1)
				dopsf=True;
			else
				dopsf=False;
			applicator.put(dopsf);
			// store modelvis and other controls #6
			applicator.put(controlRecord);
			//tell the weighting scheme
			applicator.put(weightRecord);
			/// Tell worker to process it 
            applicator.apply(*cmc);

        }
        // Wait for all outstanding processes to return
        rank = casa::applicator.nextProcessDone(*cmc, allDone);
        while (!allDone) {
			Int serialbug;
			//if(casa::applicator.isSerial())
				//casa::applicator.get(serialbug);// get that extra put
            Bool status;
            casa::applicator.get(status);
            if(status)
                cerr << "remainder rank " << rank << " successful " << endl;
            else
                cerr << "remainder rank " << rank << " failed " << endl;

            rank = casa::applicator.nextProcessDone(*cmc, allDone);
			if(casa::applicator.isSerial())
				allDone=true;
        }


	//sleep(10);
    }
  
  
	}  
  
	//sleep(10);
  
  
  
  
 exit(0); 
}
