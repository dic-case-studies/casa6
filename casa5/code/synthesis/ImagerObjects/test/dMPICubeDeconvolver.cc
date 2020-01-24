/*
 * dMPICubeDeconvolver.cc: demonstrator of Cube deconvolution using MPI
//# Copyright (C) 2020
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
 *  Created on: Jan 5, 2020
 *      Author: kgolap
 */
#include <synthesis/ImagerObjects/CubeMinorCycleAlgorithm.h>
#include <synthesis/ImagerObjects/SynthesisDeconvolver.h>
#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>

extern casa::Applicator casa::applicator;
int main(int argc, char **argv)
{
  using namespace casa;
  using namespace casacore;
  if (argc<2) {
    cout <<"Usage: dMPICubeDeconvolver rootname "<<endl;
    exit(1);
  }
  String rootname(argv[1]);
  CubeMinorCycleAlgorithm *cmc=new CubeMinorCycleAlgorithm();
  //casa::applicator.defineAlgorithm(cmc);
  ///argv and argc are needed just to callthe right overloaded init
  casa::applicator.init(argc, argv);
  Int numprocs = applicator.numProcs(); 
  cerr << "Number of procs: " << numprocs << endl;
	
  cerr << "#####CONTROLLER ?" <<  applicator.isController() <<  " worker? " <<  applicator.isWorker() <<  endl;
  Record dminfo;
  if(casa::applicator.isController())
  {
    String imagename="mpiDeconv_procs_"+String::toString(numprocs);
    {
    Table tb(rootname+".psf");
    tb.deepCopy(imagename+".psf", dminfo, Table::New, True); 
    Table tb2(rootname+".residual");
    tb2.deepCopy(imagename+".residual", dminfo, Table::New, True);
    ////sumwt is needed as disk based imagestore check for it ...though in the
    ///case of deconvolution it is never used.
    Table tb3(rootname+".sumwt");
    tb3.deepCopy(imagename+".sumwt", dminfo, Table::New, True);
    }
  
	
      SynthesisParamsDeconv decpars;
      
      {
	decpars.imageName = rootname;
	decpars.algorithm="hogbom";
      }
      Record iterBotRec;
      {
	iterBotRec.define("cycleniter", Int(100));
	iterBotRec.define("cyclethreshold", Float(0.0));
	iterBotRec.define("thresholdreached", False);
	iterBotRec.define("loopgain", Float(0.1));
	iterBotRec.define("nsigma", Float(0.0));
	
	
      }
      Record decParsRec = decpars.toRecord();
      CountedPtr<SIImageStore> imstor = new SIImageStore(imagename, False);
      // checking that psf,  residual available
      cerr << "shapes " <<  imstor->residual()->shape() <<  " " <<  imstor->psf()->shape() <<  " " <<  imstor->sumwt()->shape() <<  endl;
      Int numchan=imstor->residual()->shape()[3];
      String psfname=imstor->psf()->name();
      Vector<ImageBeamSet> chanBeams(numchan);
      ImageBeamSet fullBeamSet=imstor->getBeamSet();
      cerr << "fullBeamSet shp " << fullBeamSet.shape() << "    " << fullBeamSet.getBeams() << endl;
      for (Int k =0 ; k <numchan; ++k){
        chanBeams[k]=imstor->getChannelBeamSet(k);
      }
      Float psl=imstor->getPSFSidelobeLevel();
      String residualname=imstor->residual()->name();
      String maskname=imstor->mask()->name();
      String modelname=imstor->model()->name();
      imstor->psf()->unlock();
      imstor->residual()->unlock();
      //Lets make a simple mask and 0 model
      imstor->mask()->set(1.0);
      imstor->model()->set(0.0);
      imstor->mask()->unlock();
      imstor->model()->unlock();
      imstor->releaseLocks();
      imstor.reset();
	
      Int rank(0);
      Bool assigned; 
      Bool allDone(false);
      Vector<Int> chanRange(2);
      Record beamsetRec;
      for (Int k=0; k < numchan; ++k) {
	assigned=casa::applicator.nextAvailProcess(*cmc, rank);
	cerr << "assigned "<< assigned << endl;
	while(!assigned) {
	  cerr << "SErial ? " << casa::applicator.isSerial() << endl;
	  rank = casa::applicator.nextProcessDone(*cmc, allDone);
	  cerr << "while rank " << rank << endl;
          Vector<Int> chanprocessed;
          casa::applicator.get(chanprocessed);
	  Record status;
	  casa::applicator.get(status);
	  if(status.nfields())
	    cerr << k << " rank " << rank << " successful " << endl;
	  else
	    cerr << k << " rank " << rank << " failed " << endl;
	  cerr <<"rank " << rank << " return rec "<< status << endl;
	  assigned = casa::applicator.nextAvailProcess(*cmc, rank);
	  
	}

            ///send process info
            // put data sel params #1
			applicator.put(decParsRec);
			// put image sel params #2
			applicator.put(iterBotRec);
			// put which channel to process #3
			chanRange.set(k);
			applicator.put(chanRange);
			// psf  #4
			applicator.put(psfname);
			// residual #5
			applicator.put(residualname);
			// model #6
			applicator.put(modelname);
			// mask #7
			applicator.put(maskname);
                        //#8 beamset
                        //need to use local variable for serial case
                        beamsetRec=chanBeams[k].toRecord();
                        //cerr << "beamsetRec " << beamsetRec << endl;
                        applicator.put(beamsetRec);
                        //#9 psf lobe level
                        applicator.put(psl);
			/// Tell worker to process it 
            applicator.apply(*cmc);

        }
        // Wait for all outstanding processes to return
        rank = casa::applicator.nextProcessDone(*cmc, allDone);
        while (!allDone) {
          Vector<Int> chanprocessed;
          casa::applicator.get(chanprocessed);
          Record status;
          casa::applicator.get(status);
            if(status.nfields() >0)
                cerr << "remainder rank " << rank << " successful " << endl;
            else
                cerr << "remainder rank " << rank << " failed " << endl;

            rank = casa::applicator.nextProcessDone(*cmc, allDone);
			if(casa::applicator.isSerial())
				allDone=true;
        }


	//sleep(10);
    }
  
  
  
  
  
 exit(0); 
}
