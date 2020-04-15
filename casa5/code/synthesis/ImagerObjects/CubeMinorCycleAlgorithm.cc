//# CubeMinorCycleAlgorithm.cc: implementation of class to deconvolve cube in parallel/serial 
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
//# $Id$
#include <casacore/lattices/Lattices/LatticeLocker.h>
#include <synthesis/ImagerObjects/CubeMinorCycleAlgorithm.h>
#include <synthesis/ImagerObjects/SynthesisDeconvolver.h>
#include <casa/Containers/Record.h>
#include <synthesis/ImagerObjects/SimpleSIImageStore.h>
#include <imageanalysis/Utilities/SpectralImageUtil.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN
extern Applicator applicator;

  CubeMinorCycleAlgorithm::CubeMinorCycleAlgorithm() : myName_p("CubeMinorCycleAlgorithm"), autoMaskOn_p(False),chanFlag_p(0), status_p(False){
	
}
CubeMinorCycleAlgorithm::~CubeMinorCycleAlgorithm() {
	
}
	
void CubeMinorCycleAlgorithm::get() {
  ///New instructions reset previous state
  reset();
	//cerr << "in get for child process " << applicator.isWorker() << endl;
	Record decParsRec;
     
	// get deconv params record #1
	applicator.get(decParsRec);
	// get iter control rec #2
	applicator.get(iterBotRec_p);
	// channel range to deconvolve #3
	applicator.get(chanRange_p);
	//get psf image name #4
	applicator.get(psfName_p);
	//get residual name #5
	applicator.get(residualName_p);
	//get model name #6
	applicator.get(modelName_p);
	//get mask name #7
	applicator.get(maskName_p);
        // get pb name #8
        applicator.get(pbName_p);
        //get beamsetrec #9
        //applicator.get(beamsetRec_p);
        //get psfsidelobelev #9
        applicator.get(psfSidelobeLevel_p);
        //get chanflag #10
        Record chanflagRec;
        applicator.get(chanflagRec);
        chanFlag_p.resize();
        chanflagRec.get("chanflag", chanFlag_p);
        statsRec_p=chanflagRec.asRecord("statsrec");
	//cerr <<"GET chanRange " << chanRange_p << endl;
	decPars_p.fromRecord(decParsRec);
	
	
	
}
void CubeMinorCycleAlgorithm::put() {
	
  ///# 1  chanrange processed 
  applicator.put(chanRange_p);
	//cerr << "in put " << status_p << endl;
  //#2 chanflag
  Record chanflagRec;
  chanflagRec.define("chanflag", chanFlag_p);
  chanflagRec.defineRecord("statsrec", statsRec_p);
  applicator.put(chanflagRec);
  ///#3 return record of deconvolver
  // cerr << "nfield " << returnRec_p.nfields() << endl;
  SIMinorCycleController::compressSummaryMinor(returnRec_p);
  //Matrix<Double> lala(returnRec_p.asArrayDouble("summaryminor"));
  //cerr << "chanRange " << chanRange_p << " summaryminor " << lala.shape()   << endl;
  //cerr << "model " << lala.row(2) << endl;
  //cerr << "thresh " << lala.row(3) << endl;
  //cerr << "imageid " << lala.row(5) << endl;
  // applicator.put(lala);
  //returnRec_p.removeField("summaryminor");
  ////TESTOO
  // cerr << "nfield " << returnRec_p.nfields() << endl;
  //for (uInt k =0; k <  returnRec_p.nfields() ; ++k){
  //  cerr << " name " << returnRec_p.name(k) << endl;

  //}

  //Record laloo;
  applicator.put(returnRec_p);	
	
}
	
void CubeMinorCycleAlgorithm::task(){
	status_p = False;
	try{
          SynthesisDeconvolver subDeconv;
          subDeconv.setupDeconvolution(decPars_p);
          std::shared_ptr<SIImageStore> subimstor=subImageStore();
          //ImageBeamSet bs=ImageBeamSet::fromRecord(beamsetRec_p);
          ImageBeamSet bs=(subimstor->psf()->imageInfo()).getBeamSet();
          subimstor->setBeamSet(bs);
          subimstor->setPSFSidelobeLevel(psfSidelobeLevel_p);
          LatticeLocker lock1 (*(subimstor->model()), FileLocker::Write);
          subDeconv.initMinorCycle(subimstor);
          if(autoMaskOn_p){
            subDeconv.setChanFlag(chanFlag_p);
	    subDeconv.setRobustStats(statsRec_p);
	    //cerr << "STATSRec " << statsRec_p << endl;
            subDeconv.setIterDone(iterBotRec_p.asInt("iterdone"));
            subDeconv.setPosMask(subimstor->tempworkimage());
            subDeconv.setAutoMask();
          }
          //subDeconv.setupMask();
          returnRec_p=subDeconv.executeCoreMinorCycle(iterBotRec_p);
          chanFlag_p.resize();
          chanFlag_p=subDeconv.getChanFlag();
	  statsRec_p=Record();
	  statsRec_p=subDeconv.getRobustStats();
          writeBackToFullImage(modelName_p, chanRange_p[0], chanRange_p[1], (subimstor->model()));
          if(autoMaskOn_p){
            writeBackToFullImage(posMaskName_p, chanRange_p[0], chanRange_p[1], (subimstor->tempworkimage()));
            writeBackToFullImage(maskName_p, chanRange_p[0], chanRange_p[1], (subimstor->mask()));
          }
        }
        catch (AipsError x) {
          cerr << "Exception: " << x.getMesg() << endl;
          returnRec_p=Record();
        }
        catch(...){
          cerr << "Unknown exception" << endl;
          returnRec_p=Record();
        }
        
       	status_p = True;
}
String&	CubeMinorCycleAlgorithm::name(){
	return myName_p;
}

std::shared_ptr<SIImageStore> CubeMinorCycleAlgorithm::subImageStore(){
  std::shared_ptr<ImageInterface<Float> >subpsf=nullptr;
  std::shared_ptr<ImageInterface<Float> >subresid=nullptr;
  std::shared_ptr<ImageInterface<Float> >submodel=nullptr;
  std::shared_ptr<ImageInterface<Float> > submask=nullptr;
  std::shared_ptr<ImageInterface<Float> > subpb=nullptr;
  std::shared_ptr<ImageInterface<Float> > subposmask=nullptr;
	Int chanBeg=0;
	Int chanEnd=0;
	chanBeg=chanRange_p[0];
	chanEnd=chanRange_p[1];
	//cerr << "chanBeg " << chanBeg << " chanEnd " << chanEnd << " imId " << imId << endl;
        
	
        makeTempImage(subpsf, psfName_p, chanBeg, chanEnd);
        makeTempImage(subresid, residualName_p, chanBeg, chanEnd);
        makeTempImage(submodel, modelName_p, chanBeg, chanEnd, True);
        makeTempImage(submask, maskName_p, chanBeg, chanEnd, True);
        //       PagedImage<Float> model(modelName_p, TableLock::UserLocking);
        //   submodel.reset(SpectralImageUtil::getChannel(model, chanBeg, chanEnd, true));
        if(!pbName_p.empty()){
           makeTempImage(subpb, pbName_p, chanBeg, chanEnd);
        }
        if(iterBotRec_p.isDefined("posmaskname") ){
            iterBotRec_p.get("posmaskname", posMaskName_p);
            if(Table::isReadable(posMaskName_p)){
                makeTempImage(subposmask, posMaskName_p, chanBeg, chanEnd, True);
                if(subposmask)
                  autoMaskOn_p=True;
              }
          }

            std::shared_ptr<SIImageStore> subimstor(new SimpleSIImageStore(submodel, subresid, subpsf, nullptr, nullptr, submask, nullptr, nullptr, subpb, nullptr, subposmask));
        
	//cerr << "subimagestor TYPE" << subimstor->getType() << endl;
	return subimstor;
}
  void CubeMinorCycleAlgorithm::makeTempImage(std::shared_ptr<ImageInterface<Float> >& outptr,  const String& imagename, const Int chanBeg, const Int chanEnd, const Bool writelock){
    PagedImage<Float> im(imagename, writelock ? TableLock::UserLocking : TableLock::UserNoReadLocking);
    SubImage<Float> *tmpptr=nullptr;
    if(writelock)
              im.lock(FileLocker::Write, 30);
    ////TESTOO
    //outptr.reset(SpectralImageUtil::getChannel(im, chanBeg, chanEnd, writelock));
    

    ///END of TESTOO
    
    tmpptr=SpectralImageUtil::getChannel(im, chanBeg, chanEnd, false);
    if(tmpptr){
      IPosition tileshape=tmpptr->shape();
      tileshape[2]=1; tileshape[3]=1;
      TiledShape tshape(tmpptr->shape(),tileshape);
      outptr.reset(new TempImage<Float>(tshape, tmpptr->coordinates()));
      outptr->copyData(*tmpptr);
      //cerr << "IMAGENAME " << imagename << " masked " << im.isMasked() << " tmptr  " << tmpptr->isMasked() << endl;
      if(tmpptr->isMasked()){
	  outptr->makeMask ("mask0", true, true, false, true);
	  outptr->pixelMask().put(tmpptr->getMask());
      }
      ImageInfo iinfo=tmpptr->imageInfo();
      outptr->setImageInfo(iinfo);
      delete tmpptr;
    }
    
    im.unlock();
  }
 void CubeMinorCycleAlgorithm::writeBackToFullImage(const String imagename, const Int chanBeg, const Int chanEnd, std::shared_ptr<ImageInterface<Float> > subimptr){
    PagedImage<Float> im(imagename, TableLock::UserLocking);
    im.lock(FileLocker::Write, 30);
    SubImage<Float> *tmpptr=nullptr; 
    tmpptr=SpectralImageUtil::getChannel(im, chanBeg, chanEnd, true);
    tmpptr->copyData(*(subimptr));
                 
    im.unlock();
    delete tmpptr;
                 
  }
void CubeMinorCycleAlgorithm::reset(){
		
  iterBotRec_p=Record();
  modelName_p="";
  residualName_p="";
  psfName_p="";
  maskName_p="";
  pbName_p="";
  posMaskName_p="";
  chanRange_p.resize();
  returnRec_p=Record();
  beamsetRec_p=Record();
  //psfSidelobeLevel_p;
  autoMaskOn_p=False;
  chanFlag_p.resize();
  statsRec_p=Record();
  status_p=False;
                
	
	
}	
	
	
} //# NAMESPACE CASA - END
