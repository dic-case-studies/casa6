//# CubeMajorCycleAlgorithm.cc: implementation of class to grid and degrid (and write model vis when necessary) in parallel/serial 
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
#include <synthesis/ImagerObjects/CubeMajorCycleAlgorithm.h>
#include <synthesis/ImagerObjects/SynthesisImagerVi2.h>
#include <casa/Containers/Record.h>
#include <synthesis/ImagerObjects/SimpleSIImageStore.h>
#include <imageanalysis/Utilities/SpectralImageUtil.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN
extern Applicator applicator;

CubeMajorCycleAlgorithm::CubeMajorCycleAlgorithm() : myName_p("CubeMajorCycleAlgorithm"),  ftmRec_p(0), iftmRec_p(0), polRep_p(0),startmodel_p(0), status_p(False){
	
}
CubeMajorCycleAlgorithm::~CubeMajorCycleAlgorithm() {
	
}
	
void CubeMajorCycleAlgorithm::get() {
	reset();
	//cerr << "in get for child process " << applicator.isWorker() << endl;
	Record vecImParsRec;
	Record vecSelParsRec;
	Record vecGridParsRec;
	// get data sel params #1
	applicator.get(vecSelParsRec);
	// get image sel params #2
	applicator.get(vecImParsRec);
	// get gridders params #3
	applicator.get(vecGridParsRec);
	// get which channel to process #4
	applicator.get(chanRange_p);
	cerr <<"GET chanRange " << chanRange_p << endl;
	// psf or residual CubeMajorCycleAlgorithm #5
	applicator.get(dopsf_p);
	// store modelvis and other controls #6
	applicator.get(controlRecord_p);
	// weight params
	applicator.get(weightParams_p);
	//Somewhere before this we have to make sure that vecSelPars has more than 0 fields)
	dataSel_p.resize(vecSelParsRec.nfields());
	/// Fill the private variables
	for (Int k=0; k < Int(dataSel_p.nelements()); ++k){
		(dataSel_p[k]).fromRecord(vecSelParsRec.asRecord(String::toString(k)));
	}
	//imsel and gridsel should be the same numbers (number of image fields)
	Int nmajorcycles=0;
	if(controlRecord_p.isDefined("nmajorcycles"))
		controlRecord_p.get("nmajorcycles",nmajorcycles);
	imSel_p.resize(vecImParsRec.nfields());
	gridSel_p.resize(vecImParsRec.nfields());
	ftmRec_p.resize(vecImParsRec.nfields());
	iftmRec_p.resize(vecImParsRec.nfields());
	polRep_p.resize(vecImParsRec.nfields());
	polRep_p.set(-1);
	startmodel_p.resize(vecImParsRec.nfields());
	startmodel_p.set(Vector<String>(0));
	for (uInt k=0; k < imSel_p.nelements(); ++k){
		Record imSelRec=vecImParsRec.asRecord(String::toString(k));
		if(imSelRec.isDefined("polrep"))
			imSelRec.get("polrep", polRep_p[k]);
		//Only first major cycle we need to reset model
		if(nmajorcycles==1)
			imSelRec.get("startmodel", startmodel_p[k]);
		//have to do that as fromRecord check does not like to have both model and startmodel on disk !
		imSelRec.define("startmodel", Vector<String>(0));
		(imSel_p[k]).fromRecord(imSelRec);
		Record recGridSel=vecGridParsRec.asRecord(String::toString(k));
		(gridSel_p[k]).fromRecord(recGridSel);
		if(!recGridSel.isDefined("ftmachine")){
			ftmRec_p.resize();
			iftmRec_p.resize();
		}
		if(ftmRec_p.nelements() >0){
			ftmRec_p[k]=recGridSel.asRecord("ftmachine");
			iftmRec_p[k]=recGridSel.asRecord("iftmachine");
		}
			
	}
	
	
}
void CubeMajorCycleAlgorithm::put() {
	
	if(applicator.isSerial()){
		serialBug_p=Applicator::DONE;
		//applicator.put(serialBug_p);
		
	}
	//cerr << "in put " << status_p << endl;
	applicator.put(status_p);	
	
}
	
void CubeMajorCycleAlgorithm::task(){
	status_p = False;
	//SynthesisImagerVi2 imgr;
	//imgr.selectData(dataSel_p);
	// We do not use chanchunking in this model
	for (uInt k=0; k < gridSel_p.nelements(); ++k)
		gridSel_p[k].chanchunks = 1;

	//imgr.defineImage(imSel_p,gridSel_p);
	// need to find how many subfields/outliers have been set
	//CountedPtr<SIImageStore> imstor =imgr.imageStore(0);
	//CountedPtr<ImageInterface<Float> > resid=imstor->residual();
	//Int nchan = resid->shape()(3);
	//std::shared_ptr<SIImageStore> subImStor=imstor->getSubImageStore(0, 1, chanId_p, nchan, 0,1);
	
	SynthesisImagerVi2 subImgr;
	for (Int k=0; k < Int(dataSel_p.nelements()); ++k){
		subImgr.selectData(dataSel_p[k]);
	}
	Vector<CountedPtr<SIImageStore> > subImStor(imSel_p.nelements());
	//Do multifield in one process only for now 
	//TODO if all fields have same nchan then partition all on all subcalls
	if(chanRange_p[0]==0){
		for (uInt k=0; k < imSel_p.nelements(); ++k){
			subImStor[k]=subImageStore(k);
			if(ftmRec_p.nelements()>0){
				subImgr.defineImage(subImStor[k], ftmRec_p[k], iftmRec_p[k]);	
			}else{
				subImgr.defineImage(subImStor[k],  gridSel_p[k].ftmachine);
			}
		}
	}else{
		subImStor.resize(1);
		subImStor[0]=subImageStore(0);
		if(ftmRec_p.nelements()>0){
				subImgr.defineImage(subImStor[0], ftmRec_p[0], iftmRec_p[0]);
		}else{
			subImgr.defineImage(subImStor[0],  gridSel_p[0].ftmachine);
		}
	}
	subImgr.setCubeGridding(False);
	// TO DO get weight param and set weight
	if(!weightParams_p.isDefined("type") || weightParams_p.asString("type")=="natural")
		subImgr.weight("natural");
	else
		subImgr.weight(weightParams_p);
	if (!dopsf_p){
		subImgr.executeMajorCycle(controlRecord_p);
		if(controlRecord_p.isDefined("dividebyweight") && controlRecord_p.asBool("dividebyweight"))
		for(uInt k=0; k < subImStor.nelements(); ++k){
			{
			LatticeLocker lock1 (*(subImStor[k]->residual()), FileLocker::Write);
			LatticeLocker lock2 (*(subImStor[k]->sumwt()), FileLocker::Read);
			subImStor[k]->divideResidualByWeight();
			subImStor[k]->residual()->flush();
			}
			subImStor[k]->residual()->unlock();
			if(subImStor[k]->hasModel())
				subImStor[k]->model()->unlock();
		}
	}
	else{
		subImgr.makePSF();
		for(uInt k=0; k < subImStor.nelements(); ++k){
		if(controlRecord_p.isDefined("dividebyweight") && controlRecord_p.asBool("dividebyweight"))
		{
			LatticeLocker lock1 (*(subImStor[k]->psf()), FileLocker::Write);
			LatticeLocker lock2 (*(subImStor[k]->sumwt()), FileLocker::Read);
			subImStor[k]->dividePSFByWeight();
			subImStor[k]->psf()->flush();
		}
		subImStor[k]->psf()->unlock();
		}
	}
	
	status_p = True;
}
String&	CubeMajorCycleAlgorithm::name(){
	return myName_p;
}

CountedPtr<SIImageStore> CubeMajorCycleAlgorithm::subImageStore(const int imId){
	//For some reason multiterm deconvolver is allowed with cubes !
	String isMTdeconv="";
	if(imId==0 && imSel_p[imId].deconvolver=="mtmfs") isMTdeconv=".tt0";
	String residname=imSel_p[imId].imageName+".residual"+isMTdeconv;
	String psfname=imSel_p[imId].imageName+".psf"+isMTdeconv;
	String sumwgtname=imSel_p[imId].imageName+".sumwt"+isMTdeconv;
	if(imId > 0 && !Table::isReadable(sumwgtname)){
		return multiTermImageStore(imId);
	}
	shared_ptr<ImageInterface<Float> >subpsf=nullptr;
	shared_ptr<ImageInterface<Float> >subresid=nullptr;
	shared_ptr<ImageInterface<Float> >submodel=nullptr;
	shared_ptr<ImageInterface<Float> > subweight=nullptr;
	cerr << imId << " sumwt name " << sumwgtname << endl;
	if(!Table::isReadable(sumwgtname))
		throw(AipsError("Programmer error: sumwt disk image is non existant")); 
	PagedImage<Float> sumwt(sumwgtname, TableLock::UserNoReadLocking);
	//Should be partitioning for main image only
	//chanRange
	Int chanBeg=0;
	Int chanEnd=0;
	if(imId==0){
		chanBeg=chanRange_p[0];
		chanEnd=chanRange_p[1];
	}
	else{
		chanBeg=0;
		chanEnd=sumwt.shape()[3]-1;
	}
	
	cerr << "chanBeg " << chanBeg << " chanEnd " << chanEnd << " imId " << imId << endl;
	if(dopsf_p){
		PagedImage<Float> psf(psfname, TableLock::UserNoReadLocking);
		subpsf.reset(SpectralImageUtil::getChannel(psf, chanBeg, chanEnd, true));
	}
	else{
		//need to loop over all fields somewhere
		PagedImage<Float> resid(residname, TableLock::UserNoReadLocking);
		subresid.reset(SpectralImageUtil::getChannel(resid, chanBeg, chanEnd, true));
		//String modelname=imSel_p.imageName+".model";
		if(controlRecord_p.isDefined("modelnames")){
			Vector<String> modelnames(controlRecord_p.asArrayString("modelnames"));
			if(imId >= int(modelnames.nelements()))
				throw(AipsError("Number of model images does not match number of image fields defined"));
			if(Table::isReadable(modelnames[imId])){
				PagedImage<Float> model(modelnames[imId], TableLock::UserNoReadLocking);
				//Darn has to lock it as writable because overlap in SIMapperCollection code 
				//wants that...though we are not really modifying it here
				Bool writeisneeded=(imSel_p.nelements()!=1 || startmodel_p[imId].nelements() >0);
				submodel.reset(SpectralImageUtil::getChannel(model, chanBeg, chanEnd, writeisneeded));
			}
			
		}
		
	}
	
	Vector<String> weightnames(controlRecord_p.asArrayString("weightnames"));
	if(imId >= int(weightnames.nelements()))
		throw(AipsError("Number of model images does not match number of image fields defined"));
	if(Table::isReadable(weightnames[imId])){
		PagedImage<Float> weight(weightnames[imId], TableLock::UserNoReadLocking);
		subweight.reset(SpectralImageUtil::getChannel(weight, chanBeg, chanEnd, true));
	}
	
	shared_ptr<ImageInterface<Float> >subsumwt(SpectralImageUtil::getChannel(sumwt, chanBeg, chanEnd, true));
	bool useweightimage=(subweight) ? true : false;
	CountedPtr<SIImageStore> subimstor=new SimpleSIImageStore(submodel, subresid, subpsf, subweight, nullptr, nullptr, subsumwt, nullptr, nullptr, nullptr, useweightimage);
	if(polRep_p[imId]< 0)
		throw(AipsError("data polarization type is not defined"));
	StokesImageUtil::PolRep polrep=(StokesImageUtil::PolRep)polRep_p[imId];
	subimstor->setDataPolFrame(polrep);
	if(startmodel_p[imId].nelements() >0){
		LatticeLocker lock1 (*(subimstor->model()), FileLocker::Write);
		subimstor->setModelImage(startmodel_p[imId]);	
	}
	//cerr << "subimagestor TYPE" << subimstor->getType() << endl;
	return subimstor;
}

CountedPtr<SIImageStore> CubeMajorCycleAlgorithm::multiTermImageStore(const Int imId){
	uInt nterms=0;
	String sumwgtname=imSel_p[imId].imageName+".sumwt.tt"+String::toString(nterms);
	while (Table::isReadable(sumwgtname)){
		++nterms;
		sumwgtname=imSel_p[imId].imageName+".sumwt.tt"+String::toString(nterms);
	}
	if(nterms==0){
		throw(AipsError("outlier "+String::toString(imId)+" field weight image is not defined"));
	}
	nterms=(nterms+1)/2;
	CountedPtr<SIImageStore> subimstor=new SIImageStoreMultiTerm(imSel_p[imId].imageName, nterms, True);
	if(polRep_p[imId]< 0)
		throw(AipsError("data polarization type is not defined"));
	StokesImageUtil::PolRep polrep=(StokesImageUtil::PolRep)polRep_p[imId];
	subimstor->setDataPolFrame(polrep);
	return subimstor;
}	
void CubeMajorCycleAlgorithm::reset(){
		
		dataSel_p.resize();
		imSel_p.resize();
		gridSel_p.resize();
		ftmRec_p.resize();
		iftmRec_p.resize();
		polRep_p.resize();
		chanRange_p.resize();
		dopsf_p=False;
		controlRecord_p=Record();
		weightParams_p=Record();
		startmodel_p.resize();
		status_p=False;
	
	
}
	
	
	
	
	
	
} //# NAMESPACE CASA - END
