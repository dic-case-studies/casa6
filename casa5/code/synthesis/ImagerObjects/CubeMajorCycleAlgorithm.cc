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

CubeMajorCycleAlgorithm::CubeMajorCycleAlgorithm() : myName_p("CubeMajorCycleAlgorithm"),  status_p(False){
	
}
CubeMajorCycleAlgorithm::~CubeMajorCycleAlgorithm() {
	
}
	
void CubeMajorCycleAlgorithm::get() {
	//cerr << "in get for child process " << applicator.isWorker() << endl;
	Record imparsRec;
	Record vecSelParsRec;
	Record gridparsRec;
	// get data sel params #1
	applicator.get(vecSelParsRec);
	// get image sel params #2
	applicator.get(imparsRec);
	// get gridders params #3
	applicator.get(gridparsRec);
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
	imSel_p.fromRecord(imparsRec);
	gridSel_p.fromRecord(gridparsRec);
	
	
}
void CubeMajorCycleAlgorithm::put() {
	
	if(applicator.isSerial()){
		serialBug_p=Applicator::DONE;
		applicator.put(serialBug_p);
		
	}
	//cerr << "in put " << status_p << endl;
	applicator.put(status_p);	
	
}
	
void CubeMajorCycleAlgorithm::task(){
	status_p = False;
	//SynthesisImagerVi2 imgr;
	//imgr.selectData(dataSel_p);
	// We do not use chanchunking in this model
	gridSel_p.chanchunks = 1;
	//imgr.defineImage(imSel_p,gridSel_p);
	// need to find how many subfields/outliers have been set
	//CountedPtr<SIImageStore> imstor =imgr.imageStore(0);
	//CountedPtr<ImageInterface<Float> > resid=imstor->residual();
	//Int nchan = resid->shape()(3);
	//std::shared_ptr<SIImageStore> subImStor=imstor->getSubImageStore(0, 1, chanId_p, nchan, 0,1);
	CountedPtr<SIImageStore> subImStor=subImageStore();
	SynthesisImagerVi2 subImgr;
	for (Int k=0; k < Int(dataSel_p.nelements()); ++k){
		subImgr.selectData(dataSel_p[k]);
	}
	subImgr.defineImage(subImStor,  gridSel_p.ftmachine);
	subImgr.setCubeGridding(False);
	// TO DO get weight param and set weight
	if(!weightParams_p.isDefined("type") || weightParams_p.asString("type")=="natural")
		subImgr.weight("natural");
	else
		subImgr.weight(weightParams_p);
	if (!dopsf_p){
		subImgr.executeMajorCycle(controlRecord_p);
		if(controlRecord_p.isDefined("dividebyweight") && controlRecord_p.asBool("dividebyweight"))
		{
			LatticeLocker lock1 (*(subImStor->residual()), FileLocker::Write);
			LatticeLocker lock2 (*(subImStor->sumwt()), FileLocker::Read);
			subImStor->divideResidualByWeight();
			subImStor->residual()->flush();
		}
		subImStor->residual()->unlock();
		if(subImStor->hasModel())
			subImStor->model()->unlock();
	}
	else{
		subImgr.makePSF();
		if(controlRecord_p.isDefined("dividebyweight") && controlRecord_p.asBool("dividebyweight"))
		{
			LatticeLocker lock1 (*(subImStor->psf()), FileLocker::Write);
			LatticeLocker lock2 (*(subImStor->sumwt()), FileLocker::Read);
			subImStor->dividePSFByWeight();
			subImStor->psf()->flush();
		}
		subImStor->psf()->unlock();
	}
	
	status_p = True;
}
String&	CubeMajorCycleAlgorithm::name(){
	return myName_p;
}

CountedPtr<SIImageStore> CubeMajorCycleAlgorithm::subImageStore(){
	String residname=imSel_p.imageName+".residual";
	String psfname=imSel_p.imageName+".psf";
	String sumwgtname=imSel_p.imageName+".sumwt";
	shared_ptr<ImageInterface<Float> >subpsf=nullptr;
	shared_ptr<ImageInterface<Float> >subresid=nullptr;
	shared_ptr<ImageInterface<Float> >submodel=nullptr;
	if(dopsf_p){
		PagedImage<Float> psf(psfname, TableLock::UserNoReadLocking);
		subpsf.reset(SpectralImageUtil::getChannel(psf, chanRange_p[0], chanRange_p[1], true));
	}
	else{
		//need to loop over all fields somewhere
		PagedImage<Float> resid(residname, TableLock::UserNoReadLocking);
		subresid.reset(SpectralImageUtil::getChannel(resid, chanRange_p[0], chanRange_p[1], true));
		//String modelname=imSel_p.imageName+".model";
		if(controlRecord_p.isDefined("modelnames")){
			Vector<String> modelnames(controlRecord_p.asArrayString("modelnames"));
			if(Table::isReadable(modelnames[0])){
				PagedImage<Float> model(modelnames[0], TableLock::UserNoReadLocking);
				submodel.reset(SpectralImageUtil::getChannel(model, chanRange_p[0], chanRange_p[1], false));
			}
		}
	}
	
	
	PagedImage<Float> sumwt(sumwgtname, TableLock::UserNoReadLocking);
	shared_ptr<ImageInterface<Float> >subsumwt(SpectralImageUtil::getChannel(sumwt, chanRange_p[0], chanRange_p[1], true));
	CountedPtr<SIImageStore> subimstor=new SimpleSIImageStore(submodel, subresid, subpsf, nullptr, nullptr, nullptr, subsumwt, nullptr, nullptr, nullptr);
	//cerr << "subimagestor TYPE" << subimstor->getType() << endl;
	return subimstor;
}

	
	
	
	
	
	
	
} //# NAMESPACE CASA - END
