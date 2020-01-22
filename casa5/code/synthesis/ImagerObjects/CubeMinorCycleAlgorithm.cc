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

CubeMinorCycleAlgorithm::CubeMinorCycleAlgorithm() : myName_p("CubeMinorCycleAlgorithm"), status_p(False){
	
}
CubeMinorCycleAlgorithm::~CubeMinorCycleAlgorithm() {
	
}
	
void CubeMinorCycleAlgorithm::get() {

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
	cerr <<"GET chanRange " << chanRange_p << endl;
	decPars_p.fromRecord(decParsRec);
	
	
	
}
void CubeMinorCycleAlgorithm::put() {
	
  ///# 1  chanrange processed 
  applicator.put(chanRange_p);
	//cerr << "in put " << status_p << endl;
  ///#2 return record of deconvolver
  applicator.put(returnRec_p);	
	
}
	
void CubeMinorCycleAlgorithm::task(){
	status_p = False;
	
	SynthesisDeconvolver subDeconv;
	subDeconv.setupDeconvolution(decPars_p);
	std::shared_ptr<SIImageStore> subimstor=subImageStore();
	subDeconv.initMinorCycle(subimstor);
	returnRec_p=subDeconv.executeMinorCycle(iterBotRec_p);
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
	Int chanBeg=0;
	Int chanEnd=0;
	chanBeg=chanRange_p[0];
	chanEnd=chanRange_p[1];
	//cerr << "chanBeg " << chanBeg << " chanEnd " << chanEnd << " imId " << imId << endl;
        
	
        PagedImage<Float> psf(psfName_p, TableLock::UserNoReadLocking);
        // darn makeBeamSet wants to write in the psf ! So accessing it writable
        subpsf.reset(SpectralImageUtil::getChannel(psf, chanBeg, chanEnd, true));
        
        
	PagedImage<Float> resid(residualName_p, TableLock::UserNoReadLocking);
	subresid.reset(SpectralImageUtil::getChannel(resid, chanBeg, chanEnd, true));
	PagedImage<Float> model(modelName_p, TableLock::UserNoReadLocking);
	submodel.reset(SpectralImageUtil::getChannel(model, chanBeg, chanEnd, true));

	PagedImage<Float> mask(maskName_p, TableLock::UserNoReadLocking);
	submask.reset(SpectralImageUtil::getChannel(mask, chanBeg, chanEnd, true));

	std::shared_ptr<SIImageStore> subimstor(new SimpleSIImageStore(submodel, subresid, subpsf, nullptr, nullptr, submask, nullptr, nullptr, nullptr, nullptr));
	//cerr << "subimagestor TYPE" << subimstor->getType() << endl;
	return subimstor;
}

	
	
	
} //# NAMESPACE CASA - END
