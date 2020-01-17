//# CubeMajorCycleAlgorithm.h: class to grid and degrid (and write model vis when necessary) in parallel/serial 
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


#ifndef SYNTHESIS_CUBEMAJORCYCLEALGORITHM_H
#define SYNTHESIS_CUBEMAJORCYCLEALGORITHM_H
#include <synthesis/Parallel/Algorithm.h>
#include <synthesis/Parallel/Applicator.h>
#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>
/*ControlRecord_p : necessary keys
 * dividebyweight : True or False
 * weightnames : Vector of Empty string or names of sensitivity images
 * Optional keys:
 * nmajorcycles : number printed in the logger 
 * weightdensity: weightdensity image for non per chan briggs weighting
 * modelnames: Vector of string for possible model images to use to grid residual visibilities 
 * 
 */
namespace casa { //# NAMESPACE CASA - BEGIN

	class SIImageStore;

	class CubeMajorCycleAlgorithm : public Algorithm {
	public:
		//Constructor/desctructor
		CubeMajorCycleAlgorithm();
		~CubeMajorCycleAlgorithm();
		//Functions that needs to be overloaded
		void get();

		// Return the results to the controller
		void put();

		// Return the name of the algorithm
		casacore::String &name();
	private:
		void task();
		void reset();
		casacore::CountedPtr<SIImageStore> subImageStore(const int whichImageId=0);
		casacore::CountedPtr<SIImageStore> multiTermImageStore(const casacore::Int imId);
		casacore::String myName_p;
		casacore::Vector<SynthesisParamsSelect> dataSel_p;
		casacore::Vector<SynthesisParamsImage> imSel_p;
		casacore::Vector<SynthesisParamsGrid> gridSel_p;
		casacore::Vector<casacore::Record> ftmRec_p;
		casacore::Vector<casacore::Record> iftmRec_p;
		casacore::Vector<casacore::Int> polRep_p;
		casacore::Vector<casacore::Int> chanRange_p;
		casacore::Bool dopsf_p;
		casacore::Record controlRecord_p;
		casacore::Record weightParams_p;
		casacore::Vector<casacore::Vector<casacore::String> > startmodel_p;
          casacore::String movingSource_p;
		casacore::Bool status_p;
		casacore::Int serialBug_p; //have to send a private variable in serial case
	};
	
	
} //# NAMESPACE CASA - END
#endif //SYNTHESIS_CUBEMAJORCYCLEALGORITHM_H


