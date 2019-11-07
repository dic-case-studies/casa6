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
		casacore::CountedPtr<SIImageStore> subImageStore();
		casacore::String myName_p;
		casacore::Vector<SynthesisParamsSelect> dataSel_p;
		SynthesisParamsImage imSel_p;
		SynthesisParamsGrid gridSel_p;
		casacore::Int chanId_p;
		casacore::Bool dopsf_p;
		casacore::Record controlRecord_p;
		
		casacore::Bool status_p;
		casacore::Int serialBug_p; //have to send a private variable in serial case
	};
	
	
} //# NAMESPACE CASA - END
#endif //SYNTHESIS_CUBEMAJORCYCLEALGORITHM_H


