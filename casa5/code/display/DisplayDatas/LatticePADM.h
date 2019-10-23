//# LatticePADisplayMethod.h: base for drawing axis-bound lattice elements
//# Copyright (C) 1996,1997,1998,1999,2000,2001
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

#ifndef TRIALDISPLAY_LATTICEPADM_H
#define TRIALDISPLAY_LATTICEPADM_H

//# aips includes:
#include <casa/aips.h>

//# display library includes:
#include <display/DisplayDatas/PrincipalAxesDM.h>

namespace casacore{

	template <class T> class Array;
	template <class T> class MaskedLattice;
	class IPosition;
}

namespace casa { //# NAMESPACE CASA - BEGIN

//# forwards:
	template <class T> class LatticePADisplayData;

// <summary>
// Partial implementation of PrincipalAxesDM for casacore::Lattice-based data.
// </summary>
//
// <synopsis>
// This class is a partial (ie. base) implementation of PrincipalAxesDM
// which adds methods particular to handling casacore::Lattice-based data.
// </synopsis>

	template <class T> class LatticePADisplayMethod : public PrincipalAxesDM {

	public:

		// Constructor
		// do I need the default constructor?
		LatticePADisplayMethod();
		LatticePADisplayMethod(const casacore::uInt xAxis, const casacore::uInt yAxis,
		                       const casacore::uInt mAxis, const casacore::IPosition fixedPos,
		                       LatticePADisplayData<T> *arDat);
		// 2d version
		LatticePADisplayMethod(const casacore::uInt xAxis, const casacore::uInt yAxis,
		                       LatticePADisplayData<T> *arDat);

		// Destructor
		virtual ~LatticePADisplayMethod();
		// Extract data from the lattice: used by draw() in PrincipalAxesDM
		// this is probably not needed in this class...
		virtual casacore::Bool dataGetSlice(casacore::Matrix<T>& datMatrix,
				casacore::Matrix<casacore::Bool>& mask,
				const casacore::IPosition& start,
				const casacore::IPosition& sliceShape,
				const casacore::IPosition& stride);
	protected:

		// Query the shape of the lattice: used by draw() in PrincipalAxesDM
		virtual casacore::IPosition dataShape();



		virtual casacore::Bool dataGetSlice(casacore::Matrix<T>& datMatrix,
		                          casacore::Matrix<casacore::Bool>& mask,
		                          const casacore::IPosition& start,
		                          const casacore::IPosition& sliceShape,
		                          const casacore::IPosition& stride,
		                          casacore::MaskedLattice<T>& latt);

	private:

	};


} //# NAMESPACE CASA - END

#ifndef AIPS_NO_TEMPLATE_SRC
#include <display/DisplayDatas/LatticePADM.tcc>
#endif //# AIPS_NO_TEMPLATE_SRC
#endif
