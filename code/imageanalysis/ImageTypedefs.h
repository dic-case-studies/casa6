//# ImageTypedefs.h
//# Copyright (C) 1998,1999,2000,2001
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

#ifndef IMAGETYPEDEFS_H_
#define IMAGETYPEDEFS_H_

#include <casa/aipstype.h>
#include <casa/BasicSL/Complexfwd.h>
#include <casa/Utilities/CountedPtr.h>

#define SPIIT SHARED_PTR<casacore::ImageInterface<T>>
#define SPCIIT SHARED_PTR<const casacore::ImageInterface<T>>

#define SPIIU SHARED_PTR<casacore::ImageInterface<U>>
#define SPCIIU SHARED_PTR<const casacore::ImageInterface<U>>

#define SPIICT SHARED_PTR<casacore::ImageInterface<ComplexType>>
#define SPIIRT SHARED_PTR<casacore::ImageInterface<RealType>>

namespace casacore{

	template<class T> class ImageInterface;
}

namespace casa {

    using SPCIIF = SHARED_PTR<const casacore::ImageInterface<casacore::Float> >;
	using SPIIF = SHARED_PTR<casacore::ImageInterface<casacore::Float> >;
	using SPCIIC = SHARED_PTR<const casacore::ImageInterface<casacore::Complex> >;
	using SPIIC = SHARED_PTR<casacore::ImageInterface<casacore::Complex> >;
	using SPCIID = SHARED_PTR<const casacore::ImageInterface<casacore::Double> >;
	using SPIID = SHARED_PTR<casacore::ImageInterface<casacore::Double> >;
	using SPCIIDC = SHARED_PTR<const casacore::ImageInterface<casacore::DComplex> >;
	using SPIIDC = SHARED_PTR<casacore::ImageInterface<casacore::DComplex> >;
	using ITUPLE = std::tuple<SPIIF, SPIIC, SPIID, SPIIDC>;
}

#endif
