//# Copyright (C) 2005
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
#ifndef CONVERTERCHANNEL_H_
#define CONVERTERCHANNEL_H_

#include <display/QtPlotter/conversion/Converter.h>

namespace casa {

//Purpose of this class is to convert from channels= pixels to the world units/pixels
//used by the spectral axis.
	class ConverterChannel : public Converter {
	public:
		ConverterChannel( const QString& oldUnits, const QString& newUnits );

		virtual double toPixel( double value, casacore::SpectralCoordinate spectralCoordinate );
		virtual casacore::Vector<double> convert( const casacore::Vector<double>& oldValues, casacore::SpectralCoordinate spectralCoordinate );
		virtual double convert ( double oldValue, casacore::SpectralCoordinate spectralCoordinate );
		virtual ~ConverterChannel();
	};

} /* namespace casa */
#endif /* CONVERTERCHANNEL_H_ */
