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

#include "ConverterFrequency.h"
#include <QDebug>
using namespace casacore;
namespace casa {

	ConverterFrequency::ConverterFrequency(const QString& oldUnits, const QString& newUnits) :
		Converter( oldUnits, newUnits ) {
	}

	double ConverterFrequency::toPixel( double value, SpectralCoordinate spectralCoordinate) {
		Double pixelValue;
		Vector<String> spectralUnits = spectralCoordinate.worldAxisUnits();
		String spectralUnit = spectralUnits[0];
		QString spectralUnitStr( spectralUnit.c_str());
		if ( spectralUnitStr != oldUnits ) {
			Vector<double> freqValues(1);
			freqValues[0] = value;
			convertFrequency( freqValues, oldUnits, spectralUnitStr, spectralCoordinate );
			value = freqValues[0];
		}
		spectralCoordinate.toPixel( pixelValue, value );
		return pixelValue;
	}

	Vector<double> ConverterFrequency::convert( const Vector<double>& oldValues, SpectralCoordinate spectralCoordinate) {
		Vector<double> resultValues( oldValues.size() );
		resultValues = oldValues;
		convertFrequency( resultValues, oldUnits, newUnits, spectralCoordinate );
		return resultValues;
	}

	void ConverterFrequency::convertFrequency( Vector<double> &resultValues,
	        QString& frequencySourceUnits, QString& frequencyDestUnits, SpectralCoordinate& coord ) {
		int sourceIndex = Converter::FREQUENCY_UNITS.indexOf( frequencySourceUnits );
		int destIndex = Converter::FREQUENCY_UNITS.indexOf( frequencyDestUnits );
		Converter::convert( resultValues, sourceIndex, destIndex, coord );
	}

	ConverterFrequency::~ConverterFrequency() {
		// TODO Auto-generated destructor stub
	}

using namespace casacore;
} /* namespace casa */
