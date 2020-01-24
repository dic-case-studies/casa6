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

#ifndef CONVERTERINTENSITY_H_
#define CONVERTERINTENSITY_H_

#include <QString>
#include <QList>
#include <casa/Arrays/Vector.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>

namespace casa {

	/**
	 * Converts intensity units:  Jy/Beam, Kelvin, Fraction of Peak, etc.
	 */

	class ConverterIntensity {
	public:
		static const QString FRACTION_OF_PEAK;
		static const QString KELVIN;
		static const QString JY_SR;
		static const QString JY_ARCSEC;
		static const QString JY_BEAM;
		static const QString JY;
		static const QString ADU;
		static const QString TIMES_PIXELS;
		static bool isSupportedUnits( const QString& yUnit );
		//Hertz values are needed corresponding to the values for Jy/Beam Kelvin conversions
		//only.  Both oldUnits and newUnits refer to the old and new units of the values
		//array.  In order to do FRACTION_OF_PEAK conversions, a maximum value with
		//corresponding maximum units must be passed in.
		static void convert( casacore::Vector<float>& values, const casacore::Vector<float> hertzValues,
		                     const QString& oldUnits, const QString& newUnits,
		                     double maxValue, const QString& maxUnits,
		                     double beamAngle, double beamArea, casacore::SpectralCoordinate coord );

		//Converts between Jy/Beam units.  For example, MJy/Beam <-> Jy/Beam
		static double convertJyBeams( const QString& sourceUnits, const QString& destUnits,
				double value, casacore::SpectralCoordinate& coord );
		static double convertJY( const QString& oldUnits, const QString& newUnits,
				double value, casacore::SpectralCoordinate& coord );
		static double convertJYSR( const QString& oldUnits,const QString& newUnits,
				double value, casacore::SpectralCoordinate& coord );
		static double convertKelvin( const QString& oldUnits,const QString& newUnits,
				double value, casacore::SpectralCoordinate& coord );
		virtual ~ConverterIntensity();

	private:
		ConverterIntensity();
		static double percentToValue( double yValue, double maxValue );
		static double valueToPercent( double yValue, double maxValue );
		static double convertQuantity( double yValue, double frequencyValue,
		                               const QString& oldUnits, const QString& newUnits,
		                               double beamSolidAngle, double beamArea );
		static void convertJansky( casacore::Vector<float>& values, const QString& oldUnits,
		                           const QString& newUnits, casacore::SpectralCoordinate& coord );
		static void convertKelvin( casacore::Vector<float>& values, const QString& oldUnits,
		                           const QString& newUnits, casacore::SpectralCoordinate& coord );
		static bool isJansky( const QString& units );
		static bool isKelvin( const QString& units );
		static double convertNonKelvinUnits( double yValue, const QString& oldUnits,
				const QString& newUnits, double beamArea );
		static QString getJanskyBaseUnits( const QString& units );
		static QString getKelvinBaseUnits( const QString& units );
		static QString stripPixels( const QString& units );
		static double beamToArcseconds( double yValue, double beamArea );
		static double arcsecondsToBeam( double yValue, double beamArea );
		static double srToArcseconds( double yValue );
		static double arcsecondsToSr( double yValue );
		static const QList<QString> BEAM_UNITS;
		static const QList<QString> JY_UNITS;
		static const QList<QString> JY_SR_UNITS;
		static const QList<QString> KELVIN_UNITS;
		static const double SPEED_LIGHT_FACTOR;
		static const double FREQUENCY_FACTOR;
		static const double ARCSECONDS_PER_STERADIAN;
	};

} /* namespace casa */
#endif /* CONVERTERINTENSITY_H_ */
