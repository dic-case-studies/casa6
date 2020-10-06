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
#ifndef GAUSSIANESTIMATEDIALOG_QO_H_1
#define GAUSSIANESTIMATEDIALOG_QO_H_1

#include <QDialog>
#include <ui/ui_GaussianEstimateDialog.h>
#include <display/QtPlotter/SearchMoleculesDialog.qo.h>
#include <display/QtPlotter/SpecFitGaussian.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>

class QScrollArea;

namespace casa {

	class GaussianEstimateWidget;
	class Converter;

	class GaussianEstimateDialog : public QDialog {
		Q_OBJECT

	public:
		GaussianEstimateDialog(QWidget *parent = 0);
		~GaussianEstimateDialog();
		void setGaussCount( int count );
		void setCurveData( const casacore::Vector<float>& xValues, const casacore::Vector<float>& yValues);
		void setRangeX( float min, float max );
		void setRangeY( float min, float max );
		void setSpecFitUnits( const QString& specUnits );
		void setDisplayYUnits( const QString& units );
		QString getDisplayYUnits() const;
		QString getUnits() const;
		void setCurveColor( QColor color );
		SpecFitGaussian getEstimate(int index );
		void setEstimates( QList<SpecFitGaussian>& estimates );
		void setSpectralCoordinate( casacore::SpectralCoordinate& coord );

	private slots:
		void unitsChanged( int index );
		void plotsCoordinatedChanged( int );
		void coordinatedValuesChanged( float );

	private:
		void resetEstimates();
		void clearPlots( int removeCount);
		void setCurveData();
		void setDisplayYUnits();
		void setCurveColor();
		void initializeLimits( const casacore::Vector<float>& values, float* const min, float* const max );
		casacore::Vector<float> translateDataUnits( const casacore::Vector<float>& xValues, Converter* converter );
		casacore::Vector<float> xVals;
		casacore::Vector<float> yVals;
		QColor curveColor;
		QString specUnitStr;
		QString displayYUnits;
		QString unitStr;
		QList<GaussianEstimateWidget*> plots;
		Ui::GaussianEstimateDialogClass ui;
		bool plotsCoordinated;

		QWidget* plotHolderWidget;
		casacore::SpectralCoordinate spectralCoordinate;
	};
}
#endif // GAUSSIANESTIMATEDIALOG_QO_H
