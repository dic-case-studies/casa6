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
#ifndef MOMENTSETTINGSWIDGETRADIO_QO_H
#define MOMENTSETTINGSWIDGETRADIO_QO_H

#include <QWidget>
#include <QMap>
#include <QThread>
#include <QProgressDialog>
#include <casa/Quanta/Quantum.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <display/QtPlotter/ProfileTaskFacilitator.h>
#include <display/QtPlotter/MomentSettingsWidgetRadio.ui.h>
#include <imageanalysis/ImageAnalysis/ImageMomentsProgressMonitor.h>

namespace casacore{

	template <class T> class ImageInterface;
}

namespace casa {

	class MomentCollapseThreadRadio;
	class ThresholdingBinPlotDialog;
	class Converter;
	template <class T> class ImageMoments;

	class CollapseResult {

	public:
		CollapseResult( const casacore::String& outputName, bool tmp, SHARED_PTR<casacore::ImageInterface<float>> img );

		casacore::String getOutputFileName() const {
			return outputFileName;
		}
		bool isTemporaryOutput() const {
			return temporary;
		}
		SHARED_PTR<casacore::ImageInterface<float> > getImage() const {
			return image;
		}

	private:
		casacore::String outputFileName;
		bool temporary;

		SHARED_PTR<casacore::ImageInterface<float> > image;
	};


//Note:  ImageMomentsProgressMonitor is an interface that provides this class
//with updates concerning the progress of the moment calculation task.

	/**
	 * Responsible for running the collapse algorithm in
	 * the background so that we don't freeze the GUI.
	 */
	class MomentCollapseThreadRadio : public QThread, public ImageMomentsProgressMonitor {
		Q_OBJECT
	public:
		MomentCollapseThreadRadio( ImageMoments<float>* imageAnalysis );
		bool isSuccess() const;
		void setChannelStr( casacore::String str );
		void setMomentNames( const casacore::Vector<QString>& momentNames );
		void setOutputFileName( QString name );
		casacore::String getErrorMessage() const;
		std::vector<CollapseResult> getResults() const;
		void setData(const casacore::Vector<int>& moments, const int axis,
		             const casacore::Vector<casacore::String>& method,
		             const casacore::Vector<int>& smoothaxes,
		             const casacore::Vector<casacore::String>& smoothtypes,
		             const casacore::Vector<casacore::Quantity>& smoothwidths,
		             const casacore::Vector<float>& includepix,
		             const casacore::Vector<float>& excludepix,
		             const double peaksnr, const double stddev,
		             const casacore::String& doppler = "RADIO", const casacore::String& baseName = "");
		void run();
		void halt();
		//Methods from the ImageMomentsProgressMonitor interface
		void setStepCount( int count );
		void setStepsCompleted( int count );
		void done();
		~MomentCollapseThreadRadio();

	signals:
		void stepCountChanged( int count );
		void stepsCompletedChanged( int count );

	private:
		bool getOutputFileName( casacore::String& outName, int moment, const casacore::String& channelStr ) const;

		ImageMoments<float>* analysis;
		casacore::Vector<int> moments;
		casacore::Vector<QString> momentNames;
		int axis;
		casacore::String channelStr;
		casacore::Vector<casacore::String> method;
		casacore::Vector<int> smoothaxes;
		casacore::Vector<casacore::String> smoothtypes;
		casacore::Vector<casacore::Quantity> smoothwidths;
		casacore::Vector<float> includepix;
		casacore::Vector<float> excludepix;
		double peaksnr;
		double stddev;
		casacore::String doppler;
		casacore::String baseName;
		QString outputFileName;
		int stepSize;
		std::vector<CollapseResult> collapseResults;
		casacore::String errorMsg;
		bool collapseError;
		bool stopImmediately;
	};

//Note: ProfileTaskFacilitator abstracts out some of the common functionality
//needed for calculating moments and spectral line fitting into a single
//base class


	class MomentSettingsWidgetRadio : public QWidget, public ProfileTaskFacilitator {
		Q_OBJECT

	public:
		MomentSettingsWidgetRadio(QWidget *parent = 0);

		void setUnits( QString units );
		void setRange( double min, double max );
		void reset();
		~MomentSettingsWidgetRadio();

	signals:
		void updateProgress(int);
		void momentsFinished();

	private slots:
		void setStepCount( int count );
		void setStepsCompleted( int count );
		void thresholdingChanged();
		void thresholdSpecified();
		void adjustTableRows( int count );
		void collapseImage();
		void setCollapsedImageFile();
		void collapseDone();
		void graphicalThreshold();
		void symmetricThresholdChanged( int checkedState );
		void thresholdTextChanged( const QString& text );
		void stopMoments();

	private:
		void _initAnalysis();
		casacore::Record _makeRegionRecord( );
		enum SummationIndex {MEAN, INTEGRATED, WEIGHTED_MEAN, DISPERSION, MEDIAN,
		                     MEDIAN_VELOCITY, STDDEV,  RMS, ABS_MEAN_DEV, MAX, MAX_VELOCITY, MIN,
		                     MIN_VELOCITY, END_INDEX
		                    };
		QMap<SummationIndex, int> momentMap;
		Ui::MomentSettingsWidgetRadio ui;
		ImageMoments<float>* imageAnalysis;
		MomentCollapseThreadRadio* collapseThread;
		ThresholdingBinPlotDialog* thresholdingBinDialog;
		QString outputFileName;
		QList<QString> momentOptions;
		QProgressDialog progressBar;
		QString m_units;

		//Progress Monitor functionality
		int momentCount;
		int cycleCount;
		int baseIncrement;
		int previousCount;

		void setTableValue(int row, int col, float val );
		void getChannelMinMax( int channelIndex, QString& minStr, QString& maxStr ) const;
		void convertChannelRanges( const QString& oldUnits, const QString& newUnits );
		void convertChannelValue( const QString& channelStr, const QString& channelIdentifier,
		                          Converter* converter, int row, int col, bool toPixels,
		                          casacore::SpectralCoordinate& coord );
		casacore::String makeChannelInterval( float startChannelIndex,float endChannelIndex ) const;
		casacore::Vector<int> populateMoments( casacore::Vector<QString>& momentNames );
		casacore::Vector<casacore::String> populateMethod() const;
		casacore::String populateChannels(casacore::uInt * nSelectedChannels, bool * channelOK);
		bool populateThresholds( casacore::Vector<float>& includeThreshold, casacore::Vector<float>& excludeThreshold );
		bool populateThreshold( casacore::Vector<float>& threshold );
	};

}

#endif // MOMENTSETTINGSWIDGETRADIO_QO_H
