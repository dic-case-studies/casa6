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

#ifndef SEARCHMOLECULESWIDGET_QO_H_1
#define SEARCHMOLECULESWIDGET_QO_H_1

#include <QWidget>
#include <QMap>
#include <QThread>
#include <QProgressDialog>
#include <display/QtPlotter/SearchMoleculesWidget.ui.h>
#include <measures/Measures/MRadialVelocity.h>
#include <measures/Measures/MDoppler.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <display/QtPlotter/conversion/Converter.h>
#include <spectrallines/Splatalogue/Searcher.h>
#include <spectrallines/Splatalogue/SearcherFactory.h>
#include <display/QtPlotter/SearchMoleculesResultDisplayer.h>
#include <QDebug>

class QTimer;

namespace casa {

	class QtCanvas;

	/**
	 * Responsible for running the search algorithm in
	 * the background so that we don't freeze the GUI.
	 */
	class SearchThread : public QThread {
	public:
		SearchThread( Searcher* searcher, int offset ) {
			this->searcher= searcher;
			this->offset = offset;
			countNeeded = true;
		}
		casacore::String getErrorMessage() const {
			return errorMsg;
		}
		casacore::String getErrorMessageCount() const {
			return errorMsgCount;
		}

		void setCountNeeded( bool needed ) {
			countNeeded = needed;
		}

		long getResultsCount() const {
			return searchResultsCount;
		}

		vector<SplatResult> getResults() const {
			return searchResults;
		}

		void stopSearch() {
			searcher->stopSearch();
		}

		void run() {
			if ( offset == 0 && countNeeded ) {
				searchResultsCount = searcher->doSearchCount( errorMsgCount );
			} else {
				searchResultsCount = 1;
			}
			if ( searchResultsCount > 0 ) {
				searchResults = searcher->doSearch( errorMsg, offset );
			}
		}
		~SearchThread() {
		}
	private:
		Searcher* searcher;
		int searchResultsCount;
		int offset;
		bool countNeeded;
		vector<SplatResult> searchResults;
		string errorMsg;
		string errorMsgCount;
	};


	class SearchMoleculesWidget : public QWidget {
		Q_OBJECT

	public:
		SearchMoleculesWidget(QWidget *parent = 0);
		void setCanvas( QtCanvas* drawCanvas );
		QString getUnit() const;
		bool isLocal() const;

		void setRange( double min, double max, QString units );
		void setSpectralCoordinate(casacore::SpectralCoordinate coord );
		void updateReferenceFrame();
		static void setInitialReferenceFrame( QString initialReferenceStr );
		void setResultDisplay( SearchMoleculesResultDisplayer* resultDisplay );
		double getRedShiftedValue( bool reverseRedshift, double value, bool* valid ) const;

		vector<SplatResult> getSearchResults() const;
		casacore::MDoppler::Types getDopplerType() const;
		casacore::MRadialVelocity::Types getReferenceFrame() const;
		casacore::MFrequency::Types getReferenceFrequency() const;
		~SearchMoleculesWidget();
		static const QString SPLATALOGUE_UNITS;
		static const QString SEARCH_DEFAULT_UNITS;

	signals:
		void searchCompleted();
		void redshiftChanged();

	private slots:
		void search();
		void searchUnitsChanged( const QString& searchUnits );
		void redshiftChanged( const QString& redshiftStr );
		void dopplerShiftChanged();
		void dopplerVelocityUnitsChanged();
		void searchFinished();
		void prevResults();
		void nextResults();
		void stopSearch();

	private:

		static QString initialReferenceStr;

		void setAstronomicalFilters( Searcher* searcher );
		void convertRangeLineEdit( QLineEdit* lineEdit, Converter* converter );
		void initializeSearchRange( QLineEdit* lineEdit, double& value, bool* valid );
		vector<string> initializeChemicalNames();
		vector<string> initializeChemicalFormulas();
		QList<QString> getSearchChemicals();
		vector<string> convertStringFormats( const QList<QString>& names );
		double redShiftToVelocity( QString velocityUnits) const;
		double velocityToRedshift( QString velocityUnits ) const;
		void startSearchThread();
		void setSearchRangeDefault();
		double setRangeValue( double value, QString units );
		double getRedShift() const;
		casacore::MDoppler getRedShiftAdjustment( bool reverseRedshift) const;

		enum AstroFilters { NONE, TOP_20, PLANETARY_ATMOSPHERE,HOT_CORES,
		                    DARK_CLOUDS,DIFFUSE_CLOUDS,COMETS, AGB_PPN_PN,EXTRAGALACTIC
		                  };

		Ui::SearchMoleculesWidget ui;

		QString unitStr;
		QString dopplerVelocityUnitStr;
		vector<SplatResult> searchResults;
		QList<QString> velocityUnitsList;
		QMap<QString, casacore::MRadialVelocity::Types> radialVelocityTypeMap;
		QMap<QString, casacore::MDoppler::Types> dopplerTypeMap;
		bool dopplerInVelocity;
		bool searchInterrupted;
		SearchThread* searchThread;
		Searcher* searcher;
		QtCanvas* canvas;
		QProgressDialog progressBar;

		//For conversion
		casacore::SpectralCoordinate coord;

		//Scrolling support
		int searchResultCount;
		int searchResultOffset;
		int searchResultLimit;
		static const double SPEED_LIGHT;
		static const QString M_PER_SEC;
		static const QString KM_PER_SEC;

		static const double SPLATALOGUE_DEFAULT_MIN;
		static const double SPLATALOGUE_DEFAULT_MAX;
		SearchMoleculesResultDisplayer* resultDisplay;
	};

}

#endif // SEARCHMOLECULESWIDGET_H
