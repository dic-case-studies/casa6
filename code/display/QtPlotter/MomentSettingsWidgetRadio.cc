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


#include "MomentSettingsWidgetRadio.qo.h"
#include <imageanalysis/ImageAnalysis/ImageMoments.h>
#include <display/Display/DisplayCoordinateSystem.h>
#include <display/QtPlotter/ProfileTaskMonitor.h>
#include <display/QtPlotter/ThresholdingBinPlotDialog.qo.h>
#include <display/QtPlotter/conversion/Converter.h>
#include <display/QtPlotter/Util.h>
#include <imageanalysis/Regions/CasacRegionManager.h>
#include <ms/MSOper/MS1ToMS2Converter.h>
#include <display/Display/Options.h>
#include <casa/Logging.h>

#include <QFileDialog>
#include <QTime>
#include <QDebug>
#include <QTemporaryFile>

using namespace casacore;
namespace casa {

	MomentCollapseThreadRadio::MomentCollapseThreadRadio( ImageMoments<Float>* imageAnalysis ):
		analysis( imageAnalysis ), stepSize( 10 ), collapseError(false) {
		imageAnalysis->setProgressMonitor( this );
	}

	CollapseResult::CollapseResult( const String& outputName, bool tmp, std::shared_ptr<ImageInterface<Float>> img ):
				outputFileName(outputName),
				temporary( tmp ),
				image(img) {}

	bool MomentCollapseThreadRadio::isSuccess() const {
		bool success = false;
		if ( collapseResults.size() > 0 && !collapseError ) {
			success = true;
		}
		return success;
	}

	void MomentCollapseThreadRadio::setStepCount( int count ) {
		//We don't want to send too many notifications.
		const int STEP_LIMIT = 100;
		stepSize = count / STEP_LIMIT;
		const int BASE_STEP = 10;
		if ( stepSize < BASE_STEP ) {
			stepSize = BASE_STEP;
		}
		int adjustedCount = count / stepSize + 1;
		emit stepCountChanged( adjustedCount );
	}

	void MomentCollapseThreadRadio::halt() {
		stopImmediately = true;
	}

	void MomentCollapseThreadRadio::setStepsCompleted( int count ) {
		if ( count % stepSize == 0 ) {
			int adjustedCount = count / stepSize;
			emit stepsCompletedChanged( adjustedCount );
		}
	}

	void MomentCollapseThreadRadio::done() {
	}

	String MomentCollapseThreadRadio::getErrorMessage() const {
		return errorMsg;
	}

	void MomentCollapseThreadRadio::setData(const Vector<Int>& mments, const Int axis,
	                                        const Vector<String>& methodVec,
	                                        const Vector<Int>& smoothaxesVec,
	                                        const Vector<String>& smoothtypesVec,
	                                        const Vector<Quantity>& smoothwidthsVec,
	                                        const Vector<Float>& includepixVec,
	                                        const Vector<Float>& excludepixVec,
	                                        const Double peaksnr, const Double stddev,
	                                        const String& doppler, const String& baseName) {
		moments.resize( mments.size());

		moments = mments;
		this->axis = axis;
		method = methodVec;
		smoothaxes = smoothaxesVec;
		smoothtypes = smoothtypesVec;
		smoothwidths = smoothwidthsVec;
		includepix = includepixVec;
		excludepix = excludepixVec;
		this->peaksnr = peaksnr;
		this->stddev = stddev;
		this->doppler = doppler;
		this->baseName = baseName;
	}

	std::vector<CollapseResult> MomentCollapseThreadRadio::getResults() const {
		return collapseResults;
	}

	void MomentCollapseThreadRadio::setChannelStr( String str ) {
		channelStr = str;
	}

	void MomentCollapseThreadRadio::setOutputFileName( QString name ) {
		outputFileName = name;
	}

	void MomentCollapseThreadRadio::setMomentNames( const Vector<QString>& momentNames ) {
		this->momentNames = momentNames;
	}

	bool MomentCollapseThreadRadio::getOutputFileName( String& outName,
	        int /*moment*/, const String& channelStr ) const {

		bool tmpFile = true;
		//Use a default base name
		if (outputFileName.isEmpty()) {
			outName = baseName;
		}
		//Use the user specified name
		else {
			outName = outputFileName.toStdString();
			tmpFile = false;
		}

		//Prepend the channel and moment used to make it descriptive.
		//outName = outName + "_" + String(momentNames[moment].toStdString());
		if ( channelStr != "") {
			outName =  outName+"_"+channelStr;
		}
		if ( tmpFile ) {
			outName = viewer::options.temporaryPath( outName );
		}

		return tmpFile;
	}


	void MomentCollapseThreadRadio::run() {
		try {
			//casa::utilj::ThreadTimes t1;

			//Output file
			String outFile;
			bool outputFileTemporary = getOutputFileName( outFile, 0, channelStr );
			if ( !analysis->setMoments(moments) ){
				errorMsg = analysis->errorMessage();
				collapseError = true;
			}
			else {
				if ( !analysis->setMomentAxis( axis ) ){
					errorMsg = analysis->errorMessage();
					collapseError = true;
				}
				else {
                    try {
					    analysis->setInExCludeRange(includepix, excludepix);
   						auto newImages = analysis->createMoments(
							   false, outFile, false );
						int newImageCount = newImages.size();
						for ( int i = 0; i < newImageCount; i++ ){
							std::shared_ptr<ImageInterface<Float>> newImage = dynamic_pointer_cast<ImageInterface<Float>> (newImages[i]);
							CollapseResult result( outFile, outputFileTemporary, newImage );
							collapseResults.push_back( result );
						}
                 }
                    catch(const AipsError& x) {
						errorMsg = x.getMesg();
						collapseError = true;
					}
				}
			}
			//casa::utilj::ThreadTimes t2;
			//casa::utilj::DeltaThreadTimes dt = t2 - t1;
			//qDebug() << "Elapsed time moment="<<moments[0]<< " elapsed="<<dt.elapsed()<<" cpu="<<dt.cpu();
		}
		catch( AipsError& error ) {
			errorMsg = error.getLastMessage();
			collapseError = true;
		}
	}

	MomentCollapseThreadRadio::~MomentCollapseThreadRadio() {
	}


	MomentSettingsWidgetRadio::MomentSettingsWidgetRadio(QWidget *parent)
		: QWidget(parent), imageAnalysis( NULL ), collapseThread( NULL ),
		  thresholdingBinDialog( NULL ), progressBar( parent ) {
		ui.setupUi(this);

		//Initialize the progress bar
		progressBar.setWindowModality( Qt::ApplicationModal );
		Qt::WindowFlags flags = Qt::Dialog;
		flags |= Qt::FramelessWindowHint;
		progressBar.setWindowFlags( flags);
		progressBar.setWindowTitle( "Collapse/Moments");
		progressBar.setLabelText( "Calculating moments...");
		connect( this, SIGNAL( updateProgress(int)), &progressBar, SLOT( setValue( int )));
		connect( this, SIGNAL( momentsFinished()), &progressBar, SLOT(cancel()));
		connect( &progressBar, SIGNAL(canceled()), this, SLOT(stopMoments()));

		momentOptions << "Mean Value, Mean Intensity" <<
		              "Integrated Value, Sum" <<
		              "Weighted Mean, Velocity Field"<<
		              "Intensity-Weighted Dispersion of Spectral Coordinate, Velocity Dispersion" <<
		              "Median Value, Median Intensity" <<
		              "Spectral Coordinate of Median, Median Velocity Field" <<
		              "Standard Deviation About Mean, Noise, Intensity Scatter" <<
		              "Root Mean Square Intensity"<<
		              "Absolute Mean Deviation" <<
		              "Maximum Intensity, MaximumValue" <<
		              "Spectral Coordinate of Maximum, Velocity of Maximum"<<
		              "Minimum Intensity, MinimumValue" <<
		              "Spectral Coordinate of Minimum, Velocity of Minimum";
		for ( int i = 0; i < static_cast<int>(END_INDEX); i++ ) {
			QListWidgetItem* listItem = new QListWidgetItem( momentOptions[i], ui.momentList);
			if ( i == static_cast<int>(INTEGRATED) ) {
				ui.momentList->setCurrentItem( listItem );
			}
		}
		int columnWidth = ui.momentList->sizeHintForColumn(0);

		ui.momentList->setMinimumWidth( 2*columnWidth/3 );
		ui.momentList->setMaximumWidth( columnWidth );

		//Right now, there is not a clear need for the moment map, but if some
		//moments are no longer used in the display, it will be needed.
		momentMap[MEAN] = MomentsBase<Float>::AVERAGE;
		momentMap[INTEGRATED] = MomentsBase<Float>::INTEGRATED;
		momentMap[WEIGHTED_MEAN] = MomentsBase<Float>::WEIGHTED_MEAN_COORDINATE;
		momentMap[DISPERSION] = MomentsBase<Float>::WEIGHTED_DISPERSION_COORDINATE;
		momentMap[MEDIAN] = MomentsBase<Float>::MEDIAN;
		momentMap[MEDIAN_VELOCITY] = MomentsBase<Float>::MEDIAN_COORDINATE;
		momentMap[STDDEV] = MomentsBase<Float>::STANDARD_DEVIATION;
		momentMap[RMS] = MomentsBase<Float>::RMS;
		momentMap[ABS_MEAN_DEV] = MomentsBase<Float>::ABS_MEAN_DEVIATION;
		momentMap[MAX] = MomentsBase<Float>::MAXIMUM;
		momentMap[MAX_VELOCITY] = MomentsBase<Float>::MAXIMUM_COORDINATE;
		momentMap[MIN] = MomentsBase<Float>::MINIMUM;
		momentMap[MIN_VELOCITY] = MomentsBase<Float>::MINIMUM_COORDINATE;

		ui.channelTable->setColumnCount( 2 );
		QStringList tableHeaders =(QStringList()<< "Min" << "Max");
		ui.channelTable->setHorizontalHeaderLabels( tableHeaders );
		ui.channelTable->setSelectionBehavior(QAbstractItemView::SelectRows);
		ui.channelTable->setSelectionMode( QAbstractItemView::SingleSelection );
		ui.channelTable->setColumnWidth(0, 125);
		ui.channelTable->setColumnWidth(1, 125);

		connect( ui.collapseButton, SIGNAL(clicked()), this, SLOT( collapseImage()));
		connect( ui.channelIntervalCountSpinBox, SIGNAL( valueChanged(int)), this, SLOT(adjustTableRows(int)));
		connect( ui.includeRadioButton, SIGNAL(clicked()), this, SLOT(thresholdingChanged()));
		connect( ui.excludeRadioButton, SIGNAL(clicked()), this, SLOT(thresholdingChanged()));
		connect( ui.noneRadioButton, SIGNAL(clicked()), this, SLOT(thresholdingChanged()));
		connect( ui.outputButton, SIGNAL(clicked()), this, SLOT( setCollapsedImageFile()));
		connect( ui.graphThresholdButton, SIGNAL(clicked()), this, SLOT( graphicalThreshold()));
		connect( ui.symmetricIntervalCheckBox, SIGNAL(stateChanged(int)), this, SLOT(symmetricThresholdChanged(int)));
		connect( ui.maxThresholdLineEdit, SIGNAL(textChanged( const QString&)), this, SLOT(thresholdTextChanged( const QString&)));

		thresholdingChanged();
		ui.channelIntervalCountSpinBox->setValue( 1 );

		//Make sure the min and max thresholds accept only doubles.
		const QDoubleValidator* validator = new QDoubleValidator( this );
		ui.minThresholdLineEdit->setValidator( validator );
		ui.maxThresholdLineEdit->setValidator( validator );
	}

	String MomentSettingsWidgetRadio::makeChannelInterval( float startChannelIndex,
	        float endChannelIndex ) const {
		String channelStr=String::toString( startChannelIndex)+"~"+String::toString(endChannelIndex);

		return channelStr;
	}

	String MomentSettingsWidgetRadio::populateChannels(uInt* nSelectedChannels, bool * ok ) {
		int channelIntervalCount = ui.channelIntervalCountSpinBox->value();
		String channelStr;
		*ok = true;
		for ( int i = 0; i < channelIntervalCount; i++ ) {
			QString startStr;
			QString endStr;
			getChannelMinMax( i, startStr, endStr );

			float startChanVal = startStr.toFloat(ok);
			if ( !(*ok) ){
				return channelStr;
			}
			float endChanVal = endStr.toFloat(ok);
			if ( !(*ok) ){
				return channelStr;
			}

			std::shared_ptr<const ImageInterface<float> > image = taskMonitor->getImage();
			if ( m_units != "Channels"){
				//Convert the units to Channels
				Bool valid = true;
				SpectralCoordinate coord = taskMonitor->getSpectralCoordinate(image, valid );
				if ( valid ){
					Converter* converter = Converter::getConverter( m_units, "" );
					startChanVal = converter->toPixel( startChanVal, coord );
					endChanVal = converter->toPixel( endChanVal, coord );
					delete converter;
				}
				else {
					*ok = false;
					return channelStr;
				}
			}

			if ( isValidChannelRangeValue( QString::number(startChanVal), "Start" ) &&
					isValidChannelRangeValue( QString::number(endChanVal), "End" ) ) {
				if ( endChanVal < startChanVal ) {
					//Switch them around - the code expects the startVal
					//to be less than the endVal;
					float tempVal = startChanVal;
					startChanVal = endChanVal;
					endChanVal = tempVal;
				}
				//Do final check that values are in the range of the spec axis.
				IPosition imShape = image->shape();
				DisplayCoordinateSystem displayCoord = image->coordinates();
				int specIndex = displayCoord.spectralAxisNumber();
				int specMax = imShape[specIndex];
				if ( startChanVal < 0 || startChanVal >= specMax ){
					startChanVal = 0;
				}
				if ( endChanVal < 0 || endChanVal >= specMax ){
					endChanVal = specMax - 1;
				}

				int startChannelIndex = startChanVal;
				int endChannelIndex = endChanVal;
				*nSelectedChannels = *nSelectedChannels + (endChannelIndex - startChannelIndex + 1);

				String channelIntervalStr = makeChannelInterval( startChannelIndex, endChannelIndex );
				if ( i > 0 ) {
					channelStr = channelStr + ",";
				}
				channelStr = channelStr + channelIntervalStr;
			}
			else {
				*ok = false;
			}
		}
		return channelStr;
	}

	bool MomentSettingsWidgetRadio::populateThresholds( Vector<Float>& includeThreshold,
	        Vector<Float>& excludeThreshold ) {
		bool validThresholds = true;
		if ( ui.includeRadioButton->isChecked() ) {
			validThresholds = populateThreshold( includeThreshold );
		} else if ( ui.excludeRadioButton->isChecked() ) {
			validThresholds = populateThreshold( excludeThreshold );
		}
		return validThresholds;
	}

	bool MomentSettingsWidgetRadio::populateThreshold( Vector<Float>& threshold ) {

		bool validThreshold = true;

		//Neither threshold should be blank
		QString minThresholdStr = ui.minThresholdLineEdit->text();
		QString maxThresholdStr = ui.maxThresholdLineEdit->text();
		const int ALL_THRESHOLD = -1;
		double minThreshold = ALL_THRESHOLD;
		double maxThreshold = ALL_THRESHOLD;
		if ( ! minThresholdStr.isEmpty() ) {
			minThreshold = minThresholdStr.toDouble();
		}
		if ( ! maxThresholdStr.isEmpty() ) {
			maxThreshold = maxThresholdStr.toDouble();
		}

		//Minimum should be less than the maximum
		if ( minThreshold > maxThreshold ) {
			validThreshold = false;
			QString msg = "Minimum threshold should be less than the maximum threshold.";
			Util::showUserMessage( msg, this );
		} else {

			//Initialize the vector.
			if ( minThreshold == ALL_THRESHOLD && maxThreshold == ALL_THRESHOLD ) {
				threshold.resize( 1 );
				threshold[0] = ALL_THRESHOLD;
			} else {
				threshold.resize( 2 );
				threshold[0] = minThreshold;
				threshold[1] = maxThreshold;
			}
		}
		return validThreshold;
	}

	Vector<Int> MomentSettingsWidgetRadio::populateMoments( Vector<QString>& momentNames ) {
		//Set up which moments we want
		QList<QListWidgetItem*> selectedItems = ui.momentList->selectedItems();
		int momentCount = selectedItems.length();
		momentNames.resize( momentCount );
		Vector<Int> whichMoments(momentCount);
		if ( momentCount == 0 ) {
			QString msg = "Please select at least one moment.";
			Util::showUserMessage( msg, this );
		} else {
			for( int i = 0; i < momentCount; i++ ) {
				QString selectedText = selectedItems[i]->text();
				momentNames[i] = selectedText;
				int index = momentOptions.indexOf( selectedText );
				if ( index >= 0 ) {
					int momentIndex = static_cast<int>(momentMap[static_cast<SummationIndex>(index)]);
					whichMoments[i] = momentIndex;
				}
			}
		}
		return whichMoments;
	}

	void MomentSettingsWidgetRadio::_initAnalysis(){
		if ( imageAnalysis == nullptr ){
			if ( taskMonitor != nullptr ){
				Record region = _makeRegionRecord( );
				if ( region.nfields() == 0 ){
					String empty("");
					std::shared_ptr<const ImageInterface<float> > image = taskMonitor->getImage();
					if ( image ){
						std::shared_ptr<const SubImage<Float> > result =
							SubImageFactory<Float>::createSubImageRO(*image, region, empty, NULL);
						ImageInterface<Float>* image2 = new SubImage<Float>( *result );
						LogOrigin log("MomentSettingsWidgetRadio", "collapseImage", WHERE);
						LogIO os(log);
						imageAnalysis = new ImageMoments<Float>(*image2, os);
					}
				}

			}

		}
	}

	Record MomentSettingsWidgetRadio::_makeRegionRecord(){

		QString fileName = taskMonitor->getImagePath();
		String infile(fileName.toStdString());

		//Initialize the channels
		uInt nSelectedChannels;
		bool channelOK = true;
		Record region;
		String channelStr = populateChannels( &nSelectedChannels, & channelOK );
		if ( channelOK ){
			std::shared_ptr<const ImageInterface<float> > image = taskMonitor->getImage();
			DisplayCoordinateSystem cSys = image -> coordinates();
			IPosition pos = image->shape();
			String regionName;
			String stokesStr = "";
			CasacRegionManager crm( cSys );
			String diagnostics;
			String pixelBox="";
			region = crm.fromBCS( diagnostics, nSelectedChannels, stokesStr,
										 NULL, regionName, channelStr, CasacRegionManager::USE_FIRST_STOKES,
										 pixelBox, pos, infile);
		}
		return region;
	}


	void MomentSettingsWidgetRadio::collapseImage() {

		// Get the spectral axis number.
		// TODO: Generalize this to any hidden axis
		std::shared_ptr<const ImageInterface<float> > image = taskMonitor->getImage();
		DisplayCoordinateSystem cSys = image -> coordinates();
		int spectralAxisNumber = cSys.spectralAxisNumber();
		if ( spectralAxisNumber < 0 ){
			spectralAxisNumber = Util::getTabularFrequencyAxisIndex( image );
		}
		Vector<String> method;

		//Note default SNRPEAK is 3.  Must be nonnegative.
		Double peaksnr = 3;
		//Note default stddev is 0. Must be nonnegative.
		Double stddev = 0;

		//Initialize the include/exclude pixels
		Vector<Float> excludepix;
		Vector<Float> includepix;
		bool validThresholds = populateThresholds( includepix, excludepix );
		if ( !validThresholds ) {
			return;
		}


		//QString fileName = taskMonitor->getImagePath();
		//String infile(fileName.toStdString());

		_initAnalysis();

		//Set up the thread that will do the work.
		if ( imageAnalysis != nullptr ){
			delete collapseThread;
			collapseThread = new MomentCollapseThreadRadio( imageAnalysis );
			connect( collapseThread, SIGNAL( finished() ), this, SLOT(collapseDone()));
			connect( collapseThread, SIGNAL(stepCountChanged(int)), this, SLOT(setStepCount(int)));
			connect( collapseThread, SIGNAL(stepsCompletedChanged(int)), this, SLOT(setStepsCompleted(int)));

			//Do a collapse image for each of the moments.
			Vector<QString> momentNames;
			Vector<Int> moments = populateMoments( momentNames );
			collapseThread-> setMomentNames( momentNames );
			momentCount = moments.size();
			Vector<Int> smoothaxes;
			Vector<String> smoothtypes;
			Vector<Quantity> smoothwidths;
			String baseName( taskMonitor->getFileName().toStdString());
			collapseThread->setData(moments, spectralAxisNumber,
								     method, smoothaxes, smoothtypes, smoothwidths,
									includepix, excludepix, peaksnr, stddev,
									"RADIO", baseName);
			if ( !outputFileName.isEmpty() ) {
				collapseThread->setOutputFileName( outputFileName );
			}

			uInt nSelectedChannels;
			bool channelOK = true;
			String channelStr = populateChannels( &nSelectedChannels, &channelOK );
			if ( channelOK ){
				collapseThread->setChannelStr( channelStr );
				previousCount = 0;
				cycleCount = 0;
				if ( moments.size() == 1 ){
					progressBar.setCancelButtonText( QString() );
				}
				else {
					progressBar.setCancelButtonText( "Cancel");
				}
				progressBar.show();
		//#warning "Revert to THREADING"
				collapseThread->start();
				//collapseThread->run();
			}
		}
		else {
			QString msg = "Unable to calculate moment(s).";
			Util::showUserMessage( msg, this );
		}
	}

	void MomentSettingsWidgetRadio::collapseDone() {
		//Update the viewer with the collapsed image.
		emit momentsFinished();
		if ( collapseThread != NULL && collapseThread->isSuccess()) {
			std::vector<CollapseResult> results = collapseThread->getResults();
			for ( int i = 0; i < static_cast<int>(results.size()); i++ ) {
				String outName = results[i].getOutputFileName();
				bool outputTemporary = results[i].isTemporaryOutput();
				std::shared_ptr<ImageInterface<Float> > newImage = results[i].getImage();
				taskMonitor->imageCollapsed(outName, "image", "raster", true, outputTemporary, newImage );
			}
			taskMonitor->setPurpose(ProfileTaskMonitor::MOMENTS_COLLAPSE );
		} else {

			QString msg( "Moment calculation failed.");
			String errorMsg = collapseThread->getErrorMessage();
			if ( ! errorMsg.empty() ) {
				msg.append( "\n");
				msg.append( errorMsg.c_str() );
			}
			Util::showUserMessage( msg, this );
		}
	}

	void MomentSettingsWidgetRadio::getChannelMinMax( int channelIndex, QString& minStr, QString& maxStr ) const {
		QTableWidgetItem* minItem = ui.channelTable->item( channelIndex, 0 );
		if ( minItem != NULL ) {
			minStr = minItem->text();
		}
		QTableWidgetItem* maxItem  = ui.channelTable->item( channelIndex, 1 );
		if ( maxItem != NULL ) {
			maxStr = maxItem->text();
		}
	}

	void MomentSettingsWidgetRadio::convertChannelValue( const QString& channelStr,
	        const QString& channelIdentifier, Converter* converter, int row, int col,
	        bool toPixels, SpectralCoordinate& coord ) {
		if ( isValidChannelRangeValue( channelStr, channelIdentifier )) {
			float chanVal = channelStr.toFloat();
			if ( ! toPixels ) {
				chanVal = converter->convert( chanVal, coord );
			} else {
				chanVal = converter->toPixel( chanVal, coord );
			}
			setTableValue( row, col, chanVal );
		}
	}

	void MomentSettingsWidgetRadio::convertChannelRanges( const QString& oldUnits, const QString& newUnits ) {
		int channelIntervalCount = ui.channelIntervalCountSpinBox->value();
		bool toPixels = false;
		if ( newUnits.isEmpty() ) {
			toPixels = true;
		}
		Converter* converter = Converter::getConverter( oldUnits, newUnits );
		if ( imageAnalysis != NULL ){
			std::shared_ptr<const ImageInterface<Float> > imagePtr = taskMonitor->getImage();
			Bool validCoord;
			SpectralCoordinate coord = taskMonitor->getSpectralCoordinate(imagePtr, validCoord );
			for ( int i = 0; i < channelIntervalCount; i++ ) {
				QString startStr;
				QString endStr;
				getChannelMinMax( i, startStr, endStr );
				convertChannelValue( startStr, "Start", converter, i, 0, toPixels, coord );
				convertChannelValue( endStr, "End", converter, i, 1, toPixels, coord );
			}
		}
		delete converter;
	}

	void MomentSettingsWidgetRadio::setUnits( QString unitStr ) {
		int bracketIndex = unitStr.indexOf( "[]");
		if ( bracketIndex > 0 ) {
			unitStr = "";
		}


		QString prevUnit = ui.channelGroupBox->title();
		int startIndex = unitStr.indexOf( "[");
		int endIndex = unitStr.indexOf( "]");
		if ( startIndex > 0 && endIndex > 0 ) {
			unitStr = unitStr.mid(startIndex, endIndex + 1 - startIndex);
		}
		ui.channelGroupBox->setTitle( unitStr );
		if ( prevUnit != "Channels" && prevUnit != unitStr ) {
			QString oldUnits = Util::stripBrackets( prevUnit );
			QString newUnits = Util::stripBrackets( unitStr );
			m_units = newUnits;
			convertChannelRanges( oldUnits, newUnits );
		}
	}

	void MomentSettingsWidgetRadio::setTableValue(int row, int col, float val ) {
		QTableWidgetItem* peakItem = new QTableWidgetItem();
		peakItem -> setText( QString::number( val ) );
		ui.channelTable->setItem( row, col, peakItem );
	}

	void MomentSettingsWidgetRadio::setRange( double min, double max ) {
		if (max < min) {
			ui.channelIntervalCountSpinBox->setValue( 0 );
			ui.channelTable->setRowCount( 0 );
		} else {
			QList<QTableWidgetSelectionRange> selectionRanges = ui.channelTable->selectedRanges();
			int selectionCount = selectionRanges.length();
			if ( selectionCount == 0 ) {
				if ( ui.channelTable->isVisible() && ui.channelIntervalCountSpinBox->value() > 0 ) {
					QString msg( "Please select a row in the channel table before specifying the estimate.");
					Util::showUserMessage( msg, this );
				}
			} else if ( selectionCount == 1 ) {
				QTableWidgetSelectionRange selectionRange = selectionRanges[0];
				int tableRow = selectionRange.bottomRow();
				setTableValue( tableRow, 0, min );
				setTableValue( tableRow, 1, max );
			}
		}
	}

	void MomentSettingsWidgetRadio::reset() {
		if ( imageAnalysis != nullptr ) {
			delete imageAnalysis;
			imageAnalysis = nullptr;
		}
		delete collapseThread;
		collapseThread = nullptr;
		if ( taskMonitor != nullptr ) {
			_initAnalysis();
		}
	}

	void MomentSettingsWidgetRadio::thresholdingChanged( ) {
		bool enabled = false;
		if ( ui.includeRadioButton->isChecked() || ui.excludeRadioButton->isChecked() ) {
			enabled = true;
		}

		ui.maxThresholdLineEdit->setEnabled( enabled );
		ui.symmetricIntervalCheckBox->setEnabled( enabled );
		ui.graphThresholdButton->setEnabled( enabled );
		if ( !ui.symmetricIntervalCheckBox->isChecked() ) {
			ui.minThresholdLineEdit->setEnabled( enabled );
		}
		if ( !enabled ) {
			ui.minThresholdLineEdit->clear();
			ui.maxThresholdLineEdit->clear();
		}
	}


	void MomentSettingsWidgetRadio::adjustTableRows( int count ) {
		ui.channelTable -> setRowCount( count );
		//Select the last row of the table
		QTableWidgetSelectionRange selectionRange(count-1,0,count-1,1);
		ui.channelTable->setRangeSelected( selectionRange, true );
	}

	void MomentSettingsWidgetRadio::setCollapsedImageFile() {
		string homedir = getenv("HOME");
		QFileDialog fd( this, tr("Specify a root file for the collapsed image(s)."),
		                QString(homedir.c_str()), "");
		fd.setFileMode( QFileDialog::AnyFile );
		if ( fd.exec() ) {
			QStringList fileNames = fd.selectedFiles();
			if ( fileNames.size() > 0 ) {
				outputFileName = fileNames[0];
			}
		}

	}

	void MomentSettingsWidgetRadio::symmetricThresholdChanged( int checkedState ) {
		if ( checkedState == Qt::Checked ) {
			ui.minThresholdLineEdit->setEnabled( false );
			//Copy the next from the max to the min.
			QString maxText = ui.maxThresholdLineEdit->text();
			thresholdTextChanged( maxText );
		} else {
			ui.minThresholdLineEdit->setEnabled( true );
		}
	}

	void MomentSettingsWidgetRadio::thresholdTextChanged( const QString& text ) {
		if ( ui.symmetricIntervalCheckBox->isChecked() ) {
			Bool validDouble = false;
			if ( !text.isEmpty() ) {
				double thresholdValue = text.toDouble( &validDouble );
				if ( validDouble ) {
					thresholdValue = -1 * thresholdValue;
					QString oppositeText = QString::number( thresholdValue );
					ui.minThresholdLineEdit->setText( oppositeText );
				} else {
					//Shouldn't get here, but just in case let the user know there is
					//a problem.
					QString msg( "Please specify a valid number for the threshold.");
					Util::showUserMessage( msg, this );
				}
			}
		}
	}

	void MomentSettingsWidgetRadio::thresholdSpecified() {
		pair<double,double> minMaxValues = thresholdingBinDialog->getInterval();
		QString maxValueStr = QString::number( minMaxValues.second);
		ui.maxThresholdLineEdit->setText( maxValueStr );
		if ( ! ui.symmetricIntervalCheckBox->isChecked() ) {
			ui.minThresholdLineEdit->setText( QString::number(minMaxValues.first) );
		} else {
			thresholdTextChanged( maxValueStr );
		}
	}



	void MomentSettingsWidgetRadio::graphicalThreshold() {
		if ( thresholdingBinDialog == NULL ) {
			QString yUnits = this->getYUnit();
			thresholdingBinDialog = new ThresholdingBinPlotDialog( yUnits, this );
			connect( thresholdingBinDialog, SIGNAL(accepted()), this, SLOT(thresholdSpecified()));
		}
		// ImageInterface<Float>* image = const_cast<ImageInterface<Float>* >(taskMonitor->getImage().get());
		std::shared_ptr<ImageInterface<Float> > image(std::const_pointer_cast<ImageInterface<Float> >(taskMonitor->getImage()));
		thresholdingBinDialog->setImage( image );
		thresholdingBinDialog->show();
		QString minValueStr = ui.minThresholdLineEdit->text();
		QString maxValueStr = ui.maxThresholdLineEdit->text();
		double minValue = minValueStr.toDouble();
		double maxValue = maxValueStr.toDouble();
		thresholdingBinDialog->setInterval( minValue, maxValue );
	}
//*************************************************************************
//       Methods from the ImageMomentsProgressMonitor interface
//*************************************************************************

//Note:  because the moments computation is run in a background thread,
//and progress updates must occur in the GUI thread, communication between
//the background thread and the progress bar is via signal/slots.

	void MomentSettingsWidgetRadio::setStepCount( int count ) {
		progressBar.setMinimum( 0 );
		progressBar.setMaximum( count );
		baseIncrement = count / momentCount;
	}

	void MomentSettingsWidgetRadio::setStepsCompleted( int count ) {
		//Cycling over again with a new moment.
		if ( count < previousCount ) {
			cycleCount++;
		}
		int taskCount = cycleCount * baseIncrement + count / momentCount;
		previousCount = count;
		emit updateProgress( taskCount );
	}

	void MomentSettingsWidgetRadio::stopMoments(){
		if ( collapseThread != NULL && collapseThread->isRunning()){
			collapseThread->halt();
		}
	}

	MomentSettingsWidgetRadio::~MomentSettingsWidgetRadio() {
		if ( imageAnalysis != NULL ) {
			delete imageAnalysis;
		}
	}



}
