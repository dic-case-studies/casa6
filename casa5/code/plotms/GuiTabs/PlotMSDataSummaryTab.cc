//# Copyright (C) 2009
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
//# $Id: $
#include <plotms/GuiTabs/PlotMSDataSummaryTab.qo.h>


#include <casaqt/QwtPlotter/QPPlotter.qo.h>
#include <casaqt/QtUtilities/QtUtilities.h>
#include <plotms/Actions/PlotMSAction.h>
#include <plotms/GuiTabs/PlotMSDataCollapsible.qo.h>
#include <plotms/Gui/PlotMSPlotter.qo.h>

#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <QDebug>

#include <QMessageBox>
#include <plotms/Data/MSCache.h>

using namespace casacore;
namespace casa {


// Constructors/Destructors//

PlotMSDataSummaryTab::PlotMSDataSummaryTab(PlotMSPlotter* parent) :
		PlotMSTab(parent){

	ui.setupUi(this);

	its_force_reload = false;
	makingPlot = false;
	rowLimit = 1;
	colLimit = 1;

	 // Add as watcher.
	parent->getParent()->getPlotManager().addWatcher(this);

	layout()->setAlignment(Qt::AlignTop);
    
    scrollWidget = new QWidget( ui.dataScroll );
    //scrollWidget->setStyleSheet("QScrollArea#scrollWidget { background-color: green; }");
    //scrollWidget->setMaximumHeight(850);

    ui.dataScroll->setWidget( scrollWidget );
    ui.dataScroll->setWidgetResizable( true );
    QVBoxLayout* scrollLayout = new QVBoxLayout();
    scrollLayout->setSpacing( 2 );
    scrollLayout->setContentsMargins(2,2,2,2);
    scrollWidget->setLayout( scrollLayout );

    scrollWidget->setAutoFillBackground( true );
    QPalette pal = scrollWidget->palette();
    QColor bgColor( Qt::blue );
    pal.setColor( QPalette::Background, bgColor );
    scrollWidget->setPalette( pal );
    bottomSpacer = new QSpacerItem( 200,20, QSizePolicy::Preferred, QSizePolicy::Expanding );

    connect( ui.addSingleButton, SIGNAL(clicked()), this, SLOT(addSinglePlot()));


    //So plot parameters are up and available for the user.
    addSinglePlot();


    // Additional slot for Plot button for shift+plot forced redraw feature.
     // All this does is note if shift was down during the click.
     // This slot should be called before the main one in synchronize action.
     connect( ui.plotButton, SIGNAL(clicked()),  SLOT(observeModKeys()) );

     // To fix zoom stack first-element bug, must be sure Zoom, Pan, etc
     // in toolbar are all unclicked.   Leave it to the Plotter aka
     // "parent" know how to do this, offering a slot
     // we can plug the Plot button into.
     // The prepareForPlotting() slot is available for other things to do
     // at this point in time (like evil diagnostic experiements).
     connect( ui.plotButton, SIGNAL(clicked()),  parent, SLOT(prepareForPlotting()) );

     // Synchronize plot button.  This makes the reload/replot happen.
     itsPlotter_->synchronizeAction(PlotMSAction::PLOT, ui.plotButton);
}

PlotMSDataSummaryTab::~PlotMSDataSummaryTab() {

}

void PlotMSDataSummaryTab::emptyLayout(){
    // Start at end of QList to avoid indexing issues...
    for ( int i=(dataList.size()-1); i > 0; --i ){
        close(dataList[i]);
    }
    dataList[0]->clearData();
}

void PlotMSDataSummaryTab::setGridSize( int rowCount, int colCount ){

    //Store the maximum number of plots we can support.
    rowLimit = rowCount;
    colLimit = colCount;


    //Tell everyone to update their grid size, deleting any whose location
    //exceeds the current limits.
    int dataCount = dataList.size();
    for ( int i = dataCount-1; i >=0; --i ){
        if (i >= rowCount*colCount) {
            close(dataList[i]);
        } else {
            dataList[i]->setGridSize( rowCount, colCount );
            if (dataList[i]->isMinimized())
                dataList[i]->maximizeDisplay();
        }
    }

    //Generate new plots.
    //plot( true );
}

void PlotMSDataSummaryTab::fillLayout(){
	QLayout* scrollLayout = scrollWidget->layout();
	for ( int i = 0; i < dataList.size(); i++ ){
		scrollLayout->addWidget( dataList[i]);
	}
}


void PlotMSDataSummaryTab::addSinglePlot( int existingIndex){
	//Minimize any open plots so the new one will be seen.

	int dataCount = dataList.size();
	for ( int i = 0; i < dataCount; i++ ){
		dataList[i]->minimizeDisplay();
	}
	insertData( existingIndex );
}


void PlotMSDataSummaryTab::insertData( int index ){
	/*****
	 * NOTE::  When we support overlays, we will have
	 * to change this method so it doesn't just insert
	 * at index, but it figures out what the index should
	 * be with overlays thrown in.
	 */
	PlotMSDataCollapsible* plotTab = NULL;
	if ( 0 <= index && index < dataList.size() ){
		plotTab = dataList[index];
	}
	else {
		makingPlot = true;
		QLayout* scrollLayout = scrollWidget->layout();
		scrollLayout->removeItem( bottomSpacer );
		plotTab = new PlotMSDataCollapsible( itsPlotter_, scrollWidget, index );
		dataList.append( plotTab );
		makingPlot = false;
		connect(  plotTab, SIGNAL( close(PlotMSDataCollapsible*)),
				this, SLOT( close(PlotMSDataCollapsible*)));
		plotTab->setGridSize( rowLimit, colLimit );

		scrollLayout->addWidget( plotTab );
		scrollLayout->addItem( bottomSpacer );
	}
}

void PlotMSDataSummaryTab::parametersHaveChanged(const PlotMSWatchedParameters& p,
        int updateFlag) {
	for ( int i = 0; i < dataList.size(); i++ ){
		dataList[i]->parametersHaveChanged( p, updateFlag );
	}
}

// Implements PlotMSPlotManagerWatcher::plotsChanged().
void PlotMSDataSummaryTab::plotsChanged(const PlotMSPlotManager& manager){
	int dataCount = dataList.size();
	int managerPlotCount = manager.numPlots();
	//Make sure scriptors have not added some plots through DBUS that we
	//don't know about while showing the GUI.
	if ( managerPlotCount > dataCount && !makingPlot ){
		for ( int i = dataCount; i < managerPlotCount; i++ ){
			this->addSinglePlot(i);
		}
	}
	for ( int i = 0; i < dataCount; i++ ){
		dataList[i]->plotsChanged( manager, i );
	}

}

void PlotMSDataSummaryTab::singleDataChanged(PlotMSDataCollapsible* collapsible){
	int singleDataIndex = dataList.indexOf( collapsible );
	emit changed( singleDataIndex );
}

void PlotMSDataSummaryTab::close( PlotMSDataCollapsible* collapsible ){
	QLayout* scrollLayout = scrollWidget->layout();
	scrollLayout->removeWidget( collapsible );
	int collapseIndex = dataList.indexOf( collapsible );
	if ( collapseIndex >= 0 ){
		dataList.removeAt( collapseIndex );
	}
	delete collapsible;
	for (int i=0; i<dataList.size(); ++i) {
		dataList[i]->maximizeDisplay();
		dataList[i]->plot(false);
	}
}

void PlotMSDataSummaryTab::refreshPageHeader(){
	auto * controllerQtDataModel = new PlotMSPageHeaderDataModel(itsParent_);
	QtPageHeaderDataModelPtr controllerDataModelPtr { new QtPageHeaderDataModel(controllerQtDataModel) };

	auto plotter = itsParent_->getPlotter();

	plotter->refreshPageHeaderDataModel(controllerDataModelPtr);

}

bool PlotMSDataSummaryTab::plot(){
	bool plotted = false;
	//If there is at least one of the graphs updating its data
	//we want to return plotted=true;
	int dataCount = dataList.size();
	for ( int i = 0; i < dataCount; i++ ){
		if ( dataList[i]->plot( its_force_reload ) ){
			plotted = true;
		}
	}
	refreshPageHeader();
	return plotted;
}

void PlotMSDataSummaryTab::completePlotting( bool success, PlotMSPlot* plot ){

	//Reset the plotIndex, which may not be reset if the person is scripting.
	//and has not clicked the 'plot' button.  Please see CAS-6813.
	int completedIndex = -1;

	int dataCount = dataList.size();
	for ( int i = 0; i < dataCount; i++ ){
		if ( dataList[i]->managesPlot( plot ) ){
			completedIndex = i;
			break;
		}
	}
	if ( completedIndex >= 0 ){
		completePlotting( success, completedIndex );

	}
	refreshPageHeader();
}

void PlotMSDataSummaryTab::completePlotting( bool success, int plotIndex ){
	dataList[plotIndex]->clearData();
	dataList[plotIndex]->completePlotting( success);
}

vector<vector<PMS::Axis> > PlotMSDataSummaryTab::selectedLoadAxes() const {
	vector<vector<PMS::Axis> > loadAxes;
	for ( int i = 0; i < dataList.size(); i++ ){
		loadAxes.push_back( dataList[i]->getSelectedLoadAxes() );
	}
	return loadAxes;
}

vector<vector<PMS::Axis> > PlotMSDataSummaryTab::selectedReleaseAxes() const {
	vector<vector<PMS::Axis> > releaseAxes;
	for ( int i = 0; i < dataList.size(); i++ ){
		releaseAxes.push_back( dataList[i]->getSelectedReleaseAxes() );
	}
	return releaseAxes;
}

vector<PlotMSPlot*> PlotMSDataSummaryTab::getCurrentPlots(){
	vector<PlotMSPlot*> currentPlots;
	for ( int i = 0; i < dataList.size(); i++ ){
		PlotMSPlot* plot = dataList[i]->getPlot();
		if ( plot != NULL ){
			currentPlots.push_back( plot );
		}
	}
	return currentPlots;
}

void PlotMSDataSummaryTab::observeModKeys()   {
	// Bitflags report if shift, etc were down during click of Plot button
	Qt::KeyboardModifiers itsModKeys = QApplication::keyboardModifiers();
	bool using_shift_key = (itsModKeys & Qt::ShiftModifier) !=0;
	bool always_replot_checked = ui.forceReplotChk->isChecked();
	its_force_reload = using_shift_key  ||  always_replot_checked;
}

vector<String> PlotMSDataSummaryTab::getFiles() const {
	vector<String> loadedFiles;
	for ( int i = 0; i < dataList.size(); i++ ){
		String fileName = dataList[i]->getFile();
		if ( fileName.length() > 0 ){
			loadedFiles.push_back( fileName );
		}
	}
	return loadedFiles;
}


void PlotMSDataSummaryTab::resizeEvent( QResizeEvent* /*event*/ ){
	QSize currentSize = size();
	int usedHeight = 0;
	for ( int i = 0; i < dataList.size(); i++ ){
		QSize widgetSize = dataList[i]->sizeHint();
		usedHeight = usedHeight + widgetSize.height();
	}

	int openIndex = -1;
	for ( int i = 0; i < dataList.size(); i++ ){
		if ( !dataList[i]->isMinimized() ){
			openIndex = i;
		}
	}
	//Pass the height increase/descrease to the open one.
	if ( openIndex >= 0 ){
		int heightDiff = currentSize.height() - usedHeight;
		dataList[openIndex]->resetHeight( heightDiff );
	}
}

}
