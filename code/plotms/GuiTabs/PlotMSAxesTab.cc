//# PlotMSAxesTab.cc: Plot tab for axes parameters.
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
#include <plotms/GuiTabs/PlotMSAxesTab.qo.h>

#include <casaqt/QtUtilities/QtUtilities.h>
#include <plotms/Gui/PlotMSAxisWidget.qo.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>

using namespace casacore;
namespace casa {

///////////////////////////////
// PLOTMSAXESTAB DEFINITIONS //
///////////////////////////////

PlotMSAxesTab::PlotMSAxesTab(PlotMSPlotTab* plotTab, PlotMSPlotter* parent) :
        PlotMSPlotSubtab(plotTab, parent) {
    setupUi(this);

    // Setup x-axis widget.
    itsXWidget_ = new PlotMSAxisWidget(PMS::DEFAULT_XAXIS, X_BOTTOM | X_TOP);
    itsXWidget_->axisLabel()->setText("X Axis:");
    itsXWidget_->insertLabelDefaults( itsLabelDefaults_ );
    QtUtilities::putInFrame(xFrame, itsXWidget_);
    
    // Initially, hide multiple y-axis support.
    setMultipleAxesYEnabled();

    // Connect widgets.
    connect(itsXWidget_, SIGNAL(axisChanged()), SIGNAL(changed()));
    connect(noneRadio, SIGNAL(toggled(bool)), SIGNAL(changed()));
    connect(atmRadio, SIGNAL(toggled(bool)), SIGNAL(changed()));
    connect(tskyRadio, SIGNAL(toggled(bool)), SIGNAL(changed()));    
    connect(addYButton, SIGNAL(clicked()), this, SLOT(addYWidget()));
    connect(removeYButton, SIGNAL(clicked()), this, SLOT(removeYWidget()));
    connect(yAxisCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(yAxisSelected(int)));

    //Add an initial y-widget.
    addYWidget();
}

void PlotMSAxesTab::setMultipleAxesYEnabled(){
	bool multipleYAxes = false;
	if ( itsYWidgets_.size() > 1 ){
		multipleYAxes = true;
	}
	yAxisLabel->setVisible( multipleYAxes );
	yAxisCombo->setVisible( multipleYAxes );
	removeYButton->setVisible( multipleYAxes );
}

void PlotMSAxesTab::yAxisSelected( int index ){
	//Remove the current y-widget from the display and
	//put in the newly selected one.
	QLayout* yAxisLayout = yAxisFrame->layout();
	if ( yAxisLayout == NULL ){
		yAxisLayout = new QVBoxLayout();
	}
	else {
		QLayoutItem* layoutItem = yAxisLayout->takeAt( 0 );
		if ( layoutItem != NULL ){
			QWidget* widget = layoutItem->widget();
			if ( widget != NULL ){
				widget->setParent( NULL );
			}
			delete layoutItem;
		}
	}
	yAxisLayout->addWidget( itsYWidgets_[index]);
	yAxisFrame->setLayout( yAxisLayout );
}

void PlotMSAxesTab::removeYWidget(){

	//Remove the item from the combo
	int removeIndex = yAxisCombo->currentIndex();
	yAxisCombo->removeItem( removeIndex );

	//Remove the widget
	PlotMSAxisWidget* removeWidget = itsYWidgets_[removeIndex];
	itsYWidgets_.removeAt( removeIndex );

	//Reset the selected y-widget.
	int currentIndex = yAxisCombo->currentIndex();
	yAxisSelected( currentIndex );

	setMultipleAxesYEnabled();
	emit yAxisIdentifierRemoved( removeIndex );
	delete removeWidget;
	emit changed();
}

void PlotMSAxesTab::setYAxisLabel( PlotMSAxisWidget* yWidget,
		int /*index*/ ){
	/*QString numberStr = "";
	if ( index > 0 ){
		numberStr = QString::number( index+1 );
	}
	QString axisLabel = "Y Axis "+ numberStr;
	yWidget->axisLabel()->setText( axisLabel + ":" );*/
	QString identifier = yWidget->getIdentifier();

	yAxisCombo->addItem( identifier );

}

void PlotMSAxesTab::addYWidget(){
	PlotMSAxisWidget* yWidget = new PlotMSAxisWidget(PMS::DEFAULT_YAXIS, Y_LEFT | Y_RIGHT );

	int index = itsYWidgets_.size();

	yWidget->insertLabelDefaults( itsLabelDefaults_ );

	itsYWidgets_.append( yWidget );
	setYAxisLabel( yWidget, index);
	yAxisCombo->setCurrentIndex( index );

	setMultipleAxesYEnabled();

	connect(yWidget, SIGNAL(axisChanged()), SIGNAL(changed()));
	connect(yWidget, SIGNAL(axisIdentifierChanged(PlotMSAxisWidget*)),
			this, SLOT(axisIdentifierChanged(PlotMSAxisWidget*)));
	//Because we may not have been able to connect to listeners when we
	//had just one y-axis, let them know the first one exists, too.
	if ( index <= 1 ){
		emit yAxisIdentifierChanged( 0, itsYWidgets_[0]->getIdentifier());
	}
	emit yAxisIdentifierChanged( index, yWidget->getIdentifier());
	emit changed();
}

PlotMSAxesTab::~PlotMSAxesTab() {
	while( ! itsYWidgets_.isEmpty() ){
		PlotMSAxisWidget* widget = itsYWidgets_.takeAt(0 );
		widget->setParent( NULL );
		delete widget;
	}
}

bool PlotMSAxesTab::isAxesValid() const {
	bool axesValid = true;
	QList<PlotMSAxisWidget*> uniqueYWidgets;
	bool found = false;
	for ( int i = 0; i < itsYWidgets_.size(); i++ ){
		found = false;
		for ( int j = 0; j < uniqueYWidgets.size(); j++ ){
			if ( uniqueYWidgets[j]->matchesData( itsYWidgets_[i])){
				found = true;
				break;
			}
		}
		if ( found ){
			break;
		}
		else {
			uniqueYWidgets.append( itsYWidgets_[i]);
		}
	}
	if ( found ){
		axesValid = false;
	}
	return axesValid;
}


void PlotMSAxesTab::getValue(PlotMSPlotParameters& params) const {
    PMS_PP_Cache* c = params.typedGroup<PMS_PP_Cache>();
    PMS_PP_Axes* a = params.typedGroup<PMS_PP_Axes>();
    PMS_PP_MSData* d = params.typedGroup<PMS_PP_MSData>();

    if(c == NULL) {
        params.setGroup<PMS_PP_Cache>();
        c = params.typedGroup<PMS_PP_Cache>();
    }
    if(a == NULL) {
        params.setGroup<PMS_PP_Axes>();
        a = params.typedGroup<PMS_PP_Axes>();
    }
    if(d == NULL) {
        params.setGroup<PMS_PP_MSData>();
        d = params.typedGroup<PMS_PP_MSData>();
    }
    

    //The cache must have exactly as many x-axes as y-axes so we duplicate
    //the x-axis properties here.
    int yAxisCount = itsYWidgets_.size();

    PMS::Axis xAxis = itsXWidget_->axis();
    c->resize( yAxisCount );
    a->resize( yAxisCount );
    for ( int i = 0; i < yAxisCount; i++ ){
    	c->setXAxis(xAxis, itsXWidget_->data(), i);
    	a->setXAxis(itsXWidget_->attachAxis(), i);
		a->setXRange(itsXWidget_->rangeCustom(), itsXWidget_->range(), i);
    }

    for ( int i = 0; i < yAxisCount; i++ ){
    	c->setYAxis(itsYWidgets_[i]->axis(), itsYWidgets_[i]->data(), i);
    	a->setYAxis(itsYWidgets_[i]->attachAxis(), i);
    	a->setYRange(itsYWidgets_[i]->rangeCustom(), itsYWidgets_[i]->range(), i);
    }

    bool showatm(atmRadio->isChecked());
    bool showtsky(tskyRadio->isChecked());
    c->setShowAtm(showatm);
    c->setShowTsky(showtsky);
    if (showatm || showtsky) {
        // add ATM yaxis "under the hood" for GUI client
        if (xAxis==PMS::CHANNEL || 
            xAxis==PMS::FREQUENCY) { 
            PMS::Axis atmAxis = (showatm ? PMS::ATM : PMS::TSKY);
            bool found(false);
            const vector<PMS::Axis> yAxes = c->yAxes();
            for (uInt i=0; i<yAxes.size(); ++i) {
                if (yAxes[i] == atmAxis) {
                        found=True;
                        break;
                }
            }
            if (!found) {
                // add axis to Cache axes
                int index = c->numXAxes();
                c->setAxes(xAxis, atmAxis, c->xDataColumn(0), PMS::DEFAULT_DATACOLUMN, index);
                // set Axes positions
                a->resize(index+1, true);  // copy values for index 0
                a->setAxes(a->xAxis(index-1), Y_RIGHT, index);
				// keep same xaxis range
    			a->setXRange(itsXWidget_->rangeCustom(), itsXWidget_->range(), index);
                // set Display symbol color
                PMS_PP_Display* disp = params.typedGroup<PMS_PP_Display>();
                PlotSymbolPtr atmSymbol = disp->unflaggedSymbol(index);
                atmSymbol->setSymbol("circle");
                atmSymbol->setColor("#FF00FF");
                disp->setUnflaggedSymbol(atmSymbol, index);
                PlotSymbolPtr flaggedSymbol = disp->flaggedSymbol();
                disp->setFlaggedSymbol(flaggedSymbol, index);
            }
        } 
    }
}

void PlotMSAxesTab::setValue(const PlotMSPlotParameters& params) {
    const PMS_PP_Cache* c = params.typedGroup<PMS_PP_Cache>();
    const PMS_PP_Axes* a = params.typedGroup<PMS_PP_Axes>();
    if(c == NULL || a == NULL) return; // shouldn't happen

    PMS::Axis cacheAxis = c->xAxis();    
    PMS::DataColumn cacheColumn =  c->xDataColumn();
    PlotAxis axesAxis = a->xAxis();
    bool axesXRangeSet = a->xRangeSet();
    std::pair<double, double> axesXRange = a->xRange();
    itsXWidget_->setValue(cacheAxis, cacheColumn, axesAxis,
            axesXRangeSet, axesXRange);
    //itsXWidget_->setValue(c->xAxis(), c->xDataColumn(), a->xAxis(),
    //        a->xRangeSet(), a->xRange());
    bool atm(c->showAtm()), tsky(c->showTsky()), overlay(atm || tsky);
    atmRadio->setChecked(atm);
    tskyRadio->setChecked(tsky);
    noneRadio->setChecked(!overlay);
    int yAxisCount = a->numYAxes();
    int yWidgetSize = itsYWidgets_.size();
    int minValue = qMin( yAxisCount, yWidgetSize );
    for ( int i = 0; i < minValue; i++ ){
    	itsYWidgets_[i]->setValue(c->yAxis(i), c->yDataColumn(i), a->yAxis(i),
            a->yRangeSet(i), a->yRange(i));
    }
}

void PlotMSAxesTab::axisIdentifierChanged(PlotMSAxisWidget* axisWidget){

	int yIndex = itsYWidgets_.indexOf( axisWidget );
	if ( yIndex >= 0 ){
		QString newIdentifier = axisWidget->getIdentifier();
		yAxisCombo->setItemText(yIndex, newIdentifier );
		emit yAxisIdentifierChanged( yIndex, newIdentifier );
	}
}

void PlotMSAxesTab::update(const PlotMSPlot& plot) {
    const PlotMSPlotParameters& params = plot.parameters();
    PlotMSPlotParameters newParams(params);
    getValue(newParams);
    
    const PMS_PP_MSData* d = params.typedGroup<PMS_PP_MSData>();
    const PMS_PP_Cache* c = params.typedGroup<PMS_PP_Cache>(),
                      *c2 = newParams.typedGroup<PMS_PP_Cache>();
    const PMS_PP_Axes* a = params.typedGroup<PMS_PP_Axes>(),
                     *a2 = newParams.typedGroup<PMS_PP_Axes>();
    
    // shouldn't happen
    if(d == NULL || c == NULL || c2 == NULL || a == NULL || a2 == NULL) return;

    // Update "in cache" for widgets.
    vector<PMS::Axis> laxes = plot.cache().loadedAxes();
    bool found = false;
    for(unsigned int i = 0; !found && i < laxes.size(); i++){
        if(laxes[i] == c2->xAxis()){
        	found = true;
        	break;
        }
    }
    itsXWidget_->setInCache(found);
   
    int yAxisCount = itsYWidgets_.size();
    for ( int j = 0; j < yAxisCount; j++ ){
    	found = false;
    	PMS::Axis yData = c2->yAxis(j);
    	for(unsigned int i = 0; !found && i < laxes.size(); i++){
    		if(laxes[i] == yData){
    			found = true;
    			break;
    		}
    	}
    	itsYWidgets_[j]->setInCache(found);
    }

    // Update labels.
    QLabel* axisXLabel = itsXWidget_->axisLabel();
    bool axisDifference=false;
    if ( c->xAxis() != c2->xAxis() ){
    	axisDifference = true;
    }
    bool dataDifference = false;
    if ( PMS::axisIsData(c->xAxis()) &&
            c->xDataColumn() != c2->xDataColumn()){
    	dataDifference = true;
    }

    highlightWidgetText(axisXLabel, d->isSet() && (axisDifference || dataDifference ));
    highlightWidgetText(itsXWidget_->dataLabel(), d->isSet() && dataDifference);
    highlightWidgetText(itsXWidget_->attachLabel(), d->isSet() &&
                a->xAxis() != a2->xAxis());
    highlightWidgetText(itsXWidget_->rangeLabel(), d->isSet() &&
                (a->xRangeSet() != a2->xRangeSet() ||
                (a->xRangeSet() && a->xRange() != a2->xRange())));

    bool overlayChanged = (c->showAtm() != atmRadio->isChecked()) || (c->showTsky() != tskyRadio->isChecked());
    highlightWidgetText(overlayLabel, overlayChanged);
    
    for ( int i = 0; i < yAxisCount; i++ ){
    	QLabel* axisLabel = itsYWidgets_[i]->axisLabel();
    	highlightWidgetText(axisLabel, d->isSet() && (c->yAxis(i) != c2->yAxis(i) ||
                (PMS::axisIsData(c->yAxis(i)) &&
                 c->yDataColumn(i) != c2->yDataColumn(i))));
    	highlightWidgetText(itsYWidgets_[i]->dataLabel(), d->isSet() &&
                c->yDataColumn(i) != c2->yDataColumn(i));
    	highlightWidgetText(itsYWidgets_[i]->attachLabel(), d->isSet() &&
                a->yAxis(i) != a2->yAxis(i));
    	highlightWidgetText(itsYWidgets_[i]->rangeLabel(), d->isSet() &&
                (a->yRangeSet(i) != a2->yRangeSet(i) ||
                (a->yRangeSet(i) && a->yRange(i) != a2->yRange(i))));
    }

    // If the user hasn't set a custom range, set defaults for time axis.
    // For time axis, get bounds from cache; else use 0 for other axes types
    pair<Double,Double> timebounds;
    bool isDate;

    if (!itsXWidget_->rangeCustom()) {
    	isDate = (c2->xAxis() == PMS::TIME);
    	if (isDate) {
		timebounds = plot.cache().getTimeBounds();
		itsXWidget_->setRange(isDate, timebounds.first, timebounds.second); 
    	}
    	else {
		itsXWidget_->setRange(isDate, 0.0, 0.0);
    	}
    }

    for (int i = 0; i < yAxisCount; i++ ){
	if (!itsYWidgets_[i]->rangeCustom()) {
    		PMS::Axis yData = c2->yAxis(i);
    		isDate = (yData == PMS::TIME);
    		if (isDate) {
			timebounds = plot.cache().getTimeBounds();
			itsYWidgets_[i]->setRange(isDate, timebounds.first, timebounds.second); 
			break;
    		}
        	else {
	    	itsYWidgets_[i]->setRange(isDate, 0.0, 0.0);
        	}
    	}
    }
}

}
