//# PlotMSDataTab.cc: Plot tab for data parameters.
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
#include <plotms/GuiTabs/PlotMSDataTab.qo.h>

#include <casaqt/QtUtilities/QtUtilities.h>
#include <plotms/Actions/PlotMSAction.h>
#include <plotms/Gui/PlotMSAveragingWidget.qo.h>
#include <plotms/Gui/PlotMSPlotter.qo.h>
#include <plotms/Gui/PlotMSSelectionWidget.qo.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>

#include <QDebug>

using namespace casacore;
namespace casa {

///////////////////////////////
// PLOTMSDATATAB DEFINITIONS //
///////////////////////////////

// Constructors/Destructors//

PlotMSDataTab::PlotMSDataTab(PlotMSPlotTab* plotTab, PlotMSPlotter* parent) :
    		PlotMSPlotSubtab(plotTab, parent) {
    ui.setupUi(this);
    // only show these checkbox for bandpass plots (checked by allowAtm)
    ui.overlayLabel->hide();
    ui.noneRadio->hide();
    ui.atmRadio->hide();
    ui.tskyRadio->hide();

    // Setup widgets
    itsFileWidget_ = new QtFileWidget(true, false);
    QtUtilities::putInFrame(ui.locationFrame, itsFileWidget_);

    itsSelectionWidget_ = new PlotMSSelectionWidget();
    QtUtilities::putInFrame(ui.selectionFrame, itsSelectionWidget_);
    QtUtilities::putInScrollArea(ui.selectionFrame);

    itsAveragingWidget_ = new PlotMSAveragingWidget();
    QtUtilities::putInFrame(ui.averagingFrame, itsAveragingWidget_);
    
    itsLabelDefaults_.insert(ui.locationLabel, ui.locationLabel->text());
    itsLabelDefaults_.insert(ui.selectionLabel, ui.selectionLabel->text());
    itsLabelDefaults_.insert(ui.averagingLabel, ui.averagingLabel->text());

    // Connect widgets
    connect(itsFileWidget_, SIGNAL(changed()), SIGNAL(changed()));
    connect(itsSelectionWidget_, SIGNAL(changed()), SIGNAL(changed()));
    connect(itsAveragingWidget_, SIGNAL(changed()), SIGNAL(changed()));
    connect(ui.noneRadio, SIGNAL(toggled(bool)), SIGNAL(changed()));
    connect(ui.atmRadio, SIGNAL(toggled(bool)), SIGNAL(changed()));
    connect(ui.tskyRadio, SIGNAL(toggled(bool)), SIGNAL(changed()));

}



PlotMSDataTab::~PlotMSDataTab() { }

String PlotMSDataTab::getFileName() const {
	return itsFileWidget_->getFile();
}

String PlotMSDataTab::getSelection() const {
	PlotMSSelection selection = itsSelectionWidget_->getValue();
	return selection.toStringShort();
}

String PlotMSDataTab::getAveraging() const {
	PlotMSAveraging averaging = itsAveragingWidget_->getValue();
	return averaging.toStringShort();
}


// Public Methods //

void PlotMSDataTab::getValue(PlotMSPlotParameters& params) const {
    PMS_PP_MSData* d = params.typedGroup<PMS_PP_MSData>();
    if(d == NULL) {
        params.setGroup<PMS_PP_MSData>();
        d = params.typedGroup<PMS_PP_MSData>();
    }
    d->setFilename(itsFileWidget_->getFile());
    d->setSelection(itsSelectionWidget_->getValue());
    d->setAveraging(itsAveragingWidget_->getValue());
    // don't have to check setting showatm since checkbox is only shown
    // when atm is valid:
    d->setShowAtm(ui.atmRadio->isChecked());
    d->setShowTsky(ui.tskyRadio->isChecked());
}



void PlotMSDataTab::setValue(const PlotMSPlotParameters& params) {
    const PMS_PP_MSData* d = params.typedGroup<PMS_PP_MSData>();
    if(d == NULL) return;
    itsFileWidget_->setFile(d->filename());
    bool atm(d->showAtm()), tsky(d->showTsky()), overlay(atm || tsky);
    ui.atmRadio->setChecked(atm);
    ui.tskyRadio->setChecked(tsky);
    ui.noneRadio->setChecked(!overlay);
    itsSelectionWidget_->setValue(d->selection());
    itsAveragingWidget_->setValue(d->averaging());
}

void PlotMSDataTab::update(const PlotMSPlot& plot) {
    const PMS_PP_MSData* d = plot.parameters().typedGroup<PMS_PP_MSData>();
    if(d == NULL) return;

    // overlay radio buttons depend on whether bandpass table
    const casacore::String filename = itsFileWidget_->getFile();
    bool atmAllowed = (!filename.empty() && d->allowAtm(filename));
    ui.overlayLabel->setVisible(atmAllowed);
    ui.noneRadio->setVisible(atmAllowed);
    ui.atmRadio->setVisible(atmAllowed);
    ui.tskyRadio->setVisible(atmAllowed);

    highlightWidgetText(ui.locationLabel, itsFileWidget_->getFile() != d->filename());
    bool overlayChanged = (d->showAtm() != ui.atmRadio->isChecked()) || (d->showTsky() != ui.tskyRadio->isChecked());
    highlightWidgetText(ui.overlayLabel, overlayChanged);
    highlightWidgetText(ui.selectionLabel,
	    itsSelectionWidget_->getValue().fieldsNotEqual(d->selection()));
    highlightWidgetText(ui.averagingLabel,
            itsAveragingWidget_->getValue() != d->averaging());
}

}
