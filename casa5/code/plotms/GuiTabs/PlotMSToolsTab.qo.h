//# PlotMSToolsTab.qo.h: Subclass of PlotMSTab for tools management.
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
#ifndef PLOTMSTOOLSTAB_QO_H_
#define PLOTMSTOOLSTAB_QO_H_

#include <ui/ui_PlotMSToolsTab.h>

#include <graphics/GenericPlotter/PlotTool.h>
#include <plotms/GuiTabs/PlotMSTab.qo.h>

namespace casa {


class PlotMSToolsTab;  /* fwd */

// Registered with all Canvases so Tracker can act upon key presses
class TrackerKeyHandler : public PlotKeyEventHandler  {
	
	public:
		TrackerKeyHandler(PlotMSToolsTab *);
				
		virtual void handleKey(const PlotKeyEvent& event);
	
	private: 
		PlotMSToolsTab *tools_tab;
};




// Subclass of PlotMSTab that handles the tools for the current plot.  Watches
// no parameters.
class PlotMSToolsTab : public PlotMSTab, Ui::ToolsTab,
                       public PlotTrackerToolNotifier {
    Q_OBJECT
    
public:
    // Constructor which takes the parent plotter, and the QtActionGroup to
    // use to synchronize tool actions with the radio buttons on the tab.
    PlotMSToolsTab(PlotMSPlotter* parent);
    
    // Destructor.
    ~PlotMSToolsTab();
    
    
    // Implements PlotMSTab::tabName().
    QString tabName() const { return "Tools"; }
    
    // Overrides PlotMSTab::toolButtons().
    QList<QToolButton*> toolButtons() const;
    
    // Implements PlotMSParametersWatcher::parametersHaveChanged.  Currently
    // does nothing.
    void parametersHaveChanged(const PlotMSWatchedParameters& params,
            int updateFlag) { (void)params,(void)updateFlag; }
    
    
    // Show/hide the iteration buttons on this tab.
    void showIterationButtons(bool show);
    
    

    
public slots:
    // Slot for when all tools are turned off, and the "None" radio button
    // should be checked.
    void toolsUnchecked();

    // Tracker "snapshot" feature. Copies value in live display
    // into multi-line text box for user to copy/paste.
    // Made a slot in case it's useful to connect to a signal, but
    // for the initial version, this is not done.
    void takeSnapshotOfTrackerValue();
    
    // Erase contents of the text box holding recorded tracker values
    void clearRecordedValues();
    
protected:
    // Implements PlotTrackerToolNotifier::notifyTrackerChanged().  Updates the
    // tracker information in the line edit, if the proper checkbox is toggled.
    void notifyTrackerChanged(PlotTrackerTool& tool);

public:
    TrackerKeyHandler *tracker_key_handler;
    
};

}

#endif /* PLOTMSTOOLSTAB_QO_H_ */
