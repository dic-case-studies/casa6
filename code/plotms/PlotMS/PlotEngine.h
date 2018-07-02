//# Copyright (C) 2008
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#

#ifndef PLOTENGINE_H_
#define PLOTENGINE_H_

#include <graphics/GenericPlotter/PlotLogger.h>
#include <graphics/GenericPlotter/PlotFactory.h>
#include <plotms/Plots/PlotMSPlotManager.h>
#include <plotms/PlotMS/PlotMSExportParam.h>
#include <plotms/Plots/PlotMSPlotParameters.h>
#include <plotms/PlotMS/PlotMSFlagging.h>

namespace casa {

class PlotMSPlot;

/**
 * Abstracts the interface DBUS uses to talk to the main PlotMS
 * application.
 */

class PlotEngine {
public:

	// Gets the logger associated with this PlotMS.
	virtual PlotLoggerPtr getLogger() = 0;

	// Returns the PlotMSPlotManager associated with this PlotMS.
	virtual PlotMSPlotManager& getPlotManager() = 0;
	virtual PlotFactoryPtr getPlotFactory() = 0;
	virtual PlotMSParameters& getParameters() = 0;
	virtual PlotMSExportParam& getExportParameters() = 0;

	//Symbols
	virtual PlotSymbolPtr createSymbol (const casacore::String& descriptor,
	    		casacore::Int size, const casacore::String& color,
	        	const casacore::String& fillPattern, bool outline ) = 0;

	virtual PlotMSPlot* addOverPlot(const PlotMSPlotParameters* p = NULL) = 0;


	//Show/hide
	virtual void allowPopups(bool allow = true) = 0;
	virtual bool popupsAllowed() = 0;
	virtual void showGUI(bool show = true) = 0;
	virtual bool guiShown() const = 0;
	virtual bool isClosed() const = 0;
	virtual void clearPlots() = 0;

	virtual bool isOperationCompleted() const = 0;
	virtual void setOperationCompleted( bool complete) = 0;
	virtual bool save(const PlotExportFormat& format) = 0;
	virtual PlotMSFlagging getFlagging() const = 0;
	virtual void setFlagging(PlotMSFlagging flag) = 0;
	virtual void quitApplication() = 0;
	virtual bool isDrawing() const = 0;

	virtual casacore::Record locateInfo( casacore::Bool& success, casacore::String& errorMessage ) = 0;

protected:
	PlotEngine(){}
	virtual ~PlotEngine(){};
};

}
#endif /* PLOTENGINE_H_ */
