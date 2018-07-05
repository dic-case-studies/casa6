//# Copyright (C) 2008
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
//#
//# $Id$

#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Plots/PlotMSPlotParameterGroups.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <plotms/test/tUtil.h>


#include <iostream>
#include <msvis/MSVis/UtilJ.h>
#include <casa/namespace.h>
#include <QApplication>

/**
 * Tests whether a simple plot can be exported.
 */
int main(int /*argc*/, char** /*argv[]*/) {

	String dataPath = tUtil::getFullPath( "pm_ngc5921.ms" );
    cout << "tExport using data from "<<dataPath.c_str()<<endl;
    String exportPath = tUtil::getExportPath();
    cout << "Writing plotfiles to " << exportPath << endl;

    // Set up plotms object.
    PlotMSApp app(false, false);


    // Set up default parameters for plot.
	// Can also set parameters for averaging, transformation,
	// iteration, display, etc. in PlotMSPlotParameterGroups
    PlotMSPlotParameters plotParams = PlotMSPlot::makeParameters(&app);

    PMS_PP_MSData* ppdata = plotParams.typedGroup<PMS_PP_MSData>();
    if (ppdata == NULL) {
        plotParams.setGroup<PMS_PP_MSData>();
        ppdata = plotParams.typedGroup<PMS_PP_MSData>();
    }
    ppdata->setFilename( dataPath );
    app.addOverPlot( &plotParams );

    String outFile = exportPath + "plotMSExportTest.jpg";
    PlotExportFormat::Type type = PlotExportFormat::JPG;
	PlotExportFormat format(type, outFile );
	format.resolution = PlotExportFormat::SCREEN;

	bool ok = app.save(format);
	cout << "tExport:: Result of save=" << ok << endl;
	bool okOutput = tUtil::checkFile( outFile, 240000, 260000, -1 );
	cout << "tExport:: Result of save file check=" << okOutput << endl;
    bool test = ok && okOutput;

    // clean up
    tUtil::clearFile(outFile);
    tUtil::clearFile(exportPath);

    bool checkGui = tUtil::exitMain( false );
    return !(test && checkGui);    
}
