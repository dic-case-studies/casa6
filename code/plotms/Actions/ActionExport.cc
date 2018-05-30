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

#include "ActionExport.h"
#include <plotms/Client/Client.h>
#include <plotms/Actions/ActionFactory.h>
#include <plotms/Threads/ThreadController.h>
#include <plotms/Threads/ExportThread.h>
#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <iomanip>
#include <fstream>

#include <QDebug>

using namespace casacore;
namespace casa {

ActionExport::ActionExport( Client* client )
: PlotMSAction( client ), format( PlotExportFormat::JPG, ""){
	itsType_= PLOT_EXPORT;
	interactive = false;
}

void ActionExport::setExportFormat( const PlotExportFormat& exportFormat ){
	format = exportFormat;
}

void ActionExport::setInteractive( bool interactive ){
	this->interactive = interactive;
}

bool ActionExport::loadParameters(){
	bool parametersLoaded = false;
	if ( client != NULL ){
		plots = client->getCurrentPlots();
		interactive = client->isInteractive();
		if ( plots.size() > 0  ){
			if ( !format.location.empty() ){
				parametersLoaded = true;
			}
		}
	}
    return parametersLoaded;
}

bool ActionExport::exportText( PlotMSApp* plotms ){
	Record rec;
	CountedPtr<PlotMSAction> action = ActionFactory::getAction( SEL_INFO, client );
	bool ok = action->doActionWithResponse(plotms, rec);
	if(rec.nfields() < 1) return ok;

	const String X_AXIS("xaxis");
	const String Y_AXIS("yaxis");

	// Write record data to file
	ofstream csv_file;
	csv_file.open(format.location.c_str());
	Record firstRecord = rec.subRecord(0);

	for(uInt n = 0; n < firstRecord.nfields(); ++n) {
		Record plotRecord = firstRecord.subRecord(n);
		Record firstPoint = plotRecord.subRecord( "0" );
		bool hasChan(firstPoint.isDefined("chan")),
			hasAnt2(firstPoint.isDefined("ant2")),
			hasFreq(firstPoint.isDefined("freq"));

		// print header
		csv_file << "# x y ";
		if (hasChan) csv_file << "chan ";
		csv_file << "scan field ant1 ";
		if (hasAnt2) csv_file << "ant2 ";
		csv_file << "ant1name ";
		if (hasAnt2) csv_file << "ant2name ";
		csv_file << "time ";
		if (hasFreq) csv_file << "freq ";
		csv_file << "spw corr obs" << endl;
		// print units
		String xunit(""), yunit("");
		if ( plotRecord.isDefined( X_AXIS )){
			xunit = plotRecord.asString(X_AXIS);
			plotRecord.removeField( X_AXIS );
		}
		if ( plotRecord.isDefined( Y_AXIS)){
			yunit = plotRecord.asString(Y_AXIS);
			plotRecord.removeField( Y_AXIS );
		}
		if ( xunit.length() > 0 || yunit.length() > 0 ){
			csv_file << "# " << xunit << " " << yunit;
			if (hasChan) csv_file << " None"; // chan
			csv_file << " None None None"; // scan field ant1
			if (hasAnt2) csv_file << " None"; // ant2
			csv_file << " None"; // ant1name
			if (hasAnt2) csv_file << " None"; // ant2name
			csv_file << " MJD(seconds)"; // time
			if (hasFreq) csv_file << " GHz"; // freq
			csv_file << " None None None"  // spw corr obs
				<< endl;
		}
		// print plot number
		csv_file << "# From plot " << n << endl;

		// one field per point
		for(uInt _field = 0; _field < plotRecord.nfields(); ++_field) {
			String field_str = String::toString(_field);
			Record fieldRecord = plotRecord.subRecord( field_str );
			Double x = fieldRecord.asDouble("x");
			Double y = fieldRecord.asDouble("y");
			Int chan, ant2;
			if (hasChan) chan = fieldRecord.asInt("chan");
			Int scan = fieldRecord.asInt("scan");
			Int field = fieldRecord.asInt("field");
			Int ant1 = fieldRecord.asInt("ant1");
			if (hasAnt2) ant2 = fieldRecord.asInt("ant2");
			String ant1name = fieldRecord.asString("ant1name");
			String ant2name;
			if (hasAnt2) ant2name = fieldRecord.asString("ant2name");
			Double time = fieldRecord.asDouble("time");
			Int spw = fieldRecord.asInt("spw");
			Double freq;
			if (hasFreq) freq = fieldRecord.asDouble("freq");
			String corr = fieldRecord.asString("corr");
			Int obsId = fieldRecord.asInt("obsid");

			int precision = csv_file.precision();
			if(xunit == "Time") {
				csv_file << std::setprecision(3) << std::fixed
						<< x << " ";
				csv_file.unsetf(ios_base::fixed);
				csv_file.precision(precision);
			} else if(xunit == "Frequency") {
				csv_file << std::setprecision(9) << std::fixed
						<< x << " ";
				csv_file.unsetf(ios_base::fixed);
				csv_file.precision(precision);
			} else {
				csv_file << x << " ";
			}

			if(yunit == "Time") {
				csv_file << std::setprecision(3) << std::fixed
						<< y << " ";
				csv_file.unsetf(ios_base::fixed);
				csv_file.precision(precision);
			} else if(yunit == "Frequency") {
				csv_file << std::setprecision(9) << std::fixed
						<< y << " ";
				csv_file.unsetf(ios_base::fixed);
				csv_file.precision(precision);
			} else {
				csv_file << y << " ";
			}
			if (hasChan) csv_file << chan << " ";
			csv_file << scan << " " << field << " " << ant1 << " ";
			if (hasAnt2) csv_file  << ant2 << " "; 
			csv_file << ant1name << " ";
			if (hasAnt2) csv_file << ant2name << " ";
			csv_file << std::setprecision(3) << std::fixed
					<< time << " ";
			if (hasFreq) csv_file << std::setprecision(9) << std::fixed
					<< freq << " ";
			csv_file.unsetf(ios_base::fixed);
			csv_file.precision(precision);
			csv_file << spw << " " << corr << " " << obsId << endl;
		}
	}
	csv_file.close();
	return ok;
}

PlotExportFormat ActionExport::adjustFormat( PlotExportFormat::Type t){
	PlotExportFormat exportFormat(t, format.location);
	exportFormat.resolution = format.resolution;

	exportFormat.dpi = format.dpi;
	if(exportFormat.dpi <= 0){
		exportFormat.dpi = -1;
	}
	exportFormat.width = format.width;
	if(exportFormat.width <= 0){
		exportFormat.width = -1;
	}
	exportFormat.height = format.height;
	if(exportFormat.height <= 0){
		exportFormat.height = -1;
	}
	return exportFormat;
}

bool ActionExport::doActionSpecific(PlotMSApp* plotms){
	bool ok = true;

	String form = PlotExportFormat::exportFormat( format.type );
	PlotExportFormat::Type t = PlotExportFormat::exportFormat(form, &ok);
	if(!ok) {
		t = PlotExportFormat::typeForExtension(format.location, &ok);
		if(!ok) {
			itsDoActionResult_ = "Invalid format extension for filename '"+
						format.location + "'!";
			return ok;
		}
	}

	if ( t == PlotExportFormat::TEXT ){
		ok = exportText(plotms);
	}
	else {

		PlotExportFormat exportFormat = adjustFormat( t );
		// TODO export fix screen resolution
		// Quick hack for screen resolution images.  Taking a screenshot without
		// drawing the items is basically impossible in the non-main (non-GUI) thread,
		// so for now just turn on high resolution so that it has to draw each
		// items.  This isn't ideal because it is slow, but for now it's better to
		// have something that works and is slow than something that doesn't work.

		if((exportFormat.type == PlotExportFormat::JPG ||
			exportFormat.type == PlotExportFormat::PNG) &&
			exportFormat.resolution == PlotExportFormat::SCREEN) {
			cout << "NOTICE: Exporting to images in screen resolution is currently"
					<< " not working.  Switching to high resolution (which is slower,"
					<< " but works)." << endl;
			exportFormat.resolution = PlotExportFormat::HIGH;
		}

		// Tell parent the format so it can be used in msg to user
		plotms->setExportFormat(exportFormat);

		ExportThread* exportThread = new ExportThread();
		exportThread->setExportFormat( exportFormat );
		exportThread->setPlots( plots );
		setUpClientCommunication( exportThread, -1 );

		setUseThreading( false );

		ok = initiateWork( exportThread );
		if ( threadController == NULL ){
			delete exportThread;
		}
	}
	return ok;
}



ActionExport::~ActionExport() {
}

using namespace casacore;
} /* namespace casa */
