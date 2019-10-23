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

#ifndef QPEXPORTER_H_
#define QPEXPORTER_H_

#include <graphics/GenericPlotter/PlotOptions.h>
#include <graphics/GenericPlotter/PlotCanvas.h>
#include <QImage>

namespace casa {

class PlotCanvas;
class QPExportCanvas;
class QPPlotter;
/**
 * Utility class for exporting plots.
 */

class QPExporter {
public:


	// Exports the given canvas to the given format.
	static bool exportCanvas(PlotCanvas* canvas, const PlotExportFormat& format);

	// Exports the given plotter to the given format.
	static bool exportPlotter(QPPlotter* plotter, const PlotExportFormat& format);

	// Exports a collection of canvases to the given format.
	static bool exportCanvases(std::vector<QPExportCanvas*>& canvases,
			const PlotExportFormat& format, PlotCanvas* grabCanvas,
			QPPlotter* grabPlotter);

	virtual ~QPExporter();

private:
	QPExporter();
	static bool exportPostscript( const PlotExportFormat& format,
			std::vector<QPExportCanvas*> &qcanvases,
			QPExportCanvas* grabCanvas, QPPlotter* grabPlotter);

	static QImage produceHighResImage(
			const PlotExportFormat& format,
			std::vector<QPExportCanvas*> &qcanvases,
			int width, int height,
			int rowIndex, int columnIndex,
			bool &wasCanceled);

	static QImage produceScreenImage(
			const PlotExportFormat& format,
			std::vector<QPExportCanvas*> &qcanvases,
			int width, int height,
			int rowCount, int colCount,
			bool &wasCanceled);

	static bool  exportToImageFile(
			const PlotExportFormat& format,
			std::vector<QPExportCanvas*> &qcanvases,
			QPExportCanvas* grabCanvas,
			QPPlotter* grabPlotter);

	static int findAxisHeight( std::vector<QPExportCanvas*> &qcanvases );
	static int findAxisWidth( std::vector<QPExportCanvas*> &qcanvases );
	static int getCanvasCount( std::vector<QPExportCanvas*> &qcanvases );
    static void getAxesCount(std::vector<QPExportCanvas*> &qcanvases,
            casacore::Int& externalX, casacore::Int& externalY);
	static void findGridProperties( QPExportCanvas* grabCanvas, QPPlotter* grabPlotter,
			casacore::Int& width, casacore::Int& height, casacore::Int& gridRows, casacore::Int& gridCols);
    static void findXAxisLocations(casacore::Int numX, casacore::Bool vertical, casacore::Bool& top, casacore::Bool& bottom);
    static void findYAxisLocations(casacore::Int numY, casacore::Bool vertical, casacore::Bool& left, casacore::Bool& right);
    static void findYAxisSecondRow(casacore::Int numY, casacore::Bool isLeftAxis, casacore::Bool& left, casacore::Bool& right);
    static void findYAxisSecondRow(casacore::Int numY, casacore::Int nCols, std::vector<QPExportCanvas*> &qcanvases,
            casacore::Bool& left, casacore::Bool& right);

	static const casacore::String CLASS_NAME;
	static const casacore::String EXPORT_NAME;
};

} /* namespace casa */
#endif /* QPEXPORTER_H_ */
