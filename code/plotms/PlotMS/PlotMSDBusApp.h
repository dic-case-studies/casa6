//# PlotMSDBusApp.h: Controller for plotms using DBus.
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
#ifndef PLOTMSDBUSAPP_H_
#define PLOTMSDBUSAPP_H_

#include <casaqt/QtUtilities/QtDBusXmlApp.qo.h>
#include <plotms/PlotMS/PlotMSParameters.h>
#include <plotms/Plots/PlotMSPlotManager.h>
#include <plotms/PlotMS/PlotEngine.h>

namespace casa {

//# Forward declarations.
//class PlotEngine;


// Subclass of QtDBusXmlApp to control plotms using DBus communication.
class PlotMSDBusApp: public QtDBusXmlApp, public PlotMSParametersWatcher,
                     public PlotMSPlotManagerWatcher {
    
    //# Friend class declarations.
    friend class PlotMSDBusAppWatcher;
    
public:

    static const QString &name( );
    QString dbusName( ) const { return QString(name( )); }

    // Static //
    
    // Constants for the casaplotms standalone executable.
    // <group>
    static const casacore::String APP_NAME;
    static const casacore::String APP_CASAPY_SWITCH;
    static const casacore::String APP_LOGFILENAME_SWITCH;
    static const casacore::String APP_LOGFILTER_SWITCH;
    // </group>
    
    
    // PARAMETERS //
    
    // Parameter names.
    // <group>
    static const casacore::String PARAM_AVERAGING; // casacore::Record (see PlotMSAveraging)
    static const casacore::String PARAM_AXIS_X; // String
    static const casacore::String PARAM_AXIS_Y; // String
    static const casacore::String PARAM_AXIS_Y_LOCATION;
    static const casacore::String PARAM_SHOWATM; // bool
    static const casacore::String PARAM_SHOWTSKY; // bool
    static const casacore::String PARAM_GRIDROWS; //int
    static const casacore::String PARAM_GRIDCOLS; //int
    static const casacore::String PARAM_CLEARSELECTIONS; // bool
    static const casacore::String PARAM_DATACOLUMN_X; // String
    static const casacore::String PARAM_DATACOLUMN_Y; // String
    static const casacore::String PARAM_FILENAME; // String
    static const casacore::String PARAM_FLAGGING; // Record
    static const casacore::String PARAM_HEIGHT; // int or uInt
    static const casacore::String PARAM_ITERATE; // casacore::Record (see PlotMSIterParam)
    static const casacore::String PARAM_PLOTINDEX; // int or uInt
    static const casacore::String PARAM_PRIORITY; // String
    static const casacore::String PARAM_SELECTION; // casacore::Record (see PlotMSSelection)
    static const casacore::String PARAM_TRANSFORMATIONS; // casacore::Record (see PlotMSTransformations)
    static const casacore::String PARAM_CALIBRATION; // casacore::Record (see PlotMSCalibration)
    static const casacore::String PARAM_PAGE_HEADER_ITEMS; // String
    static const casacore::String PARAM_UPDATEIMMEDIATELY; // bool
    static const casacore::String PARAM_WIDTH; // int or uInt
    static const casacore::String PARAM_EXPORT_FILENAME; // String
    static const casacore::String PARAM_EXPORT_FORMAT; //String
    static const casacore::String PARAM_EXPORT_VERBOSE; // bool
    static const casacore::String PARAM_EXPORT_RANGE; //String
    static const casacore::String PARAM_EXPORT_HIGHRES; // bool
    static const casacore::String PARAM_EXPORT_DPI; // int
    static const casacore::String PARAM_EXPORT_WIDTH; // int
    static const casacore::String PARAM_EXPORT_HEIGHT; // int
    static const casacore::String PARAM_EXPORT_INTERACTIVE; // bool
    static const casacore::String PARAM_EXPORT_ASYNC;  // bool
    static const casacore::String PARAM_COLORIZE;      // bool
    static const casacore::String PARAM_COLORAXIS;     // string
    static const casacore::String PARAM_CANVASTITLE;    // string
    static const casacore::String PARAM_CANVASTITLEFONT;  // int
    static const casacore::String PARAM_DATA_INDEX;    //int
    static const casacore::String PARAM_XAXISLABEL;    // string
    static const casacore::String PARAM_YAXISLABEL;    // string
    static const casacore::String PARAM_XAXISFONT;    // int
    static const casacore::String PARAM_YAXISFONT;    // int
    static const casacore::String PARAM_SHOWMAJORGRID;  // bool
    static const casacore::String PARAM_SHOWMINORGRID;  // bool
    static const casacore::String PARAM_MAJORCOLOR;    // string
    static const casacore::String PARAM_MINORCOLOR;    // string
    static const casacore::String PARAM_MAJORSTYLE;    // string 
    static const casacore::String PARAM_MINORSTYLE;    // string
    static const casacore::String PARAM_MAJORWIDTH;    // int 
    static const casacore::String PARAM_MINORWIDTH;    // int 
    static const casacore::String PARAM_SHOWLEGEND;    //bool
    static const casacore::String PARAM_LEGENDPOSITION; //string
    static const casacore::String PARAM_XAUTORANGE;    // bool
    static const casacore::String PARAM_XMIN;          // double
    static const casacore::String PARAM_XMAX;          // double
    static const casacore::String PARAM_YAUTORANGE;    // bool
    static const casacore::String PARAM_YMIN;          // double
    static const casacore::String PARAM_YMAX;          // double
    static const casacore::String PARAM_SYMBOL;        // bool
    static const casacore::String PARAM_SYMBOLSHAPE;   // string
    static const casacore::String PARAM_SYMBOLSIZE;    // int
    static const casacore::String PARAM_SYMBOLCOLOR;   // string
    static const casacore::String PARAM_SYMBOLFILL;    // string
    static const casacore::String PARAM_SYMBOLOUTLINE; // bool
    static const casacore::String PARAM_FLAGGEDSYMBOL;        // bool
    static const casacore::String PARAM_FLAGGEDSYMBOLSHAPE;   // string
    static const casacore::String PARAM_FLAGGEDSYMBOLSIZE;    // int
    static const casacore::String PARAM_FLAGGEDSYMBOLCOLOR;   // string
    static const casacore::String PARAM_FLAGGEDSYMBOLFILL;    // string
    static const casacore::String PARAM_FLAGGEDSYMBOLOUTLINE; // bool
    static const casacore::String PARAM_XCONNECTOR;    // string
    static const casacore::String PARAM_TIMECONNECTOR; // bool
    
    
    // </group>

    // DBus method name for getting the log parameters, including: the sink
    // filename (PARAM_FILENAME) and the filter priority (PARAM_PRIORITY).
    // PARAMETERS: none.
    // RETURNS: value (casacore::Record), unless invalid or run asynchronously.
    static const casacore::String METHOD_GETLOGPARAMS;
    
    // DBus method name for setting the log parameters, using a casacore::Record with
    // zero or more of the parameters set (see METHOD_GETLOGPARAMS).
    // PARAMETERS: value (casacore::Record).
    // RETURNS: none.
    static const casacore::String METHOD_SETLOGPARAMS;
    
    // DBus method name for getting the plotms parameters, including: the
    // "clear selections on axes change" flag (PARAM_CLEARSELECTIONS), and the
    // cached image width (PARAM_WIDTH) and height (PARAM_HEIGHT).
    // PARAMETERS: none.
    // RETURNS: value (casacore::Record), unless invalid or run asynchronously.
    static const casacore::String METHOD_GETPLOTMSPARAMS;
    
    // DBus method name for setting the plotms parameters, using a casacore::Record with
    // zero or more of the parameters set (see METHOD_GETPLOTMSPARAMS).
    // PARAMETERS: value (casacore::Record).
    // RETURNS: none.
    static const casacore::String METHOD_SETPLOTMSPARAMS;
    
    // DBus method name for setting the cached image size to the current screen
    // resolution.
    // PARAMETERS: none.
    // RETURNS: none.
    static const casacore::String METHOD_SETCACHEDIMAGESIZETOSCREENRES;
    
    // DBus method name for getting the plot parameters at the given index
    // (PARAM_PLOTINDEX), including: the casacore::MS filename (PARAM_FILENAME), the x
    // axis (PARAM_AXIS_X) and data column (PARAM_DATACOLUMN_X), the y axis
    // (PARAM_AXIS_Y) and data column (PARAM_DATACOLUMN_Y), averaging
    // (PARAM_AVERAGING), selection (PARAM_SELECTION), transformations
    // (PARAM_TRANSFORMATIONS), and calibration (PARAM_CALIBRATION)
    // PARAMETERS: plot index.
    // RETURNS: value (casacore::Record), unless invalid or run asynchronously.
    static const casacore::String METHOD_GETPLOTPARAMS;
    
    // DBus method name for setting the plot parameters at the given index
    // (PARAM_PLOTINDEX), using a casacore::Record with zero or more of the parameters
    // set (see METHOD_GETPLOTPARAMS).
    // PARAMETERS: value (casacore::Record).
    // RETURNS: none.
    static const casacore::String METHOD_SETPLOTPARAMS;
    
    //Sets the export parameters.
    // PARAMETERS: value (casacore::Record).
    // RETURNS: none.
    static const casacore::String METHOD_SETEXPORTPARAMS;

    // DBus method name for getting the flag extension parameters
    // (PARAM_FLAGGING).
    // PARAMETERS: none.
    // RETURNS: value (casacore::Record), unless invalid or run asynchronously.
    static const casacore::String METHOD_GETFLAGGING;
    
    // DBus method name for setting the flag extension parameters.
    // PARAMETERS: flagging value.
    // RETURNS: none.
    static const casacore::String METHOD_SETFLAGGING;
    
    // DBus method names for showing/hiding the window.  Does NOT quit the
    // entire application.
    // PARAMETERS: none.
    // RETURNS: none.
    // <group>
    static const casacore::String METHOD_SHOW;
    static const casacore::String METHOD_HIDE;
    // </group>
    
    // DBus method name for updating the running PlotMS with any attributes
    // that were set with updateImmediately = false.
    // PARAMETERS: none.
    // RETURNS: none.
    static const casacore::String METHOD_UPDATE;
   
    //Existing plots should be removed.
    //PARAMETERS: none.
    //RETURNS: none.
    static const casacore::String METHOD_CLEARPLOTS;

    // DBus method name for quitting the entire application.
    // PARAMETERS: none.
    // RETURNS: none.
    static const casacore::String METHOD_QUIT;
    
    //DBus method name for exporting plot file.

    static const casacore::String METHOD_SAVE;

    //DBus method name for determining if a plot is being drawn
    static const casacore::String METHOD_ISDRAWING;
    
    //is the top widget still shown
    static const casacore::String METHOD_ISCLOSED;

    // DBus method name for locating points in a specified region
    // PARAMETERS: upper left and lower right bounding box coordinates
    // RETURNS: meta data of located points (casacore::Record)
    static const casacore::String METHOD_LOCATEINFO;

    // Returns the name that the plotms in the process with the given ID is (or
    // would be) registered with in the CASA DBus server.
    static casacore::String dbusName(pid_t pid);
    
    // Non-Static //
    
    // Constructor which takes PlotMS parent object.
    PlotMSDBusApp(PlotEngine& plotms);
    
    // Destructor.
    ~PlotMSDBusApp();
    
    // Connects to the DBus server using the dbusName() method with the current
    // process ID.  Returns whether the connection succeeded or not.
    bool connectToDBus( const QString &dbus_name="" );
    
    // Implements PlotMSParametersWatcher::parametersHaveChanged().
    void parametersHaveChanged(
    	const PlotMSWatchedParameters& params, int updateFlag
    );
    
    // Implements PlotMSPlotManagerWatcher::plotsChanged().
    void plotsChanged(const PlotMSPlotManager& manager);
    
protected:
    // Implements QtDBusXmlApp::dbusRunXmlMethod().
    void dbusRunXmlMethod(
    	const casacore::String& methodName, const casacore::Record& parameters,
        casacore::Record& retValue, const casacore::String& callerName, bool isAsync
    );
    
    // Overrides QtDBusXmlApp::dbusXmlReceived() to print the message to the log
    // as needed.
    void dbusXmlReceived(const QtDBusXML& xml);
    
private:
    // Parent PlotMS.
    //PlotMSApp& itsPlotms_;
    PlotEngine& itsPlotms_;
    
    // Set PlotMS parameters that haven't yet been transferred to the current
    // PlotMS.
    PlotMSParameters itsParams_;
    
    // Set PlotMSSinglePlot parameters that haven't yet been transfered to the
    // current PlotMS.
    vector<PlotMSPlotParameters> itsPlotParams_;
    
    // Flag for whether to call update() during show() or not.  This will be
    // true if the user updates something while the GUI is hidden.
    bool itsUpdateFlag_;

    // Helper methods for posting log messages.
    void log(const casacore::String& message);
    void logWarn(const casacore::String& message);
    
    // Adjusts the given plot index to be an acceptable, and returns whether
    // the parameters were resized or not.
    bool plotParameters(int& plotIndex) const;
    
    // Helper for updating.
    bool update();

    // helper for saving
    bool _savePlot(const casacore::Record& parameters);

    // helper for locate
    casacore::Record _locateInfo(const casacore::Record& parameters);

    //Make sure users don't set the plot index to an invalid value.
    bool checkPlotIndex( int index );

};

}

#endif /* PLOTMSDBUSAPP_H_ */
