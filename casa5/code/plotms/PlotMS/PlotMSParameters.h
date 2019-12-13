//# PlotMSParameters.h: Parameter classes for plotms.
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
//#
//# $Id:  $
#ifndef PLOTMSPARAMETERS_H_
#define PLOTMSPARAMETERS_H_

#include <plotms/PlotMS/PlotMSConstants.h>
#include <plotms/PlotMS/PlotMSWatchedParameters.h>

namespace casa {

// Subclass of PlotMSWatchedParameters that hold parameters for the whole
// plotter.  These parameters include:
// * log file name
// * log events flag
// * log minimum priority filter
// * whether to clear any selections when axes are changed or not
// * width and height for the cached image
class PlotMSParameters : public PlotMSWatchedParameters {
public:
    // Static //
    
    // Update flags.
    // <group>
    static const int UPDATE_LOG;
    static const int UPDATE_PLOTMS_OPTIONS;
    // </group>
    
    // Gets/Sets the file chooser history limit.  (See QtFileDialog.)  Static
    // parameter.
    // <group>
    static int chooserHistoryLimit();
    static void setChooserListoryLimit(int histLimit);
    // </group>
    
    
    // Non-Static //
    
    // Constructor, with default values for parameters.
    PlotMSParameters(const casacore::String& logFilename = PMS::DEFAULT_LOG_FILENAME,
            int logEvents = PMS::DEFAULT_LOG_EVENTS,
            casacore::LogMessage::Priority logPriority = PMS::DEFAULT_LOG_PRIORITY,
            bool clearSelections = PMS::DEFAULT_CLEAR_SELECTIONS,
            int cachedImageWidth = PMS::DEFAULT_CACHED_IMAGE_WIDTH,
            int cachedImageHeight = PMS::DEFAULT_CACHED_IMAGE_HEIGHT,
            int rowCount = PMS::DEFAULT_GRID_ROWS,
            int colCount = PMS::DEFAULT_GRID_COLS);
    
    // Copy constructor.  See operator=().
    PlotMSParameters(const PlotMSParameters& copy);
    
    // Destructor.
    ~PlotMSParameters();

    
    // Gets/Sets the log sink location/filename.
    // <group>
    casacore::String logFilename() const;
    void setLogFilename(const casacore::String& filename);
    // </group>
    
    // Returns the current log events.
    int logEvents() const;
    
    // Returns the current log minimum priority.
    casacore::LogMessage::Priority logPriority() const;
    
    // Sets the current log filter.
    void setLogFilter(int logEvents, casacore::LogMessage::Priority priority);
    
    // Gets/Sets whether any selections are cleared when plot axes are changed
    // or not.
    // <group>
    bool clearSelectionsOnAxesChange() const;
    void setClearSelectionsOnAxesChange(bool flag);
    // </group>
    
    // Gets/Sets the cached image size.  See
    // PlotCanvas::cachedAxesStackImageSize().
    // <group>
    std::pair<int, int> cachedImageSize() const;
    void setCachedImageSize(int width, int height);
    // </group>
    
    // Sets the cached image size to the current screen resolution.
    void setCachedImageSizeToResolution();
    


    //<group>
    bool setGridSize( int rows, int cols );
    int getRowCount() const;
    int getColCount() const;
    void setRowCount( int rowCount );
    void setColCount( int colCount );
    //</group>

    // Copy operator.
    PlotMSParameters& operator=(const PlotMSParameters& copy);
    
    
    // Implements PlotMSWatchedParameters::equals().  Will return false if the
    // other parameters are not of type PlotMSParameters.
    bool equals(const PlotMSWatchedParameters& other, int updateFlags) const;
    
private:
    // Log filename.
    casacore::String itsLogFilename_;
    
    // Log events flag.
    int itsLogEvents_;
    
    // Log minimum priority.
    casacore::LogMessage::Priority itsLogPriority_;
    
    // Clear selections on axes change flag.
    bool itsClearSelectionsOnAxesChange_;
    
    // Cached image sizes.
    int itsCachedImageWidth_, itsCachedImageHeight_;

    int rowCount;
    int colCount;


};

//Removal of compile warnings for unused variables.
class DummyClass {
private:
	static const int dummyDraw;
	static const int dummyData;
	static const int dummyCache;
	static const int dummyAxes;
	static const int dummyCanvas;
	static const int dummyDisplay;
	static const int dummyIter;
};

class DirectionAxisParams {
public:
	DirectionAxisParams(
			PMS::CoordSystem coordSystem=PMS::DEFAULT_COORDSYSTEM,
			PMS::InterpMethod interpMethod=PMS::DEFAULT_INTERPMETHOD
			);
	PMS::CoordSystem getCoordSystem() const;
	PMS::InterpMethod getInterpMethod() const;
	friend bool operator<(const DirectionAxisParams & p1, const DirectionAxisParams & p2);
	friend bool operator!=(const DirectionAxisParams & p1, const DirectionAxisParams & p2);
private:
	PMS::CoordSystem coordSystem_;
	PMS::InterpMethod interpMethod_;
};

}

#endif /* PLOTMSPARAMETERS_H_ */
