//# QPCanvas.qo.h: Qwt implementation of generic PlotCanvas class.
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
//# $Id: $
#ifndef QPCANVAS_QO_H_
#define QPCANVAS_QO_H_

#ifdef AIPS_HAS_QWT

#include <casaqt/QwtConfig.h>
#include <graphics/GenericPlotter/PlotOptions.h>
#include <graphics/GenericPlotter/PlotLogger.h>
#include <graphics/GenericPlotter/Plotter.h>
#include <casaqt/QwtPlotter/QPImageCache.h>
#include <casaqt/QwtPlotter/QPLayeredCanvas.qo.h>
#include <casaqt/QwtPlotter/QPExportCanvas.h>
#include <casaqt/QwtPlotter/QPOptions.h>
#include <casaqt/QwtPlotter/QPPlotItem.qo.h>

#include <qwt_plot_picker.h>

#include <QtGui>

#include <vector>

namespace casa {

//# Forward declarations
class QPPlotter;
class AxisListener;

// Implementation of PlotCanvas for the Qwt plotter.  Mainly consists of
// wrappers and containers around a QwtPlot object.
class QPCanvas : public QFrame, public QPExportCanvas {
    Q_OBJECT
    
    friend class QPAxesCache;
    friend class QPDrawThread;
    friend class QPLayeredCanvas;
    friend class QPPlotItem;
    friend class QPPlotter;
    
public:
    // Static //
    
    // Keep a z-order variable to increment for subsequent items.
    static double zOrder;
    
    // Convenient access to class name (QPCanvas).
       static const casacore::String CLASS_NAME;
    
    // Convenient access to "origin" names for logging.
    // <group>
    static const casacore::String DRAW_NAME;

    // </group>
    
    

    // Non-Static //
    
    // Constructor which takes (optional) parent QPPlotter.
    QPCanvas(QPPlotter* parent = NULL);

    // Destructor.
    ~QPCanvas();
    
    
    // Include overloaded methods.
    using PlotCanvas::setBackground;
    using PlotCanvas::setSelectLine;
    using PlotCanvas::setTitleFont;
    using PlotCanvas::showAxes;
    using PlotCanvas::showCartesianAxis;
    using PlotCanvas::setAxisFont;
    using PlotCanvas::setAxisRange;
    using PlotCanvas::setAxesRanges;
    using PlotCanvas::showGrid;
    using PlotCanvas::setGridMajorLine;
    using PlotCanvas::setGridMinorLine;
    using PlotCanvas::setLegendLine;
    using PlotCanvas::setLegendFill;
    using PlotCanvas::setLegendFont;



    // PlotCanvas Methods //
    
    // Implements PlotCanvas::size().
    std::pair<int, int> size() const;
    virtual void setMinimumSize( int width, int height ){
    	QFrame::setMinimumSize( width, height );
    }

    virtual void show(){
    	QFrame::show();
    }

    virtual void hide(){
    	QFrame::setVisible( false );
    }

    // Implements PlotCanvas::title().
    casacore::String title() const;

    // Implements PlotCanvas::setTitle().
    void setTitle(const casacore::String& title);
    
    // Implements PlotCanvas::titleFont().
    PlotFontPtr titleFont() const;
    
    // Implements PlotCanvas::setTitleFont().
    void setTitleFont(const PlotFont& font);
    
    // Implements PlotCanvas::background().
    PlotAreaFillPtr background() const;
    PlotAreaFillPtr defaultBackground() const;
    
    // Implements PlotCanvas::setBackground().
    void setBackground(const PlotAreaFill& areaFill);

    // Implements PlotCanvas::cursor().
    PlotCursor cursor() const;
    
    // Implements PlotCanvas::setCursor().
    void setCursor(PlotCursor cursor);
    
    // Implements PlotCanvas::refresh().
    // <group>
    void refresh();
    void refresh(int drawLayersFlag);
    // </group>
    
    // Implements PlotCanvas::isQWidget().
    bool isQWidget() const { return true; }

    
    // Implements PlotCanvas::shownAxes().
    // Returns a bit set (really an int) of bits defined by PlotAxis enum or'd together
    PlotAxisBitset shownAxes() const;

    // Implements PlotCanvas::showAxes().
    // Takes a bit set, as an int, of bits defined by PlotAxis enum or'd together
    void showAxes(PlotAxisBitset axes);
    
    // Implements PlotCanvas::axisScale().
    PlotAxisScale axisScale(PlotAxis axis) const;

    // Implements PlotCanvas::setAxisScale().
    void setAxisScale(PlotAxis axis, PlotAxisScale scale, casacore::uInt base=10);

    // Implements PlotCanvas::setAxisScaleSortDirection().
    bool setAxisScaleSortDirection(PlotAxis axis, SortDirection direction);
    std::pair<bool,SortDirection> axisScaleSortDirection(PlotAxis axis) const;

    // Sets/gets the angle format of the scale for the given axis
    void setAxisScaleAngleFormat(PlotAxis axis, AngleFormat format);
    AngleFormat axisScaleAngleFormat(PlotAxis axis) const;

    // Implements PlotCanvas::axisReferenceValueSet().
    bool axisReferenceValueSet(PlotAxis axis) const;
    
    // Implements PlotCanvas::axisReferenceValueValue().
    double axisReferenceValue(PlotAxis axis) const;
    
    // Implements PlotCanvas::setAxisReferenceValue().
    void setAxisReferenceValue(PlotAxis axis, bool on, double value = 0);
    
    // Implements PlotCanvas::cartesianAxisShown().
    bool cartesianAxisShown(PlotAxis axis) const;

    // Implements PlotCanvas::showCartesianAxis().
    void showCartesianAxis(PlotAxis mirrorAxis, PlotAxis secondaryAxis,
            bool show = true, bool hideNormalAxis = true);
    
    // Implements PlotCanvas::axisLabel().
    casacore::String axisLabel(PlotAxis axis) const;

    // Implements PlotCanvas::setAxisLabel().
    void setAxisLabel(PlotAxis axis, const casacore::String& title);

    // Implements PlotCanvas::axisFont().
    PlotFontPtr axisFont(PlotAxis axis) const;
    
    // Implements PlotCanvas::setAxisFont().
    void setAxisFont(PlotAxis axis, const PlotFont& font);

    // Implements PlotCanvas::colorBarShown().
    bool colorBarShown(PlotAxis axis = Y_RIGHT) const;

    // Implements PlotCanvas::showColorBar().
    void showColorBar(bool show = true, PlotAxis axis = Y_RIGHT);

   
    // Implements PlotCanvas::axisRange().
    prange_t axisRange(PlotAxis axis) const;

    // Implements PlotCanvas::setAxisRange().
    void setAxisRange(PlotAxis axis, double from, double to);
    
    // Implements PlotCanvas::invertAxis().
    void invertAxis(PlotAxis axis);

    // Overrides PlotCanvas::setAxesRanges().
    void setAxesRanges(PlotAxis xAxis, double xFrom, double xTo,
                       PlotAxis yAxis, double yFrom, double yTo);

    // Implements PlotCanvas::axesAutoRescale().
    bool axesAutoRescale() const;

    // Implements PlotCanvas::setAxesAutoRescale().
    void setAxesAutoRescale(bool autoRescale = true);

    // Implements PlotCanvas::rescaleAxes().
    void rescaleAxes();

    // Implements PlotCanvas::axesRatioLocked().
    bool axesRatioLocked() const;
    
    // Implements PlotCanvas::setAxesRatioLocked().
    void setAxesRatioLocked(bool locked = true);
    
       
    // Implements PlotCanvas::cachedAxesStackSizeLimit().
    int cachedAxesStackSizeLimit() const;

    // Implements PlotCanvas::setCachedAxesStackSizeLimit().
    void setCachedAxesStackSizeLimit(int sizeInKilobytes);
    
    // Overrides PlotCanvas::cachedAxesStackImageSize().
    std::pair<int, int> cachedAxesStackImageSize() const;
    
    // Overrides PlotCanvas::setCachedAxesStackImageSize().
    void setCachedAxesStackImageSize(int width, int height);


    // Implements PlotCanvas::plotItem().  If the given items is NOT an
    // instance of a QPPlotItem, a copy of the given items is made.  The
    // original item is NOT kept by the canvas, so any subsequent changes to
    // the original items will not be reflected on the canvas.
    bool plotItem(PlotItemPtr item, PlotCanvasLayer layer = MAIN);

    // Implements PlotCanvas::allPlotItems().
    std::vector<PlotItemPtr> allPlotItems() const;
    
    // Implements PlotCanvas::layerPlotItems().
    std::vector<PlotItemPtr> layerPlotItems(PlotCanvasLayer layer) const;

    // Overrides PlotCanvas::numPlotItems().
    unsigned int numPlotItems() const;

    // Overrides PlotCanvas::numLayerPlotItems().
    unsigned int numLayerPlotItems(PlotCanvasLayer layer) const;
    
    // Implements PlotCanvas::removePlotItems().
    void removePlotItems(const std::vector<PlotItemPtr>& items);
    
    // Overrides PlotCanvas::clearPlotItems().
    void clearPlotItems();
    
    // Overrides PlotCanvas::clearPlots().
    void clearPlots();
    
    // Overrides PlotCanvas::clearLayer().
    void clearLayer(PlotCanvasLayer layer);

    
    // Implements PlotCanvas::holdDrawing().
    void holdDrawing();
    
    // Implements PlotCanvas::releaseDrawing().
    void releaseDrawing();
    
    // Implements PlotCanvas::drawingIsHeld().
    bool drawingIsHeld() const;

        
    // Implements PlotCanvas::setSelectLineShown().
    void setSelectLineShown(bool shown = true);
    
    // Implements PlotCanvas::selectLine().
    PlotLinePtr selectLine() const;
    
    // Implements PlotCanvas::setSelectLine().
    void setSelectLine(const PlotLine& line);

    
    // Implements PlotCanvas::gridShown().
    bool gridShown(bool* xMajor = NULL, bool* xMinor = NULL,
            bool* yMajor = NULL, bool* yMinor = NULL) const;
    
    // Implements PlotCanvas::showGrid().
    void showGrid(bool xMajor, bool xMinor, bool yMajor,bool yMinor);
    
    // Implements PlotCanvas::gridMajorLine().
    PlotLinePtr gridMajorLine() const;

    // Implements PlotCanvas::setGridMajorLine().
    void setGridMajorLine(const PlotLine& line);

    // Implements PlotCanvas::gridMinorLine().
    PlotLinePtr gridMinorLine() const;

    // Implements PlotCanvas::setGridMinorLine().
    void setGridMinorLine(const PlotLine& line);

    
    // Implements PlotCanvas::legendShown().
    bool legendShown() const;
    
    // Implements PlotCanvas::showLegend().
    void showLegend(bool on = true, LegendPosition position = EXT_BOTTOM);
    
    // Implements PlotCanvas::legendPosition().
    LegendPosition legendPosition() const;
    
    // Implements PlotCanvas::setLegendPosition().
    void setLegendPosition(LegendPosition position);
    
    // Implements PlotCanvas::legendLine().
    PlotLinePtr legendLine() const;
    
    // Implements PlotCanvas::setLegendLine().
    void setLegendLine(const PlotLine& line);
    
    // Implements PlotCanvas::legendFill().
    PlotAreaFillPtr legendFill() const;
    
    // Implements PlotCanvas::setLegendFill().
    void setLegendFill(const PlotAreaFill& area);
    
    // Implements PlotCanvas::legendFont().
    PlotFontPtr legendFont() const;
    
    // Implements PlotCanvas::setLegendFont().
    void setLegendFont(const PlotFont& font);


    // Implements PlotCanvas::autoIncrementColors().
    bool autoIncrementColors() const;

    // Implements PlotCanvas::setAutoIncrementColors().
    void setAutoIncrementColors(bool autoInc = true);

    // Implements PlotCanvas::exportToFile().
    bool exportToFile(const PlotExportFormat& format);

    // Implements PlotCanvas::fileChooserDialog().
    casacore::String fileChooserDialog(const casacore::String& title = "File Chooser",
            const casacore::String& directory = "");
    
    // Implements PlotCanvas::dateFormat().
    const casacore::String& dateFormat() const;
    
    // Implements PlotCanvas::setDateFormat().
    void setDateFormat(const casacore::String& dateFormat);
    
    // Implements PlotCanvas::relativeDateFormat().
    const casacore::String& relativeDateFormat() const;
    
    // Implements PlotCanvas::setRelativeDateFormat().
    void setRelativeDateFormat(const casacore::String& dateFormat);

    // Implements PlotCanvas::convertCoordinate().
    PlotCoordinate convertCoordinate(const PlotCoordinate& coord,
           PlotCoordinate::System newSystem = PlotCoordinate::WORLD) const;

    // Implements PlotCanvas::textWidthHeightDescent().
    std::vector<double> textWidthHeightDescent(const casacore::String& text,
            PlotFontPtr font) const;
    
    // Implements PlotCanvas::implementation().
    int implementation() const { return Plotter::QWT; }
    
    // Implements PlotCanvas::implementationFactory().
    PlotFactory* implementationFactory() const;
    

    // QPCanvas Methods //
    
    // Provides access to the underlying QPLayeredCanvas.
    // <group>
    QPLayeredCanvas& asQwtPlot();
    const QPLayeredCanvas& asQwtPlot() const;
    // </group>
    
    // Provides access to the QwtPlotPicker used for selection events.
    QwtPlotPicker& getSelecter();
    
    // Reinstalls the tracker filter (in case another QwtPlotPicker is added to
    // the QwtPlotCanvas).
    void reinstallTrackerFilter();
    
    // Overrides QWidget::sizeHint() to return an invalid size.
    QSize sizeHint() const;
    
    // Overrides QWidget::minimumSizeHint() to return an invalid size.
    QSize minimumSizeHint() const;

    //Overriden to take into account the screen size available to this canvas
    //based on how many other canvases, etc are being displayed.  Without this
    //method, the draw thread hangs when the size available to the canvas is
    //less that the minimum size hint.
    virtual void setMinimumSizeHint( int width, int height );

    virtual void setCommonAxes( bool commonX, bool commonY );
    void addAxisListener( AxisListener* listener );
    void clearAxisListeners();

    virtual bool isDrawing();
protected:
    // Sets the parent QPPlotter to the given.  This MUST be done when a canvas
    // is added to the plotter so that it can use the plotter's logger if
    // needed.
    void setQPPlotter(QPPlotter* parent);
    
    // Returns the parent's logger.
    virtual PlotLoggerPtr logger() const;
    
    // See QPPlotter::logObject().  If called before setQPPlotter() is called,
    // creates a queue that is then posted when setQPPlotter() is called.
    void logObject(const casacore::String& className, void* address, bool creation,
            const casacore::String& message = casacore::String());
    
    // See QPPlotter::logMethod().  Does NOT queue messages if called before
    // setQPPlotter() is called.
    void logMethod(const casacore::String& className, const casacore::String& methodName,
            bool entering, const casacore::String& message = casacore::String());
    
    // Provides access to the cached axes stack.
    // <group>
    QPAxesCache& axesCache();
    const QPAxesCache& axesCache() const;
    // </group>
    
    // For catching Qt press events.
    void mousePressEvent(QMouseEvent* event);
    
    // For catching Qt click and release events.
    void mouseReleaseEvent(QMouseEvent* event);
    
    // For catching Qt double-click events.
    void mouseDoubleClickEvent(QMouseEvent* event);
    
    // For catching Qt key events.
    void keyReleaseEvent(QKeyEvent* event);
    
    // For catching Qt scroll wheel events.
    void wheelEvent(QWheelEvent* event);
    
    // For catching Qt resize events.
    void resizeEvent(QResizeEvent* event);

    bool isThreading() const;


private:
    // Parent QPPlotter.
    QPPlotter* m_parent;
    
    // Queued log messages before parent is set.
    std::vector<PlotLogObject> m_queuedLogs;
    
    // Main QwtPlot object.
    QPLayeredCanvas m_canvas;

    // Main-layer plot items.
    std::vector<std::pair<PlotItemPtr, QPPlotItem*> > m_plotItems;
    
    // Annotation-layer plot items.
    std::vector<std::pair<PlotItemPtr, QPPlotItem*> > m_layeredItems;
    
    // Scale draws.
    QPScaleDraw* m_scaleDraws[QwtPlot::axisCnt];
    
    // Whether the axes ratio is locked or not.
    bool m_axesRatioLocked;

    bool isCommonAxis( PlotAxis axis ) const;
    bool commonX;
    bool commonY;
    
    // Used for recalculating axes ranges if the ratio is locked.
    std::vector<double> m_axesRatios;
    
    // Cached axes stack.
    QPAxesCache m_stackCache;

    // Whether auto-increment colors is turned on or not.
    bool m_autoIncColors;
    
    // Used auto-incremented colors.
    std::vector<int> m_usedColors;
    
    // Picker used for select events.
    QwtPlotPicker m_picker;
    
    // Filter used for mouse move events.  Has to initialize after the picker
    // to be first in the filter.
    QPMouseFilter m_mouseFilter;
    
    // Legend, and properties.
    // <group>
    QPLegendHolder* m_legend;
    QPFont m_legendFont;
    bool m_legendFontSet;
    // </group>
    
    // Flag for whether we're in mouse dragging mode or not.
    bool m_inDraggingMode;
    
    QList<AxisListener*> axisListeners;

    /*
    // For catching single vs. double clicks.
    // <group>
    bool m_ignoreNextRelease;
    QTimer m_timer;
    QMouseEvent* m_clickEvent;
    // </group>
     */
    
    // Date formats.
    // <group>
    casacore::String m_dateFormat;
    casacore::String m_relativeDateFormat;
    // </group>
    
       
    // Converts the given Qt global position to a pixel PlotCoordinate.
    // <group>
    PlotCoordinate globalPosToPixelCoord(int x, int y);
    PlotCoordinate globalPosToPixelCoord(QMouseEvent* event) {
        return globalPosToPixelCoord(event->globalX(), event->globalY()); }
    PlotCoordinate globalPosToPixelCoord(QWheelEvent* event) {
        return globalPosToPixelCoord(event->globalX(), event->globalY()); }
    // </group>
    

    virtual bool print( QPrinter& printer );
    virtual bool print(  QPainter* painter, PlotAreaFillPtr paf, double widthWidth,
    		double widgetHeight, int externalAxisWidth, int externalAxisHeight,
    		int rowIndex, int colIndex, QRect imageRect );
    virtual bool printRect( QPainter* painter, QRect rect);

    virtual int canvasWidth() const{
    	return width();
    }
    virtual int canvasHeight() const {
    	return height();
    }
    virtual const QPalette& palette() const {
    	return asQwtPlot().palette();
    }
    virtual QPalette::ColorRole backgroundRole() const{
    	return asQwtPlot().backgroundRole();
    }


    QImage  grabImageFromCanvas( const PlotExportFormat& format );

    
    // Converts between axes bitset flags (1,2,4,8 in PlotAxis and vector indices (0-3).
    // (Does not deal with bitsets for combinations of axes, only single axis)
    // <group>
    static unsigned int axisIndex(PlotAxis a);   
    static PlotAxis axisIndex(unsigned int i);
    // </group>
    
    const QwtScaleDiv* getAxisScaleDiv(int axisId) const;
	// set time scale to even hh:mm
	void setTimeScaleDiv(PlotAxis axis, double from, double to);
	
    QSize minSizeHint;

    // default background
    PlotAreaFillPtr defaultBackground_;


private slots:    
    // For when the selecter has selected a region; emit a PlotSelectEvent.
    // Include Qwt5 and Qwt6 definitions for moc;
    // signal comes from QwtPlotPicker
    void regionSelected(const QwtDoubleRect&);
    void regionSelected2(const QRectF&);
    
    // For catching single vs. double clicks.
    // void timeout();
    
    // For catching mouse move events from the filter.
    void trackerMouseEvent(QMouseEvent* event);

    void enableAxis( QwtPlot::Axis axis, bool enable );
};

}

#endif

#endif /*QPCANVAS_QO_H_*/
