//# QPCanvasHelpers.qo.h: Helper classes for QPCanvas.
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
#ifndef QPCANVASHELPERS_QO_H_
#define QPCANVASHELPERS_QO_H_

#ifdef AIPS_HAS_QWT

#include <casaqt/QwtPlotter/QPOptions.h>
#include <casaqt/QwtPlotter/QPPlotItem.qo.h>
#include <graphics/GenericPlotter/PlotCanvas.h>

#include <casacore/casa/Quanta/MVAngle.h>

#include <qwt/qwt_legend.h>
#if QWT_VERSION >= 0x060000
#include <qwt/qwt_plot_legenditem.h>
#endif
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_scale_draw.h>

#include <QObject>
#include <QMouseEvent>
#include <QSpacerItem>

namespace casa {

//# Forward Declarations
class QPCanvas;


// Miscellaneous Classes //

// Filter to install on the QwtPlotCanvas to catch mouse move events (which
// are otherwise stolen by other filters in QwtPicker classes).
class QPMouseFilter : public QObject {
    Q_OBJECT
    
public:
    // Constructor which takes the canvas to install itself on.
    QPMouseFilter(QwtPlotCanvas* canvas);
    QPMouseFilter(QWidget* widget);

    // Destructor.
    ~QPMouseFilter();
    
    // Turns mouse tracking on or off on the underlying canvas.
    void turnTracking(bool on);
    
signals:
    // This signal is emitted if tracking is turned on, and when a
    // QMouseEvent is received for a mouse move on the underlying canvas.
    void mouseMoveEvent(QMouseEvent* e);
    
protected:
    // Filter the events on the given object.  Checks only for mouse events
    // on the underlying canvas.
    bool eventFilter(QObject* obj, QEvent* ev);
    
private:
    // Canvas.
    QwtPlotCanvas* m_canvas;
};


// Subclass of QwtScaleDraw to implement additional functionality.  Keeps track
// of which PlotAxisScale is set, and can use reference values instead of
// absolute values (see PlotCanvas::setAxisReferenceValue()).  Additionally,
// for date values it can convert a double in either modified julian seconds or
// modified julian days into a casacore::String representation of the date, using a
// format that can be set.
class QPScaleDraw : public QwtScaleDraw {
public:	
	// Constructor that takes the parent and axis to attach to.
    QPScaleDraw(QwtPlot* parent, QwtPlot::Axis axis);
    
    // Destructor.
    ~QPScaleDraw();
    
    // Gets/Sets the scale used.
    // <group>
    PlotAxisScale scale() const;
    void setScale(PlotAxisScale scale, casacore::uInt base=10);
    // </group>
    
    // Gets/Sets the format used for date values.  See Plotter::dateFormat().
    // <group>
    const casacore::String& dateFormat() const;
    void setDateFormat(const casacore::String& newFormat);
    // </group>
    
    // Gets/Sets the format used for relative date values.  See
    // Plotter::relativeDateFormat().
    // <group>
    const casacore::String& relativeDateFormat() const;
    void setRelativeDateFormat(const casacore::String& newFormat);
    // </group>
    
    // Gets/Sets the reference value.  See PlotCanvas::setAxisReferenceValue().
    // <group>
    bool referenceValueSet() const;    
    double referenceValue() const;    
    void setReferenceValue(bool on, double value = 0);
    // </group>
    
    // Gets/Sets the format used for angle values
    AngleFormat angleFormat() const;
    void setAngleFormat(AngleFormat newFormat);

    // Overrides QwtScaleDraw::label().
    QwtText label(double value) const;
    
    virtual void draw(QPainter* painter, const QPalette& palette) const;

#if QWT_VERSION < 0x060000
    virtual int extent( const QPen& pen, const QFont& font ) const;
#endif

private:
	// finds step between tick values for precision
	int getTickPrecision() const;

	// Parent.
	QwtPlot* m_parent;
	
	// Axis.
	QwtPlot::Axis m_axis;
	
	// Scale.
	PlotAxisScale m_scale;
	
	// Date formats.
	// <group>
	casacore::String m_dateFormat;
	casacore::String m_relativeDateFormat;
	// </group>
	
	// Reference value.
	// <group>
	bool m_referenceSet;
	double m_referenceValue;
	// </group>

	// Angle Format
	AngleFormat m_angleFormat;
	casacore::MVAngle::Format m_mvAngleFormat;
	static const casacore::uInt m_angleFormatDefaultPrecision = 6;
};

// Legend Classes //

// Subclass of QwtLegend to handle outline and background, and be able to draw
// itself using a QPainter.
class QPLegend : public QwtLegend {
public:
    // Constructor which takes optional parent widget.
    QPLegend(QWidget* parent = NULL);
    
    // Destructor.
    ~QPLegend();
    
    
    // Overrides QwtLegend::sizeHint() to take border into account.
    QSize sizeHint() const;
    
    // Gets/Sets the outline.
    // <group>
    const QPLine& line() const;
    const QPen& pen() const;
    void setLine(const PlotLine& line);
    void setLine(const QPen& pen) { setLine(QPLine(pen)); }
    // </group>
    
    // Gets/Sets the background.
    // <group>
    const QPAreaFill& areaFill() const;
    const QBrush& brush() const;
    //The update flag allows an export operation to set the fill without affecting
    //what is shown on the GUI.
    void setAreaFill(const PlotAreaFill& fill, bool update=true);
    void setAreaFill(const QBrush& brush) { setAreaFill(QPAreaFill(brush)); }
    // </group>
    
    // Draws the legend's outline and bacgkround using the given painter on the
    // given rect.  If useQwtPainter is true, QwtPainter is used (which
    // preserves any set QwtMetricsMap).
    void drawOutlineAndBackground(QPainter* painter, const QRect& rect,
            bool useQwtPainter = false);
    
protected:
    // Overrides QWidget::paintEvent() to draw outline and background if
    // needed.
    void paintEvent(QPaintEvent* event);
    
private:
    // Outline.
    QPLine m_line;
    
    // Background.
    QPAreaFill m_areaFill;
};


// Holder for QPLegend that is responsible for its placement.
class QPLegendHolder : public QWidget {
    Q_OBJECT
    
public:
    // Default padding for internal legends.
    static const int DEFAULT_INTERNAL_PADDING;
    
    // Converts from our legend position to Qwt's.
    static QwtPlot::LegendPosition legendPosition(
            PlotCanvas::LegendPosition pos);
    
    
    // Constructor which takes the canvas and initial position.
    QPLegendHolder(QPCanvas* canvas, PlotCanvas::LegendPosition position,
            int padding = DEFAULT_INTERNAL_PADDING);
    
    // Destructor.
    ~QPLegendHolder();

    QPLegend* legend() { return m_legend; }
    // Shows/hides the legend.
    // <group>
    bool legendShown() const;
    void showLegend(bool show = true);
    // </group>
    
    // Gets/Sets the legend position.
    // <group>
    bool isInternal() const;
    PlotCanvas::LegendPosition position() const;
    void setPosition(PlotCanvas::LegendPosition position);
    // </group>
    
    // Gets/Sets the outline.
    // <group>
    const QPLine& line() const { return m_legend->line(); }
    const QPen& pen() const { return m_legend->pen(); }
    void setLine(const PlotLine& line) { m_legend->setLine(line); }
    void setLine(const QPen& pen) { m_legend->setLine(pen); }
    // </group>
    
    // Gets/Sets the background.
    // <group>
    const QPAreaFill& areaFill() const { return m_legend->areaFill(); }
    const QBrush& brush() const { return m_legend->brush(); }
    //The update flag allows an export operation to set the fill without
    //affecting what is shown on the GUI.
    void setAreaFill(const PlotAreaFill& fill, bool update = true ) {
        if (m_legend) m_legend->setAreaFill(fill, update); }
    void setAreaFill(const QBrush& brush) { 
        if (m_legend) m_legend->setAreaFill(brush); }
    // </group>
    
    // Returns the rect for the internal legend, given a canvas rect.
    QRect internalLegendRect(const QRect& canvasRect) const;
    
    // See QPLegend::drawOutlineAndBackground.
    void drawOutlineAndBackground(QPainter* painter, const QRect& rect,
            bool useQwtPainter = false) {
        m_legend->drawOutlineAndBackground(painter, rect, useQwtPainter); }
    
private:
    // Canvas.
    QPCanvas* m_canvas;
    
    // Actual legend.
    QPLegend* m_legend;
    
#if QWT_VERSION >= 0x060000
    QwtPlotLegendItem* m_legendItem;
    void setupLegendItem();
    void setAlignment(PlotCanvas::LegendPosition pos);
#endif

    // Current position.
    PlotCanvas::LegendPosition m_position;
    
    // Spacers for internal legends.
    QSpacerItem* m_spaceTop, *m_spaceLeft, *m_spaceRight, *m_spaceBottom;
    
    // Padding for internal legends.
    int m_padding;
    
    
    // Update the spacers for the position, padding, and line width.
    void updateSpacers();
};


// "Base" Item Classes //

// Abstract superclass for all "base" item classes.  Base items are a special
// type of plot item that are not PlotItems but *are* QwtPlotItems that go
// below the normal items.
class QPBaseItem : public QPLayerItem {
    friend class QPLayeredCanvas;
    
public:
    // Z indexes for known base items.
    // <group>
    static const double BASE_Z_CARTAXIS;
    static const double BASE_Z_GRID;
    // </group>    
    
    
    // Constructor.
    QPBaseItem();
    
    // Destructor.
    virtual ~QPBaseItem();
    
    
    // Implements QPLayerItem::itemChanged().  Calls default definition.
    virtual void itemChanged() { QwtPlotItem::itemChanged(); }
    
    // Implements QPLayerItem::itemDrawCount().
    virtual unsigned int itemDrawCount() const { return 1; }
    
protected:
    // Attached canvas, or NULL for none.
    QPCanvas* m_canvas;
    
    
    // Attaches this item to the given canvas.
    void qpAttach(QPCanvas* canvas);
    
    // Detaches this item from its canvas.
    void qpDetach();
};


// Subclass of QPBaseItem for drawing grids.  Basically just a wrapper for
// QwtPlotGrid.
class QPGrid : public QPBaseItem, public QwtPlotGrid {
public:
    // Constructor.
    QPGrid();
    
    // Destructor.
    ~QPGrid();
    

    // Overrides QwtPlotItem::itemChanged() to use QPBaseItem's definition.
    void itemChanged() { QPBaseItem::itemChanged(); }
    
    // Implements QPLayerItem::shouldDraw().
    bool shouldDraw() const;
    
    // Implements QPLayerItem::itemTitle().
    casacore::String itemTitle() const { return "grid"; }
    
    // Overrides QwtPlotItem::boundingRect() to use QwtPlotGrid's definition.
    QwtDoubleRect boundingRect() const { return QwtPlotGrid::boundingRect(); }
    
    // Overrides QwtPlotItem::updateScaleDiv() to use QwtPlotGrid's definition.
    void updateScaleDiv(const QwtScaleDiv& xDiv, const QwtScaleDiv& yDiv) {
        QwtPlotGrid::updateScaleDiv(xDiv, yDiv); }
    
protected:
    // Implements QPLayerItem::draw_().  Ignores draw index and count.
#if QWT_VERSION >= 0x060000
    void draw_(QPainter* p, const QwtScaleMap& xMap,
                const QwtScaleMap& yMap, const QRectF& drawRect,
                unsigned int drawIndex, unsigned int drawCount) const
#else
    void draw_(QPainter* p, const QwtScaleMap& xMap,
                const QwtScaleMap& yMap, const QRect& drawRect,
                unsigned int drawIndex, unsigned int drawCount) const
#endif
        {
		(void)drawIndex; (void)drawCount;
        QwtPlotGrid::draw(p, xMap, yMap, drawRect); }
};


// Subclass of QPBaseItem for drawing cartesian axes.  See
// http://pyqwt.sourceforge.net/examples/CartesianDemo.py.html .
class QPCartesianAxis : public QPBaseItem {
public:
    // "Master" is the one being drawn; "slave" is the one that the master
    // is attached to.
    QPCartesianAxis(QwtPlot::Axis master, QwtPlot::Axis slave);
    
    // Destructor.
    ~QPCartesianAxis();

    
    // Implements QPLayerItem::shouldDraw().
    bool shouldDraw() const { return true; }
    
    // Implements QPLayerItem::itemTitle().
    casacore::String itemTitle() const { return "cartesian axis"; }
    
protected:
    // Implements QPLayerItem::draw_().  Ignores draw index and count.
#if QWT_VERSION >= 0x060000
    void draw_(QPainter* p, const QwtScaleMap& xMap,
                const QwtScaleMap& yMap, const QRectF& drawRect,
                unsigned int drawIndex, unsigned int drawCount) const;
#else
    void draw_(QPainter* p, const QwtScaleMap& xMap,
                const QwtScaleMap& yMap, const QRect& drawRect,
                unsigned int drawIndex, unsigned int drawCount) const;
#endif

private:
    // Master axis.
    QwtPlot::Axis m_axis;
    
    // Scale draw.
    QwtScaleDraw m_scaleDraw;
};

} // end namespace casa

#endif

#endif /* QPCANVASHELPERS_QO_H_ */
