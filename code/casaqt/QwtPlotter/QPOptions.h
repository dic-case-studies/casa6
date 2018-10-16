//# QPOptions.h: Qwt implementation of generic PlotOption classes.
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
#ifndef QPOPTIONS_QO_H_
#define QPOPTIONS_QO_H_

#ifdef AIPS_HAS_QWT

#include <graphics/GenericPlotter/PlotOperation.h>
#include <graphics/GenericPlotter/PlotOptions.h>

#include <QColor>
#include <QFont>
#include <QMutex>

#include <qwt_color_map.h>
#include <qwt_plot.h>
#include <qwt_symbol.h>

namespace casa {

// Implementation of PlotColor for Qwt plotter, using QColor.
class QPColor : public PlotColor {
public:
    // Defaults to QColor::QColor().
    QPColor();

    // Constructor for the given hex value or name.
    QPColor(const casacore::String& name);
    
    // Constructor using a QColor.
    QPColor(const QColor& color);
    
    // Copy constructor for generic PlotColor.
    // <group>
    QPColor(const PlotColor& color);
    QPColor(const PlotColorPtr color);
    // </group>
    
    // Copy constructor.
    QPColor(const QPColor& color);
    
    // Destructor.
    ~QPColor();
    
    
    // Include overloaded methods.
    using PlotColor::operator=;
    
    
    // PlotColor Methods //
    
    // Implements PlotColor::asHexadecimal().
    casacore::String asHexadecimal() const;
       
    // Implements PlotColor::asName().
    casacore::String asName() const;    

    // Implements PlotColor::asHexadecimalOrName().
    void setAsHexadecimalOrName(const casacore::String& str);
    
    // Implements PlotColor::alpha().
    double alpha() const;
    
    // Implements PlotColor::setAlpha().
    void setAlpha(double a);
    
    
    // QPColor Methods //
    
    // Provides access to the underlying QColor.
    // <group>
    const QColor& asQColor() const;
    void setAsQColor(const QColor& color);
    // </group>
    
private:
    QColor m_color; // Color
};


// Implementation of PlotFont for Qwt plotter, using QFont and QPColor.
class QPFont : public PlotFont {
public:
    // Defaults to QFont::QFont() and QPColor::QPColor().
    QPFont();
    
    // Constructor that takes a QFont and QColor.
    QPFont(const QFont& font, const QColor& color);
    
    // Copy constructor for generic PlotFont.
    // <group>
    QPFont(const PlotFont& font);
    QPFont(const PlotFontPtr font);
    // </group>
    
    // Copy constructor.
    QPFont(const QPFont& font);
    
    // Destructor.
    ~QPFont();
    
    
    // Include overloaded methods.
    using PlotFont::operator=;
    using PlotFont::setColor;
    
    
    // PlotFont Methods //
    
    // Implements PlotFont::pointSize().
    double pointSize() const;

    // Implements PlotFont::setPointSize().
    void setPointSize(double size);
    
    // Implements PlotFont::pixelSize().
    int pixelSize() const;
    
    // Implements PlotFont::setPixelSize().
    void setPixelSize(int size);
    
    // Implements PlotFont::fontFamily().
    casacore::String fontFamily() const;
    
    // Implements PlotFont::setFontFamily().
    void setFontFamily(const casacore::String& font);
    
    // Implements PlotFont::color().
    PlotColorPtr color() const;
    
    // Implements PlotFont::setColor().
    void setColor(const PlotColor& color);
    
    // Implements PlotFont::italics().
    bool italics() const;
    
    // Implements PlotFont::setItalics().
    void setItalics(bool i = true);
    
    // Implements PlotFont::bold().
    bool bold() const;
    
    // Implements PlotFont::setBold().
    void setBold(bool b = true);
    
    // Implements PlotFont::underline().
    bool underline() const;
    
    // Implements PlotFont::setUnderline().
    void setUnderline(bool u = true);
    
    
    // QPFont Methods //
    
    // Provides access to the underlying QFont.
    // <group>
    const QFont& asQFont() const;
    void setAsQFont(const QFont& font);
    // </group>
    
    // Provides access to the underlying QColor through QPColor.
    // <group>
    const QColor& asQColor() const;
    void setAsQColor(const QColor& color);
    // </group>
    
private:
    QFont m_font;    // Font
    QPColor m_color; // Color
};


// Implementation of PlotLine for Qwt plotter, using QPen.
class QPLine : public PlotLine {
public:
    // Defaults to QPen::QPen().
    QPLine();
    
    // Constructor that takes a QPen.
    QPLine(const QPen& pen);
    
    // Copy constructor for generic PlotLine.
    // <group>
    QPLine(const PlotLine& copy);
    QPLine(const PlotLinePtr copy);
    // </group>
    
    // Destructor.
    ~QPLine();
    
    
    // Include overloaded methods.
    using PlotLine::operator=;
    using PlotLine::setColor;
    
    
    // PlotLine Methods //
    
    // Implements PlotLine::width().
    double width() const;
    
    // Implements PlotLine::setWidth().
    void setWidth(double width);
    
    // Implements PlotLine::style().
    Style style() const;
    
    // Implements PlotLine::setStyle().
    void setStyle(Style style);
    
    // Implements PlotLine::color().
    PlotColorPtr color() const;
    
    // Implements PlotLine::setColor().
    void setColor(const PlotColor& color);
    
    
    // QPLine Methods //
    
    // Provides access to the underlying QPen.
    // <group>
    const QPen& asQPen() const;
    void setAsQPen(const QPen& pen);
    // </group>
    
private:
    QPen m_pen;      // Pen
    QPColor m_color; // Color (for convenience)
};


// Implementation of PlotAreaFill for Qwt plotter, using QBrush.
class QPAreaFill : public PlotAreaFill {
public:
    // Defaults to QBrush::QBrush().
    QPAreaFill();
    
    // Constructor that takes a QBrush.
    QPAreaFill(const QBrush& brush);
    
    // Copy constructor for generic PlotAreaFill.
    // <group>
    QPAreaFill(const PlotAreaFill& copy);
    QPAreaFill(const PlotAreaFillPtr copy);
    // </group>
    
    // Destructor.
    ~QPAreaFill();
    
    
    // Include overloaded methods.
    using PlotAreaFill::operator=;
    using PlotAreaFill::setColor;
    
    
    // PlotAreaFill Methods //
    
    // Implements PlotAreaFill::color().
    PlotColorPtr color() const;
    
    // Implements PlotAreaFill::setColor().
    void setColor(const PlotColor& color);
    
    // Implements PlotAreaFill::pattern().
    Pattern pattern() const;
    
    // Implements PlotAreaFill::setPattern().
    void setPattern(Pattern pattern);
    
    
    // QPAreaFill Methods //
    
    // Provides access to the underlying QBrush.
    // <group>
    const QBrush& asQBrush() const;    
    void setAsQBrush(const QBrush& brush);
    // </group>
    
private:
    QBrush m_brush;  // Brush.
    QPColor m_color; // color (for convenience).
};


// Implementation of QPSymbol for Qwt plotter, using QwtSymbol.
class QPSymbol : public PlotSymbol, public QwtSymbol {
public:
    // Defaults to QwtSymbol::QwtSymbol().
    QPSymbol();
   
#if QWT_VERSION >= 0x060000
    // Implement QwtSymbol constructor
    QPSymbol(QwtSymbol::Style style, const QBrush & brush, 
	const QPen & pen, const QSize & sz);
#else
    // Copy constructor for QwtSymbol (private in Qwt6)
    QPSymbol(const QwtSymbol& s);
#endif
    
    // Copy constructor for generic PlotSymbol.
    // <group>
    QPSymbol(const PlotSymbol& copy);
    QPSymbol(const PlotSymbolPtr copy);
    // </group>
    
    // Destructor.
    ~QPSymbol();
    
    
    // Include overloaded methods.
    using PlotSymbol::operator=;
    using PlotSymbol::setSize;
    using PlotSymbol::setAreaFill;
    using PlotSymbol::setLine;
    
    
    // PlotSymbol Methods //
    
    // Implements PlotSymbol::size().
    psize_t size() const;
    
    // Implements PlotSymbol::setSize().  Note that the underlying QwtSymbol
    // only supports int sizes.
    void setSize(double width, double height, bool heightIsPixel = true);
    
    // Implements PlotSymbol::heightIsPixel().
    bool heightIsPixel() const;
    
    // Implements PlotSymbol::setHeightIsPixel().
    void setHeightIsPixel(bool pixel = true);
    
    // Implements PlotSymbol::symbol().
    Symbol symbol() const;
    
    // Implements PlotSymbol::symbolChar().
    char symbolChar() const;
    
    // Implements PlotSymbol::symbolUChar().
    unsigned short symbolUChar() const;
    
    // Implements PlotSymbol::setSymbol().
    // <group>
    void setSymbol(Symbol symbol);
    void setSymbol(char c);
    // </group>
    
    // Implements PlotSymbol::setUSymbol().
    void setUSymbol(unsigned short unicode);
    
    // Implements PlotSymbol::line().
    PlotLinePtr line() const;
    
    // Implements PlotSymbol::setLine().
    void setLine(const PlotLine& line);
    
    // Implements PlotSymbol::areaFill().
    PlotAreaFillPtr areaFill() const;
    
    // Implements PlotSymbol::setAreaFill().
    void setAreaFill(const PlotAreaFill& fill);
    
    
    // QwtSymbol Methods //
    
    // Overrides QwtSymbol::draw().  Assumes pen and brush have been set
    // beforehand.
    void draw(QPainter* p, const QRect& r) const;
    
#if QWT_VERSION >= 0x060000
    inline void drawSymbols(QPainter* p, const QPointF* points, int numPoints) const
        { QwtSymbol::drawSymbols(p, points, numPoints); };
#endif
    // Overrides QwtSymbol::clone().
    //QwtSymbol* clone() const;
    
    
    // QPSymbol Methods //
    
    // Provides access to the underlying QPen and QBrush in a safe way,
    // since QPainter::drawPoint works differently than the other draw
    // methods.  WARNING: the pixel pen is NOT updated if you call
    // QwtSymbol::setPen() AFTER a QPSymbol::setSymbol() or a
    // QPSymbol::setLine() or a QPSymbol::setAreaFill().
    // <group>
    const QPen& drawPen() const;
    const QBrush& drawBrush() const;
    // </group>
    
private:
    Symbol m_style; // Style
    QChar m_char;   // Character
    
    bool m_heightIsPixel; // Flag for whether height is in pixels or points.
    
    QPen m_pixelPen;     // Pen used to draw pixels.
    QBrush m_pixelBrush; // Brush used to draw pixels (empty).
    
    // Points to either QwtSymbol::pen() or the pixel pen.
    const QPen* m_drawPen;
    
    // Points to either QwtSymbol::brush() or the pixel brush.
    const QBrush* m_drawBrush;
};


// Implementation of PlotMutex for the qwt plotter.  Very thin layer on top of
// QMutex.
class QPMutex : public PlotMutex {
public:
    // Constructor.
    QPMutex();
    
    // Destructor.
    ~QPMutex();
    
    
    // Implements PlotMutex::lock().
    void lock();
    
    // Implements PlotMutex::unlock().
    void unlock();
    
    // Implements PlotMutex::tryLock().
    bool tryLock();
    
private:
    // Mutex.
    QMutex m_mutex;
};


// Class containing useful conversions between the generic plotting stuff and
// the Qt/Qwt plotting stuff.
class QPOptions {
public:
    // Converts between Qwt's axis enum and CASA's.
    // <group>
    static QwtPlot::Axis axis(PlotAxis a);
    static PlotAxis axis(QwtPlot::Axis a);
    static bool isAxisX( PlotAxis a );

    // </group>
    
    // Converts between Qt's cursors and CASA's.
    // <group>
    static Qt::CursorShape cursor(PlotCursor c) {
        switch(c) {
        case NORMAL_CURSOR: return Qt::ArrowCursor;
        case HAND_OPEN:     return Qt::OpenHandCursor;
        case HAND_CLOSED:   return Qt::ClosedHandCursor;
        case CROSSHAIR:     return Qt::CrossCursor;
        case WAIT:          return Qt::WaitCursor;
        case TEXT:          return Qt::IBeamCursor;
        
        default: return Qt::ArrowCursor;
        }
    }    
    static PlotCursor cursor(Qt::CursorShape c) {
        switch(c) {
        case Qt::ArrowCursor:                     return NORMAL_CURSOR;
        case Qt::OpenHandCursor:                  return HAND_OPEN;
        case Qt::ClosedHandCursor:                return HAND_CLOSED;
        case Qt::CrossCursor:                     return CROSSHAIR;
        case Qt::WaitCursor: case Qt::BusyCursor: return WAIT;
        case Qt::IBeamCursor:                     return TEXT;
        
        default:                                  return NORMAL_CURSOR;
        }
    }
    // </group>
    
    // Converts between Qt's line style and CASA's.
    // <group>
    static Qt::PenStyle linestyle(PlotLine::Style s) {
        switch(s) {
        // supported
        case PlotLine::SOLID: return Qt::SolidLine;
        case PlotLine::DASHED: return Qt::DashLine;
        case PlotLine::DOTTED: return Qt::DotLine;
        
        // unsupported
        default: return Qt::NoPen;
        }
    }    
    static PlotLine::Style linestyle(Qt::PenStyle s) {
        switch(s) {
        // supported
        case Qt::SolidLine: return PlotLine::SOLID;
        case Qt::DashLine: return PlotLine::DASHED;
        case Qt::DotLine: return PlotLine::DOTTED;
        
        // unsupported
        default: return PlotLine::NOLINE;
        }
    }
    // </group>
    
    // Converts between Qt's brush style and CASA's.
    // <group>
    static Qt::BrushStyle areapattern(PlotAreaFill::Pattern p) {
        switch(p) {
        // supported
        case PlotAreaFill::FILL: return Qt::SolidPattern;
        case PlotAreaFill::MESH1: return Qt::Dense2Pattern;
        case PlotAreaFill::MESH2: return Qt::Dense4Pattern;
        case PlotAreaFill::MESH3: return Qt::Dense6Pattern;
        
        // unsupported
        default: return Qt::NoBrush;
        }
    }    
    static PlotAreaFill::Pattern areapattern(Qt::BrushStyle p) {
        switch(p) {
        // supported
        case Qt::SolidPattern: return PlotAreaFill::FILL;
        case Qt::Dense1Pattern: case Qt::Dense2Pattern:
            return PlotAreaFill::MESH1;
        case Qt::Dense3Pattern: case Qt::Dense4Pattern:
            return PlotAreaFill::MESH2;
        case Qt::Dense5Pattern: case Qt::Dense6Pattern: case Qt::Dense7Pattern:
            return PlotAreaFill::MESH3;
        
        // unsupported
        default: return PlotAreaFill::NOFILL;
        }
    }
    // </group>
    
    // Converts between Qwt's symbol style and CASA's.
    // <group>
    static QwtSymbol::Style symbol(PlotSymbol::Symbol s) {
        switch(s) {
        // supported
        case PlotSymbol::AUTOSCALING:
        case PlotSymbol::CIRCLE: return QwtSymbol::Ellipse;
        case PlotSymbol::SQUARE: return QwtSymbol::Rect;
        case PlotSymbol::DIAMOND: return QwtSymbol::Diamond;
        
        // unsupported
        default: return QwtSymbol::NoSymbol;
        }
    }    
    static PlotSymbol::Symbol symbol(QwtSymbol::Style s) {
        switch(s) {
        // supported
        case QwtSymbol::Ellipse: return PlotSymbol::CIRCLE;
        case QwtSymbol::Rect: return PlotSymbol::SQUARE;
        case QwtSymbol::Diamond: return PlotSymbol::DIAMOND;
        
        // unsupported
        default: return PlotSymbol::NOSYMBOL;
        }
    }
    // </group>
    
    // Returns a string to be used in stylesheet css for the given style.
    static casacore::String cssLineStyle(PlotLine::Style style) {
        switch(style) {
        case PlotLine::SOLID:  return "solid";
        case PlotLine::DASHED: return "dashed";
        case PlotLine::DOTTED: return "dotted";
        
        default: return "";
        }
    }
    
    // Rounds the given double to the nearest int.
    static int round(double d) {
        if(d >= 0) return (int)(d + 0.5);
        else return (int)(d - 0.5);
    }
    
    // Returns color maps for standard CASA objects.
    // <group>
    static QwtLinearColorMap* standardSpectrogramMap();    
    static QwtLinearColorMap* standardRasterMap();   
    static QwtLinearColorMap* rasterMap(const vector<double>& values);
    // </group>
};

}

#endif

#endif /*QPOPTIONS_QO_H_*/
