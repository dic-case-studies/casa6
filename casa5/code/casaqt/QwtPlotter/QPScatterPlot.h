//# QPScatterPlot.h: Qwt implementation of generic ScatterPlot class.
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
#ifndef QPSCATTERPLOT_H_
#define QPSCATTERPLOT_H_

#ifdef AIPS_HAS_QWT

#include <graphics/GenericPlotter/Plot.h>
#include <casaqt/QwtPlotter/QPOptions.h>
#include <casaqt/QwtPlotter/QPPlotItem.qo.h>

#include <QObject>

namespace casa {

// Implementation of MaskedPlot, ErrorPlot, and ColoredPlot for Qwt plotter.
class QPScatterPlot : public QPPlotItem, public MaskedScatterPlot,
                      public ErrorPlot, public ColoredPlot {
public:
    // Static //
    
    // Convenient access to class name (QPScatterPlot).
    const static casacore::String CLASS_NAME;
        
    // Non-Static //
    
    // Constructor which takes the data (and determines its type) and an
    // optional title.
    QPScatterPlot(PlotPointDataPtr data, const casacore::String& title = "Scatter Plot");
    
    // Copy constructor for generic ScatterPlot.
    QPScatterPlot(const ScatterPlot& copy);
    
    // Destructor.
    ~QPScatterPlot();
    
    
    // Include overloaded methods.
    using Plot::setLine;
    using ScatterPlot::setSymbol;
    using MaskedScatterPlot::setMaskedLine;
    using MaskedScatterPlot::setMaskedSymbol;
    using ErrorPlot::setErrorLine;
    
    
    // PlotItem Methods //
    
    // Implements PlotItem::isValid().
    bool isValid() const;
    
    
    // QPPlotItem Methods //
    
    // Overrides QPPlotItem::shouldDraw().
    bool shouldDraw() const;
    
    // Overrides QwtPlotItem::boundingRect();
    QwtDoubleRect boundingRect() const;

#if QWT_VERSION >= 0x060000
    // implements QwtPlotItem::legendIcon
    QwtGraphic legendIcon(int index, const QSizeF& size) const;
#else 
    // Overrides QwtPlotItem::legendItem().
    QWidget* legendItem() const;
#endif
    
    
    // Plot Methods //
    
    // Implements Plot::dataChanged().
    void dataChanged() { itemChanged(); }
    
    // Implements Plot::linesShown().
    bool linesShown() const;
    
    // Implements Plot::setLinesShown().
    void setLinesShown(bool linesShown = true);
    
    // Implements Plot::line().
    PlotLinePtr line() const;
    
    // Implements Plot::setLine().
    void setLine(const PlotLine& line);
    
    // Implements ScatterPlot::linesStep() and setLinesStep().
    inline bool linesStep() const { return m_step; }
    inline void setLinesStep(bool linesStep = true) { m_step = linesStep; }
    
    // ScatterPlot Methods //
    
    // Implements ScatterPlot::pointData().
    PlotPointDataPtr pointData() const;
    
    // Implements ScatterPlot::symbolsShown().
    bool symbolsShown() const;
    
    // Implements ScatterPlot::setSymbolsShown().
    void setSymbolsShown(bool symbolsShown = true);
    
    // Implements ScatterPlot::symbol().
    PlotSymbolPtr symbol() const;
    
    // Implements ScatterPlot::setSymbol().
    void setSymbol(const PlotSymbol& symbol);
    
    
    // MaskedScatterPlot Methods //
    
    // Implements MaskedScatterPlot::maskedData().
    PlotMaskedPointDataPtr maskedData() const;

    // When underlying data is deleted
    void clearData();

    // Implements MaskedScatterPlot::maskedLinesShown().
    bool maskedLinesShown() const;

    // Implements MaskedScatterPlot::setMaskedLinesShown().
    void setMaskedLinesShown(bool linesShown = true);

    // Implements MaskedScatterPlot::maskedLine().
    PlotLinePtr maskedLine() const;

    // Implements MaskedScatterPlot::setMaskedLine().
    void setMaskedLine(const PlotLine& line);

    // Implements MaskedScatterPlot::maskedLinesStep() and setMaskedLinesStep().
    inline bool maskedLinesStep() const { return m_maskedStep; }
    inline void setMaskedLinesStep(bool linesStep = true) { m_maskedStep = linesStep; }

    // Implements MaskedScatterPlot::maskedSymbolsShown().
    bool maskedSymbolsShown() const;

    // Implements MaskedScatterPlot::setMaskedSymbolsShown().
    void setMaskedSymbolsShown(bool symbolsShown = true);

    // Implements MaskedScatterPlot::maskedSymbol().
    PlotSymbolPtr maskedSymbol() const;

    // Implements MaskedScatterPlot::setMaskedSymbol().
    void setMaskedSymbol(const PlotSymbol& symbol);
    
    
    // ErrorPlot Methods //
    
    // Implements ErrorPlot::errorData().
    PlotErrorDataPtr errorData() const;
    
    // Implements ErrorPlot::errorLineShown().
    bool errorLineShown() const;
    
    // Implements ErrorPlot::setErrorLineShown().
    void setErrorLineShown(bool show = true);
    
    // Implements ErrorPlot::errorLine().
    PlotLinePtr errorLine() const;
    
    // Implements ErrorPlot::setErrorLine().
    void setErrorLine(const PlotLine& line);
    
    // Implements ErrorPlot::errorCapSize().
    unsigned int errorCapSize() const;
    
    // Implements ErrorPlot::setErrorCapSize().
    void setErrorCapSize(unsigned int capSize);
    
    
    // ColoredPlot Methods //
    
    // Implements ColoredPlot::binnedColorData().
    PlotBinnedDataPtr binnedColorData() const;
    
    // Implements ColoredPlot::colorForBin().
    PlotColorPtr colorForBin(unsigned int bin) const;
    
    // Implements ColoredPlot::setColorForBin().
    void setColorForBin(unsigned int bin, const PlotColorPtr color);
    
protected:
    // Implements QPPlotItem::className().
    const casacore::String& className() const { return CLASS_NAME; }
    
    // Implements QPLayerItem::draw_().
#if QWT_VERSION >= 0x060000
    void draw_(QPainter* painter, const QwtScaleMap& xMap,
              const QwtScaleMap& yMap, const QRectF& canvasRect,
              unsigned int drawIndex, unsigned int drawCount) const;
#else
    void draw_(QPainter* painter, const QwtScaleMap& xMap,
              const QwtScaleMap& yMap, const QRect& canvasRect,
              unsigned int drawIndex, unsigned int drawCount) const;
#endif
    
private:
    // casacore::Data pointers.
    // <group>
    PlotPointDataPtr m_data;
    PlotMaskedPointDataPtr m_maskedData;
    PlotErrorDataPtr m_errorData;
    PlotBinnedDataPtr m_coloredData;
    // </group>
    
    // Customization objects.
    // <group>
    QPSymbol m_symbol;
    QPLine m_line;
    bool m_step;
    QPSymbol m_maskedSymbol;
    QPLine m_maskedLine;
    QPLine m_errorLine;
    unsigned int m_errorCap;
    bool m_maskedStep;
    // </group>
    
    // Binned colors.
    // <group>
    QList<QPColor*> m_colors;
    QList<QBrush> m_coloredBrushes;
    // </group>
    
#if QWT_VERSION >= 0x060000
    // Make non-const symbol from m_symbol for colorized plots
    QPSymbol* coloredSymbol(const QColor& color) const;
#endif

    // Updates the binned color brushes.
    void updateBrushes();
};

}

#endif

#endif /*QPSCATTERPLOT_H_*/
