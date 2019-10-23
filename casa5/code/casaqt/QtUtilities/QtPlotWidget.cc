//# PlotWidget.cc: Classes for GUI editing of plot customization objects.
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
#include <casaqt/QtUtilities/QtPlotWidget.qo.h>

#include <casaqt/QtUtilities/QtUtilities.h>

#include <QColorDialog>
#include <QDebug>

using namespace std;

using namespace casacore;
namespace casa {

//////////////////////////////
// QTPLOTWIDGET DEFINITIONS //
//////////////////////////////

QtPlotWidget::QtPlotWidget(PlotFactoryPtr factory, QWidget* parent) :
        QtEditingWidget(parent), itsFactory_(factory) { }
QtPlotWidget::~QtPlotWidget() { }


/////////////////////////////////
// PLOTCOLORWIDGET DEFINITIONS //
/////////////////////////////////

PlotColorWidget::PlotColorWidget(PlotFactoryPtr factory, bool showAlpha,
        QWidget* parent) : QtPlotWidget(factory, parent) {
    setupUi(this);
    if(!showAlpha) {
        alphaLabel->setVisible(false);
        alpha->setVisible(false);
    }
    setColor(factory->color("000000"));
    
    connect(colorEdit, SIGNAL(editingFinished()), SLOT(colorChanged()));
    connect(choose, SIGNAL(clicked()), SLOT(colorChoose()));
    if(showAlpha)
        connect(alpha, SIGNAL(valueChanged(double)), SLOT(colorChanged()));
}

PlotColorWidget::~PlotColorWidget() { }

PlotColorPtr PlotColorWidget::getColor() const {
    PlotColorPtr color = itsFactory_->color(colorEdit->text().toStdString());
    if(alpha->isVisible()) color->setAlpha(alpha->value());
    return color;
}

void PlotColorWidget::setColor(PlotColorPtr color) {
    if(!color.null()) {
        blockSignals(true);
        bool changed = itsColor_.null() || *itsColor_ != *color;
        
        itsColor_ = itsFactory_->color(*color);
        colorEdit->setText(color->asHexadecimal().c_str());
        if(alpha->isVisible()) alpha->setValue(color->alpha());
        
        blockSignals(false);
        if(changed){

        	emit this->changed();
        }
    }
}

String PlotColorWidget::getColorString() const {
    return colorEdit->text().toStdString(); }

void PlotColorWidget::setColor(const String& color) {
    setColor(itsFactory_->color(color)); }

void PlotColorWidget::colorChoose() {
    QColor color = QColorDialog::getColor(QColor("#" + colorEdit->text()));
    if(color.isValid()){
    	QString colorStr = color.name().replace( '#', "");
    	//setColor( colorStr.toStdString());
    	colorEdit->setText( colorStr );
    	colorChanged();
    }
}

void PlotColorWidget::colorChanged() {
    emit changed();
    PlotColorPtr currColor = getColor();
    if(*currColor != *itsColor_) emit differentFromSet();
}


////////////////////////////////
// PLOTFILLWIDGET DEFINITIONS //
////////////////////////////////

PlotFillWidget::PlotFillWidget(PlotFactoryPtr factory, bool showAlpha,
		String fillColor, QWidget* parent) :
		QtPlotWidget(factory, parent) {
    setupUi(this);
    itsColorWidget_ = new PlotColorWidget(factory, showAlpha);
    QtUtilities::putInFrame(colorFrame, itsColorWidget_);
    
    setFill(factory->areaFill(fillColor));
    
    connect(itsColorWidget_, SIGNAL(changed()), SLOT(fillChanged()));
    connect(fillChooser, SIGNAL(currentIndexChanged(int)), SLOT(fillChanged()));
}

PlotFillWidget::~PlotFillWidget() { }

PlotAreaFillPtr PlotFillWidget::getFill() const {
    PlotAreaFillPtr fill = itsFactory_->areaFill("");
    fill->setColor(itsColorWidget_->getColor());
    
    int index = fillChooser->currentIndex();
    PlotAreaFill::Pattern pattern = PlotAreaFill::NOFILL;
    if(index == 0)      pattern = PlotAreaFill::FILL;
    else if(index == 1) pattern = PlotAreaFill::MESH1;
    else if(index == 2) pattern = PlotAreaFill::MESH2;
    else if(index == 3) pattern = PlotAreaFill::MESH3;
    fill->setPattern(pattern);
    
    return fill;
}

void PlotFillWidget::setFill(PlotAreaFillPtr fill) {
    if(!fill.null()) {
        blockSignals(true);
        bool changed = itsFill_.null() || *itsFill_ != *fill;
        
        itsFill_ = itsFactory_->areaFill(*fill);
        itsColorWidget_->setColor(fill->color());
        
        PlotAreaFill::Pattern pattern = fill->pattern();
        int index = 4;
        if(pattern == PlotAreaFill::FILL)       index = 0;
        else if(pattern == PlotAreaFill::MESH1) index = 1;
        else if(pattern == PlotAreaFill::MESH2) index = 2;
        else if(pattern == PlotAreaFill::MESH3) index = 3;
        fillChooser->setCurrentIndex(index);
        
        blockSignals(false);
        if(changed) emit this->changed();
    }
}

void PlotFillWidget::fillChanged() {
    emit changed();
    PlotAreaFillPtr currFill = getFill();
    if(*currFill != *itsFill_) emit differentFromSet();
}




////////////////////////////////
// PLOTLINEWIDGET DEFINITIONS //
////////////////////////////////

PlotLineWidget::PlotLineWidget(PlotFactoryPtr factory, bool useCompact,
        bool showAlpha, QWidget* parent) : QtPlotWidget(factory, parent) {
    setupUi(this);
    stackedWidget->setCurrentIndex(useCompact ? 0 : 1);
    
    itsColorWidget_ = new PlotColorWidget(factory, showAlpha);
    QtUtilities::putInFrame(useCompact ? cColorFrame : nColorFrame,
                            itsColorWidget_);
    
    setLine(factory->line("black"));
    
    if(useCompact) {
        cWidth->setValidator(new QIntValidator(1, 99, cWidth));
        connect(cWidth, SIGNAL(textChanged(const QString&)),
                SLOT(lineChanged()));
        connect(cStyle, SIGNAL(currentIndexChanged(int)), SLOT(lineChanged()));
    } else {
        connect(nWidth, SIGNAL(valueChanged(int)), SLOT(lineChanged()));
        connect(nStyle, SIGNAL(currentIndexChanged(int)), SLOT(lineChanged()));
    }
    connect(itsColorWidget_, SIGNAL(changed()), SLOT(lineChanged()));
}

PlotLineWidget::~PlotLineWidget() { }

PlotLinePtr PlotLineWidget::getLine() const {
    PlotLinePtr line = itsFactory_->line("");
    int width = cWidth->text().toInt();
#ifdef Q_WS_MAC
	// cannot see dotted grid on Mac, make bigger
	if (width==1 && lineStyle() == PlotLine::DOTTED) {
		width += 1;
	}
#endif
    line->setWidth((stackedWidget->currentIndex() == 0) ?
                   width : nWidth->value());
    line->setStyle(lineStyle());
    line->setColor(itsColorWidget_->getColor());
    return line;
}

void PlotLineWidget::setLine(PlotLinePtr line) {
    if(!line.null()) {
        blockSignals(true);
        bool changed = itsLine_.null() || *itsLine_ != *line;
        
        itsLine_ = itsFactory_->line(*line);
        int width = (int)(line->width() + 0.5);
        if(stackedWidget->currentIndex() == 0)
            cWidth->setText(QString::number(width));
        else nWidth->setValue(width);
        setLineStyle(line->style());
        itsColorWidget_->setColor(line->color());
        
        blockSignals(false);
        if(changed) emit this->changed();
    }
}

PlotLine::Style PlotLineWidget::lineStyle() const {
    int index = (stackedWidget->currentIndex() == 0) ?
                cStyle->currentIndex() : nStyle->currentIndex();
    if(index == 0) return PlotLine::SOLID;
    else if(index == 1) return PlotLine::DASHED;
    else if(index == 2) return PlotLine::DOTTED;
    else return PlotLine::NOLINE;
}

void PlotLineWidget::setLineStyle(PlotLine::Style style) {
    int index;
    if(style == PlotLine::SOLID)       index = 0;
    else if(style == PlotLine::DASHED) index = 1;
    else if(style == PlotLine::DOTTED) index = 2;
    else                               index = 3;
    if(stackedWidget->currentIndex() == 0) cStyle->setCurrentIndex(index);
    else                                   nStyle->setCurrentIndex(index);
}

void PlotLineWidget::lineChanged() {
    emit changed();
    PlotLinePtr currLine = getLine();
    if(*currLine != *itsLine_) emit differentFromSet();
}


//////////////////////////////////
// PLOTSYMBOLWIDGET DEFINITIONS //
//////////////////////////////////

PlotSymbolWidget::PlotSymbolWidget(PlotFactoryPtr factory,
        PlotSymbolPtr defaultSymbol, bool showAlphaFill, bool showCustom,
        bool showAlphaLine, bool showCharacter, QWidget* parent) :
        QtPlotWidget(factory, parent) {
    setupUi(this);

    //Line
    itsLineWidget_ = new PlotLineWidget(factory, false, showAlphaLine);
    QtUtilities::putInFrame(outlineFrame, itsLineWidget_);
    outlineFrame->setEnabled(false);
    
    if(!showCustom) outlineCustomFrame->setVisible(false);
    
    //Symbol
    if(defaultSymbol.null()){
        itsDefault_ = itsFactory_->symbol(PlotSymbol::AUTOSCALING);
    }
    else {
        itsDefault_ = itsFactory_->symbol(*defaultSymbol);
    }
    setSymbol(itsDefault_);
    
    //Fill
    String symbolColor = itsDefault_->getColor();
    itsFillWidget_ = new PlotFillWidget(factory, showAlphaFill, symbolColor );
    QtUtilities::putInFrame(fillFrame, itsFillWidget_);

    if(!showCharacter && itsDefault_->symbol() != PlotSymbol::CHARACTER) {
        SymbolWidget::style->removeItem(4);
        charEdit->hide();
    }
    
    // only emit change for radio buttons turned on
    connect(noneButton, SIGNAL(toggled(bool)), SLOT(symbolChanged(bool)));
    connect(defaultButton, SIGNAL(toggled(bool)), SLOT(symbolChanged(bool)));
    connect(customButton, SIGNAL(toggled(bool)), SLOT(symbolChanged(bool)));
    connect(outlineNone, SIGNAL(toggled(bool)), SLOT(symbolChanged(bool)));
    connect(outlineDefault, SIGNAL(toggled(bool)), SLOT(symbolChanged(bool)));
    connect(outlineCustom, SIGNAL(toggled(bool)), SLOT(symbolChanged(bool)));
    
    connect(SymbolWidget::size, SIGNAL(valueChanged(int)),
            SLOT(symbolChanged()));
    connect(SymbolWidget::style, SIGNAL(currentIndexChanged(int)),
            SLOT(symbolChanged()));
    connect(charEdit, SIGNAL(editingFinished()), SLOT(symbolChanged()));
    connect(itsFillWidget_, SIGNAL(changed()), SLOT(symbolChanged()));
    connect(itsLineWidget_, SIGNAL(changed()), SLOT(symbolChanged()));
}

PlotSymbolWidget::~PlotSymbolWidget() { }

PlotSymbolPtr PlotSymbolWidget::getSymbol() const {
    if(defaultButton->isChecked()) return itsFactory_->symbol(*itsDefault_);
    else {
        PlotSymbol::Symbol s = PlotSymbol::NOSYMBOL;
        
        if(customButton->isChecked()) {
            int i = SymbolWidget::style->currentIndex();
            if(i == 0)      s = PlotSymbol::CIRCLE;
            else if(i == 1) s = PlotSymbol::SQUARE;
            else if(i == 2) s = PlotSymbol::DIAMOND;
            else if(i == 3) s = PlotSymbol::PIXEL;
            else if(i == 4) s = PlotSymbol::AUTOSCALING;
        }
        
        PlotSymbolPtr sym = itsFactory_->symbol(s);
        QString text = charEdit->text();
        if(s == PlotSymbol::CHARACTER && text.size() >= 1)
            sym->setUSymbol(text[0].unicode());
        
        int i = SymbolWidget::size->value();
        sym->setSize(i, i);

        PlotAreaFillPtr areaFill = itsFillWidget_->getFill();
        sym->setAreaFill(areaFill);
        sym->setColor( areaFill->color());
        
        if(outlineNone->isChecked())
            sym->setLine(itsFactory_->line("black", PlotLine::NOLINE, 1));
        else if(outlineDefault->isChecked())
            sym->setLine(itsFactory_->line("black", PlotLine::SOLID, 1));
        else
            sym->setLine(itsLineWidget_->getLine());
        
        return sym;
    }
}

void PlotSymbolWidget::setSymbol(PlotSymbolPtr symbol) {    
    if(symbol.null()) itsSymbol_ = itsFactory_->symbol(PlotSymbol::NOSYMBOL);
    else              itsSymbol_ = itsFactory_->symbol(*symbol);

    blockSignals(true);
    bool changed = itsSymbol_.null() || *itsSymbol_ != *symbol;

    if(itsSymbol_->symbol() == PlotSymbol::NOSYMBOL)
        noneButton->setChecked(true);
    else if(*itsSymbol_ == *itsDefault_) defaultButton->setChecked(true);
    else customButton->setChecked(true);

    if(itsSymbol_->symbol() == PlotSymbol::PIXEL)
        itsSymbol_->setSize(1, 1);

    if(itsMinSizes_.find(itsSymbol_->symbol()) != itsMinSizes_.end())
        SymbolWidget::size->setMinimum(itsMinSizes_[itsSymbol_->symbol()]);
    else SymbolWidget::size->setMinimum(1);

    SymbolWidget::size->setValue((int)(itsSymbol_->size().first + 0.5));
    SymbolWidget::size->setEnabled(
        itsSymbol_->symbol() != PlotSymbol::PIXEL &&
        itsSymbol_->symbol() != PlotSymbol::AUTOSCALING);

    PlotSymbol::Symbol s = itsSymbol_->symbol();
    int index = 0;
    if(s == PlotSymbol::SQUARE)    index = 1;
    else if(s == PlotSymbol::DIAMOND)   index = 2;
    else if(s == PlotSymbol::PIXEL)     index = 3;
    //else if(s == PlotSymbol::CHARACTER) index = 4;
    else if(s == PlotSymbol::AUTOSCALING) index = 4;
    SymbolWidget::style->setCurrentIndex(index);
    charEdit->setEnabled(s == PlotSymbol::CHARACTER);
    if(s == PlotSymbol::CHARACTER)
        charEdit->setText(QString(itsSymbol_->symbolUChar()));

    //itsFillWidget_->setFill(itsSymbol_->areaFill());
    PlotLinePtr line = itsSymbol_->line();
    if(line->style() == PlotLine::NOLINE) outlineNone->setChecked(true);
    else if(*line == *itsFactory_->line("black", PlotLine::SOLID, 1))
        outlineDefault->setChecked(true);
    else {
        outlineCustom->setChecked(true);
        itsLineWidget_->setLine(line);
    }
    
    blockSignals(false);
    if(changed) emit this->changed();
}

void PlotSymbolWidget::addRadioButtonsToGroup(QButtonGroup* group) const {
    if(group == NULL) return;
    group->addButton(noneButton);
    group->addButton(defaultButton);
    group->addButton(customButton);
}

void PlotSymbolWidget::setMinimumSizes(const map<PlotSymbol::Symbol, int>& m) {
    for(map<PlotSymbol::Symbol, int>::const_iterator iter = m.begin();
        iter != m.end(); iter++)
        itsMinSizes_[iter->first] = iter->second;
    
    PlotSymbolPtr currSymbol = getSymbol();
    
    if(itsMinSizes_.find(currSymbol->symbol()) != itsMinSizes_.end())
        SymbolWidget::size->setMinimum(itsMinSizes_[currSymbol->symbol()]);
    else SymbolWidget::size->setMinimum(1);
}

void PlotSymbolWidget::symbolChanged(bool check) {
    if(!check) return;
    
    PlotSymbolPtr currSymbol = getSymbol();
    
    if(itsMinSizes_.find(currSymbol->symbol()) != itsMinSizes_.end())
        SymbolWidget::size->setMinimum(itsMinSizes_[currSymbol->symbol()]);
    else SymbolWidget::size->setMinimum(1);

    if(currSymbol->symbol() == PlotSymbol::PIXEL) {
        currSymbol->setSize(1, 1);
        SymbolWidget::size->setValue(1);
    } else if(currSymbol->symbol() == PlotSymbol::AUTOSCALING) {
        currSymbol->setSize(2, 2);
        SymbolWidget::size->setValue(2);
    }

    charEdit->setEnabled(currSymbol->symbol() == PlotSymbol::CHARACTER);
    SymbolWidget::size->setEnabled(
        currSymbol->symbol() != PlotSymbol::PIXEL &&
        currSymbol->symbol() != PlotSymbol::AUTOSCALING);
    emit changed();
    if(*currSymbol != *itsSymbol_) emit differentFromSet();
}


////////////////////////////////
// PLOTFONTWIDGET DEFINITIONS //
////////////////////////////////

PlotFontWidget::PlotFontWidget(PlotFactoryPtr factory, bool showAlpha,
        QWidget* parent) : QtPlotWidget(factory, parent) {
    setupUi(this);
    itsColorWidget_ = new PlotColorWidget(factory, showAlpha);
    QtUtilities::putInFrame(colorFrame, itsColorWidget_);
        
    setFont(itsFactory_->font());
    
    // only emit change for radio buttons turned on
    connect(family, SIGNAL(currentIndexChanged(int)), SLOT(fontChanged()));
    connect(FontWidget::size, SIGNAL(valueChanged(int)), SLOT(fontChanged()));
    connect(sizeUnit, SIGNAL(currentIndexChanged(int)), SLOT(fontChanged()));
    connect(itsColorWidget_, SIGNAL(changed()), SLOT(fontChanged()));
    connect(bold, SIGNAL(toggled(bool)), SLOT(fontChanged()));
    connect(italic, SIGNAL(toggled(bool)), SLOT(fontChanged()));
    connect(underline, SIGNAL(toggled(bool)), SLOT(fontChanged()));
}

PlotFontWidget::~PlotFontWidget() { }

PlotFontPtr PlotFontWidget::getFont() const {
    PlotFontPtr font = itsFactory_->font(family->currentText().toStdString());
    font->setPointSize(sizeUnit->currentIndex() == 0 ?
                       FontWidget::size->value() : -1);
    font->setPixelSize(sizeUnit->currentIndex() == 1 ?
                       FontWidget::size->value() : -1);
    font->setColor(itsColorWidget_->getColor());
    font->setBold(bold->isChecked());
    font->setItalics(italic->isChecked());
    font->setUnderline(underline->isChecked());
    return font;
}

void PlotFontWidget::setFont(PlotFontPtr font) {
    if(!font.null()) {
        blockSignals(true);
        bool changed = itsFont_.null() || *itsFont_ != *font;
        
        itsFont_ = itsFactory_->font(*font);
        family->setCurrentFont(QFont(itsFont_->fontFamily().c_str()));
        if(itsFont_->pointSize() >= 0) {
            FontWidget::size->setValue((int)(itsFont_->pointSize() + 0.5));
            sizeUnit->setCurrentIndex(0);
        } else {
            FontWidget::size->setValue((int)(itsFont_->pixelSize() + 0.5));
            sizeUnit->setCurrentIndex(1);
        }
        itsColorWidget_->setColor(itsFont_->color());
        bold->setChecked(itsFont_->bold());
        italic->setChecked(itsFont_->italics());
        underline->setChecked(itsFont_->underline());
        
        blockSignals(false);
        if(changed) emit this->changed();
    }
}

void PlotFontWidget::fontChanged() {
    emit changed();
    PlotFontPtr currFont = getFont();
    if(*currFont != *itsFont_) emit differentFromSet();
}

}
