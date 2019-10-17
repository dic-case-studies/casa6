//# TBPlotter.cc: Widget to collect plot parameters and plot on the canvas.
//# Copyright (C) 2005
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

#include <casaqt/QtBrowser/TBPlotter.qo.h>
#include <casaqt/QtBrowser/TBBrowser.qo.h>
#include <casaqt/QtBrowser/TBTableTabs.qo.h>
#include <casaqt/QtBrowser/TBTable.h>
#include <casaqt/QtBrowser/TBField.h>
#include <casaqt/QtBrowser/TBFilterRules.qo.h>
#include <casaqt/QtBrowser/TBConstants.h>
#include <casaqt/QtBrowser/QProgressPanel.qo.h>
#include <casaqt/QtBrowser/TBRowLocate.qo.h>

#include <algorithm>

using namespace std;

using namespace casacore;
namespace casa {

////////////////////////////
// PLOTSLICER DEFINITIONS //
////////////////////////////

// Constructors/Destructors //

PlotSlicer::PlotSlicer(): QHBoxLayout(), complex(false) {
    setSpacing(3);
    setContentsMargins(0, 0, 0, 0);
    complexChooser = new QComboBox();
    complexChooser->addItem("Real");
    complexChooser->addItem("Imaginary");

    axisLabel = new QLabel("Plot axis");
    plotAxisChooser = new QSpinBox();
    plotAxisChooser->setWrapping(true);
    connect(plotAxisChooser, SIGNAL(valueChanged(int)),
            this, SLOT(axisChosen(int)));
}

PlotSlicer::~PlotSlicer() {
    delete complexChooser;
    delete plotAxisChooser;
    delete axisLabel;
    for(unsigned int i = 0; i < spinners.size(); i++)
        delete spinners.at(i);
}

// Public Methods //

bool PlotSlicer::setDimension(vector<int>* d, bool c) {
  return setDimension(d,c,false);
}

bool PlotSlicer::setDimension(vector<int>* d, bool c, bool index) {
    bool valid = true;
    
    // Clear out old widgets
    while(count() > 0) {
        removeItem(itemAt(0));
    }

    for(unsigned int i = 0; i < spinners.size(); i++)
        delete spinners.at(i);
    spinners.clear();

    // Add spinners
    if(d != NULL && d->size() > 0) {
        for(unsigned int i = 0; i < d->size(); i++) {
            QSpinBox* s = new QSpinBox();
            s->setWrapping(true);
            int n = d->at(i);
            if(n > 0)
                s->setMaximum(n - 1);
            else {
                s->setEnabled(false);
                valid = false;
            }
            addWidget(s);
            spinners.push_back(s);
        }
    }

    // ADD axis Slicer
    bool indexPlot = false;
    if (index && d != NULL && d->size() > 0){
      addWidget(axisLabel);
      plotAxisChooser->setMaximum(d->size() - 1);
      addWidget(plotAxisChooser);
      plotAxisChooser->setValue(0);
      axisChosen(0);
      indexPlot = true;
    }
    axisLabel ->setVisible(indexPlot);    
    plotAxisChooser->setVisible(indexPlot);


    // Add complexChooser if necessary
    complex = c;
    complexChooser->setCurrentIndex(0);
    if(c) addWidget(complexChooser);
    complexChooser->setVisible(c);

    return valid;
}

void PlotSlicer::getDimension(vector<int>& d, bool& c, bool& a) {
    for(unsigned int i = 0; i < spinners.size(); i++) {
        d.push_back(spinners.at(i)->value());
    }
    c = complex;
    if(c)
        a = complexChooser->currentIndex() == 0;
}

void PlotSlicer::getDimension(vector<int>& d, bool& c, bool& a, int& axis) {
  getDimension(d,c,a);
  if (plotAxisChooser->isVisible()){
    axis = plotAxisChooser->value();
  } else {
    axis = -1;
  }
}

void PlotSlicer::axisChosen(int axis){
  if (axis >= 0 && axis < int(spinners.size())){
    for (unsigned int i = 0; i < spinners.size(); i++){
      spinners.at(i)->setEnabled(true);
    }
    spinners.at(axis)->setEnabled(false);
  }
}


///////////////////////////
// TBPLOTTER DEFINITIONS //
///////////////////////////

// Constructors/Destructors //

TBPlotter::TBPlotter(TBBrowser* b, PlotFactoryPtr f): QMainWindow(),
        browser(b), factory(f), plotCanvas(NULL), isIndexPlot(false) {
    setupUi(this);
    
    if(f.null()) {
        String message = "Error: given a null factory, which probably means ";
        message += "there's a problem with your plotting implementation.";
#ifndef AIPS_HAS_QWT
        if(TBConstants::defaultPlotterImplementation == casa::Plotter::QWT) {
            message += "\nThe browser is currently set to use QWT, but your";
            message += " makedefs say that you don't have it installed.  ";
            message += "Please install QWT and rebuild.";
        }
#endif
        QMessageBox::critical(this, "Plot Factory Error", message.c_str());
        close();
        deleteLater();
        return;
    }

    vector<String> v = PlotExportFormat::supportedImageFormatStrings();
    for(unsigned int i = 0; i < v.size(); i++)
        exportChooser->addItem(v[i].c_str());
    
    plotCanvas = new TBPlotCanvas(factory);
    TBConstants::insert(canvasFrame, plotCanvas);
    TBConstants::insert(xSliceFrame, &xSlice);
    TBConstants::insert(ySliceFrame, &ySlice);
    
    selectBox->setVisible(false);
    connect(selectLocateButton, SIGNAL(clicked()), this, SLOT(selectLocate()));
    connect(selectClearButton, SIGNAL(clicked()), this,SLOT(clearSelection()));
    
    // Set up tables
    vector<String> tables = browser->openedTableNames();
    for(unsigned int i = 0; i < tables.size(); i++)
        tableChooser->addItem(tables.at(i).c_str());

    if(tables.size() <= 1)
        tableChooser->setEnabled(false);

    // Need to set up model to disable [index] selection
    xChooserModel = new QStandardItemModel();
    yChooserModel = new QStandardItemModel();
    xChooser->setModel(xChooserModel);
    yChooser->setModel(yChooserModel);

    tableChooser->setCurrentIndex(b->getTabWidget()->currentIndex());
    tableChosen(tableChooser->currentText());

    connect(xChooser, SIGNAL(currentIndexChanged(int)),
            this, SLOT(xChosen(int)));
    connect(yChooser, SIGNAL(currentIndexChanged(int)),
            this, SLOT(yChosen(int)));
    xChosen(0);
    yChosen(0);

    // Format
    
    PlotLinePtr line = TBConstants::defaultPlotLine(factory);
    PlotSymbolPtr symbol = TBConstants::defaultPlotSymbol(factory);
    
    lineWidthSpinner->setValue((int)(line->width() + 0.5));
    String color = line->color()->asHexadecimal();
    if(color.size() > 0 && color[0] != '#') color = "#" + color;
    lineColorEdit->setText(color.c_str());
    
    pointSize1->setValue((int)(symbol->size().first + 0.5));
    pointSize2->setValue((int)(symbol->size().second + 0.5));
    color = symbol->areaFill()->color()->asHexadecimal();
    if(color.size() > 0 && color[0] != '#') color = "#" + color;
    symbolColorEdit->setText(color.c_str());
    
    connect(lineColorChoose, SIGNAL(clicked()), this, SLOT(setLineColor()));
    connect(symbolColorChoose, SIGNAL(clicked()), this,SLOT(setSymbolColor()));

    // Connect signals
    connect(tableChooser, SIGNAL(currentIndexChanged(QString)),
            this, SLOT(tableChosen(QString)));
    connect(plotButton, SIGNAL(clicked()), this, SLOT(plot()));
    connect(overplotButton, SIGNAL(clicked()), this, SLOT(overplot()));
    connect(openButton, SIGNAL(clicked()), this, SLOT(openNewPlotter()));
    connect(browser, SIGNAL(tableOpened(casacore::String, casacore::String)),
            this, SLOT(tableOpened(casacore::String, casacore::String)));
    connect(browser, SIGNAL(tableClosed(casacore::String)),
            this, SLOT(tableClosed(casacore::String)));
    connect(clearButton, SIGNAL(clicked()),
            plotCanvas, SLOT(clearAndHideAxes()));
    connect(allRowsButton, SIGNAL(clicked()), this, SLOT(allRows()));
    connect(exportButton, SIGNAL(clicked()), this, SLOT(exportImage()));
    connect(b, SIGNAL(filterRuleAvailable(int)),
            this, SLOT(filterRuleEntered(int)));
    connect(b, SIGNAL(filterRuleCleared(int)),
                this, SLOT(filterRuleCleared(int)));
    connect(gridXmaj, SIGNAL(clicked()), this, SLOT(gridChanged()));
    connect(gridXmin, SIGNAL(clicked()), this, SLOT(gridChanged()));
    connect(gridYmaj, SIGNAL(clicked()), this, SLOT(gridChanged()));
    connect(gridYmin, SIGNAL(clicked()), this, SLOT(gridChanged()));
    connect(plotCanvas, SIGNAL(regionSelected(bool)),
            this, SLOT(regionSelected(bool)));
}

TBPlotter::~TBPlotter() {
    map<String, vector<vector<int>*> >::iterator iter = dimensions.begin();

    for(; iter != dimensions.end(); iter++) {
        vector<vector<int>*> d = iter->second;
        for(unsigned int i = 0; i < d.size(); i++)
            delete d.at(i);
    }

    if(plotCanvas != NULL) delete plotCanvas;

    delete xChooserModel;
    delete yChooserModel;
}

// Public Methods //

QProgressPanel* TBPlotter::addProgressPanel(String label, bool hideable,
                                            bool cancelable) {
    // disable GUI components
    plotCanvas->setEnabled(false);
    paramsBox->setEnabled(false);
    formatBox->setEnabled(false);
    toolsBox->setEnabled(false);
    
    QProgressPanel* qpp = new QProgressPanel(label, hideable, cancelable);
    gridLayout->addWidget(qpp, 0, 0, -1, -1);
    qpp->show();

    // force a repaint
    QCoreApplication::processEvents();

    return qpp;
}

void TBPlotter::removeProgressPanel(QProgressPanel* panel) {
    gridLayout->removeWidget(panel);
    plotCanvas->setEnabled(true);
    paramsBox->setEnabled(true);
    formatBox->setEnabled(true);
    toolsBox->setEnabled(true);
}


// Protected Methods //

void TBPlotter::closeEvent(QCloseEvent* event) {
    if(dockWidget->isFloating())
        dockWidget->close();

    event->accept();
}


// Private Methods //

void TBPlotter::doPlot(bool overplot) {
    if(!xValid || !yValid) return;
    
    TBTableTabs* tt = browser->table(qPrintable(tableChooser->currentText()));
    if(tt == NULL) return;

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

    // Collect plot parameters
    String xName = qPrintable(xChooser->currentText());
    int xi = xChooser->currentIndex();
    String xType = "";
    if(xi > 1) xType = types.at(xi - 2);
    String yName = qPrintable(yChooser->currentText());
    int yi = yChooser->currentIndex();
    String yType = "";
    if(yi > 1) yType = types.at(yi - 2);
        
    int rf = rowFrom->value();
    int rt = rowTo->value();
    int ri = rowInterval->value();

    // Make sure isIndexPlot is set properly
    if ((isIndexPlot && xi != 1 && yi != 1) || 
	(!isIndexPlot && (xi ==1 || yi == 1)) || (xi == 1 && yi == 1)) {
        QMessageBox::critical(this, "Plot Error",
			      "Invalid data selection for index plotting.");
        QApplication::restoreOverrideCursor();
        return;
      }

    if(isIndexPlot){
      rt = rf;
      ri = 1;
    } else if (rf > rt) {
        int temp = rf;
        rf = rt;
        rt = temp;
        rowFrom->setValue(rf);
        rowTo->setValue(rt);
    }
        
    if((rf == rt && ri != 1) || ri > (rt - rf + 1)) {
        QMessageBox::critical(this, "Plot Error", "Invalid row selection.");
        QApplication::restoreOverrideCursor();
        return;
    }

    bool xComplex(false), yComplex(false), xAmp(false), yAmp(false);
    vector<int> xs, ys;
    int xAxis, yAxis;
    //xSlice.getDimension(xs, xComplex, xAmp);
    //ySlice.getDimension(ys, yComplex, yAmp);
    xSlice.getDimension(xs, xComplex, xAmp, xAxis);
    ySlice.getDimension(ys, yComplex, yAmp, yAxis);
    
    stringstream ss;
    ss << "For table " << tt->getName() << ", plotting:" << endl;
    ss << "\tRows " << rf << " : " << ri << " : " << rt << endl;
    ss << "\tX: " << xName;
    if(xs.size() > 0) ss << " [array]";
    if(xComplex) {
        ss << " [complex - ";
        if(xAmp) ss << "real";
        else ss << "imaginary";
        ss << "]";
    }
    ss << endl;
        
    ss << "\tY: " << yName;
    if(ys.size() > 0) ss << " [array]";
    if(yComplex) {
        ss << " [complex - ";
        if(yAmp) ss << "real";
        else ss << "imaginary";
        ss << "]";
    }
    TBConstants::dprint(TBConstants::DEBUG_HIGH, ss.str());

    int adjustedX = 0;
    if(xi > 1) adjustedX = adjustedIndices.at(xi - 2);
    int adjustedY = 0;
    if(yi > 1) adjustedY = adjustedIndices.at(yi - 2);

    TBPlotData* data = NULL;
    try {
        PlotParams xp;
        xp.rowNumbers = xi == 0;
        xp.complex = xComplex;
        xp.complexAmp = xAmp;
        xp.colIndex = adjustedX;
        xp.slice = xs;

        PlotParams yp;
        yp.rowNumbers = yi == 0;
        yp.complex = yComplex;
        yp.complexAmp = yAmp;
        yp.colIndex = adjustedY;
        yp.slice = ys;

        //ph.step();
        
        TBFilterRuleSequence* rule = NULL;
        if(filterBox->isEnabled() && filterBox->isChecked()) {
            rule = browser->filterAt(tableChooser->currentIndex());
        }
        
        // Get data
	if (isIndexPlot) {
	  bool isX = (xi == 1);
	  data = tt->getTable()->plotIndices((isX ? yp : xp),
					     (isX ? yAxis : xAxis), 
					     isX, rf, rule);
	} else {
	  data = tt->getTable()->plotRows(xp, yp, rf, rt, ri, rule);
	}
    } catch(const char* s) {
        QMessageBox::critical(this, "Plot Error",
                              "Error loading data from table.");
        String mesg = "Plot error: ";
        mesg += s;
        TBConstants::dprint(TBConstants::DEBUG_HIGH, mesg);
        QApplication::restoreOverrideCursor();
        return;
    }
    
    // If data is valid, plot it
    if(data != NULL && !data->data.null() && data->data->isValid()) {
        int n = data->data->size();
        data->table = tt;

        TBConstants::dprint(TBConstants::DEBUG_LOW, "Plotting:");
        for(int i = 0; i < n; i++) {
            TBConstants::dprint(TBConstants::DEBUG_LOW, "(" +
                                TBConstants::dtoa(data->data->xAt(i)) + "," +
                                TBConstants::dtoa(data->data->yAt(i)) + ")", 1);
        }

        if(xs.size() > 0) {
            xName += " [";
            for(unsigned int i = 0; i < xs.size(); i++) {
	        if (i != static_cast<unsigned int>(xAxis)) {
		  xName += TBConstants::itoa(xs.at(i));
		} else {
		  xName += "*";
		}
                if(i < xs.size() - 1) xName += " ";
            }
            xName += "]";
        }
        if(xComplex) {
            if(xAmp) xName += " - real";
            else xName += " - imaginary";
        }
        if(ys.size() > 0) {
            yName += " [";
            for(unsigned int i = 0; i < ys.size(); i++) {
	        if (i != static_cast<unsigned int>(yAxis)) {
		  yName += TBConstants::itoa(ys.at(i));
		} else {
		  yName += "*";
		}
                if(i < ys.size() - 1) yName += " ";
            }
            yName += "]";
        }
        if(yComplex) {
            if(yAmp) yName += " - real";
            else yName += " - imaginary";
        }

        // set up plot format
        TBPlotFormat format(factory);
        format.setCurveStyle(lineStyleChooser->currentText());
        format.line->setWidth(lineWidthSpinner->value());
        
        QColor qcolor(lineColorEdit->text());
        String scolor;
        if(qcolor.isValid()) {
            scolor = qPrintable(qcolor.name());
            if(scolor.size() > 0 && scolor[0] == '#') scolor.erase(0,1);
            format.line->setColor(scolor);
        }
        
        format.setPointStyle(pointStyleChooser->currentText());
        format.symbol->setSize(pointSize1->value(), pointSize2->value());
        qcolor.setNamedColor(symbolColorEdit->text());
        if(qcolor.isValid()) {
            scolor = qPrintable(qcolor.name());
            if(scolor.size() > 0 && scolor[0] == '#') scolor.erase(0,1);
            PlotAreaFillPtr area = format.symbol->areaFill();
            area->setColor(scolor);
            format.symbol->setAreaFill(area);
        }
        if(symbolOutlineBox->isChecked()) format.symbol->setLine("000000");
            
        // plot on canvas
        plotCanvas->setXAxisTitle(xName);
        plotCanvas->setYAxisTitle(yName);
        data->title = xName + " vs. " + yName;
        
        plotCanvas->setXAxisDate(xType == TBConstants::TYPE_DATE);
        plotCanvas->setYAxisDate(yType == TBConstants::TYPE_DATE);
        plotCanvas->plot(data, format, overplot);
    } else {
        QMessageBox::critical(this, "Plot Error",
                              "Error loading data from table.");
    }

    QApplication::restoreOverrideCursor();
}

// Private Slots //

void TBPlotter::tableChosen(QString n) {
    String name = qPrintable(n);
    TBTableTabs* tt = browser->table(name);
    TBTable* table;

    if(tt != NULL && ((table = tt->getTable()) != NULL)) {
        map<String, vector<vector<int>*> >::iterator iter =
                                                dimensions.find(name);

        if(iter == dimensions.end()) {
            vector<vector<int>*> dims;

            TBConstants::dprint(TBConstants::DEBUG_MED,
                                "Fetching dimensions for table " +
                                table->getName() + "...");
            for(int i = 0; i < table->getNumFields(); i++) {
                String type = table->field(i)->getType();
                if(TBConstants::typeIsPlottable(type)) {
                    vector<int>* d = NULL;
                    if(TBConstants::typeIsArray(type)) {
                        try {
                            d = new vector<int>(table->dataDimensionsAt(i));
                            stringstream ss;
                            ss << table->getFields()->at(i)->getName();
                            ss << ": [ ";
                            for(unsigned int j = 0; j < d->size(); j++)
                                ss << d->at(j) << ' ';
                            ss << ']';
                            TBConstants::dprint(TBConstants::DEBUG_MED,
                                                ss.str(), 1);
                        } catch(const char* s) {
                            d = NULL;
                            stringstream ss;
                            ss << "Cannot fetch dimensions of field ";
                            ss << table->getFields()->at(i)->getName();
                            ss << ".";
                            String msg = ss.str();
                            TBConstants::dprint(TBConstants::DEBUG_HIGH, s);
                            QMessageBox::critical(this, "Plot Error",
                                                    msg.c_str());
                        }
                    }
                    dims.push_back(d);
                }
            }
            dimensions[name] = dims;
        }

        rowFrom->setMaximum(table->getTotalRows() - 1);
        rowTo->setMaximum(table->getTotalRows() - 1);
        rowInterval->setMaximum(table->getTotalRows());
        if((int)TBConstants::DEFAULT_SELECT_NUM <= table->getTotalRows())
            rowTo->setValue(TBConstants::DEFAULT_SELECT_NUM);
        else
            rowTo->setValue(table->getTotalRows());
        rowInterval->setValue(TBConstants::DEFAULT_ROW_INTERVAL);

        update = false;
        adjustedIndices.clear();
        xChooser->clear();
        yChooser->clear();
        types.clear();
        vector<TBField*>* fields = table->getFields();
        xChooser->addItem("[row #]");
        yChooser->addItem("[row #]");
        xChooser->addItem("[index]");
        yChooser->addItem("[index]");
        for(unsigned int i = 0; i < fields->size(); i++) {
            TBField* field = fields->at(i);
            String type = field->getType();
            
            if(TBConstants::typeIsPlottable(type)) {
                xChooser->addItem(fields->at(i)->getName().c_str());
                yChooser->addItem(fields->at(i)->getName().c_str());
                adjustedIndices.push_back(i);
                types.push_back(type);
            }
        }
        update = true;
        xChosen(0);
        yChosen(0);

        filterBox->setEnabled(browser->filterAvailable(browser->indexOf(tt)));
        plotCanvas->setTable(name);
    }
}

void TBPlotter::xChosen(int x) {
    chosen(true, x);
}

void TBPlotter::yChosen(int y) {
    chosen(false, y);
}

void TBPlotter::chosen(bool x, int i) {
    if(!update) return;

    (x ? xSliceLabel : ySliceLabel)->setEnabled(false);

    // Check for the selection of the other axis
    int j = (x ? yChooser : xChooser)->currentIndex();
    if (isIndexPlot && i != 1 && j != 1){
      // index plot is canceled by this selection
      indexReleased(x, j);
    }

    if(i == 0) {
        (x ? xSlice : ySlice).setDimension(NULL, false);
        (x ? xValid : yValid) = true;
        return;
    }
    
    if(i == 1) {
        //index is selected
        indexChosen(x);
        (x ? xValid : yValid) = true;
        return;
    }
        
    map<String, vector<vector<int>*> >::iterator iter =
        dimensions.find(qPrintable(tableChooser->currentText()));
    
    if(iter != dimensions.end()) {
        vector<vector<int>*> dims = iter->second;

        vector<int>* d = dims.at(i - 2);
        String type = types.at(i - 2);
        bool c = TBConstants::typeIsComplex(type);
        
        if((d != NULL && d->size() > 0) || c)
            (x ? xSliceLabel : ySliceLabel)->setEnabled(true);
        //(x ? xValid : yValid) = (x ? xSlice : ySlice).setDimension(d, c);
        (x ? xValid : yValid) = (x ? xSlice : ySlice).setDimension(d, c, (j == 1));
        
        if(!(x ? xValid : yValid)) {
            stringstream ss;
            ss << "Field " << qPrintable((x?xChooser:yChooser)->itemText(i));
            ss << " contains no data or is otherwise invalid.";
            String msg = ss.str();
            QMessageBox::warning(this, "Plot Error", msg.c_str());
        }
    }
}

void TBPlotter::indexChosen(bool x) {
    // No spinbox for index axis
    (x ? xValid : yValid) = (x ? xSlice : ySlice).setDimension(NULL,false);

    // Disable "[index]" of the other axis
    (x ? yChooserModel : xChooserModel)->item((x ? yChooser : xChooser)->findText("[index]"))->setEnabled(false);

    // Reset spinbox of the other axis
    int otherId = (x ? yChooser : xChooser)->currentIndex();
    if (otherId > 2){
      map<String, vector<vector<int>*> >::iterator iter =
        dimensions.find(qPrintable(tableChooser->currentText()));
    
      if(iter != dimensions.end()) {
        vector<vector<int>*> dims = iter->second;

	vector<int>* od = dims.at(otherId - 2);
        String type = types.at(otherId - 2);
        bool c = TBConstants::typeIsComplex(type);

        if((od != NULL && od->size() > 0) || c)
            (x ? xSliceLabel : ySliceLabel)->setEnabled(true);
	(x ? yValid : xValid) = (x ? ySlice : xSlice).setDimension(od, c, true);
      }
    } else { // the other axis is row number
      (x ? yValid : xValid) = (x ? ySlice : xSlice).setDimension(NULL, false, true);
    }
    // disable Row iteration
    enableRowIteration(false);

    isIndexPlot = true;
}

void TBPlotter::indexReleased(bool x, int i) {
    // this should be set first
    isIndexPlot = false;

    // re-configure spinbox of the other axis
    (x ? yChosen(i) : xChosen(i));
    // enable "[index]" of the other axis
    (x ? yChooserModel : xChooserModel)->item((x ? yChooser : xChooser)->findText("[index]"))->setEnabled(true);
    // enable Row iteration
    enableRowIteration(true);
}

void TBPlotter::enableRowIteration(bool visible) {
    rowTo->setVisible(visible);
    rowToLabel->setVisible(visible);
    rowInterval->setVisible(visible);
    rowIntervalLabel->setVisible(visible);
    allRowsButton->setVisible(visible);
}

void TBPlotter::plot() {
    doPlot(false);
}

void TBPlotter::overplot() {
    doPlot(true);
}

void TBPlotter::openNewPlotter() {
    TBPlotter* plotter = new TBPlotter(browser, factory);
    plotter->show();
}

void TBPlotter::tableOpened(String table, String /*fullpath*/) {
    tableChooser->addItem(table.c_str());
    tableChooser->setEnabled(true);
}

void TBPlotter::tableClosed(String table) {
    for(int i = 0; i < tableChooser->count(); i++) {
        String t = qPrintable(tableChooser->itemText(i));
        if(t == table) {
            tableChooser->removeItem(i);
            if(tableChooser->count() <= 1)
                tableChooser->setEnabled(false);
        }
    }
}

void TBPlotter::clear() {   
    plotCanvas->clearAndHideAxes();
}

void TBPlotter::allRows() {
    rowFrom->setValue(0);
    rowTo->setValue(rowTo->maximum());
    rowInterval->setValue(1);
}

void TBPlotter::setColor(QLineEdit* lineEdit) {
    QColor c = QColorDialog::getColor(QColor(lineEdit->text()), this);
    if(c.isValid())
        lineEdit->setText(c.name());
}

void TBPlotter::exportImage() {
    if(exportChooser->count() < 1) return; // shouldn't happen, but make sure
    
    if(plotCanvas->getNumPlots() < 1) {
        QMessageBox::warning(this, "Export Error",
                             "There are no plots to export!");
        return;
    }

    QString format = exportChooser->currentText();
    QString file = QFileDialog::getSaveFileName(this, "Save as " + format,
                                 "", format + " Image (*." + format.toLower() +
                                 ")");

    if(!file.isEmpty()) {
        String sformat = qPrintable(format);
        String sfile = qPrintable(file);
        
        Result r = plotCanvas->exportToImage(sformat, sfile);
        if(!r.valid)
            QMessageBox::critical(this, "Export Error", r.result.c_str());
    }
}

void TBPlotter::filterRuleEntered(int i) {
    if(i == tableChooser->currentIndex()) {
        filterBox->setEnabled(true);
    }
}

void TBPlotter::filterRuleCleared(int i) {
    if(i == tableChooser->currentIndex()) {
        filterBox->setChecked(false);
        filterBox->setEnabled(false);
    }
}

void TBPlotter::gridChanged() {
    plotCanvas->setShownGrids(gridXmaj->isChecked(), gridXmin->isChecked(),
                              gridYmaj->isChecked(), gridYmin->isChecked());
}

void TBPlotter::regionSelected(bool selected) {
    selectBox->setVisible(selected);
    selectLocateButton->setEnabled(selected);
    selectClearButton->setEnabled(selected);
}

void TBPlotter::clearSelection() {
    plotCanvas->clearSelectedRectangle();
    selectBox->setVisible(false);
    selectLocateButton->setEnabled(false);
    selectClearButton->setEnabled(false);
}

void TBPlotter::selectLocate() {
    vector<TBPlotData*> data = plotCanvas->allData();
    
    map<TBTableTabs*, pair<vector<int>*, bool> > rows;
    
    PlotRegion region = plotCanvas->currentSelection();
    double x1 = region.upperLeft().x(), x2 = region.lowerRight().x();
    double y1 = region.lowerRight().y(), y2 = region.upperLeft().y();
    vector<int>* r;
    int index;
    unsigned int n;
    double x, y;
    for(unsigned int i = 0; i < data.size(); i++) {
        TBPlotData* d = data[i];
        TBTableTabs* t = d->table;
        
        
        // check that t is still open
        index = browser->indexOf(t);
        if(index < 0) continue;
        
        map<TBTableTabs*, pair<vector<int>*, bool> >::iterator it=rows.find(t);
        if(it == rows.end()) {
            r = new vector<int>();
            rows[t] = pair<vector<int>*, bool>(r, false);
        } else {
            r = it->second.first;
            it->second.second = true;
        }
        
        n = d->data->size();
        for(unsigned int j = 0; j < n; j++) {
            x = d->data->xAt(j), y = d->data->yAt(j);
            if(x >= x1 && x <= x2 && y >= y1 && y <= y2)
                r->push_back(d->rows[j]);
        }
    }
    
    // make sure the map has at least one entry in it
    // (in case all the selected tables were closed)
    map<TBTableTabs*, pair<vector<int>*, bool> >::iterator it = rows.begin();
    if(it == rows.end()) return;
    
    bool found;
    int ir;
    for(map<TBTableTabs*,pair<vector<int>*,bool> >::iterator iter=rows.begin();
        iter != rows.end(); iter++) {
        if(iter->second.second) {
            // multiple plots from one table, so we have to merge them
            r = new vector<int>();
            for(unsigned int i = 0; i < iter->second.first->size(); i++) {
                ir = iter->second.first->at(i);
                found = false;
                for(unsigned int j = 0; j < r->size(); j++) {
                    if(ir == r->at(j)) {
                        found = true;
                        break;
                    }
                }
                if(!found) r->push_back(ir);
            }
            std::sort(r->begin(), r->end());
            delete iter->second.first;
            iter->second.first = r;
        }
    }
    
    TBLocatedRows* lr = new TBLocatedRows();
    for(map<TBTableTabs*,pair<vector<int>*,bool> >::iterator iter=rows.begin();
    iter != rows.end(); iter++)
        lr->put(iter->first, iter->second.first);

    TBRowLocate* rl = new TBRowLocate(lr);
    rl->setVisible(true);
}

}
