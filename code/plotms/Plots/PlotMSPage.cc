//# PlotMSPage.cc: Layout of PlotCanvases on a single "page".
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
#include <plotms/Plots/PlotMSPage.h>
#include <plotms/PlotMS/PlotMS.h>
#include <plotms/Plots/PlotMSPlot.h>
#include <QDebug>

using namespace casacore;
namespace casa {

////////////////////////////
// PLOTMSPAGE DEFINITIONS //
////////////////////////////

// Public Methods //

PlotMSPage::PlotMSPage(const PlotMSPage& copy) {
    operator=(copy);
}

PlotMSPage::~PlotMSPage() { }


unsigned int PlotMSPage::canvasRows() const {
	return itsCanvases_.size();
}
unsigned int PlotMSPage::canvasCols() const {
    if(itsCanvases_.size() > 0) return itsCanvases_[0].size();
    else                        return 0;
}


PlotMSPage& PlotMSPage::operator=(const PlotMSPage& copy) {
    itsParent_ = copy.itsParent_;
    itsCanvases_ = copy.itsCanvases_;
    itsCanvasOwners_ = copy.itsCanvasOwners_;
    return *this;
}


// Private Methods //

PlotMSPage::PlotMSPage(PlotMSPages& parent) :
        itsParent_(&parent){

}


void PlotMSPage::resize(unsigned int nrows, unsigned int ncols) {

    unsigned int oldrows = canvasRows();
    unsigned int oldcols = canvasCols();
    if(nrows == oldrows && ncols == oldcols) return;

    // we're shrinking in rows or columns, so disown those that will be
    // deleted
    for(unsigned int r = 0; r < itsCanvases_.size(); r++){
    	for(unsigned int c = 0; c < itsCanvases_[r].size(); c++){
    		disown(r, c);
    	}
    }

    itsCanvases_.clear();
    itsCanvasOwners_.clear();
    itsCanvases_.resize(nrows);
    itsCanvasOwners_.resize(nrows);
    
    PlotMSApp* plotms = itsParent_->itsManager_->parent();
    PlotFactoryPtr factory = plotms->getPlotFactory();
    PlotStandardMouseToolGroupPtr tools;
    pair<int, int> cimg = plotms->getParameters().cachedImageSize();
    PlotCanvasPtr canvas;
    for(unsigned int r = 0; r < nrows; r++) {
        itsCanvases_[r].resize(ncols);
        itsCanvasOwners_[r].resize(ncols);
        for(unsigned int c = 0; c < ncols; c++) {
        	canvas = factory->canvas();
        	itsCanvases_[r][c] = canvas;
        	plotms->canvasAdded( canvas );

        	// set cached image size
        	canvas->setCachedAxesStackImageSize(cimg.first, cimg.second);
        }
    }
}

PlotCanvasPtr PlotMSPage::canvas(unsigned int row, unsigned int col) {
	PlotCanvasPtr canvas = PlotCanvasPtr();
    if(row < itsCanvases_.size() && col < itsCanvases_[row].size()){
        canvas = itsCanvases_[row][col];
    }
    return canvas;
}

QList<PlotMSPlot*> PlotMSPage::owner(unsigned int row, unsigned int col) const {
	QList<PlotMSPlot*> owners;
    if(row < itsCanvasOwners_.size() && col < itsCanvasOwners_[row].size() ){
        owners = itsCanvasOwners_[row][col];
    }
    return owners;
}

bool PlotMSPage::isOwned(unsigned int row, unsigned int col) {
	bool canvasOwned = true;
	QList<PlotMSPlot*> owners = owner( row, col );
	if ( owners.isEmpty() ){
		canvasOwned = false;
	}
    return canvasOwned;
}

bool PlotMSPage::isOwner( int rowIndex, int colIndex, PlotMSPlot* plot ) const {
	bool canvasOwner = false;
	QList<PlotMSPlot*> canvasOwners = owner( rowIndex, colIndex );
	if ( canvasOwners.contains( plot )){
		canvasOwner = true;
	}
	return canvasOwner;
}

bool PlotMSPage::isSpot( int rowIndex, int colIndex, PlotMSPlot* /*plot*/ ) const {
	bool availableSpace = true;
	int canvasCount = itsCanvases_.size();
	if ( rowIndex >= canvasCount){
		availableSpace = false;
	}
	else if ( colIndex >= static_cast<int>(itsCanvases_[rowIndex].size())){
		availableSpace = false;
	}
	/*else {
		if ( !itsCanvasOwners_[rowIndex][colIndex] .isEmpty() &&
				!itsCanvasOwners_[rowIndex][colIndex].contains(plot)){
			availableSpace = false;
		}
	}*/
	return availableSpace;
}

pair<int,int> PlotMSPage::findEmptySpot() const {
	pair<int,int> location(-1,1);
	int ownerRowCount = itsCanvasOwners_.size();
	bool locationFound = false;
	for ( int i = 0; i < ownerRowCount; i++ ){
		int ownerColCount = itsCanvasOwners_[0].size();
		for ( int j = 0; j < ownerColCount; j++ ){
			if ( itsCanvasOwners_[i][j].isEmpty() ){
				location.first = i;
				location.second = j;
				locationFound = true;
				break;
			}
		}
		if ( locationFound ){
			break;
		}
	}
	return location;
}

bool PlotMSPage::setOwner(unsigned int row, unsigned int col, PlotMSPlot* plot) {
	bool ownerSet = true;
    if(row >= itsCanvasOwners_.size() ||
    		col >= itsCanvasOwners_[row].size() ||
    		plot == NULL ||
    		( !itsCanvasOwners_[row][col].isEmpty() && itsCanvasOwners_[row][col].contains(plot))){
    	ownerSet = false;
    }
    else {
    	itsCanvasOwners_[row][col].append(plot);
    }
    return ownerSet;
}

void PlotMSPage::disown( PlotMSPlot* plot ){
	int rowCount = itsCanvasOwners_.size();
	for ( int i = 0; i < rowCount; i++ ){
		int colCount = itsCanvasOwners_[i].size();
		for ( int j = 0; j < colCount; j++ ){
			if ( itsCanvasOwners_[i][j].contains( plot ) ){
				disown( i, j, plot );
			}
		}
	}
}

bool PlotMSPage::disown( int row, int col ){
	bool disowned = false;
	int rowCount = itsCanvasOwners_.size();
	if ( row < rowCount ){
		int colCount = itsCanvasOwners_[row].size();
		if ( col < colCount ){
			for ( auto iter = itsCanvasOwners_[row][col].begin( ); iter != itsCanvasOwners_[row][col].end( ); ++iter ) {
				disown( row, col, *iter);
			}
		}
		disowned = true;
	}
	return disowned;
}



bool PlotMSPage::disown(unsigned int row, unsigned int col, PlotMSPlot* plot ) {
    if(row >= itsCanvases_.size() || col >= itsCanvases_[row].size() ){
    	return false;
    }

    if ( row < itsCanvasOwners_.size() && col < itsCanvasOwners_[row].size() ){
    	if ( !itsCanvasOwners_[row][col].isEmpty() && itsCanvasOwners_[row][col].contains(plot) ){
    		/*int ownerCount = itsCanvasOwners_[row][col].size();
    		for ( int i = 0; i < ownerCount; i++ ){
    			itsCanvasOwners_[row][col][i]->canvasWasDisowned(itsCanvases_[row][col]);
    		}*/
    		plot->canvasWasDisowned( itsCanvases_[row][col]);
    		itsCanvasOwners_[row][col].removeOne( plot );
    	}
    }
    return true;
}

void PlotMSPage::clearCanvases(){
	int rowCount = itsCanvases_.size();
	for ( int i = 0; i < rowCount; i++ ){
		int colCount = itsCanvases_[i].size();
		for ( int j = 0; j < colCount; j++ ){
			clearCanvas( i, j );
		}
	}
}

void PlotMSPage::clearCanvas( int row, int col ){
	bool canvasDisowned = disown( row, col );
	if ( canvasDisowned ){
		PlotCanvasPtr canvas = itsCanvases_[row][col];
		if(canvas.null()){
			return;
		}
		canvas->showAllAxes( false );
		canvas->setTitle( "" );
		canvas->setCommonAxes( false, false );
	}
}



void PlotMSPage::setupPage() {

	// Clear out old canvases.
    PlotterPtr plotter = itsParent_->itsManager_->plotter();
    plotter->setCanvasLayout(PlotCanvasLayoutPtr());
    
    // Create a new grid layout.
    PlotMSParameters params = itsParent_->getPageParameters();
    unsigned int rows = params.getRowCount();
    unsigned int cols = params.getColCount();
    resize( rows, cols );

    PlotLayoutGrid* grid = new PlotLayoutGrid(rows, cols);
    PlotGridCoordinate coord(0, 0);
    for(unsigned int r = 0; r < rows; r++) {

        for(unsigned int c = 0; c < cols; c++) {
            coord.row = r;
            coord.col = c;
            grid->setCanvasAt(coord, itsCanvases_[r][c]);
        }
    }
    plotter->setCanvasLayout(grid);
}
}



