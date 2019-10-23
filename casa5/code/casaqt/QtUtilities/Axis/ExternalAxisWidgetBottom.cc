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

#include "ExternalAxisWidgetBottom.h"
#include <QDebug>
#include <QPainter>
#include <qwt_plot.h>
#include <qwt_plot_canvas.h>
#include <qwt_scale_div.h>
#include <casaqt/QwtPlotter/QPCanvasHelpers.qo.h>

using namespace casacore;
namespace casa {

ExternalAxisWidgetBottom::ExternalAxisWidgetBottom( QWidget* parent, QwtPlot* plot,
		bool leftAxisInternal, bool bottomAxisInternal,
		bool rightAxisInternal ):
	ExternalAxisWidgetHorizontal( parent, plot, leftAxisInternal, bottomAxisInternal,
			rightAxisInternal ){
	setMinimumHeight( AXIS_SMALL_SIDE / 2 );
	scaleDraw = new QPScaleDraw( NULL, QwtPlot::xBottom );
	setSizePolicy( QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding );
	plotAxis = QwtPlot::xBottom;
}



void ExternalAxisWidgetBottom::defineAxis( QLine& axisLine ){
	int left = getStartX();
	int right = getEndX();
	int top = 1;
	int bottom = 1;

	QPoint firstPt( left, top );
	QPoint secondPt( right, bottom );
	axisLine.setP1( firstPt );
	axisLine.setP2( secondPt );
}

void ExternalAxisWidgetBottom::drawTick( QPainter* painter, double xPixel, double value,
		int tickLength ){
	int yEnd = 0;
	int xPosition = static_cast<int>(xPixel);
	painter->drawLine( xPosition, yEnd, xPosition, yEnd + tickLength  );
	QFont font = painter->font();

	QString tickLabel = formatLabel( value );
	QRect fontBoundingRect = QFontMetrics(font).boundingRect( tickLabel );
	int letterWidth = fontBoundingRect.width();
	int letterHeight = fontBoundingRect.height();
	xPosition = static_cast<int>( xPixel - letterWidth / 2 );
	painter->drawText( xPosition, yEnd + tickLength + letterHeight, tickLabel );
}


void ExternalAxisWidgetBottom::drawAxisLabel( QPainter* painter ){
	 QFont font = painter->font();
	 QRect fontBoundingRect = QFontMetrics(font).boundingRect( axisLabel );
	 int yPosition = height() / 2;
	 int xPosition = width()/2 - fontBoundingRect.width() / 2;
	 painter->drawText( xPosition, yPosition, fontBoundingRect.width(), height(),
					  Qt::AlignHCenter|Qt::AlignTop, axisLabel);
}

ExternalAxisWidgetBottom::~ExternalAxisWidgetBottom() {
    delete scaleDraw;
}

using namespace casacore;
} /* namespace casa */
