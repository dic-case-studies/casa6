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

#include "ExternalAxisControlLeft.h"
#include <QDebug>
#include <QPainter>
#include <qwt_plot.h>
#include <qwt_plot_canvas.h>
#include <qwt_scale_div.h>

namespace casa {

ExternalAxisControlLeft::ExternalAxisControlLeft(QWidget* parent):
	ExternalAxisControl( parent ){
	setSizePolicy( QSizePolicy::Fixed, QSizePolicy::MinimumExpanding );
	setFixedWidth( AXIS_SMALL_SIDE );
}



void ExternalAxisControlLeft::defineAxis( QLine& line ){
	const int MARGIN = 1;
	int x = width() - MARGIN;;
	int top = getStartY();
	int bottom = height() - 2*MARGIN;
	QPoint firstPt( x, top );
	QPoint secondPt( x, bottom );
	line.setP1( firstPt );
	line.setP2( secondPt );
}

void ExternalAxisControlLeft::drawTick( QPainter* painter, float yPixel, double value,
		int tickLength ){

	//Draw the tick
	int xEnd = width();
	int tickStart = xEnd - tickLength;
	int yValue = static_cast<int>( yPixel );
	painter->drawLine( tickStart, yValue, xEnd, yValue  );

	//Draw the tick label
	QString numberStr = QString::number( value );
	QFont font = painter->font();
	QRect fontBoundingRect = QFontMetrics(font).boundingRect( numberStr );
	int labelStart = tickStart - fontBoundingRect.width() - 5;
	int letterHeight = fontBoundingRect.height();
	int position = static_cast<int>(yPixel + letterHeight/3);
	painter->drawText( labelStart, position, numberStr);
}

int ExternalAxisControlLeft::getStartY() const {
	int canvasHeight = plot->canvas()->height();
	int heightDiff =  height() - canvasHeight;
	if ( canvasHeight < MIN_START_Y ){
		heightDiff = MIN_START_Y;
	}
	return heightDiff;
}

void ExternalAxisControlLeft::drawTicks( QPainter* painter, int tickLength ){

	//Figure out how far out to start drawing ticks.
	double startPixelY = getTickStartPixel(QwtPlot::yLeft);

	//We don't want to draw too many ticks so adjust the number
	//of ticks we draw accordingly.
#if QWT_VERSION >= 0x060000
	QwtScaleDiv scaleDiv = plot->axisScaleDiv( QwtPlot::yLeft );
	const QList<double> axisTicks = scaleDiv.ticks(QwtPlot::yLeft);
#else
	QwtScaleDiv* scaleDiv = plot->axisScaleDiv( QwtPlot::yLeft );
	const QList<double> axisTicks = scaleDiv->ticks(QwtPlot::yLeft);
#endif
	int originalTickCount = axisTicks.size();
	int tickIncrement = getTickIncrement( originalTickCount );

	//Now figure out the yIncrement, how far apart the ticks should be.
	double tickDistance = getTickDistance( QwtPlot::yLeft);
	double yIncrement = getTickIncrement( tickDistance, QwtPlot::yLeft );
	for ( int i = 0; i < originalTickCount; i = i + tickIncrement ){
		//Sometimes the automatic tick system puts uneven number of ticks in.
		//Definitely weird - but that is why the incrementCount;
		int incrementCount = qRound((axisTicks[i] - axisTicks[0]) / tickDistance);
		double tickPosition = startPixelY + incrementCount * yIncrement;
		int tickIndex = originalTickCount - i - 1;
		drawTick( painter, tickPosition, axisTicks[tickIndex], tickLength);
	}
}

void ExternalAxisControlLeft::drawAxisLabel( QPainter* painter ){
	  QFont font = painter->font();
	  QString mainLabel = axisLabel.trimmed();

	  painter->rotate(-90);

	  //Draw the main label
	  QRect fontBoundingRect = QFontMetrics(font).boundingRect( mainLabel );
	  int yPosition = fontBoundingRect.height() / 2;
	  int xPosition = -height() + (height() - fontBoundingRect.width())/2;
	  painter->drawText( xPosition, yPosition, fontBoundingRect.width(),
			  fontBoundingRect.height(), Qt::AlignHCenter|Qt::AlignTop, mainLabel);

	  painter->rotate(90);
}

ExternalAxisControlLeft::~ExternalAxisControlLeft() {
}

} /* namespace casa */
