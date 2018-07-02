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

#include <QScrollArea>
#include "DisplayOptionsDialog.h"
#include <display/QtViewer/QtDisplayDataGui.qo.h>
#include <display/QtViewer/QtDisplayData.qo.h>

using namespace casacore;
namespace casa {

DisplayOptionsDialog::DisplayOptionsDialog(QWidget* parent): QDialog(parent),
		displayedDataOptions( NULL ) {
}

void DisplayOptionsDialog::showOptions( QtDisplayData* displayData ){

	QLayout* dialogLayout = layout();
	if ( dialogLayout == NULL ){
		dialogLayout = new QHBoxLayout();
	}
	else {
		while( !dialogLayout->isEmpty()){
			dialogLayout->takeAt( 0 );
		}
		delete displayedDataOptions;
		displayedDataOptions = NULL;
	}

	displayedDataOptions= new QtDisplayDataGui( displayData );
	QScrollArea* scrollArea = new QScrollArea();
	scrollArea->setWidget(displayedDataOptions);
	dialogLayout->addWidget( scrollArea );
	setLayout( dialogLayout );
	setWindowTitle( displayData->name().c_str() );
	resize(485,600);
}

DisplayOptionsDialog::~DisplayOptionsDialog() {
	delete displayedDataOptions;
}

using namespace casacore;
} /* namespace casa */
