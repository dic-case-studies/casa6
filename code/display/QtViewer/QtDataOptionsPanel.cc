//# QtDataOptionsPanel.cc: Qt implementation DD options adjustment window
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
//# $Id$

#include <QScrollArea>
#include <casa/BasicSL/String.h>
#include <display/QtViewer/QtDataOptionsPanel.qo.h>
#include <display/QtViewer/QtDisplayPanelGui.qo.h>
#include <display/QtViewer/QtDisplayDataGui.qo.h>


using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN


	QtDataOptionsPanel::QtDataOptionsPanel(QtDisplayPanelGui* panel,
			QWidget* parent) :
		QWidget(parent), panel_(panel)  {

		setWindowTitle("Viewer Data Options Panel");

		setupUi(this);
		connect( globalColorCheckBox, SIGNAL(toggled(bool)), this, SIGNAL(globalColorSettingsChanged(bool)));
		connect( auto_apply, SIGNAL(clicked(bool)), SLOT(auto_apply_state_change(bool)) );

		//#dk tabs_->setFixedWidth(585);
		//#dk tabs_->setMinimumHeight(600);
		//setMinimumWidth(570);		//#dk
		//setMinimumHeight(600);	//#dk
		// setMinimumWidth(200);			//#dk
		// setMinimumHeight(200);		//#dk
		tabs_->setMinimumWidth(350);	//#dk
		tabs_->setMinimumHeight(150);	//#dk
		connect( tabs_, SIGNAL(currentChanged(int)), this, SLOT(tabChanged(int)));
		//#dk with auto-app:    resize(595, 705);		//#dk
		// resize(510, 705);			//#dk
		resize(485, 665);			//#dk

		//if(tabs_->layout()==0) new QVBoxLayout(tabs_);	//#dk
		//tabs_->layout()->setSizeConstraint(QLayout::SetFixedSize);	//#dk
		// does no good because you cant add tabs_'s children to the layout).

		// if(layout()!=0) layout()->setSizeConstraint(QLayout::SetFixedSize); //#dk


		while(tabs_->count()>0) tabs_->removeTab(0);
		// (Qt designer insists on putting some tabs in my
		// TabWidget, whether I want them there yet or not...).

		//List<QtDisplayData*> ddlist(panel_->dds());
		//for (ListIter<QtDisplayData*> qdds(ddlist); !qdds.atEnd(); qdds++) {
		//createDDTab_(qdds.getRight());
		DisplayDataHolder::DisplayDataIterator iter = panel_->beginDD();
		while ( iter != panel_->endDD() ) {
			createDDTab_(*iter, false, -1);
			iter++;
		}
		// These are the tabs I _really_ want....


		connect( panel_, SIGNAL(ddCreated(QtDisplayData*, bool, int, bool)),
		         SLOT(createDDTab_(QtDisplayData*, bool, int)) );

		/*connect( panel_, SIGNAL(ddRemoved(QtDisplayData*)),
		         SLOT(removeDDTab_(QtDisplayData*)) );*/
	}


	QtDataOptionsPanel::~QtDataOptionsPanel() { }
	// (elaboration probably unnecessary because of Qt parenting...)


// Slots for QDD creation/destruction.

	void QtDataOptionsPanel::createDDTab_(QtDisplayData* qdd,
			Bool /*autoregister*/,
			int insertPosition ) {

		QtDisplayDataGui* qddg = new QtDisplayDataGui(qdd);

		connect( this, SIGNAL(setAutoApply(bool)), qddg, SLOT(autoApplyState(bool)) );

		//cerr<<"QDO:crT:nT:"<<tabs_->count()<<" nm:"<<qdd->nameChrs();  //#dg

		//cerr<<"QDO:crTb:tbs.szHnt:"<<tabs_->sizeHint().width()<<","
		//<<tabs_->sizeHint().height()<<endl;	//#dg

		QScrollArea* sca = new QScrollArea;
		sca->setWidget(qddg);
		if ( insertPosition < 0  || insertPosition>=tabs_->count() ){
			tabs_->addTab(sca, qdd->nameChrs());
		}
		else {
			tabs_->insertTab( insertPosition, sca, qdd->nameChrs());
		}

		if(tabs_->count()==1 && panel_->autoDDOptionsShow) panel_->showDataOptionsPanel();
		// Users want to see this window automatically once there
		// are DDs to tweak.  (Apps like interactive clean can turn
		// panel_->autoDDOptionsShow off to prevent this, if desired).

		tabs_->setTabToolTip(tabs_->indexOf(sca), qdd->nameChrs());
		tabs_->show();

		// set state for data loaded after clicking auto/manual apply...
		emit setAutoApply(auto_apply->text( ) == "auto apply" ? true : false);

	}



	void QtDataOptionsPanel::removeDDTab_(QtDisplayData* qdd) {

		for(Int i=0; i<tabs_->count(); i++) {

			if(tabs_->tabText(i) == qdd->nameChrs()) {

				QScrollArea* sca = dynamic_cast<QScrollArea*>(tabs_->widget(i));

				tabs_->removeTab(i);	// (NB: does not delete sca).

				if(tabs_->count()==0) panel_->hideDataOptionsPanel();
				// (An empty options panel is just confusing clutter).


				if(sca!=0) delete static_cast<QtDisplayDataGui*>(sca->widget());
				// This may be unnecessary, since sca is the parent, but I want
				// to assure that the QDDG is fully deleted (~QWigdet is not
				// virtual... hmm...).

				delete sca;

				break;
			}
		}
	}


	void QtDataOptionsPanel::tabChanged( int index ) {
		QString tabName = tabs_->tabText( index );
		emit dataOptionsTabChanged( tabName );
	}

	void QtDataOptionsPanel::auto_apply_state_change(bool) {
		//  green: rgb(20, 227, 67), red: rgb(255, 53, 43)
		if ( auto_apply->text( ) == "auto apply" ) {
			auto_apply->setText("manual apply");
			auto_apply->setStyleSheet( "QPushButton { color: rgb(255, 53, 43) }" );
			emit setAutoApply(false);
		} else {
			auto_apply->setText("auto apply");
			auto_apply->setStyleSheet( "QPushButton { color: rgb(20, 227, 67) }" );
			emit setAutoApply(true);
		}
	}


//void QtDataOptionsPanel::paintEvent ( QPaintEvent * event ) {
//  resize(sizeHint());  }
	//#dk hye-type trick (a very suspect action for a paintEvent....
	//#dk I have seen it cause some (but not infinite) recursion).

	void QtDataOptionsPanel::resizeEvent (QResizeEvent* /*ev*/) {	//#dg
//  cerr<<"DOPrsz x:"<<ev->size().width()<<" y:"<<ev->size().height()<<endl;
	}				//#dg -- to show size.



} //# NAMESPACE CASA - END

