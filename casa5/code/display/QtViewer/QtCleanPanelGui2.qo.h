//# QtCleanPanelGui2.qo.h:  interactive clean display panel
//# Copyright (C) 2005,2009
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

#ifndef QTCLEANPANELGUI2_H_
#define QTCLEANPANELGUI2_H_

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Arrays/Vector.h>
#include <display/Display/DisplayCoordinateSystem.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>

#include <graphics/X11/X_enter.h>
#  include <QObject>
#  include <QGroupBox>
#  include <QRadioButton>
#  include <QPushButton>
#  include <QHBoxLayout>
#  include <QVBoxLayout>
#  include <QLineEdit>
#include <graphics/X11/X_exit.h>
#include <display/QtViewer/QtDisplayPanelGui.qo.h>


namespace casacore{

	class ImageRegion;
}

namespace casa {

	class QtViewer;
	class QtDisplayPanel;
	class QtDisplayData;
	class WorldCanvasHolder;

	// <synopsis>
	// Demo class to encapsulate 'serial' running of qtviewer into callable
	// methods of a class; this example also applies it to the task of
	// interactive selection of CLEAN boxes.
	// </synopsis>
	class QtCleanPanelGui2: public QtDisplayPanelGui {

		Q_OBJECT	//# Allows slot/signal definition.  Must only occur in
		//# implement/.../*.h files; also, makefile must include
		//# name of this file in 'mocs' section.


	protected:
		friend class QtViewer;
		QtCleanPanelGui2( QtViewer* v, QWidget* parent=0,
		                 const std::list<std::string> &args = std::list<std::string>( ) );

	public:
		~QtCleanPanelGui2();

		bool supports( SCRIPTING_OPTION option ) const;
		QVariant start_interact( const QVariant &input, int id );
		QVariant setoptions( const QMap<QString,QVariant> &input, int id);

		// the QtDBusViewerAdaptor can handle loading & registering data itself,
		// but to connect up extra functionality, the upper-level QtDisplayPanelGui
		// (or in the current case, the derived QtCleanPanelGui2) would have to be
		// notified that data has been added. This will allow it to set up the
		// callbacks for drawing regions...
		void addedData( QString type, QtDisplayData * );

	protected slots:

		virtual void exitStop();
		virtual void exitDone();
		virtual void exitNoMore();

		// Connected to the rectangle region mouse tools new rectangle signal.
		// Accumulates [/ displays] selected boxes.
		virtual void newMouseRegion(casacore::Record mouseRegion, WorldCanvasHolder* wch);

		virtual void changeMaskAxis(casacore::String, casacore::String, casacore::String, std::vector<int> );
		virtual void changeImageAxis(casacore::String, casacore::String, casacore::String, std::vector<int> );
		virtual void changeMaskSelectionText( casacore::String x, casacore::String y, casacore::String z );

	signals:
		void interact( QVariant );

	protected:

		// scripted (via dbus) panels should override the closeEvent( ) and hide the gui
		// instead of deleting it when it was created via a dbus script...
		void closeEvent(QCloseEvent *event);

		std::list<QWidget*> disabled_widgets;

		QRadioButton* addRB_;
		QRadioButton* eraseRB_;
		QRadioButton* allChanRB_;
		QRadioButton* thisPlaneRB_;
		QRadioButton* allHiddenRB_;
		QRadioButton* thisHiddenRB_;
		QPushButton* maskNoMorePB_;
		QPushButton* maskDonePB_;
		QPushButton* stopPB_;
		QLineEdit* niterED_;
		QLineEdit* cycleniterED_;
		QLineEdit* threshED_;
		QLineEdit* cyclethreshED_;
		casacore::Record buttonState_;

		// standard palette...
		QPalette default_palette;
		// palette used when input is expected for clean...
		QPalette input_palette;

	private:
		bool in_interact_mode;
		int interact_id;

		void writeRegionText(const casacore::ImageRegion& imageReg, const casacore::String& filename, float value);

		QtDisplayData* imagedd_;
		QtDisplayData* maskdd_;			// later: to display clean region.
		std::string axis_change;

		DisplayCoordinateSystem csys_p;
		casacore::DirectionCoordinate dirCoord_p;

	};


} //# NAMESPACE CASA - END

#endif


