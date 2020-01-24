//# QtDBusViewerAdaptor.qo.h: provides viewer services via dbus
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
//# $Id$

#ifndef QTDBUSVIEWERADAPTOR_QO_H_
#define QTDBUSVIEWERADAPTOR_QO_H_

#include <map>
#include <limits>
#include <casaqt/QtUtilities/QtDBusXmlApp.qo.h>
#include <casa/namespace.h>
#include <QVariantMap>
#include <display/QtViewer/QtViewer.qo.h>

namespace casa {

	class QtDisplayData;
	class QtDisplayPanel;
	class QtDisplayPanelGui;

	class QtDBusViewerAdaptor : public QDBusAbstractAdaptor, public QtDBusApp {
		Q_OBJECT
		Q_CLASSINFO("D-Bus Interface", "edu.nrao.casa.viewer")

	public:

		QString dbusName( ) const {
			return QString(QtViewer::name( ));
		}
		bool connectToDBus( const QString &dbus_name="" ) {
			return QtDBusApp::connectToDBus( parent(), dbus_name );
		}

		// Constructor which takes the application.
		QtDBusViewerAdaptor( QtViewer * );

		// Destructor.
		~QtDBusViewerAdaptor();

	public slots:

		QDBusVariant start_interact( const QDBusVariant &input, int panel );
		QDBusVariant datarange( const QList<double> &range, int data=0 );
		QDBusVariant contourlevels( const QList<double> &levels,
		                            double baselevel=2147483648.0,
		                            double unitlevel=2147483648.0, int panel_or_data=0 );
		QDBusVariant contourcolor( const QString &color, int panel_or_data=0 );
		QDBusVariant contourthickness( double thickness, int panel_or_data=0 );
		QDBusVariant axes( const QString &x = "", const QString &y = "", const QString &z = "", int panel=0 );
		QDBusVariant colormap( const QString &map, int panel_or_data=0 );
		QDBusVariant colorwedge( bool show=true, int panel_or_data=0 );

		QDBusVariant load( const QString &path, const QString &displaytype = "raster", int panel=0, double scaling=0 );
		QDBusVariant reload( int panel_or_data );
		QDBusVariant unload( int data );

		QDBusVariant restore( const QString &path, int panel=0 );

		QDBusVariant panel( const QString &type="viewer", bool hidden=false  );
		QDBusVariant hide( int panel=0 );
		QDBusVariant show( int panel=0 );
		QDBusVariant close( int panel=0 );
		QDBusVariant popup( const QString &what, int panel=0 );

		QDBusVariant cwd( const QString &new_path = "" );

		QDBusVariant freeze( int panel=0 );
		QDBusVariant unfreeze( int panel=0 );

		QDBusVariant channel( int num=-1, int panel=0 );
		QDBusVariant zoom( int level, int panel=0 );
		QDBusVariant zoom( const QList<double> &blc, const QList<double> &trc, const QString &coordinates="pixel", int panel=0 );

		// like "close()", but leaves the closing up to the user if the window is not hidden
		QDBusVariant release( int panel );
		//To change certain things so the panel
		QDBusVariant setoptions( const QDBusVariant &input, int panel);
		// device:	file name or printer name
		// devicetype:	"file", "printer", or "ghostscript"
		// format:	"jpg", "pdf", "eps", "ps", "png", "xbm", "xpm", "ppm"
		//		[only used with devicetype == "file"] [extension on device, e.g. "outfile.pdf" overrides "format"]
		// scale:	size scale factor for raster output (e.g. jpg etc.)
		// dpi:		resolution of PS or EPS images
		// orientation: "portrait", "landscape"
		// media:	"letter" or "a4"
		bool output( const QString &device, const QString &devicetype = "file", int panel = 0,
		             double scale = 1.0, int dpi = 300, const QString &format = "jpg",
		             const QString &orientation = "portrait", const QString &media = "letter" );


		QDBusVariant fileinfo( const QString &path );
		QStringList keyinfo( int );

		bool done( );

	signals:
		void interact( QDBusVariant );
		void exiting( QDBusVariant );

	protected slots:

		void handle_interact( QVariant );
		void handle_destroyed_panel(QObject*);

	private:

		inline QDBusVariant error( const QString &message ) {
			QMap<QString,QVariant> err;
			err["*error*"] = message;
			return QDBusVariant(QVariant(err));
		}

		class data_desc {
		public:
			data_desc( int idx, const QString &pathx, const QString &typex,
			           QtDisplayData *ddx, QtDisplayPanelGui *dpx ) :
				id_(idx), path_(pathx), type_(typex), dd_(ddx), dp_(dpx) { }

			data_desc( int idx ) : id_(idx), dd_(0) { }
			data_desc( ) : id_(0), dd_(0) { }

			int &id( ) {
				return id_;
			}
			int id( ) const {
				return id_;
			}
			QString &path( ) {
				return path_;
			}
			const QString &path( ) const {
				return path_;
			}
			QString &type( ) {
				return type_;
			}
			const QString &type( ) const {
				return type_;
			}
			QtDisplayData *&data( ) {
				return dd_;
			}
			const QtDisplayData *data( ) const {
				return dd_;
			}
			QtDisplayPanelGui *&panel( ) {
				return dp_;
			}
			const QtDisplayPanelGui *panel( ) const {
				return dp_;
			}

		private:
			int id_;
			QString path_;
			QString type_;
			QtDisplayData *dd_;
			QtDisplayPanelGui *dp_;

			// QtDisplayData does not have a copy constructor...
			// wonder if we'll need to copy our descriptor...
			data_desc( const data_desc &other);
			data_desc &operator=( const data_desc &);
		};


		class panel_desc {
		public:

			panel_desc(QtDisplayPanelGui*p) : panel_(p) { }

			std::list<int> &data( ) {
				return data_;
			}
			const std::list<int> &data( ) const {
				return data_;
			}
			QtDisplayPanelGui *&panel( ) {
				return panel_;
			}
			const QtDisplayPanelGui *panel( ) const {
				return panel_;
			}

		private:
			std::list<int> data_;
			QtDisplayPanelGui *panel_;
		};


		QtViewer *viewer_;

		bool printps( QtDisplayPanel *panel, const QString &type, const QString &file, int dpi,
		              const QString &orientation, const QString &media );
		bool printraster( QtDisplayPanel *panel, const QString &type, const QString &file, double scale );
		void  adjusteps( const char *from, const char *to, const QSize &wcmax, const QRect &viewport );

		typedef std::map<int,panel_desc*> panelmap;
		typedef std::map<int,data_desc*> datamap;
		datamap managed_datas;
		panelmap managed_panels;

		int get_id( QtDisplayPanelGui *, QtDisplayData *, const QString &path, const QString &type );
		int get_id( QtDisplayPanelGui * );
		QtDisplayPanelGui *findpanel( int key, bool create=true );
		QtDisplayData *finddata( int key );

	protected:
		QtDisplayPanelGui *create_panel( );
		void load_data( QtDisplayPanelGui *panel, int index );
		void unload_data( QtDisplayPanelGui *panel, int index, bool erase=true );
		QtDisplayPanelGui *erase_panel( int panel );
		void erase_data( int );

	};

}

#endif
