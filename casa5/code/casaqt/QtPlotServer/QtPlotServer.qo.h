//# QtPlotServer.qo.h: main class for plotting; panels are created from here
//# with surrounding Gui functionality
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
//# $Id: QtPlotServer.qo.h,v 1.7 2006/10/10 21:42:05 dking Exp $

#ifndef QTDBUSPLOTSERVER_QO_H_
#define QTDBUSPLOTSERVER_QO_H_

#include <QObject>
#include <memory>
#if defined(CASATOOLS)
#include <queue>
#include <mutex>
#endif

namespace casa {

#if defined(CASATOOLS)
	class grpcPlotServerState;
#else
    class QtDBusPlotSvrAdaptor;
#endif
    class QtPlotSvrPanel;

    class QtPlotServer : public QObject {
    Q_OBJECT
	public:

	    ~QtPlotServer( );
#if defined(CASATOOLS)
	    QtPlotServer( const char *server_string, const char *event_sink, const char *logfile );
#else
	    QtPlotServer( const char *dbus_name=0 );
	    // name used to initialize connection to dbus
	    static const QString &name( );
#endif

	    QtPlotSvrPanel *panel( const QString &title, const QString &xlabel="", const QString &ylabel="", const QString &window_title="",
				   const QList<int> &size=QList<int>( ), const QString &legend="bottom", const QString &zoom="bottom",
				   QtPlotSvrPanel *with_panel=0, bool new_row=false  );

#if defined(CASATOOLS)
	public slots:
		// Exits Qt loop.  (Note that the loop can be restarted (and is, in
		// interactive clean, e.g.), with existing widgets intact.  This
		// does not in itself delete objects or exit the process, although
		// the driver program might do that).  Also, some of the panels may
		// have WA_DeleteOnClose set, which would cause their deletion (see,
		// e.g., createDPG()).
		virtual void quit();

		void grpc_handle_op( );

	private:
		std::shared_ptr<grpcPlotServerState> grpc_;

	public:
		std::mutex grpc_queue_mutex;
		std::queue<std::function<void()>> grpc_queue;
#else
	private:
	    static QString name_;
	    QString dbus_name_;
	    QtDBusPlotSvrAdaptor* dbus_;
#endif

  };

}

#endif
