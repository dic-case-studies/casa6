//# grpcPlotServer.qo.h: provides viewer services via dbus
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

#ifndef QTDBUSPLOTSVRADAPTOR_QO_H_
#define QTDBUSPLOTSVRADAPTOR_QO_H_

#include <map>
#include <list>
#include <QVariantMap>
#include <QString>
#include <casaqt/QtPlotServer/QtPlotServer.qo.h>
#include <qwt/qwt_plot_spectrogram.h>
#include <casagrpc/protos/shutdown.grpc.pb.h>
#include <casagrpc/protos/ping.grpc.pb.h>
#include <casagrpc/protos/plotserver.grpc.pb.h>
#include <casagrpc/protos/plotserver_events.grpc.pb.h>
#include <grpc++/grpc++.h>

class QwtPlotItem;
class QDockWidget;

namespace casa {

#if defined(CASATOOLS)
    class QtPlotServer;
#endif

    class QtPlotSvrPanel;

    namespace plotserver {

        class grpcPing : public ::casatools::rpc::Ping::Service {
        public:
            grpcPing( ) { }
            ::grpc::Status now(::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty*);
        };

        class grpcShutdown : public QObject, public ::casatools::rpc::Shutdown::Service {

            Q_OBJECT    //# Allows slot/signal definition.  Must only occur in
                        //# implement/.../*.h files
        public:
            grpcShutdown( QtPlotServer *qtv );
            ::grpc::Status now(::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty*);

        signals:
            void exit_now( );

        private:
            QtPlotServer *plot_;

        };

    }

    class grpcPlotServer : public QObject, public ::rpc::gui::plotserver::Service {

        Q_OBJECT

    protected:
        void qtGO( std::function<void()> );

    public:

         grpcPlotServer( std::string response_uri, QtPlotServer * );
    
        // duplicates the shutdown::now rpc call...
        // duplicates grpcShutdown::now(...)
        ::grpc::Status done( ::grpc::ServerContext*,
                             const ::google::protobuf::Empty*,
                             ::google::protobuf::Empty* );

        ::grpc::Status panel( ::grpc::ServerContext *context,
                              const ::rpc::gui::NewPanel *req,
                              ::rpc::gui::Id *reply );

        ::grpc::Status colors( ::grpc::ClientContext* context,
                               const ::google::protobuf::Empty* request,
                               ::rpc::gui::Colors* response );

         ::grpc::Status colormaps( ::grpc::ClientContext* context,
                                   const ::google::protobuf::Empty* request,
                                   ::rpc::gui::ColorMaps* response );

         ::grpc::Status symbols( ::grpc::ClientContext* context,
                                 const ::google::protobuf::Empty* request,
                                 ::rpc::gui::Symbols* response);

         ::grpc::Status line( ::grpc::ClientContext* context, 
                              const ::rpc::gui::NewLine* request, 
                              ::rpc::gui::Id* response );

         ::grpc::Status scatter( ::grpc::ServerContext* context, 
                                 const ::rpc::gui::NewScatter* request,
                                 ::rpc::gui::Id* response );

         ::grpc::Status histogram( ::grpc::ServerContext* context,
                                   const ::rpc::gui::NewHistogram* request,
                                   ::rpc::gui::Id* response );

         ::grpc::Status raster( ::grpc::ServerContext* context,
                                const ::rpc::gui::NewRaster* request,
                                ::rpc::gui::Id* response );


         ::grpc::Status setlabel( ::grpc::ServerContext* context,
                                  const ::rpc::gui::Label* request,
                                  ::google::protobuf::Empty* response );

         ::grpc::Status erase( ::grpc::ServerContext* context,
                               const ::rpc::gui::Id* request,
                               ::google::protobuf::Empty* response );

         ::grpc::Status close( ::grpc::ServerContext* context,
                               const ::rpc::gui::Id* request,
                               ::google::protobuf::Empty* response );

         ::grpc::Status release( ::grpc::ServerContext* context,
                                 const ::rpc::gui::Id* request,
                                 ::google::protobuf::Empty* response );

         ::grpc::Status hide( ::grpc::ServerContext* context,
                              const ::rpc::gui::Id* request,
                              ::google::protobuf::Empty* response );

         ::grpc::Status show( ::grpc::ServerContext* context,
                              const ::rpc::gui::Id* request,
                              ::google::protobuf::Empty* response );

         ::grpc::Status loaddock( ::grpc::ServerContext* context,
                                  const ::rpc::gui::DockSpec* request,
                                  ::rpc::gui::Id* response);

	class panel_desc {
	public:

	    panel_desc(QtPlotSvrPanel*p) : panel_(p) { }

	    std::list<int> &data( ) { return data_; }
	    const std::list<int> &data( ) const { return data_; }
	    QtPlotSvrPanel *&panel( ) { return panel_; }
	    const QtPlotSvrPanel *panel( ) const { return panel_; }

	    ~panel_desc( ) { }

	private:
	    std::list<int> data_;
	    QtPlotSvrPanel *panel_;
	};


	class data_desc {
	public:

	    data_desc( int index, QtPlotSvrPanel *panel, QwtPlotItem *data ) : id_(index), data_(data), panel_(panel) { }

	    int &id( ) { return id_; }
	    int id( ) const { return id_; }
	    QwtPlotItem *&data( ) { return data_; }
	    const QwtPlotItem *data( ) const { return data_; }
	    QtPlotSvrPanel *&panel( ) { return panel_; }
	    const QtPlotSvrPanel *panel( ) const { return panel_; }
	    ~data_desc( ) { delete data_; }

	private:
	    int id_;
	    QwtPlotItem *data_;
	    QtPlotSvrPanel *panel_;

	    // QtDisplayData does not have a copy constructor...
	    // wonder if we'll need to copy our descriptor...
	    data_desc( const data_desc &other);
	    data_desc &operator=( const data_desc &);
	};

	int get_id( QtPlotSvrPanel *panel );
	int get_id( QtPlotSvrPanel *panel, QwtPlotItem *data );

	void close_everything( );
	void release_everything( );
	QtPlotServer *server;

	typedef std::map<int,panel_desc*> panelmap;
	typedef std::map<int,data_desc*> datamap;
	typedef std::map<int,QDockWidget*> dockmap;

    std::unique_ptr<rpc::gui::plotserver_events::Stub> response_stub;

    std::list<QtPlotSvrPanel*> top_panels;
	panelmap managed_panels;
	datamap managed_datas;
	dockmap managed_docks;

        void release( panelmap::iterator & );

    protected slots:

        void emit_button( QtPlotSvrPanel*, QString name );
        void emit_check( QtPlotSvrPanel*, QString name, int );
        void emit_radio( QtPlotSvrPanel*, QString name, bool );
        void emit_linetext( QtPlotSvrPanel*, QString name, const QString &value );
        void emit_slidevalue( QtPlotSvrPanel*, QString name, int );
        void emit_closing( QtPlotSvrPanel*, bool gone );

    signals:

        void exit_now( );
        void new_op( );

    };

    class grpcPlotServerState {
    public:
        grpcPlotServerState( std::string response_uri, QtPlotServer *v );

        std::string uri;

        std::unique_ptr<grpc::Server> server;

        std::unique_ptr<grpcPlotServer> plot_service;
        std::unique_ptr<plotserver::grpcShutdown> shutdown_service;
        std::unique_ptr<plotserver::grpcPing> ping_service;

        ~grpcPlotServerState( ) {
            if (getenv("GRPC_DEBUG")) {
                fprintf(stdout, "stopping grpc server...\n");
                fflush(stdout);
            }
            if ( server ) server->Shutdown( );
        }

    };

}
#endif
