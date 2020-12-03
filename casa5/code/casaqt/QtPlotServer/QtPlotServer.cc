//# QtPlotServer.cc: Qt implementation of main 2D plot server display window.
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

#include <thread>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <QCoreApplication>
#include <casaqt/QtPlotServer/QtPlotServer.qo.h>
#include <casaqt/QtPlotServer/QtPlotSvrPanel.qo.h>
#if defined(CASATOOLS)
#include <casaqt/QtPlotServer/grpcPlotServer.qo.h>
#include <casagrpc/protos/registrar.grpc.pb.h>
#else
#include <casaqt/QtPlotServer/QtDBusPlotSvrAdaptor.qo.h>
using namespace casacore;
#endif

namespace casa {

#if ! defined(CASATOOLS)
    QString QtPlotServer::name_;

    const QString &QtPlotServer::name( ) {
	static bool initialized = false;
	if ( ! initialized ) {
	    name_ = "plot_server";
	}
	return name_;
    }
#endif

#if defined(CASATOOLS)
    QtPlotServer::QtPlotServer( const char *server_string, const char *event_sink, const char *logfile ) {
        static const auto debug = getenv("GRPC_DEBUG");

        grpcPlotServerState *state = new grpcPlotServerState( event_sink, this );

        //***
        //*** set up a default address (grpc picks port) and address buffers
        //***
        char address_buf[100];
        constexpr char address_template[] = "0.0.0.0:%d";
        snprintf(address_buf,sizeof(address_buf),address_template,0);
        std::string server_address(address_buf);
        int selected_port = 0;

        //***
        //*** build grpc service
        //***
        grpc::ServerBuilder builder;
        // Listen on the given address without any authentication mechanism.
        builder.AddListeningPort(server_address, grpc::InsecureServerCredentials(), &selected_port);
        // Register "service" as the instance through which we'll communicate with
        // clients. In this case it corresponds to an *synchronous* service.
        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // image viewer service (and currently interactive clean service though this needs
        // to eventually move to a seperate grpc service description which could e.g. be
        // shared with carta
        auto plot_svc = state->plot_service.get( );
        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // all gui operations must happen in the "gui thread" because Qt is not
        // thread-safe... so we need create & result signals and slots
        connect( plot_svc, SIGNAL(new_op( )),
                 this, SLOT(grpc_handle_op( )),
                 Qt::BlockingQueuedConnection );
        connect( plot_svc, SIGNAL(exit_now( )),
                 this, SLOT(quit( )) );
        builder.RegisterService(plot_svc);
        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // shutdown service is used by casatools etc. to notify gui services
        // when the system is shutting down...
        auto shutdown_svc = state->shutdown_service.get( );
        builder.RegisterService(shutdown_svc);
        connect( shutdown_svc, SIGNAL(exit_now( )),
                 this, SLOT(quit( )) );

        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // ping service is used by casatools etc. to check to see if gRPC
        // server is still running...
        auto ping_svc = state->ping_service.get( );
        builder.RegisterService(ping_svc);

        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
        // Launch server...
        state->server = builder.BuildAndStart( );
        if ( selected_port > 0 ) {
            // if an available port can be found, selected_port is set to a value greater than zero
            snprintf(address_buf,sizeof(address_buf),address_template,selected_port);
            state->uri = address_buf;
            if ( debug ) {
                std::cout << "plotserver available at " << state->uri << std::endl;
                fflush(stdout);
            }
            if ( server_string ) {
                if ( access( server_string, F_OK ) == -1 ) {
                    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
                    // create the connection to the registrar for service registration using the uri
                    // provided on the command line...
                    std::unique_ptr<casatools::rpc::Registrar::Stub> proxy =
                        casatools::rpc::Registrar::NewStub(grpc::CreateChannel(server_string, grpc::InsecureChannelCredentials( )));
                    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
                    // register our "shutdown, "image-view" and "interactive-clean" services with
                    // the registrar...
                    casatools::rpc::ServiceId sid;
                    sid.set_uri(state->uri);
                    sid.add_types("ping");
                    sid.add_types("shutdown");
                    sid.add_types("plotserver");
                    grpc::ClientContext context;
                    casatools::rpc::ServiceId accepted_sid;
                    if ( debug ) {
                        std::cout << "registering services with registrar (at " << server_string << ")" << std::endl;
                        fflush(stdout);
                    }
                    ::grpc::Status status = proxy->add(&context,sid,&accepted_sid);
                    if ( ! status.ok( ) ) {
                        // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
                        // if registration was not successful, we exit...
                        std::cerr << "registration failed, exiting..." << std::endl;
                        fflush(stderr);
                        state->server->Shutdown( );
                        QCoreApplication::exit(1);
                        exit(1);
                    }
                    if ( debug ) {
                        std::cout << "accepted service id: ( " << accepted_sid.id( ) << ", " << accepted_sid.uri( ) << ", ";
                        for ( auto i=accepted_sid.types( ).begin( ); i != accepted_sid.types( ).end( ); ++i )
                            std::cout << "'" << (*i) << "' ";
                        std::cout << ")" << std::endl;
                        fflush(stdout);
                    }
                } else {
                    if ( debug ) {
                        std::cout << "writing to file " << server_string << " (thread " <<
                            std::this_thread::get_id() << ")" << std::endl;
                        fflush(stdout);
                    }
                    int fd;
                    if ( (fd = open(server_string, O_WRONLY)) != -1 ) {
                        char uri[strlen(state->uri.c_str( ))+2];
                        sprintf(uri,"%s\n",state->uri.c_str( ));
                        if ( (size_t) write( fd, uri, strlen(uri)) != strlen(uri) ) {
                            qWarning("server failed to write gRPC URI to named pipe...");
                            qFatal("exiting...");
                            exit(1);
                        }
                        close(fd);
                    } else {
                        qWarning("server failed to open gRPC URI named pipe...");
                        qFatal("exiting...");
                        exit(1);
                    }
                }

                // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
                // complete startup
                grpc_.reset(state);
            } else {
                grpc_.reset( );
                if ( debug ) {
                    std::cout << "no registrar provided... skipped registering." << std::endl;
                    fflush(stdout);
                }
            }

        } else {
            if ( debug ) {
                std::cout << "clearing registrar>> " << state->uri << std::endl;
                fflush(stdout);
            }
            grpc_.reset( );
        }

#else
    QtPlotServer::QtPlotServer( const char *dbus_name ) {
	dbus_name_ = (dbus_name ? dbus_name : 0);
	dbus_ = new QtDBusPlotSvrAdaptor(this);
	dbus_->connectToDBus(dbus_name_);
#endif
    }

    QtPlotServer::~QtPlotServer( ) {
#if ! defined(CASATOOLS)
	delete dbus_;
#endif
    }

    QtPlotSvrPanel *QtPlotServer::panel( const QString &title, const QString &xlabel, const QString &ylabel, const QString &window_title,
					 const QList<int> &size, const QString &legend, const QString &zoom, QtPlotSvrPanel *with_panel,
					 bool new_row ) {
	QtPlotSvrPanel *result = new QtPlotSvrPanel(title,xlabel,ylabel,window_title,size,legend,zoom,with_panel,new_row);
	return result;
    }

#if defined(CASATOOLS)
	void QtPlotServer::quit() {
		static const auto debug = getenv("GRPC_DEBUG");
		if ( grpc_ && grpc_->server ) {
			if ( debug ) {
				std::cout << "entering QtPlotServer::quit( )..." << std::endl;
				std::cout << "		  ...shutting down grpc server..." << std::endl;
				fflush(stdout);
			}
			grpc_->server->Shutdown( );
		}

		if ( debug && grpc_->server ) {
			std::cout << "		  ...shutting down qt..." << std::endl;
			fflush(stdout);
		}

		if ( grpc_ && grpc_->server ) {
			QCoreApplication::exit( );
			// calling the system exit( ) here causes immediate
			// shutdown, but does not allow global cleanup...
		}
	}

	void QtPlotServer::grpc_handle_op( ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if ( debug ) {
            std::cerr << "QtPlotServer::grpc_handle_op in progress (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
		std::lock_guard<std::mutex> exc(grpc_queue_mutex);
		if ( ! grpc_queue.empty( ) ) {
			std::function<void()> f = grpc_queue.front( );
			grpc_queue.pop( );
			f( );
		}
	}

#endif

}
