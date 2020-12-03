//# grpcPlotServer.cc: provide plotting services via DBus
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
#include <future>
#include <unistd.h>

#include <casaqt/QtPlotServer/grpcPlotServer.qo.h>
#include <casaqt/QtPlotServer/QtPlotSvrPanel.qo.h>
#include <casaqt/QtUtilities/QtId.h>
#include <display/QtViewer/QtApp.h>

using namespace casacore;

namespace casa {

    namespace plotserver {

        /**********************************************************************************************************
        **********************************************************************************************************
        *****  grpcPing                                                                                      *****
        **********************************************************************************************************
        **********************************************************************************************************/
        ::grpc::Status grpcPing::now( ::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty* ) {
            static const auto debug = getenv("GRPC_DEBUG");
            if ( debug ) {
                std::cerr << "received ping event..." << std::endl;
                fflush(stderr);
            }
            return grpc::Status::OK;
        }

        /**********************************************************************************************************
        **********************************************************************************************************
        *****  grpcShutdown                                                                                  *****
        **********************************************************************************************************
        **********************************************************************************************************/
        grpcShutdown::grpcShutdown( QtPlotServer *ps ) : plot_(ps) {  }

        ::grpc::Status grpcShutdown::now( ::grpc::ServerContext*, const ::google::protobuf::Empty*, ::google::protobuf::Empty* ) {
            if (getenv("GRPC_DEBUG")) {
                std::cerr << "received shutdown notification..." << std::endl;
                fflush(stderr);
            }
            static auto bye_bye = std::async( std::launch::async, [&]( ) { sleep(2); emit exit_now( ); } );
            return grpc::Status::OK;
        }
    }

    grpcPlotServerState::grpcPlotServerState( std::string response_uri, QtPlotServer *v ) : 
                plot_service(new grpcPlotServer( response_uri, v )),
                shutdown_service(new plotserver::grpcShutdown(v)),
                ping_service(new plotserver::grpcPing( )) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "sending plotserver events to " << response_uri << " (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::done( ::grpc::ServerContext*,
                                         const ::google::protobuf::Empty*,
                                         ::google::protobuf::Empty* ) {

        static const auto debug = getenv("GRPC_DEBUG");

        if (debug) {
            std::cerr << "received grpc done( ) event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        static auto bye_bye = std::async( std::launch::async, [&]( ) { sleep(2); emit exit_now( ); } );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::panel( ::grpc::ServerContext *context,
                                           const ::rpc::gui::NewPanel *req,
                                           ::rpc::gui::Id *reply) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc panel event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        QString title = QString::fromStdString(req->title( ));
        QString xlabel = QString::fromStdString(req->xlabel( ));
        QString ylabel = QString::fromStdString(req->ylabel( ));
        QString window_title = QString::fromStdString(req->window_title( ));
        QString legend = QString::fromStdString(req->legend( ));
        QString zoom = QString::fromStdString(req->zoom( ));
        bool new_row = req->new_row( );
        bool hidden = req->hidden( );

        QList<int> size;
        if ( req->size( ).size( ) > 0 ) {
            for ( auto iter = req->size( ).begin( ); iter != req->size( ).end( ); ++iter ) {
                size.append(*iter);
            }
        }

        int with_panel = req->with_panel( );
        QtPlotSvrPanel *companion = 0;

        if ( with_panel != 0 ) {
            if ( managed_panels.find( with_panel ) == managed_panels.end( ) ) {
                char buf[250];
                sprintf( buf, "companion panel '%d' not found", with_panel );
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
            } else {
                companion = managed_panels.find( with_panel )->second->panel();
            }
        }

        QtPlotSvrPanel *panel = 0;
        std::promise<int> prom;

        qtGO( [&]( ) {
                if (debug) {
                    std::cerr << "creating qt panel... (thread " <<
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }

                panel = server->panel( title, xlabel, ylabel, window_title, size, legend, zoom, companion, new_row );
                if ( hidden ) panel->hide( );
                else panel->show( );

                connect( panel, SIGNAL(button(QtPlotSvrPanel*,QString)), SLOT(emit_button(QtPlotSvrPanel*,QString)) );
                connect( panel, SIGNAL(check(QtPlotSvrPanel*,QString,int)), SLOT(emit_check(QtPlotSvrPanel*,QString,int)) );
                connect( panel, SIGNAL(radio(QtPlotSvrPanel*,QString,bool)), SLOT(emit_radio(QtPlotSvrPanel*,QString,bool)) );
                connect( panel, SIGNAL(linetext(QtPlotSvrPanel*,QString,const QString &)), SLOT(emit_linetext(QtPlotSvrPanel*,QString,const QString &)) );
                connect( panel, SIGNAL(slidevalue(QtPlotSvrPanel*,QString,int)), SLOT(emit_slidevalue(QtPlotSvrPanel*,QString,int)) );
                connect( panel, SIGNAL(closing(QtPlotSvrPanel*,bool)), SLOT(emit_closing(QtPlotSvrPanel*,bool)) );

                prom.set_value(get_id(panel));

            } );

        if (debug) {
            std::cerr << "waiting for grpc panel result... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        auto fut = prom.get_future( );
        fut.wait( );
        auto ret = fut.get( );
        if (debug) {
            std::cerr << "returning grpc panel result: " <<
                ret << " (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        reply->set_id(ret);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::colors( ::grpc::ClientContext* context,
                                           const ::google::protobuf::Empty* request,
                                           ::rpc::gui::Colors* response ) {
        QStringList qcolors = QtPlotSvrPanel::colors( );
        for ( int i = 0; i < qcolors.size(); ++i )
            response->add_names( qcolors.at(i).toUtf8().constData() );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::colormaps( ::grpc::ClientContext* context,
                                              const ::google::protobuf::Empty* request,
                                              ::rpc::gui::ColorMaps* response ) {
        QStringList maps = QtPlotSvrPanel::colormaps( );
        for ( int i = 0; i < maps.size(); ++i )
            response->add_names( maps.at(i).toUtf8().constData() );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::symbols( ::grpc::ClientContext* context,
                                            const ::google::protobuf::Empty* request,
                                            ::rpc::gui::Symbols* response ) {
        QStringList syms = QtPlotSvrPanel::symbols( );
        for ( int i = 0; i < syms.size(); ++i )
            response->add_names( syms.at(i).toUtf8().constData() );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::line( ::grpc::ClientContext* context,
                                         const ::rpc::gui::NewLine* req,
                                         ::rpc::gui::Id* reply ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc line event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        int panel_id = req->panel( ).id( );
        QString color = QString::fromStdString(req->color( ));
        QString label = QString::fromStdString(req->label( ));
        auto xin = req->x( );
        auto yin = req->y( );

        QList<double> x;
        x.reserve(xin.size());
        std::copy(xin.begin(), xin.end(), std::back_inserter(x));
        QList<double> y;
        y.reserve(yin.size());
        std::copy(yin.begin(), yin.end(), std::back_inserter(y));
        
        if ( panel_id != 0 && managed_panels.find( panel_id ) == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel_id );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        panel_desc *paneldesc = 0;
        std::promise<int> prom;

        qtGO( [&]( ) {
                if (debug) {
                    std::cerr << "creating qwt line... (thread " <<
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }

                if ( panel_id == 0 ) {
                    if ( managed_panels.size( ) == 0 ) {
                        QtPlotSvrPanel *panel = server->panel( "", "bottom" );
                        panel_id = get_id(panel);					// adds it to the map of managed panels
                        paneldesc = managed_panels.find( panel_id )->second;
                    } else {
                        paneldesc = managed_panels.begin( )->second;
                    }
                } else {
                    paneldesc = managed_panels.find( panel_id )->second;
                }

                QwtPlotCurve *plot = paneldesc->panel( )->line( x, y, color, label );
                prom.set_value(get_id(paneldesc->panel( ),plot));
            } );

        auto fut = prom.get_future( );
        fut.wait( );
        auto data_id = fut.get( );
        if (debug) {
            std::cerr << "returning grpc line result: " <<
                data_id << " (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        paneldesc->data( ).push_back(data_id);
        reply->set_id(data_id);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::scatter( ::grpc::ServerContext* context, 
                                            const ::rpc::gui::NewScatter* req,
                                            ::rpc::gui::Id* reply ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc scatter event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        int panel_id = req->panel( ).id( );
        QString color = QString::fromStdString(req->color( ));
        QString label = QString::fromStdString(req->label( ));
        QString symbol = QString::fromStdString(req->symbol( ));
        int symbol_size = req->symbol_size( );
        int dot_size = req->dot_size( );

        auto xin = req->x( );
        auto yin = req->y( );

        QList<double> x;
        x.reserve(xin.size());
        std::copy(xin.begin(), xin.end(), std::back_inserter(x));
        QList<double> y;
        y.reserve(yin.size());
        std::copy(yin.begin(), yin.end(), std::back_inserter(y));
        
        if ( panel_id != 0 && managed_panels.find( panel_id ) == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel_id );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        panel_desc *paneldesc = 0;
        std::promise<int> prom;

        qtGO( [&]( ) {

                if (debug) {
                    std::cerr << "creating qwt scatter... (thread " <<
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }

                if ( panel_id == 0 ) {
                    if ( managed_panels.size( ) == 0 ) {
                        QtPlotSvrPanel *panel = server->panel( "", "bottom" );
                        panel_id = get_id(panel);					// adds it to the map of managed panels
                        paneldesc = managed_panels.find( panel_id )->second;
                    } else {
                        paneldesc = managed_panels.begin( )->second;
                    }
                } else {
                    paneldesc = managed_panels.find( panel_id )->second;
                }

                QwtPlotCurve *plot = paneldesc->panel( )->scatter( x, y, color, label, symbol, symbol_size, dot_size );
                prom.set_value(get_id(paneldesc->panel( ),plot));
            } );

        auto fut = prom.get_future( );
        fut.wait( );
        auto data_id = fut.get( );
        if (debug) {
            std::cerr << "returning grpc scatter result: " <<
                data_id << " (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        paneldesc->data( ).push_back(data_id);
        reply->set_id(data_id);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::histogram( ::grpc::ServerContext* context,
                                              const ::rpc::gui::NewHistogram* req,
                                              ::rpc::gui::Id* reply ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc histogram event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        int panel_id = req->panel( ).id( );
        QString color = QString::fromStdString(req->color( ));
        QString label = QString::fromStdString(req->label( ));
        int bins = req->bins( );

        auto valuesin = req->values( );

        QList<double> values;
        values.reserve(valuesin.size());
        std::copy(valuesin.begin(), valuesin.end(), std::back_inserter(values));
        
        if ( panel_id != 0 && managed_panels.find( panel_id ) == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel_id );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        panel_desc *paneldesc = 0;
        std::promise<int> prom;

        qtGO( [&]( ) {
                if (debug) {
                    std::cerr << "creating qwt histogram... (thread " <<
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }

                if ( panel_id == 0 ) {
                    if ( managed_panels.size( ) == 0 ) {
                        QtPlotSvrPanel *panel = server->panel( "", "bottom" );
                        panel_id = get_id(panel);					// adds it to the map of managed panels
                        paneldesc = managed_panels.find( panel_id )->second;
                    } else {
                        paneldesc = managed_panels.begin( )->second;
                    }
                } else {
                    paneldesc = managed_panels.find( panel_id )->second;
                }

                QwtPlotItem *plot = paneldesc->panel( )->histogram( values, bins, color, label );
                prom.set_value(get_id(paneldesc->panel( ),plot));
            } );

        auto fut = prom.get_future( );
        fut.wait( );
        auto data_id = fut.get( );
        if (debug) {
            std::cerr << "returning grpc histogram result: " <<
                data_id << " (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        paneldesc->data( ).push_back(data_id);
        reply->set_id(data_id);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::raster( ::grpc::ServerContext* context,
                                           const ::rpc::gui::NewRaster* req,
                                           ::rpc::gui::Id* reply ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc raster event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        int panel_id = req->panel( ).id( );
        QString colormap = QString::fromStdString(req->colormap( ));
        int sizex = req->sizex( );
        int sizey = req->sizey( );

        auto matrixin = req->matrix( );

        QList<double> matrix;
        matrix.reserve(matrixin.size());
        std::copy(matrixin.begin(), matrixin.end(), std::back_inserter(matrix));
        
        if ( panel_id != 0 && managed_panels.find( panel_id ) == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel_id );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        panel_desc *paneldesc = 0;
        std::promise<int> prom;

        qtGO( [&]( ) {
                if (debug) {
                    std::cerr << "creating qwt raster... (thread " <<
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }

                if ( panel_id == 0 ) {
                    if ( managed_panels.size( ) == 0 ) {
                        QtPlotSvrPanel *panel = server->panel( "", "bottom" );
                        panel_id = get_id(panel);					// adds it to the map of managed panels
                        paneldesc = managed_panels.find( panel_id )->second;
                    } else {
                        paneldesc = managed_panels.begin( )->second;
                    }
                } else {
                    paneldesc = managed_panels.find( panel_id )->second;
                }

                QwtPlotSpectrogram *spect = paneldesc->panel( )->raster(matrix, sizex, sizey, colormap);
                prom.set_value(get_id(paneldesc->panel( ),spect));
            } );

        auto fut = prom.get_future( );
        fut.wait( );
        auto data_id = fut.get( );
        if (debug) {
            std::cerr << "returning grpc raster result: " <<
                data_id << " (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        paneldesc->data( ).push_back(data_id);
        reply->set_id(data_id);
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::setlabel( ::grpc::ServerContext* context,
                                             const ::rpc::gui::Label* req,
                                             ::google::protobuf::Empty* reply ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc setlabel event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        int panel_id = req->panel( ).id( );
        QString xlabel = QString::fromStdString(req->xlabel( ));
        QString ylabel = QString::fromStdString(req->ylabel( ));
        QString title = QString::fromStdString(req->title( ));
        
        if ( panel_id != 0 && managed_panels.find( panel_id ) == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel_id );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        panel_desc *paneldesc = 0;
        if ( panel_id == 0 ) {
            if ( managed_panels.size( ) == 0 ) {
                return grpc::Status(grpc::StatusCode::NOT_FOUND, "no panels have been created");
            } else {
                paneldesc = managed_panels.begin( )->second;
            }
        } else {
            paneldesc = managed_panels.find( panel_id )->second;
        }

        if ( paneldesc == 0 ) return grpc::Status(grpc::StatusCode::NOT_FOUND, "no such panel id");

        std::promise<bool> prom;
        qtGO( [&]( ) {
                if (debug) {
                    std::cerr << "setting label... (thread " <<
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }

                if ( xlabel != "" ) { paneldesc->panel( )->setxlabel( xlabel ); }
                if ( ylabel != "" ) { paneldesc->panel( )->setylabel( ylabel ); }
                if ( title != "" ) { paneldesc->panel( )->settitle( title ); }

                prom.set_value(true);
            } );

        auto fut = prom.get_future( );
        fut.wait( );
        auto status = fut.get( );
        if (debug) {
            std::cerr << "returning after setting labels" <<
                " (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::erase( ::grpc::ServerContext* context,
                                          const ::rpc::gui::Id* req,
                                          ::google::protobuf::Empty* response ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc erase event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        int data = req->id( );
        if ( data == 0 ) {
            for ( datamap::iterator iter = managed_datas.begin();
                  iter != managed_datas.end( ); ++iter ) {
                iter->second->data( )->detach( );
                delete iter->second;
            }
            managed_datas.erase( managed_datas.begin( ), managed_datas.end( ) );
            return grpc::Status::OK;
        }

        datamap::iterator dataiter = managed_datas.find( data );
        if ( dataiter == managed_datas.end( ) ) {
            // see if the id we have maches any panels...
            panelmap::iterator paneliter = managed_panels.find( data );
            if ( paneliter == managed_panels.end( ) ) {
                char buf[250];
                sprintf( buf, "data (or panel) '%d' not found", data );
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
            } else {
                // fish through the data assigned to this panel and erase these plots...
                std::list<int> &datas = paneliter->second->data( );
                for ( std::list<int>::iterator iter = datas.begin( );
                      iter != datas.end( ); ++iter ) {
                    datamap::iterator data = managed_datas.find(*iter);
                    if ( data != managed_datas.end( ) ) {
                        data->second->data( )->detach( );
                        delete data->second;
                        managed_datas.erase(data);
                    }
                }
                std::promise<bool> prom;
                qtGO( [&]( ) {
                        paneliter->second->panel( )->replot( );
                        prom.set_value(true);
                    } );
                
                auto fut = prom.get_future( );
                fut.wait( );
                auto status = fut.get( );
                return grpc::Status::OK;
            }
        }

        for ( panelmap::iterator pi = managed_panels.begin( ); pi != managed_panels.end(); ++pi ) {
            if ( pi->second->panel() == dataiter->second->panel( ) ) {
                std::list<int> &pd = pi->second->data();
                for ( std::list<int>::iterator pdi=pd.begin( ); pdi != pd.end( ); ++pdi ) {
                    if ( *pdi == dataiter->second->id( ) ) {
                        pd.erase(pdi);
                        break;
                    }
                }
                break;
            }
        }

        // erase the one curve that matches
        dataiter->second->data( )->detach( );

        std::promise<bool> prom;
        qtGO( [&]( ) {
                dataiter->second->panel( )->replot( );
                prom.set_value(true);
            } );
                
        auto fut = prom.get_future( );
        fut.wait( );
        auto status = fut.get( );

        delete dataiter->second;
        managed_datas.erase( dataiter );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::close( ::grpc::ServerContext* context,
                                          const ::rpc::gui::Id* req,
                                          ::google::protobuf::Empty* reply ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc close event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        int panel = req->id( );

        if ( panel == 0 ) {
            close_everything( );
            return grpc::Status::OK;
        }
        
        panelmap::iterator iter = managed_panels.find( panel );
        if ( iter == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        // fish through the data assigned to this panel and remove
        // these from our cache...
        std::list<int> &datas = iter->second->data( );
        for ( std::list<int>::iterator iter = datas.begin( );
              iter != datas.end( ); ++iter ) {
            datamap::iterator data = managed_datas.find(*iter);
            if ( data != managed_datas.end( ) ) {
                delete data->second;
                managed_datas.erase(data);
            }
        }

        // now close the panel
        QtPlotSvrPanel *pp = iter->second->panel();
        delete iter->second;
        managed_panels.erase(iter);

        std::promise<bool> prom;
        qtGO( [&]( ) {
                pp->closeMainPanel( );
                prom.set_value(true);
            } );
                
        auto fut = prom.get_future( );
        fut.wait( );
        auto status = fut.get( );

        if (debug) {
            std::cerr << "returning after closing panel" <<
                " (thread " << std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::release( ::grpc::ServerContext* context,
                                            const ::rpc::gui::Id* req,
                                            ::google::protobuf::Empty* reply ) {
        int panel = req->id( );

        if ( panel == 0 ) {
            release_everything( );
            return grpc::Status::OK;
        }

        panelmap::iterator iter = managed_panels.find( panel );
        if ( iter == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        release( iter );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::hide( ::grpc::ServerContext* context,
                                         const ::rpc::gui::Id* req,
                                         ::google::protobuf::Empty* response ) {
        int panel = req->id( );
        if ( panel == 0 ) {
            std::promise<bool> prom;
            qtGO( [&]( ) {
                    for ( panelmap::iterator iter = managed_panels.begin();
                          iter != managed_panels.end(); ++iter )
                        iter->second->panel( )->hide( );
                    prom.set_value(true);
                } );
            auto fut = prom.get_future( );
            fut.wait( );
            return grpc::Status::OK;
        }

        panelmap::iterator iter = managed_panels.find( panel );
        if ( iter == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        std::promise<bool> prom;
        qtGO( [&]( ) {
                iter->second->panel( )->hide( );
                prom.set_value(true);
            } );
        auto fut = prom.get_future( );
        fut.wait( );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::show( ::grpc::ServerContext* context,
                                         const ::rpc::gui::Id* req,
                                         ::google::protobuf::Empty* response ) {
        int panel = req->id( );

        if ( panel == 0 ) {
            std::promise<bool> prom;
            qtGO( [&]( ) {
                    for ( panelmap::iterator iter = managed_panels.begin();
                          iter != managed_panels.end(); ++iter )
                        iter->second->panel( )->show( );
                    prom.set_value(true);
                } );

            auto fut = prom.get_future( );
            fut.wait( );
            return grpc::Status::OK;
        }

        panelmap::iterator iter = managed_panels.find( panel );
        if ( iter == managed_panels.end( ) ) {
            char buf[250];
            sprintf( buf, "panel '%d' not found", panel );
            return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, buf);
        }

        std::promise<bool> prom;
        qtGO( [&]( ) {
                iter->second->panel( )->show( );
                prom.set_value(true);
            } );
        auto fut = prom.get_future( );
        fut.wait( );
        return grpc::Status::OK;
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    ::grpc::Status grpcPlotServer::loaddock( ::grpc::ServerContext* context,
                                             const ::rpc::gui::DockSpec* req,
                                             ::rpc::gui::Id* reply) {
        static const auto debug = getenv("GRPC_DEBUG");
        if (debug) {
            std::cerr << "received grpc loaddock event... (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        int panel = req->panel( ).id( );
        QString file_or_xml = QString::fromStdString(req->file_or_xml( ));
        QString loc = QString::fromStdString(req->loc( ));
        auto dockablein = req->dockable( );
        QStringList dockable;
        dockable.reserve(dockablein.size( ));

        for ( auto iter = dockablein.begin(); iter != dockablein.end(); ++iter ) {
            dockable.append(QString::fromStdString(*iter));
        }

        if ( panel == 0 ) {
            if ( managed_panels.size( ) == 1 ) {
                std::promise<std::pair<QDockWidget*,QString>> prom;
                qtGO( [&]( ) {
                        if (debug) {
                            std::cerr << "creating qwt dock... (thread " <<
                                std::this_thread::get_id() << ")" << std::endl;
                            fflush(stderr);
                        }
                        prom.set_value(managed_panels.begin()->second->panel( )->loaddock( file_or_xml, loc, dockable ));
                    } );

                auto fut = prom.get_future( );
                fut.wait( );
                auto result = fut.get( );
                if ( result.first == 0 ) {
                    return grpc::Status( grpc::StatusCode::INVALID_ARGUMENT, 
                                         result.second == "" ? 
                                         "dock widget creation failure" : 
                                         result.second.toUtf8( ).constData() );
                } else {
                    int id = QtId::get_id( );
                    managed_docks.insert(dockmap::value_type(id,result.first));
                    reply->set_id(id);
                    return grpc::Status::OK;                    
                }
            } else {
                return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, "must specify a panel when multiple panels exist");
            }
        }

        panelmap::iterator iter = managed_panels.find( panel );
        if ( iter == managed_panels.end( ) ) {
            return grpc::Status( grpc::StatusCode::INVALID_ARGUMENT, "could now find requested panel" );
        }

        std::promise<std::pair<QDockWidget*,QString>> prom;
        qtGO( [&]( ) {
                if (debug) {
                    std::cerr << "creating qwt dock... (thread " <<
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }
                prom.set_value(iter->second->panel( )->loaddock( file_or_xml, loc, dockable ));
            } );
        
        auto fut = prom.get_future( );
        fut.wait( );
        auto result = fut.get( );

        if ( result.first == 0 ) {
            return grpc::Status( grpc::StatusCode::INVALID_ARGUMENT, 
                                 result.second == "" ? 
                                         "dock widget creation failure" : 
                                         result.second.toUtf8( ).constData() );
        } else {
            int id = QtId::get_id( );
            managed_docks.insert(dockmap::value_type(id,result.first));
            reply->set_id(id);
            return grpc::Status::OK;
        }
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::qtGO( std::function<void()> func ) {
        { std::lock_guard<std::mutex> exc(server->grpc_queue_mutex);
          server->grpc_queue.push(func);
        }
        emit new_op( );
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::close_everything( ) {
        // Close all open panels, which will exit Qt loop. This does not in
        // itself delete objects or exit the process, although the driver
        // program might do that. Also, some of the panels may have
        // WA_DeleteOnClose set, which would cause their deletion (see, e.g.,
        // QtPlotServer::createDPG()).
        std::promise<bool> prom;
        qtGO( [&]( ) {
                for ( panelmap::iterator iter = managed_panels.begin();
                      iter != managed_panels.end(); ++iter ) {
                    iter->second->panel()->closeMainPanel( );
                }
                prom.set_value(true);
            } );
        auto fut = prom.get_future( );
        fut.wait( );
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::release( panelmap::iterator &iter ) {

        if ( iter == managed_panels.end( ) ) { return; }

        // now close the panel
        QtPlotSvrPanel *pp = iter->second->panel();
        if ( pp->isVisible( ) == false ) {
            // releasing this panel will result in closing it...
            // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
            // fish through the data assigned to this panel and remove
            // these from our cache...
            std::list<int> &datas = iter->second->data( );
            for ( std::list<int>::iterator dataiter = datas.begin( );
                  dataiter != datas.end( ); ++dataiter ) {
                datamap::iterator data = managed_datas.find(*dataiter);
                if ( data != managed_datas.end( ) ) {
                    delete data->second;
                    managed_datas.erase(data);
                }
            }

            delete iter->second;
            managed_panels.erase(iter);
        }
        std::promise<bool> prom;
        qtGO( [&]( ) {
                pp->releaseMainPanel( );
                pp->deleteLater( );
                prom.set_value(true);
            } );
        auto fut = prom.get_future( );
        fut.wait( );
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::release_everything( ) {
        // Close all open panels, which will exit Qt loop. This does not in
        // itself delete objects or exit the process, although the driver
        // program might do that. Also, some of the panels may have
        // WA_DeleteOnClose set, which would cause their deletion (see, e.g.,
        // QtPlotServer::createDPG()).
        for ( panelmap::iterator iter = managed_panels.begin();
              iter != managed_panels.end(); ++iter ) {
            release( iter );
        }
    }

    grpcPlotServer::grpcPlotServer( std::string response_uri, QtPlotServer *s ) : server(s),
            response_stub(rpc::gui::plotserver_events::NewStub( grpc::CreateChannel( response_uri,
                                                                                     grpc::InsecureChannelCredentials( ) ) )) {
        static const auto debug = getenv("GRPC_DEBUG");
        if ( debug ) {
            std::cerr << "grpcPlotServer::grpcPlotServer( ) (thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::emit_closing(QtPlotSvrPanel *panel, bool gone ) {
        static const auto debug = getenv("GRPC_DEBUG");
        grpc::ClientContext context;
        ::rpc::gui::ClosingEvent ce;
        ::google::protobuf::Empty resp;
        ce.mutable_panel( )->set_id(get_id(panel));
        ce.set_gone(gone);
        if ( debug ) {
            std::cerr << "plotserver generating closing event " <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        ::grpc::Status st = response_stub->closing(&context,ce,&resp);
        if ( debug ) {
            std::cerr << "plotserver completed closing event " <<
                ( st.ok( ) ? "[SUCCESS]" : "[FAIL]" ) <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::emit_button(QtPlotSvrPanel *panel, QString name ) {
        static const auto debug = getenv("GRPC_DEBUG");
        grpc::ClientContext context;
        ::rpc::gui::ButtonEvent be;
        ::google::protobuf::Empty resp;
        be.mutable_panel( )->set_id(get_id(panel));
        be.set_name(name.toUtf8().constData());
        if ( debug ) {
            std::cerr << "plotserver generating button event " <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        ::grpc::Status st = response_stub->button(&context,be,&resp);
        if ( debug ) {
            std::cerr << "plotserver completed button event " <<
                ( st.ok( ) ? "[SUCCESS]" : "[FAIL]" ) <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }


    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::emit_radio(QtPlotSvrPanel *panel, QString name, bool state ) {
        static const auto debug = getenv("GRPC_DEBUG");
        grpc::ClientContext context;
        ::rpc::gui::RadioEvent re;
        ::google::protobuf::Empty resp;
        re.mutable_panel( )->set_id(get_id(panel));
        re.set_name(name.toUtf8().constData());
        re.set_state(state);
        if ( debug ) {
            std::cerr << "plotserver generating radio event " <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        ::grpc::Status st = response_stub->radio(&context,re,&resp);
        if ( debug ) {
            std::cerr << "plotserver completed radio event " <<
                ( st.ok( ) ? "[SUCCESS]" : "[FAIL]" ) <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::emit_check(QtPlotSvrPanel *panel, QString name, int state ) {
        static const auto debug = getenv("GRPC_DEBUG");
        grpc::ClientContext context;
        ::rpc::gui::CheckEvent ce;
        ::google::protobuf::Empty resp;
        ce.mutable_panel( )->set_id(get_id(panel));
        ce.set_name(name.toUtf8().constData());
        ce.set_state(state);
        if ( debug ) {
            std::cerr << "plotserver generating check event " <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        ::grpc::Status st = response_stub->check(&context,ce,&resp);
        if ( debug ) {
            std::cerr << "plotserver completed check event " <<
                ( st.ok( ) ? "[SUCCESS]" : "[FAIL]" ) <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::emit_linetext(QtPlotSvrPanel *panel, QString name, const QString &text ) {
        static const auto debug = getenv("GRPC_DEBUG");
        grpc::ClientContext context;
        ::rpc::gui::LineTextEvent lte;
        ::google::protobuf::Empty resp;
        lte.mutable_panel( )->set_id(get_id(panel));
        lte.set_name(name.toUtf8().constData());
        lte.set_text(text.toUtf8().constData());
        if ( debug ) {
            std::cerr << "plotserver generating linetext event " <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        ::grpc::Status st = response_stub->linetext(&context,lte,&resp);
        if ( debug ) {
            std::cerr << "plotserver completed linetext event " <<
                ( st.ok( ) ? "[SUCCESS]" : "[FAIL]" ) <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }

    // -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----
    void grpcPlotServer::emit_slidevalue(QtPlotSvrPanel *panel, QString name, int value ) {
        static const auto debug = getenv("GRPC_DEBUG");
        grpc::ClientContext context;
        ::rpc::gui::SlideValueEvent sve;
        ::google::protobuf::Empty resp;
        sve.mutable_panel( )->set_id(get_id(panel));
        sve.set_name(name.toUtf8().constData());
        sve.set_value(value);
        if ( debug ) {
            std::cerr << "plotserver generating slidevalue event " <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        ::grpc::Status st = response_stub->slidevalue(&context,sve,&resp);
        if ( debug ) {
            std::cerr << "plotserver completed slidevalue event " <<
                ( st.ok( ) ? "[SUCCESS]" : "[FAIL]" ) <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }

    int grpcPlotServer::get_id( QtPlotSvrPanel *panel ) {

	for ( panelmap::iterator iter = managed_panels.begin(); iter != managed_panels.end(); ++iter ) {
	    if ( iter->second->panel() == panel )
		return iter->first;
	}

	int index = QtId::get_id( );
	managed_panels.insert(panelmap::value_type(index, new panel_desc(panel)));
	return index;
    }

    int grpcPlotServer::get_id( QtPlotSvrPanel *panel, QwtPlotItem *data ) {

	for ( datamap::iterator iter = managed_datas.begin(); iter != managed_datas.end(); ++iter ) {
	    if ( iter->second->data() == data )
		return iter->second->id();
	}

	int index = QtId::get_id( );
	data_desc *dd = new data_desc(index, panel, data );
	managed_datas.insert(datamap::value_type(index, dd));
	return index;
    }

}
