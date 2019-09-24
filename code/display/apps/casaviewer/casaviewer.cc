//# qtviewer.cc:  main program for standalone Qt viewer
//# Copyright (C) 2005,2009,2010
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

#include <ios>
#include <iostream>
#include <casa/aips.h>
#include <casa/Inputs/Input.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Exceptions/Error.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <display/Display/StandAloneDisplayApp.h>
// (Configures pgplot for stand-alone Display Library apps).

#include <casa/Logging/StreamLogSink.h>
#include <casa/Logging/LogSink.h>

#include <display/QtViewer/QtDisplayData.qo.h>
#include <display/QtViewer/QtDisplayPanelGui.qo.h>
#include <display/QtViewer/QtViewer.qo.h>

#include <display/Utilities/Lowlevel.h>

#include <unistd.h>
#ifndef NO_CRASH_REPORTER
#include <stdcasa/StdCasa/CrashReporter.h>
#endif
#include <sys/stat.h>

/*
#include <graphics/X11/X_enter.h>
#include   <QApplication>
#include <graphics/X11/X_exit.h>
*/

#if ! defined(WITHOUT_DBUS)
#include <casadbus/utilities/Diagnostic.h>
#endif
#include <display/DisplayErrors.h>
#include <casacore/casa/System/AppState.h>
#include <algorithm>

#if defined(__APPLE__)
// for executable_path( )
#include <mach-o/dyld.h>
#else
// for executable_path( )
std::string read_link( const std::string &path ) {
    int buffer_size = 128;
    char *buffer = new char[buffer_size+1];
    int nchars = readlink( path.c_str( ), buffer, buffer_size );
    while ( nchars == buffer_size ) {
        buffer_size *= 2;
        delete [] buffer;
        buffer = new char[buffer_size+1];
        nchars = readlink( path.c_str( ), buffer, buffer_size );
    }
    std::string result;
    if ( nchars > 0 ) {
        buffer[nchars] = '\0';
        char *exe = realpath(buffer,NULL);
        result = exe;
        free(exe);
    }
    delete [] buffer;
    return result;
}
#endif

#include <casa/namespace.h>
using namespace casa;

static pid_t manager_root_pid = 0;
static pid_t manager_xvfb_pid = 0;
static bool sigterm_received = false;
static void preprocess_args( int argc, const char *argv[], int &numargs, char **&args,
                             char *&server_string, bool &do_dbus, bool &inital_run,
                             bool &server_startup, bool &daemon,
                             bool &without_gui, bool &persistent, bool &casapy_start,
                             char *&logfile_path );
static void start_manager_root( const char *origname, int numargs, char **args,
                                const char *dbusname, bool without_gui, pid_t root_pid );
static void launch_server( const char *origname, int numargs, char **args,
                           const char *dbusname, bool without_gui,
                           bool persistent, bool casapy_start );
static char *find_xvfb( const char *paths );
static pid_t launch_xvfb( const char *name, pid_t pid, char *&display, char *&authority );

static void exiting_server( int /*sig*/ ) {
	if ( manager_xvfb_pid ) kill( manager_xvfb_pid, SIGKILL );
	exit(0);
}
static void signal_manager_root( int sig ) {
	if ( manager_root_pid && ! sigterm_received ) {
		killpg( manager_root_pid, sig );
		sigterm_received = true;
	}
	signal( sig, signal_manager_root );
	exit(0);
}

static std::string executable_path( ) {
#if defined(__APPLE__)
    uint32_t size = PATH_MAX;
    char *buffer = (char *) malloc(sizeof(char)*size);
    if ( _NSGetExecutablePath(buffer, &size) == -1 ) {
        ++size;
        buffer = (char *) realloc(buffer,sizeof(char)*size);
        if ( _NSGetExecutablePath(buffer, &size) != 0 ) {
            free(buffer);
            fprintf( stderr, "cannot discover path to executable...\n" );
            return "";
        }
    }
    char *exepath = realpath(buffer,NULL);
    std::string result(exepath);
    free(buffer);
    free(exepath);
    return result;
#else
    char buffer[256];
    sprintf( buffer, "/proc/%d/exe", getpid( ) );
    struct stat statbuf;
    if ( lstat( buffer, &statbuf ) == 0 && S_ISLNK(statbuf.st_mode) ) {
        return read_link(buffer);
    }
    return "";
#endif
}


class ViewerApp : public QApplication {
public:
	ViewerApp( int &argc, char **argv, bool gui_enabled ) : QApplication(argc, argv, gui_enabled), viewer_(0) {
#if QT_VERSION >= 0x050000
		local_argc_ = argc;
		local_argv_ = new char*[local_argc_+1];
		for ( int i=0; i < local_argc_; ++i )
			local_argv_[i] = strdup(argv[i]);
		local_argv_[local_argc_] = 0;
#endif
    }
	bool notify( QObject *receiver, QEvent *e );
	void subscribe( QtViewer *v ) {
		viewer_ = v;
	}
private:
	QtViewer *viewer_;

#if QT_VERSION >= 0x050000
	int local_argc_;
	char **local_argv_;
public:
	int argc( ) { return local_argc_; }
	char **argv( ) { return local_argv_; }
#endif
};

bool ViewerApp::notify( QObject *receiver, QEvent *e ) {
	// cap qt event recursion limit at 500 events deep...
	static unsigned long recursion_count = 0;
	if ( recursion_count > 500 ) {
		qWarning( ) << "qt event recursion limit reached...";
		return false;
	}

	try {
		recursion_count++;

		// notify QtViewer when application is activated/deactivated, e.g. when OSX switches to
		// "mission control", allowing the viewer to recognize that the mouse cursor has left...
		if ( e->type() == QEvent::ApplicationActivate || e->type() == QEvent::ApplicationDeactivate )
			if ( viewer_ != 0 ) viewer_->activate( e->type() == QEvent::ApplicationActivate ? true : false );

		bool result = QApplication::notify(receiver,e);
		if ( recursion_count != 0 ) recursion_count--;
		return result;
	} catch ( AipsError e ) {
		qWarning( ) << "unhandled exception:" << e.getMesg().c_str();
		if ( recursion_count != 0 ) recursion_count--;
		return false;
	} catch ( viewer::internal_error e ) {
		qWarning( ) << "internal exception:" << e.what( );
		if ( recursion_count != 0 ) recursion_count--;
		return false;
	} catch ( std::exception e ) {
		qWarning( ) << "std exception:" << e.what( );
		if ( recursion_count != 0 ) recursion_count--;
		return false;
	} catch ( ... ) {
		qWarning("unhandled exception... ");
		qDebug() << "from" << receiver->objectName() << "from event type" << e->type();
		qFatal("exiting...");
		exit(1);
	}
	return true;
}


class ViewerDataState: public casacore::AppState {
public:

    ViewerDataState(const std::list<std::string> &path ) : data_path(path) { }
    virtual bool initialized( ) const { return true; }
    virtual std::list<std::string> dataPath( ) const { return data_path; }
private:
    std::list<std::string> data_path;
};


int main( int argc, const char *argv[] ) {

    std::string exepath(executable_path( ));

#ifndef NO_CRASH_REPORTER
    CrashReporter::initializeFromApplication(argv[0]);
#endif
#if ! defined(WITHOUT_DBUS)
	casa::dbus::diagnostic.argv( argc, argv );
#endif

	bool server_startup = false;
    bool daemon = false;
	bool without_gui = false;
	bool persistent = false;
	bool casapy_start = false;
	char *server_string = 0;
	bool with_dbus = false;
	char *logfile_path = 0;
	bool initial_run = false;

	char **args;
	int numargs;

	signal( SIGTERM, exiting_server );

	// On Mac OS X by default, Qt swaps the Control and Meta (Command) keys (i.e., whenever
	// Control is pressed, Qt sends Meta, and whenever Meta is pressed Control is sent).
	// When this attribute is true, Qt will not do the flip. QKeySequence::StandardShortcuts
	// will also flip accordingly (i.e., QKeySequence::Copy will be Command+C on the keyboard
	// regardless of the value set, though what is output for
	// QKeySequence::toString(QKeySequence::PortableText) will be different).
	//
	// This avoids this swapping, making the behavior more consistent when running on an X11
	// system displaying on OSX or just when using it with OSX.
	//
	// This option makes *NO*DIFFERENCE* for the context menu on the "Point Button Tool".
	// Thu Aug 30 10:32:15 EDT 2012 <drs>
	//
//     QCoreApplication::setAttribute(Qt::AA_MacDontSwapCtrlAndMeta);

	preprocess_args( argc, argv, numargs, args, server_string, with_dbus,
	                 initial_run, server_startup, daemon, without_gui,
                     persistent, casapy_start, logfile_path );

	//
	// configure datapath for casacore and colormaps...
	//
	auto ends_with = []( const std::string& str, const std::string& ending ) {
		return ( str.size( ) >= ending.size( ) ) && equal( ending.rbegin( ), ending.rend( ), str.rbegin( ) );
	};
	//
	// on linux argv[0] will be "casaviewer"
	// on OSX with the viewer packaged with CASA argv[0] will be "CASAViewer"
	// on OSX with the viewer packaged separately argv[0] will be "CASAviewer"
	if ( ends_with(exepath, "Contents/MacOS/CASAviewer") ||
         ends_with(exepath, "Contents/MacOS/casaviewer") ) {
		// initialize CASAviewer app data...
		if ( ! casacore::AppStateSource::fetch( ).initialized( ) ) {
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// Mac OSX  --  path is specific to package format
			// -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - 
			// initialize CASAviewer app data...
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// generate path to data...
			std::string datapath(exepath);
			datapath.erase( datapath.end( ) -  16, datapath.end( ) );
			std::string pgplotpath = datapath;             // save for later...
			std::string pluginpath = datapath;             // save for later...
			datapath += "Resources/casa-data";
			// initialize casacore...
			std::list<std::string> datadirs;
			datadirs.push_back(datapath);
			casacore::AppStateSource::initialize(new ViewerDataState(datadirs));
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// initialize CASAviewer app data...
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// configure pgpplot...
			std::string rgbpath = std::string("PGPLOT_RGB=") + pgplotpath + "Resources/pgplot/rgb.txt";
			std::string fontpath = std::string("PGPLOT_FONT=") + pgplotpath + "Resources/pgplot/grfont.dat";
			putenv(strdup(rgbpath.c_str( )));
			putenv(strdup(fontpath.c_str( )));
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			// set up Qt Plugin Path
			// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
			pluginpath += "Plugins";
			QCoreApplication::addLibraryPath(QString(pluginpath.c_str( )));
		}

	} else if ( ends_with(exepath, "/AppRun") ||
                ends_with(exepath, "/CASAviewer.app/usr/bin/CASAviewer") ||
                ends_with(exepath, "/casaviewer.app/usr/bin/casaviewer") ) {

		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		// linux  --  path is specific to package format
		//
		//    .../AppRun implies AppImage bash script startup, e.g. from an unpacked AppImage
		//    .../CASAviewer.app/usr/bin/CASAviewer implies debugging or running from the
		//                                          build tree before it has been packaged
		//
		// -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - 
		// initialize CASAviewer app data...
		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		// generate path to data...
		bool packed_app = ends_with(exepath, "/AppRun");
		std::string datapath(exepath);
		//     packed_app -> .../AppRun
		// not packed_app -> .../CASAviewer.app/usr/bin/CASAviewer
		datapath.erase( datapath.end( ) -  (packed_app ? 6 : 18), datapath.end( ) );
		std::string pgplotpath = datapath;			   // save for later...
		std::string pluginpath = datapath;			   // save for later...
		datapath += "data";
		// initialize casacore...
		std::list<std::string> datadirs;
		datadirs.push_back(datapath);
		casacore::AppStateSource::initialize(new ViewerDataState(datadirs));
		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		// initialize CASAviewer app data...
		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		// configure pgpplot...
		std::string rgbpath = std::string("PGPLOT_RGB=") + pgplotpath + "usr/lib/pgplot/rgb.txt";
		std::string fontpath = std::string("PGPLOT_FONT=") + pgplotpath + "usr/lib/pgplot/grfont.dat";
		putenv(strdup(rgbpath.c_str( )));
		putenv(strdup(fontpath.c_str( )));
		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		// set up Qt Plugin Path
		// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
		pluginpath += "usr/lib/plugins";
		QCoreApplication::addLibraryPath(QString(pluginpath.c_str( )));
	}

	//
	// setup casa logging's global sink, if the user supplied a path...
	//
	if ( logfile_path ) {
		FILE *file = fopen(logfile_path,"a");
		if ( file ) {
			fclose(file);
			casacore::LogSinkInterface *sink = new casacore::StreamLogSink(new ofstream(logfile_path,std::ios_base::app));
			casacore::LogSink::globalSink(sink);
		}
	}

	if ( (server_startup || without_gui) && initial_run ) {
        if ( daemon ) {
            launch_server( argv[0], numargs, args, server_string, without_gui,
                           persistent, casapy_start );
        } else {
            start_manager_root( argv[0], numargs, args, server_string, without_gui, getpid( ) );
        }
		exit(0);
	}

	INITIALIZE_PGPLOT
	try {

		ViewerApp qapp(numargs, args, true);

		// if it's a server, stick around even if all windows are closed...
		if ( server_startup ) {
			qapp.setQuitOnLastWindowClosed(false);
		}

		String	   filename    = "",
		           displaytype = "",
		           datatype    = "",
		           arg2        = "",
		           arg3        = "";


		Int narg;

#ifndef AIPS_DARWIN
		narg = qapp.argc();
		if(narg>1) filename = qapp.argv()[1];
		if(narg>2) arg2     = qapp.argv()[2];
		if(narg>3) arg3     = qapp.argv()[3];
#else
		narg = numargs;
		if(narg>1) filename = args[1];
		if(narg>2) arg2     = args[2];
		if(narg>3) arg3     = args[3];
#endif

		// Workaround for python task's "empty parameter" disability....
		if(filename==".") filename="";

		if ( filename != "" && filename[0] == '-' && filename[1] == '-' ) {
			struct stat statbuf;
			if ( stat( filename.c_str( ), &statbuf ) == -1 ) {
				filename = "";
			}
		}

		// Pass along the remaining arguments to QtViewer...
		// instead of littering the ctor arguments...
		std::list<std::string> stdargs;
		for ( int arg_index=0; args[arg_index]; ++arg_index )
			stdargs.push_back(args[arg_index]);

		QtViewer* v = new QtViewer( stdargs, server_startup || with_dbus, server_string );
		qapp.subscribe(v);

		if ( ! server_startup ) {

			// define the panel
			QtDisplayPanelGui* dpg;

			dpg = v->createDPG( );

			QtDisplayData* qdd = 0;

			// Data files are now typed automatically (see v_->filetype(filename),
			// below; e.g.: "image" or "ms").  arg2 need be used only to specify a
			// displaytype, and then only when it is not the default displaytype
			// for the datatype (e.g.  viewer "my.im", "contour" ).
			//
			// The user can enter an lel expression in place of filename, but such
			// an expression _cannot_ be automatically typed.  In this case the user
			// must have "lel" in arg2 (or in arg3: the only case where arg3 is vaguely
			// useful is something like:
			//
			//   casaviewer "'my.im'-'other.im'"  contour  lel
			//
			// arg3 is not even offered in the viewer casapy task).
			//
			// The logic below allows displaytypes or datatypes to be entered in
			// either order, and for old datatypes to be used other than "lel" (these
			// are simply ignored).  This allows old (deprecated) parameter usage in
			// scripts (such as viewer("my.ms", "ms")) to continue to be understood.
			//
			// However, the canonical 'allowed' parameter set (per user documentation)
			// is now just:
			//
			//   viewer [filename [displaytype]]
			//

			if(filename!="") {

				Bool tryDDcreate = true;

				if(arg3=="lel" || arg2=="lel") {

					// (this means that first ('filename') parameter is supposed to
					// contain a valid lel (image expression) string; this is advanced
					// (and undocumented) parameter usage).

					datatype = "lel";
					displaytype = (arg3=="lel")? arg2 : arg3;
					v->dataDisplaysAs(datatype, displaytype);
				} else {

					datatype = v->filetype(filename);


					if(datatype=="restore") {

						// filename is a restore file.

						tryDDcreate = false;

						dpg->restorePanelState(filename);
					}

					else {

						if(datatype=="nonexistent") {
							cerr << "***Can't find  " << filename << "***" << endl;
							tryDDcreate = false;
						}

						if(datatype=="unknown") {
							cerr << "***Unknown file type for  " << filename << "***" << endl;
							tryDDcreate = false;
						}

						// filename names a normal data file.  If user has passed a valid
						// displaytype in either arg2 or arg3, use it; otherwise, the
						// default displaytype for datatype will be inserted.

						displaytype = arg2;
						if(!v->dataDisplaysAs(datatype, displaytype)) {
							displaytype = arg3;
							v->dataDisplaysAs(datatype, displaytype);
						}
					}
				}


				if(tryDDcreate) {

					qdd = dpg->createDD(filename, datatype, displaytype, true,
										-1, false, false, false/*, ddo*/ );
					dpg->addedData( QString::fromStdString(displaytype), qdd );


					if(qdd==0)  cerr << dpg->errMsg() << endl;
				}
			}

			dpg->show();

			if( dpg->isEmptyDD() ) dpg->showDataManager();

		}

		Int stat = qapp.exec();

		//delete dpg;		// Used to lead to crash (double-deletion
		// of MWCTools); should work now.

		delete v;

		// cerr<<"Normal exit -- status: "<<stat<<endl;	//#diag

		return stat;
	}

	catch (const casacore::AipsError& err) {
		cerr<<"**"<<err.getMesg()<<endl;
	} catch (...) {
		cerr<<"**non-AipsError exception**"<<endl;
	}

}

// processes argv into args stripping out the the viewer setup flags... numargs has the number
// of args, and the last arg (not included in numargs count) is null (for execvp)
static void preprocess_args( int argc, const char *argv[], int &numargs, char **&args,
                             char *&server_string, bool &with_dbus, bool &initial_run,
                             bool &server_startup, bool &daemon, bool &without_gui, bool &persistent,
                             bool &casapy_start, char *&logfile_path ) {

	without_gui = false;
	persistent = false;
	casapy_start = false;
	server_string = 0;

	initial_run = (isdigit(argv[0][0]) ? false : true);

	for ( int x = 0; x < argc; ++x ) {
		if ( ! strncmp(argv[x],"--nogui",7) ) {
			without_gui = true;
			if ( argv[x][7] == '=' ) {
				char *name = strdup( &argv[x][8] );
				if ( strlen(name) <= 0 ) {
					free( name );
#if defined(WITHOUT_DBUS)
					qWarning("no gRPC registry provided with '--server=...'");
					qFatal("exiting...");
					exit(1);
#endif
				} else {
					server_string = name;
				}
			}
		} else if ( ! strncmp(argv[x],"--server",8) ) {
			server_startup = true;
			if ( argv[x][8] == '=' ) {
				char *name = strdup( &argv[x][9] );
				if ( strlen(name) <= 0 ) {
					free( name );
#if defined(WITHOUT_DBUS)
					qWarning("no gRPC registry provided with '--server=...'");
					qFatal("exiting...");
					exit(1);
#endif
				} else {
					server_string = name;
				}
			}
#if ! defined(WITHOUT_DBUS)
        } else if ( ! strncmp(argv[x],"--dbusname",10) ) {
			if ( argv[x][10] == '=' ) {
				char *name = strdup( &argv[x][11] );
				if ( strlen(name) <= 0 ) {
					free( name );
				} else {
					server_string = name;
				}
			} else if ( x + 1 < argc ) {
				server_string = strdup(argv[++x]);
			}
		} else if ( ! strcmp(argv[x],"--daemon") ) {
			daemon = true;
		} else if ( ! strcmp(argv[x],"--persist") ) {
			persistent = true;
		} else if ( ! strcmp(argv[x],"--casapy") ) {
			casapy_start = true;
		} else if ( ! strcmp(argv[x],"--dbus") ) {
			with_dbus = true;
#endif
		}  else if ( ! strncmp(argv[x],"--rcdir",7) ) {
			if ( argv[x][7] == '=' ) {
				viewer::setrcDir(&argv[x][8]);
			} else if ( x + 1 < argc ) {
				viewer::setrcDir(argv[++x]);
			}
		} else if ( ! strncmp(argv[x],"--casalogfile",13) ) {
			if ( argv[x][13] == '=' ) {
				char *file = strdup( &argv[x][14] );
				if ( strlen(file) <= 0 ) {
					free( file );
				} else {
					logfile_path = file;
				}
			} else if ( x + 1 < argc ) {
				logfile_path = strdup(argv[++x]);
			}
		} else if ( ! strncmp(argv[x],"--xvfb-pid=",11) ) {
			sscanf( argv[x], "--xvfb-pid=%u", &manager_xvfb_pid );
		}
	}

	char *orig_name = strdup(argv[0]);
	char *name = orig_name + strlen(orig_name) - 1;
	while ( name > orig_name && *name != '/' ) --name;
	if ( *name == '/' ) ++name;

	// pre-process the command line arguments (for now), it looks like
	// the current scheme is to hard-wire a limited set of command line
	// parameters... should be reworked in the future.
	args = (char**) malloc(sizeof(char*)*(argc+2));

	numargs = 0;
	if ( name != orig_name ) {
		args[numargs++] = strdup(name);
		free(orig_name);
	} else {
		args[numargs++] = orig_name;
	}


	if ( initial_run && (without_gui || server_startup) ) {
		// --nogui imples --server and for consistency, --server also forks child
		char *arg;
		if ( server_string ) {
			arg = (char*) malloc( sizeof(char)*(strlen(server_string)+15) );
			sprintf( arg, "--server=%s", server_string );
		} else {
#if defined(WITHOUT_DBUS)
			qWarning("no gRPC registry provided with '--server=...'");
			qFatal("exiting...");
			exit(1);
#else
			arg = strdup( "--server" );
#endif
        }
		args[numargs++] = arg;
	}

	for ( int x = 1; x < argc; ++x ) {
		if ( ! strcmp(argv[x], "--dbusname") ) {
			++x;
		} else if ( strncmp( argv[x], "--server", 8 ) &&
		            strncmp( argv[x], "--nogui", 7 ) &&
		            strcmp( argv[x], "--persist" ) &&
		            strcmp( argv[x], "--casapy" ) &&
		            strcmp( argv[x], "--eso3d" ) &&
		            strncmp( argv[x],"--dbusname",10) )
			args[numargs++] = strdup(argv[x]);
	}
	args[numargs] = 0;

	std::list<char*> files;
	std::list<char*> others;
	std::list<char*> flags;
	struct stat statbuf;
	for ( int x = 1; x < numargs; ++x ) {
		if ( ! stat( args[x], &statbuf ) ) {
			files.push_back(args[x]);
		} else if ( args[x][0] == '-' && args[x][1] == '-' ) {
			flags.push_back(args[x]);
		} else {
			others.push_back(args[x]);
		}
	}

	int offset = 1;
	for ( std::list<char*>::iterator iter = files.begin( );
	        iter != files.end( ); ++iter ) {
		args[offset++] = *iter;
	}
	for ( std::list<char*>::iterator iter = others.begin( );
	        iter != others.end( ); ++iter ) {
		args[offset++] = *iter;
	}
	for ( std::list<char*>::iterator iter = flags.begin( );
	        iter != flags.end( ); ++iter ) {
		args[offset++] = *iter;
	}
}

void start_manager_root( const char *origname, int numargs, char **args, const char */*dbusname*/,
                         bool without_gui, pid_t root_pid ) {

	char *display = 0;
	char *authority = 0;
	if ( without_gui ) {
		manager_xvfb_pid = launch_xvfb( args[0], root_pid, display, authority );
		sleep(2);
	}

	char *name = args[0];

	if ( display ) {
		putenv(display);
	}
	if ( authority ) {
		putenv( authority );
	}

	char *new_name = (char*) malloc( sizeof(char)*(strlen(name)+25) );
	sprintf( new_name, "%d-%s-svr", root_pid, name );
	args[0] = new_name;

	free(name);

	char **newargs = 0;
	if ( manager_xvfb_pid == 0 )
		newargs = args;
	else {
		newargs = (char**) malloc( sizeof(char*)*(numargs+2) );
		for (int i=0; i < numargs; ++i)
			newargs[i] = args[i];
		char xvfbb[124];
		sprintf( xvfbb, "--xvfb-pid=%u", manager_xvfb_pid );
		newargs[numargs] = xvfbb;
		newargs[numargs+1] = 0;
	}

	execvp( origname, newargs );
}

void launch_server( const char *origname, int numargs, char **args,
                    const char *dbusname, bool without_gui,
                    bool persistent, bool casapy_start ) {

	// ______________________________________________________________________________________________________
	//
	// root (receives sigterm & is part of casapy process group)
	//   |
	//   +-- child-root (sets new process group, forks off child to be independent, passes back pid of child, exits)
	//         |
	//         +-- manager-root (sets up xvfb if necessary, starts lowly_viewer, waits for lowly_viewer to exit or SIGTERM from root)
	//                |
	//                +-- xvfb (without_gui only)
	//                |
	//                +-- actual lowly_viewer (perhaps with DISPLAY set)
	//
	// ______________________________________________________________________________________________________

	pid_t root_pid = getpid();

	// don't let ^C [from casapy] kill the viewer...
	signal(SIGINT,SIG_IGN);
	// sent when STDIN/STDOUT is closed...
	// e.g. when running regressions...
	//      or when using expect...
	signal(SIGHUP,SIG_IGN);

	int vio[2];
	if ( pipe(vio) < 0 ) {
		perror( "child-root pipe" );
		exit(1);
	}

	pid_t child_root = fork( );
	if ( child_root == 0 ) {
		signal(SIGINT,SIG_IGN);
		close(vio[0]);

#ifdef SIGTTOU
		signal(SIGTTOU, SIG_IGN);
#endif
#ifdef SIGTTIN
		signal(SIGTTIN, SIG_IGN);
#endif
#ifdef SIGTSTP
		signal(SIGTSTP, SIG_IGN);
#endif

		if ( fork( ) != 0 ) {
			exit(0);
		}

		if ( setpgrp( ) == -1 ) {
			fprintf( stderr, "trying to do the setpgrp...%d...\n", getpid() );
			perror( "setpgrp issue" );
		}
// 	if ( setsid( ) == -1 ) {
// 	    fprintf( stderr, "trying to do the setsid...%d...\n", getpid() );
// 	    perror( "setsid issue" );
// 	}

		char buffer[50];
		sprintf( buffer,"%d", getpid() );
		write( vio[1], buffer, strlen(buffer) + 1 );
		close(vio[1]);

		start_manager_root( origname, numargs, args, dbusname, without_gui, root_pid );
		fprintf( stderr, "start_manager_root( ) should not have returned...%d...\n", getpid() );
		exit(1);
	}

	char buffer[50];
	close(vio[1]);
	read( vio[0], buffer, 50 );
	manager_root_pid = (pid_t) atoi(buffer);
	close(vio[0]);

	int child_root_status;
	waitpid( child_root, &child_root_status, 0 );
	if ( persistent || ! casapy_start ) {
		exit(0);
	}

	// catch the exit signal, and signal the manager_root
	signal(SIGTERM, signal_manager_root);

	while ( ! sigterm_received ) {
		sleep(5);
	}
}

char *find_xvfb( const char *paths ) {

	if ( !paths ) return 0;

	const char *pathptr = paths;
	char buffer[PATH_MAX];		// security: too long paths could result
	// in exec'ing an unintended binary
	char *result = 0;
	const char *p = paths;
	do {
		if ( *p == ':' || *p == '\0' ) {
			if ( p == pathptr ) {
				pathptr = p+1;
				continue;
			}
			if ((p - pathptr + 6) > (int)sizeof(buffer) ) {
				pathptr = p+1;
				continue;
			}
			memcpy(buffer, pathptr, p - pathptr);
			memcpy(&buffer[p-pathptr], "/Xvfb", 6);
			if ( access( buffer, X_OK ) == 0 ) {
				result = strdup(buffer);
				break;
			}
			pathptr = p+1;
		}
	} while ( *p++ );

	return result;
}

pid_t launch_xvfb( const char *name, pid_t pid, char *&display, char *&authority ) {

	pid_t child_xvfb = 0;

#ifdef Q_WS_X11

	char *home = getenv("HOME");

	display = 0;
	authority = 0;

	char *xvfb_path = find_xvfb( getenv("PATH") );
	if ( xvfb_path == 0 ) {
		xvfb_path = find_xvfb( "/usr/X11R6/bin:/usr/X11/bin:/usr/bin" );
		if ( xvfb_path == 0 ) {
			fprintf( stderr, "casaviewer: could not start Xvfb\n" );
			exit(1);
		}
	}

#if defined(__APPLE___)
	srandomdev( );
#else
	union {
		void *foo;
		unsigned bar;
	};
	foo = &home;
	srandom(bar);
#endif

	char xauth[33];
	for ( int x=0; x < 32; x += 8 ) {
		sprintf( &xauth[x], "%08x", (int) random( ) );
	}
	xauth[32] = '\0';


	if ( ! home ) {
		fprintf( stderr, "HOME is not defined in your environment\n" );
		exit(1);
	}

	char *xvfb_name = (char*) malloc( sizeof(char)*(strlen(name)+50) );
	sprintf( xvfb_name, "%d-%s-xvfb", pid, name );

	authority = (char*) malloc( sizeof(char)*(strlen(home)+160) );
	sprintf( authority, "%s/.casa", home );
	viewer::create_dir( authority );
	sprintf( authority, "%s/.casa/xauthority", home );

	const int display_start=6;
	int display_num;
	for ( display_num=display_start; display_num < (display_start+80); ++display_num ) {

		int io[2];
		if ( pipe(io) < 0 ) {
			perror( "xvfb pipe" );
			exit(1);
		}

		child_xvfb = fork( );
		if ( ! child_xvfb ) {

			close(io[0]);
			// don't let ^C [from casapy] kill the viewer...
// 	    signal(SIGINT,SIG_IGN);

// 	    close( fileno(stdin) );
// 	    dup2( io[1], fileno(stdout) );
			dup2( io[1], fileno(stderr) );

			char screen_arg[50];
			sprintf( screen_arg, ":%d", display_num );

			execl( xvfb_path, xvfb_name, screen_arg, "-screen", "0","2048x2048x24+32",
			       "-auth", authority, (char*) 0 );
			fprintf( stderr, "exec of Xvfb should not have returned...%d...\n", getpid() );
			perror( "virtual frame buffer" );
			exit(1);
		}

		close(io[1]);
		bool failure = false;
		fd_set read_set;
		struct timeval wait_time = { 3, 0 };
		while ( 1 ) {
			FD_ZERO( &read_set );
			FD_SET( io[0], &read_set );
			int result = select( io[0] + 1, &read_set, 0, 0, &wait_time );
			if ( result > 0 ) {
				char buffer[1024];
				int nbytes = read( io[0], buffer, sizeof(buffer)-1 );
				if ( nbytes > 0 ) {
					buffer[nbytes] = '\0';
					// HANDLE: 'FreeFontPath: FPE "unix/:7100" refcount is 2, should be 1; fixing.'
					if ( strstr( buffer, "already active" ) ||
					        strstr( buffer, "server already running" ) ||
					        strstr( buffer, "SocketCreateListener() failed") ) {
						failure = true;
						int status;
						waitpid(child_xvfb, &status, 0);
						break;
					}
				} else if ( nbytes == 0 ) {
					fprintf( stderr, "OK, no bytes read while starting xvfb, timeout?...%d...\n", getpid() );
					failure = true;
					break;
				} else if ( nbytes < 0 ) {
					// some error has occurred
					fprintf( stderr, "oops something strange occurred when starting up xvfb...%d...\n", getpid() );
					perror( "xvfb launch" );
					exit(1);
				}
			} else if ( result == 0 ) {
				// timeout, things must be OK...
				break;
			} else if ( result < 0 ) {
				fprintf( stderr, "Oops, some error occurred... try again...%d...\n", getpid() );
				perror( "what was the problem... " );
				continue;
			}

		}

		if ( failure == false )
			break;

	}

	char cmd[1024];
	sprintf( cmd, "xauth -f %s add :%d . %s", authority, display_num, xauth );
	int status = system( cmd );
	if ( status != 0 ) {
		perror( "xauth setup" );
	}

	display = (char*) malloc( sizeof(char) * 50 );
	sprintf( display, "DISPLAY=:%d.0", display_num );
	sprintf( authority, "XAUTHORITY=%s/.casa/xauthority", home );

#endif
	return child_xvfb;
}
