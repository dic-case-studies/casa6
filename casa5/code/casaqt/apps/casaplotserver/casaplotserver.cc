#include <math.h>
#include <unistd.h>
#include <signal.h>
#include <QApplication>
#include <casaqt/QtPlotServer/QtPlotServer.qo.h>
#include <casaqt/QtPlotServer/QtPlotSvrPanel.qo.h>

#if defined(CASATOOLS)
static void preprocess_args( int argc, const char *argv[], int &numargs, char **&args,
                             char *&server_string, char *&event_sink, char *&logfile );
#else
static void preprocess_args( int argc, const char *argv[], int &numargs, char **&args,
			     char *&dbus_name, bool &casapy_start );
#endif


int main( int argc, const char *argv[] ) {

    // don't let ^C [from casapy] kill the plot server...

    char **args;
    int numargs;

#if defined(CASATOOLS)
    char *server_string;
    char *event_sink;
    char *logfile;
    preprocess_args( argc, argv, numargs, args, server_string, event_sink, logfile );
	signal(SIGINT,SIG_IGN);

//  QApplication qapp(qargc, qargs);
    QApplication qapp(numargs, args);
    casa::QtPlotServer plot_server(server_string,event_sink,logfile);

#else
    char *dbus_name = 0;
    bool casapy_start = false;
    preprocess_args( argc, argv, numargs, args, dbus_name, casapy_start );

    if ( casapy_start ) {
	signal(SIGINT,SIG_IGN);
    }

    QApplication qapp(numargs, args); 

    casa::QtPlotServer plot_server(dbus_name);
#endif

    return qapp.exec();

}

// processes argv into args stripping out the the viewer setup flags... numargs has the number
// of args, and the last arg (not included in numargs count) is null (for execvp)
#if defined(CASATOOLS)
static void preprocess_args( int argc, const char *argv[], int &numargs, char **&args,
                             char *&server_string, char *&event_sink, char *&logfile ) {
#else
static void preprocess_args( int argc, const char *argv[], int &numargs, char **&args,
			     char *&dbus_name, bool &casapy_start ) {
    casapy_start = false;
    dbus_name = 0;
#endif

    // pre-process the command line arguments (for now), it looks like
    // the current scheme is to hard-wire a limited set of command line
    // parameters... should be reworked in the future.
    args = (char**) malloc(sizeof(char*)*(argc+2));
    numargs = 0;

    args[numargs++] = strdup(argv[0]);
    for ( int x = 1; x < argc; ++x ) {
#if defined(CASATOOLS)
        if ( ! strncmp(argv[x], "--server=", 9) ) {
            char *sink = strdup( &argv[x][9] );
            if ( strlen(sink) <= 0 ) {
                free( sink );
            } else {
                server_string = sink;
            }
        } else if ( ! strncmp(argv[x], "--event-uri=", 12) ) {
            char *uri = strdup( &argv[x][12] );
            if ( strlen(uri) <= 0 ) {
                free( uri );
            } else {
                event_sink = uri;
            }
        } else if ( ! strncmp(argv[x], "--logfile=", 10) ) {
            char *lf = strdup( &argv[x][10] );
            if ( strlen(lf) <= 0 ) {
                free( lf );
            } else {
                logfile = lf;
            }
#else
	if ( ! strncmp(argv[x],"--dbusname",10) ) {
	    if ( argv[x][10] == '=' ) {
		char *name = strdup( &argv[x][11] );
		if ( strlen(name) <= 0 ) {
		    free( name );
		} else {
		    dbus_name = name;
		}
	    } else if ( x + 1 < argc ) {
		dbus_name = strdup(argv[++x]);
	    }
	} else if ( ! strcmp(argv[x],"--casapy") ) {
	    casapy_start = true;
#endif
	} else {
	    args[numargs++] = strdup(argv[x]);
	}
    }

    args[numargs] = 0;
}
