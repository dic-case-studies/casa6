/***
 * Framework independent implementation file for imager...
 *
 * Implement the imager component here.
 * 
 * // TODO: WRITE YOUR DESCRIPTION HERE! 
 
 * @author Wes Young
 * @version 
 ***/

#include <iostream>
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Record.h>
#include <casa/Utilities/Assert.h>
#include <ms/MeasurementSets.h>
#include <ms/MeasurementSets/MSHistoryHandler.h>
#include <casa/Logging/LogIO.h>

#include <synthesis/ImagerObjects/grpcInteractiveClean.h>

#include <iterbotsink_cmpt.h>

using namespace std;
using namespace casacore;
using namespace casa;

     
using namespace casacore;
namespace casac {

    iterbotsink::iterbotsink( ): state(grpcInteractiveClean::getManager( )) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::iterbotsink( )\n" );
    }

    iterbotsink::~iterbotsink( ) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::~iterbotsink( )\n" );
    }

    casac::record* iterbotsink::setupiteration(const casac::record& iterpars) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::setupiteration( )\n" );
        state.setIterationDetails( iterpars );
        return getiterationdetails();
    }


    casac::record* iterbotsink::getiterationdetails( ) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::getiterationdetails( )\n" );
        return fromRecord(state.getDetailsRecord( ));
    }

    casac::record* iterbotsink::pauseforinteraction( ) {
        
        fprintf( stderr, ">>>>>--------->> iterbotsink::pauseforinteraction( )\n" );
        return fromRecord( state.pauseForUserInteraction() );
    }

    casac::record* iterbotsink::getiterationsummary( ) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::getiterationsummary( )\n" );
        return fromRecord(state.getSummaryRecord( ));
    }

    int iterbotsink::cleanComplete(const bool lastcyclecheck) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::cleanComplete(const bool lastcyclecheck)\n" );
        return state.cleanComplete(lastcyclecheck);
    }

    bool iterbotsink::endmajorcycle( ) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::endmajorcycle( )\n" );
        state.incrementMajorCycleCount( );
        state.addSummaryMajor( );
        return false;
    }

    bool iterbotsink::resetminorcycleinfo( ) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::resetminorcycleinfo( )\n" );
        state.resetMinorCycleInitInfo( );
        return false;
    }

    casac::record* iterbotsink::getminorcyclecontrols( ) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::getminorcyclecontrols( )\n" );
        return fromRecord(state.getMinorCycleControls( ));
    }  

    bool iterbotsink::mergeinitrecord(const casac::record& initrecord) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::mergeinitrecord(const casac::record& initrecord)\n" );
        // ****memory leak here???****
        casacore::Record recpars = *toRecord( initrecord );
        state.mergeCycleInitializationRecord(recpars);
        return false;
    }

    bool iterbotsink::mergeexecrecord(const casac::record& execrecord) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::mergeexecrecord(const casac::record& execrecord)\n" );
        // ****memory leak here???****
        casacore::Record recpars = *toRecord(execrecord);
        state.mergeCycleExecutionRecord(recpars);
        return false;
    }

    bool iterbotsink::changestopflag(const bool stopflag) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::changestopflag(const bool stopflag)\n" );
        state.changeStopFlag(stopflag);
        return true;
    }

    bool iterbotsink::done( ) {
        fprintf( stderr, ">>>>>--------->> iterbotsink::done( )\n" );
        state.closePanel( );
        return false;
    }

} // casac namespace
