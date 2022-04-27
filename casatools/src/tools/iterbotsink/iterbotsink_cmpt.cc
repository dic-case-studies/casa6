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
    }

    iterbotsink::~iterbotsink( ) {
    }

    casac::record* iterbotsink::setupiteration(const casac::record& iterpars) {
        state.setIterationDetails( iterpars );
        return getiterationdetails();
    }


    casac::record* iterbotsink::getiterationdetails( ) {
        return fromRecord(state.getDetailsRecord( ));
    }

    casac::record* iterbotsink::pauseforinteraction( ) {
        static const auto debug = getenv("GRPC_DEBUG");
        auto ret = state.pauseForUserInteraction( );
        if ( debug ) {
            std::cerr << "-------------------------------------------" << std::endl;
            std::cerr << "--- pauseforinteraction result:         ---" << std::endl;
            std::cerr << "-------------------------------------------" << std::endl;
            std::cerr << ret << std::endl;
            std::cerr << "-------------------------------------------" << std::endl;
        }
        return fromRecord( ret );
    }

    casac::record* iterbotsink::getiterationsummary( ) {
        return fromRecord(state.getSummaryRecord( ));
    }

    long iterbotsink::cleanComplete(const bool lastcyclecheck, const bool reachedMajorLimit) {
        return state.cleanComplete(lastcyclecheck, reachedMajorLimit);
    }

    bool iterbotsink::endmajorcycle( ) {
        state.incrementMajorCycleCount( );
        state.addSummaryMajor( );
        return false;
    }

    bool iterbotsink::resetminorcycleinfo( ) {
        state.resetMinorCycleInitInfo( );
        return false;
    }

    casac::record* iterbotsink::getminorcyclecontrols( ) {
        return fromRecord(state.getMinorCycleControls( ));
    }  

    bool iterbotsink::mergeinitrecord(const casac::record& initrecord) {
        const std::unique_ptr<casacore::Record> recpars(toRecord(initrecord));
        state.mergeCycleInitializationRecord( *recpars );
        return false;
    }

    bool iterbotsink::mergeexecrecord(const casac::record& execrecord, const long immod) {
        const std::unique_ptr<casacore::Record> recpars(toRecord(execrecord));
        state.mergeCycleExecutionRecord( *recpars, immod );
        return false;
    }

    bool iterbotsink::changestopflag(const bool stopflag) {
        state.changeStopFlag(stopflag);
        return true;
    }

    bool iterbotsink::done( ) {
        state.closePanel( );
        return false;
    }

} // casac namespace
