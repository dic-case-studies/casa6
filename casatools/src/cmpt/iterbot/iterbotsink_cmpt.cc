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

#include <synthesis/ImagerObjects/SIIterBot.h>

#include <iterbotsink_cmpt.h>

using namespace std;
using namespace casacore;
using namespace casa;

     
using namespace casacore;
namespace casac {

    iterbotsink::iterbotsink( ): state(nullptr) {
        cb.reset(new SIIterBot_callback( ));
        state = new SIIterBot_state( cb );
    }

    iterbotsink::~iterbotsink( ) {
        if ( state ) {
            delete state;
        }
    }

    casac::record* iterbotsink::setupiteration(const casac::record& iterpars) {
        const std::unique_ptr<const casacore::Record> recpars(toRecord(iterpars));
        state->setControlsFromRecord(*recpars);
        return getiterationdetails();
    }


    casac::record* iterbotsink::getiterationdetails( ) {
        return fromRecord(state->getDetailsRecord( ));
    }

    casac::record* iterbotsink::pauseforinteraction( ) {
        return getiterationdetails( );
    }

    casac::record* iterbotsink::getiterationsummary( ) {
        return fromRecord(state->getSummaryRecord( ));
    }

    int iterbotsink::cleanComplete(const bool lastcyclecheck) {
        return state->cleanComplete(lastcyclecheck);
    }

    bool iterbotsink::endmajorcycle( ) {
        state->incrementMajorCycleCount( );
        state->addSummaryMajor( );
        return false;
    }

    bool iterbotsink::resetminorcycleinfo( ) {
        state->resetMinorCycleInitInfo( );
        return false;
    }

    casac::record* iterbotsink::getminorcyclecontrols( ) {
        return fromRecord(state->getMinorCycleControls( ));
    }  

    bool iterbotsink::mergeinitrecord(const casac::record& initrecord) {
        const std::unique_ptr<const casacore::Record> recpars(toRecord(initrecord));
        state->mergeCycleInitializationRecord(*recpars);
        return false;
    }

    bool iterbotsink::mergeexecrecord(const casac::record& execrecord) {
        const std::unique_ptr<const casacore::Record> recpars(toRecord(execrecord));
        state->mergeCycleExecutionRecord(*recpars);
        return false;
    }

    bool iterbotsink::changestopflag(const bool stopflag) {
        state->changeStopFlag(stopflag);
        return true;
    }

    bool iterbotsink::done( ) {
	  delete state;
      state = 0;
      return false;
    }

} // casac namespace
