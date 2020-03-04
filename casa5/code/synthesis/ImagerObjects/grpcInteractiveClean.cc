#include <synthesis/ImagerObjects/grpcInteractiveClean.h>
#include <casatools/Config/State.h>
#include <casacore/casa/Logging/LogIO.h>
#include <images/Images/PagedImage.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <array>
#include <regex>
#include <string>

#include <algorithm> 
#include <cctype>
#include <locale>

#include <grpc++/grpc++.h>
#include "shutdown.grpc.pb.h"
#include "img.grpc.pb.h"
#include "ping.grpc.pb.h"

#include <stdcasa/StdCasa/CasacSupport.h>

using namespace casacore;

// https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
// C++ is so ridiculous... trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

namespace casa { //# NAMESPACE CASA - BEGIN

    grpcInteractiveCleanManager &grpcInteractiveClean::getManager( ) {
        static grpcInteractiveCleanManager mgr;
        return mgr;
    }

	void grpcInteractiveCleanManager::pushDetails() {
#ifdef INTERACTIVE_ITERATION
        /*FIXME    detailUpdate(fromRecord(getDetailsRecord())); */
#endif
	}

    grpcInteractiveCleanState::grpcInteractiveCleanState( ) : SummaryMinor(casacore::IPosition(2,6,0)),
                                                              SummaryMajor(casacore::IPosition(1,0)) {
		LogIO os( LogOrigin("grpcInteractiveCleanState",__FUNCTION__,WHERE) );
        reset( );
    }

    void grpcInteractiveCleanState::reset( ) {
        Niter = 0;
        CycleNiter = 0;
        InteractiveNiter = 0;
        Threshold = 0;
        CycleThreshold = 0;
        InteractiveThreshold = 0.0;
        IsCycleThresholdAuto = true;
        IsThresholdAuto = false;
        CycleFactor = 1.0;
        LoopGain = 0.1;
        StopFlag = false;
        PauseFlag = false;
        InteractiveMode = false;
        UpdatedModelFlag = false;
        IterDone = 0;
        StopCode = 0;
        Nsigma = 0.0;
        MaxPsfSidelobe = 0.0;
        MinPsfFraction = 0.05;
        MaxPsfFraction = 0.8;
        PeakResidual = 0.0;
        MinorCyclePeakResidual = 0.0;
        PrevPeakResidual = -1.0;
        NsigmaThreshold = 0.0;
        PrevMajorCycleCount = 0;
        PeakResidualNoMask = 0.0;
        PrevPeakResidualNoMask = -1.0;
        MinPeakResidualNoMask = 1e+9;
        MinPeakResidual = 1e+9;
        MaskSum = -1.0;
        MadRMS = 0.0;
        NSummaryFields = 6;
    }

    void grpcInteractiveCleanManager::setControls( int niter, int ncycle, float threshold ) {
        LogIO os( LogOrigin("grpcInteractiveCleanManager", __FUNCTION__, WHERE) );
        static const auto debug = getenv("GRPC_DEBUG");
        if ( debug ) std::cerr << "setting clean controls:";
        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {

			   state.Niter = niter;
			   if ( debug ) std::cerr << " niter=" << state.Niter;
			   state.CycleNiter = ncycle;
			   if ( debug ) std::cerr << " cycleniter=" << state.CycleNiter;
			   state.Threshold = threshold;
			   if ( debug ) std::cerr << " threshold=" << state.Threshold;
                           return dummy;

                       } ) );

        if ( debug ) {
            std::cerr << " (process " << getpid( ) << ", thread " <<
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

    }
	void grpcInteractiveCleanManager::setControlsFromRecord(const casac::record &iterpars) {
		LogIO os( LogOrigin("grpcInteractiveCleanManager", __FUNCTION__, WHERE) );
        static const auto debug = getenv("GRPC_DEBUG");

        if ( debug ) std::cerr << "initializing clean controls:";

        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {

                           /* Note it is important that niter get set first as we catch
                              negative values in the cycleniter, and set it equal to niter */
                           auto niter = iterpars.find("niter");
                           if ( niter != iterpars.end( ) ) {
                               state.Niter = niter->second.toInt( );
                               if ( debug ) std::cerr << " niter=" << state.Niter;
                           }

                           auto cycleniter = iterpars.find("cycleniter");
                           if ( cycleniter != iterpars.end( ) ) {
                               int val = cycleniter->second.toInt( );
                               if ( val <= 0 )
                                   state.CycleNiter = state.Niter;
                               else
                                   state.CycleNiter = val;
                               if ( debug ) std::cerr << " cycleniter=" << state.CycleNiter;
                           }

                           auto interactiveniter = iterpars.find("interactiveniter");
                           if ( interactiveniter != iterpars.end( ) ) {
                               state.InteractiveNiter = interactiveniter->second.toInt( );
                               if ( debug ) std::cerr << " interactiveniter=" << state.InteractiveNiter;
                           }

                           auto threshold = iterpars.find("threshold");
                           if ( threshold != iterpars.end( ) ) {
                               auto val = threshold->second.toDouble( );
                               if ( val == -1.0 ) {
                                   state.Threshold = 0.0;
                                   state.IsThresholdAuto = true;
                               } else {
                                   state.Threshold = (float) val; 
                                   state.IsThresholdAuto = false;
                               }
                               if ( debug ) {
                                   std::cerr << " threshold=" << state.Threshold;
                                   std::cerr << " isthresholdauto=" <<
                                       (state.IsThresholdAuto ? "true" : "false");
                               }
                           }
                               
                           auto cyclethreshold = iterpars.find("cyclethreshold");
                           if ( cyclethreshold != iterpars.end( ) ) {
                               state.CycleThreshold = (float) cyclethreshold->second.toDouble( );
                               state.IsCycleThresholdAuto = false;
                               if ( debug ) {
                                   std::cerr << " cyclethreshold=" << state.CycleThreshold;
                                   std::cerr << " iscyclethresholdauto=" <<
                                       (state.IsCycleThresholdAuto ? "true" : "false");
				   fflush(stderr);
                               }
                           }

                           auto interactivethreshold = iterpars.find("interactivethreshold");
                           if ( interactivethreshold != iterpars.end( ) ) {
                               state.InteractiveThreshold = (float) interactivethreshold->second.toDouble( );
                               if ( debug ) std::cerr << " interactivethreshold=" << state.InteractiveThreshold;
                           }

                           auto loopgain = iterpars.find("loopgain");
                           if ( loopgain != iterpars.end( ) ) {
                               state.LoopGain = (float) loopgain->second.toDouble( );
                               if ( debug ) std::cerr << " loopgain=" << state.LoopGain;
                           }
                           auto cyclefactor = iterpars.find("cyclefactor");
                           if ( cyclefactor != iterpars.end( ) ) {
                               state.CycleFactor = (float) cyclefactor->second.toDouble( );
                               if ( debug ) std::cerr << " cyclefactor=" << state.CycleFactor;
                           }

                           auto interactivemode = iterpars.find("interactive");
                           if ( interactivemode != iterpars.end( ) ) {
                               state.InteractiveMode = interactivemode->second.toBool( );
                               if ( debug ) std::cerr << " interactive=" <<
                                                (state.InteractiveMode ? "true" : "false");
                           }

                           auto minpsffraction = iterpars.find("minpsffraction");
                           if ( minpsffraction != iterpars.end( ) ) {
                               state.MinPsfFraction = (float) minpsffraction->second.toDouble( );
                               if ( debug ) std::cerr << " minpsffraction=" << state.MinPsfFraction;
                           }

                           auto maxpsffraction = iterpars.find("maxpsffraction");
                           if ( maxpsffraction != iterpars.end( ) ) {
                               state.MaxPsfFraction = (float) maxpsffraction->second.toDouble( );
                               if ( debug ) std::cerr << " maxpsffraction=" << state.MaxPsfFraction;
                           }

                           auto nsigma = iterpars.find("nsigma");
                           if ( nsigma != iterpars.end( ) ) {
                               state.Nsigma = (float) nsigma->second.toDouble( );
                               if ( debug ) std::cerr << " nsigma=" << state.Nsigma;
                           }
                           return dummy;

                       } ) );

        if ( debug ) {
            std::cerr << " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
    }
/*
        Float SIIterBot_state::readThreshold( Record recordIn, String id )  {
		LogIO os( LogOrigin("SIIterBot_state",__FUNCTION__,WHERE) );
		std::lock_guard<std::recursive_mutex> guard(recordMutex);    
		// Threshold can be a variant, either Float or String(with units).
		Float fthresh=0.0;
		// If a number, treat it as a number in units of Jy.
		if( recordIn.dataType(id) == TpFloat || 
		    recordIn.dataType(id) == TpDouble || 
		    recordIn.dataType(id) == TpInt )
		  { fthresh = recordIn.asFloat( RecordFieldId(id)); }
		// If a string, try to convert to a Quantity
		else if( recordIn.dataType(id) == TpString )
		  {
		    Quantity thresh; 
		    // If it cannot be converted to a Quantity.... complain, and use zero.
		    if( ! casacore::Quantity::read( thresh, recordIn.asString( RecordFieldId(id) ) ) )
		      {os << LogIO::WARN << "Cannot parse threshold value. Setting to zero." << LogIO::POST;  
			fthresh=0.0;}
		    // If converted to Quantity, get value in Jy. 
		    // ( Note : This does not check for wrong units, e.g. if the user says '100m' ! )
		    else { fthresh = thresh.getValue(Unit("Jy")); }
		  }
		// If neither valid datatype, print a warning and use zero.
		else {os << LogIO::WARN << id << " is neither a number nor a string Quantity. Setting to zero." << LogIO::POST;
		  fthresh=0.0; }

		return fthresh;
        }
*/
	void grpcInteractiveCleanManager::setIterationDetails(const casac::record &iterpars) {
		LogIO os( LogOrigin("grpcInteractiveCleanManager",__FUNCTION__,WHERE) );
        static const auto debug = getenv("GRPC_DEBUG");

        // Setup interactive masking : list of image names.
        if ( clean_images.size( ) == 0 ) {
            try {
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                ///// START : code to get a list of image names for interactive masking

                auto allimages = iterpars.find("allimages");
                if ( allimages != iterpars.end( ) ) {
                    auto rec = allimages->second.getRecord( );
                    for ( auto it = rec.begin( ); it != rec.end( ); ++it ) {
                        auto oneimg = it->second.getRecord( );
                        auto img_name = oneimg.find("imagename");
                        auto img_multiterm = oneimg.find("multiterm");
                        if ( img_name != oneimg.end( ) && img_multiterm != oneimg.end( ) ) {
                            clean_images.push_back( std::make_tuple( img_name->second.getString(),
                                                                     img_multiterm->second.toBool( ), false) );
                        }
                    }
                } else {
                    throw( AipsError("Need image names and nterms in iteration parameter list") );
                }
                if ( clean_images.size( ) <= 0 ) {
                    throw( AipsError("Need image names for iteration") );
                }

                if ( debug ) {
                    std::cerr << "clean images specified: ";
                    for ( auto it = clean_images.begin( ); it != clean_images.end( ); ++it ) {
                        if ( it != clean_images.begin( ) ) std::cerr << ", ";
                        std::cerr << std::get<0>(*it) << " [" << (std::get<1>(*it) ? "true" : "false") << "]";
                    }
                    std::cerr << " (process " << getpid( ) << ", thread " << 
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }

                ///// END : code to get a list of image names for interactive masking
                ////////////////////////////////////////////////////////////////////////////////////////////////////////

                setControlsFromRecord( iterpars );

            } catch( AipsError &x ) {
                throw( AipsError("Error in updating iteration parameters : " + x.getMesg()) );
            }
        }
	}

	void grpcInteractiveCleanManager::updateCycleThreshold( grpcInteractiveCleanState &state ) {

		Float psffraction = state.MaxPsfSidelobe * state.CycleFactor;

		psffraction = max(psffraction, state.MinPsfFraction);
		psffraction = min(psffraction, state.MaxPsfFraction);

		state.CycleThreshold = state.PeakResidual * psffraction;
		pushDetails();
	}

	void grpcInteractiveCleanManager::addSummaryMajor( ) {
        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {
                           IPosition shp = state.SummaryMajor.shape();
                           if( shp.nelements() != 1 )
                               throw(AipsError("Internal error in shape of major-cycle summary record"));

                           state.SummaryMajor.resize( IPosition( 1, shp[0]+1 ) , true );
                           state.SummaryMajor( IPosition(1, shp[0] ) ) = state.IterDone;
                           return dummy; } ) );
	}

    casacore::Record grpcInteractiveCleanManager::getDetailsRecord( bool includeSummary ) {
		LogIO os( LogOrigin("grpcInteractiveCleanManager",__FUNCTION__,WHERE) );

		Record returnRecord;

        Record result = access( returnRecord,
                       std::function< casacore::Record ( casacore::Record, grpcInteractiveCleanState & )>(
                       [&]( casacore::Record rec, grpcInteractiveCleanState &state )->casacore::Record {
                           //*** Control Variables **************************************************
                           rec.define( RecordFieldId("niter"), state.Niter );
                           rec.define( RecordFieldId("cycleniter"), state.CycleNiter );
                           rec.define( RecordFieldId("interactiveniter"), state.InteractiveNiter );

                           rec.define( RecordFieldId("threshold"),  state.Threshold );
                           rec.define( RecordFieldId("nsigma"),  state.Nsigma );

                           if( state.IsCycleThresholdAuto == true ) updateCycleThreshold(state);
                           state.IsCycleThresholdAuto = true;        // Reset this, for the next round

                           rec.define( RecordFieldId("cyclethreshold"), state.CycleThreshold );
                           rec.define( RecordFieldId("interactivethreshold"), state.InteractiveThreshold );

                           rec.define( RecordFieldId("loopgain"), state.LoopGain );
                           rec.define( RecordFieldId("cyclefactor"), state.CycleFactor );

                           //*** Status Reporting Variables *****************************************
                           rec.define( RecordFieldId("iterdone"),  state.IterDone );
                           rec.define( RecordFieldId("cycleiterdone"), state.MaxCycleIterDone );
                           rec.define( RecordFieldId("interactiveiterdone"),
                                       state.InteractiveIterDone + state.MaxCycleIterDone);

                           rec.define( RecordFieldId("nmajordone"), state.MajorDone );
                           rec.define( RecordFieldId("maxpsfsidelobe"), state.MaxPsfSidelobe );
                           rec.define( RecordFieldId("maxpsffraction"), state.MaxPsfFraction );
                           rec.define( RecordFieldId("minpsffraction"), state.MinPsfFraction );
                           rec.define( RecordFieldId("interactivemode"), state.InteractiveMode );

                           rec.define( RecordFieldId("stopcode"), state.StopCode );

                           //*** report clean's state ***********************************************
                           rec.define( RecordFieldId("cleanstate"),
                                       state.StopFlag ? "stopped" : state.PauseFlag ? "paused" : "running" );

                           if ( includeSummary ) {
                               rec.define( RecordFieldId("summaryminor"), state.SummaryMinor );
                               rec.define( RecordFieldId("summarymajor"), state.SummaryMajor );
                           }

                           return rec; }) );

        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        std::cout << result << std::endl;
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        return result;


/*      return access( returnRecord,
                       std::function< casacore::Record ( casacore::Record, grpcInteractiveCleanState & )>(
                       [&]( casacore::Record rec, grpcInteractiveCleanState &state )->casacore::Record {
                           //*** Control Variables **************************************************
                           rec.define( RecordFieldId("niter"), state.Niter );
                           rec.define( RecordFieldId("cycleniter"), state.CycleNiter );
                           rec.define( RecordFieldId("interactiveniter"), state.InteractiveNiter );

                           rec.define( RecordFieldId("threshold"),  state.Threshold );
                           rec.define( RecordFieldId("nsigma"),  state.Nsigma );
                           if( state.IsCycleThresholdAuto == true ) updateCycleThreshold(state);
                           state.IsCycleThresholdAuto = true;        // Reset this, for the next round

                           rec.define( RecordFieldId("cyclethreshold"), state.CycleThreshold );
                           rec.define( RecordFieldId("interactivethreshold"), state.InteractiveThreshold );

                           rec.define( RecordFieldId("loopgain"), state.LoopGain );
                           rec.define( RecordFieldId("cyclefactor"), state.CycleFactor );

                           //*** Status Reporting Variables *****************************************
                           rec.define( RecordFieldId("iterdone"),  state.IterDone );
                           rec.define( RecordFieldId("cycleiterdone"), state.MaxCycleIterDone );
                           rec.define( RecordFieldId("interactiveiterdone"),
                                       state.InteractiveIterDone + state.MaxCycleIterDone);

                           rec.define( RecordFieldId("nmajordone"), state.MajorDone );
                           rec.define( RecordFieldId("maxpsfsidelobe"), state.MaxPsfSidelobe );
                           rec.define( RecordFieldId("maxpsffraction"), state.MaxPsfFraction );
                           rec.define( RecordFieldId("minpsffraction"), state.MinPsfFraction );
                           rec.define( RecordFieldId("interactivemode"), state.InteractiveMode );

                           rec.define( RecordFieldId("stopcode"), state.StopCode );

                           //*** report clean's state ***********************************************
                           rec.define( RecordFieldId("cleanstate"),
                                       state.StopFlag ? "stopped" : state.PauseFlag ? "paused" : "running" );

                           if ( includeSummary ) {
                               rec.define( RecordFieldId("summaryminor"), state.SummaryMinor );
                               rec.define( RecordFieldId("summarymajor"), state.SummaryMajor );
                           }

                           return rec; }) );
*/
    }


	Record grpcInteractiveCleanManager::getMinorCycleControls( ){
		LogIO os( LogOrigin("grpcInteractiveCleanManager",__FUNCTION__,WHERE) );

		/* This returns a record suitable for initializing the minor cycle controls. */
		Record returnRecord;

        return access( returnRecord,
                       std::function< casacore::Record ( casacore::Record, grpcInteractiveCleanState & )>(
                       [&]( casacore::Record rec, grpcInteractiveCleanState &state )->casacore::Record {

                           /* If autocalc, compute cyclethresh from peak res, cyclefactor and psf sidelobe
                              Otherwise, the user has explicitly set it (interactively) for this minor cycle */
                           if( state.IsCycleThresholdAuto == true ) { updateCycleThreshold(state); }
                           state.IsCycleThresholdAuto = true; /* Reset this, for the next round */

                           /* The minor cycle will stop based on the cycle parameters. */
                           int maxCycleIterations = state.CycleNiter;
                           float cycleThreshold     = state.CycleThreshold;
                           maxCycleIterations = min(maxCycleIterations, state.Niter - state.IterDone);
                           cycleThreshold = max(cycleThreshold, state.Threshold);
                           bool thresholdReached = (cycleThreshold==state.Threshold)? True : False;

                           rec.define( RecordFieldId("cycleniter"),  maxCycleIterations);
                           rec.define( RecordFieldId("cyclethreshold"), cycleThreshold);
                           rec.define( RecordFieldId("loopgain"), state.LoopGain);
                           rec.define( RecordFieldId("thresholdreached"), thresholdReached);
                           rec.define( RecordFieldId("nsigma"), state.Nsigma);

                           return rec; }) );
    }

	int grpcInteractiveCleanManager::cleanComplete( bool lastcyclecheck ){
        LogIO os( LogOrigin("grpcInteractiveCleanManager",__FUNCTION__,WHERE) );

		int stopCode=0;

        return access( stopCode,
                       std::function< int ( int, grpcInteractiveCleanState & )>(
                       [&]( int stop_code, grpcInteractiveCleanState &state ) -> int {

                           float usePeakRes;

                           if( lastcyclecheck==True ) { usePeakRes = state.MinorCyclePeakResidual; }
                           else { usePeakRes = state.PeakResidual; }

                           // for debugging, remove it later
                           os<<LogIO::DEBUG1<<"cleanComplete-- CycleThreshold without Threshold limit="<<state.CycleThreshold<<LogIO::POST;

                           if( state.PeakResidual > 0 && state.PrevPeakResidual>0 &&
                               fabs( state.PeakResidual - state.PrevPeakResidual)/fabs(state.PrevPeakResidual) > 2.0 ) {
                               os << "[WARN] Peak residual (within the mask) increased from " << state.PrevPeakResidual << " to " << state.PeakResidual << LogIO::POST;
                           }
                           // for debugging, remove it later
                           os <<LogIO::DEBUG1<<"Threshold="<<state.Threshold<<" itsNsigmaThreshold===="<<state.NsigmaThreshold<<LogIO::POST;
                           os <<LogIO::DEBUG1<<"usePeakRes="<<usePeakRes<<" itsPeakResidual="<<state.PeakResidual<<" itsPrevPeakRes="<<state.PrevPeakResidual<<LogIO::POST;
                           os <<LogIO::DEBUG1<<"itsIterDone="<<state.IterDone<<" itsNiter="<<state.Niter<<LogIO::POST;

                           /// This may interfere with some other criterion... check.
                           float tol = 0.01; // threshold test torelance (CAS-11278)
                           if ( state.MajorDone==0 && state.IterDone==0 ) {
                               if (state.MaskSum==0.0) { stopCode=7; } // if zero mask is detected it should exit right away
                               else { stopCode=0; }
                           } else if ( state.IterDone >= state.Niter ||
                                       state.PeakResidual <= state.Threshold ||
                                       state.PeakResidual <= state.NsigmaThreshold ||
                                       abs(state.PeakResidual - state.Threshold)/state.Threshold < tol ||
                                       abs(state.PeakResidual - state.NsigmaThreshold)/state.NsigmaThreshold < tol ||
                                       state.StopFlag ) {
                               //		    os << "Reached global stopping criteria : ";

                               if ( state.IterDone >= state.Niter ) { stopCode=1; }
                               //os << "Numer of iterations. "; // (" << state.IterDone << ") >= limit (" << state.Niter << ")" ;
                               if( usePeakRes <= state.Threshold || (usePeakRes-state.Threshold)/state.Threshold < tol) {stopCode=2; }
                               else if ( usePeakRes <= state.NsigmaThreshold || (state.PeakResidual - state.NsigmaThreshold)/state.NsigmaThreshold < tol ) {
                                   if (state.NsigmaThreshold!=0.0) { stopCode=8; } // for nsigma=0.0 this mode is turned off
                               }

                               //os << "Peak residual (" << state.PeakResidual << ") <= threshold(" << state.Threshold << ")";
                               if( state.StopFlag ) {stopCode=3;}
                               //os << "Forced stop. ";
                               //		    os << LogIO::POST;

                               //return true;

                           } else { // not converged yet... but....if nothing has changed in this round... also stop

                               if (state.MaskSum==0.0) {
                                   //cout << "(7) Mask is all zero.Stopping" << endl;
                                   stopCode = 7;
                               }
                               // Nothing has changed across the last set of minor cycle iterations and major cycle.
                               else if( state.IterDone>0 && (state.MajorDone>state.PrevMajorCycleCount) &&
                                        fabs(state.PrevPeakResidual - state.PeakResidual)<1e-10)
                                   {stopCode = 4;}

                               // another non-convergent condition: diverging (relative increase is more than 3 times across one major cycle)
                               else if ( state.IterDone > 0 &&
                                         fabs(state.PeakResidualNoMask-state.PrevPeakResidualNoMask)/fabs(state.PrevPeakResidualNoMask)  > 3.0) {
                                   //cout << "(5) Peak res (no mask) : " << state.PeakResidualNoMask
                                   //     << "  Dev from prev peak res " << state.PrevPeakResidualNoMask << endl;
                                   stopCode = 5;}

                               // divergence check, 3 times increase from the minimum peak residual so far (across all previous major cycles).
                               else if ( state.IterDone > 0 &&
			      (fabs(state.PeakResidualNoMask)-state.MinPeakResidualNoMask)/state.MinPeakResidualNoMask  > 3.0 )
                      {
			//cout << "(6) Peak res (no mask): " << state.PeakResidualNoMask
			//    <<  "    Dev from min peak res " << state.MinPeakResidualNoMask << endl;
			stopCode = 6;
		      }

		  }

		/*
		if( lastcyclecheck==False)
		  {
		    cout << "*****" << endl;
		    cout << "Peak residual : " << state.PeakResidual << "  No Mask : " << state.PeakResidualNoMask << endl;
		    cout << "Prev Peak residual : " << state.PrevPeakResidual << "  No Mask : " << state.PrevPeakResidualNoMask << endl;
		    cout << "Min Peak residual : " << state.MinPeakResidual << "  No Mask : " << state.MinPeakResidualNoMask << endl;
		  }
		*/

		//		os << "Peak residual : " << state.PeakResidual << " and " << state.IterDone << " iterations."<< LogIO::POST;
		//cout << "cleancomp : stopcode : " << stopCode << endl;

		//cout << "peak res : " << state.PeakResidual << "   state.minPR : " << state.MinPeakResidual << endl;

		if( lastcyclecheck==False)
		  {
		    if( fabs(state.PeakResidual) < state.MinPeakResidual )
		      {state.MinPeakResidual = fabs(state.PeakResidual);}

		    state.PrevPeakResidual = state.PeakResidual;


		    if( fabs(state.PeakResidualNoMask) < state.MinPeakResidualNoMask )
		      {state.MinPeakResidualNoMask = fabs(state.PeakResidualNoMask);}

		    state.PrevPeakResidualNoMask = state.PeakResidualNoMask;

		    state.PrevMajorCycleCount = state.MajorDone;

		  }

		state.StopCode=stopCode;
		return stopCode; } ) );
	}

	void grpcInteractiveCleanManager::resetMinorCycleInitInfo( grpcInteractiveCleanState &state ) {
		/* Get ready to do the minor cycle */
		state.PeakResidual = 0;
		state.PeakResidualNoMask = 0;
		state.MaxPsfSidelobe = 0;
		state.MaxCycleIterDone = 0;
		state.MaskSum = -1.0;
	}

    void grpcInteractiveCleanManager::resetMinorCycleInitInfo( ) {
        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {
                           resetMinorCycleInitInfo(state);
                           return dummy; } ) );
    }

	void grpcInteractiveCleanManager::incrementMajorCycleCount( ) {

        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {
                           state.PrevMajorCycleCount = state.MajorDone;
                           state.MajorDone++;

                           /* Interactive iteractions update */
                           state.InteractiveIterDone += state.MaxCycleIterDone;

                           resetMinorCycleInitInfo(state);
                           return dummy; } ) );
	}

	void grpcInteractiveCleanManager::mergeCycleInitializationRecord( Record &initRecord ){
        LogIO os( LogOrigin("grpcInteractiveCleanManager",__FUNCTION__,WHERE) );

        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {

                           state.PeakResidual = max(state.PeakResidual, initRecord.asFloat(RecordFieldId("peakresidual")));
                           state.MaxPsfSidelobe =  max(state.MaxPsfSidelobe, initRecord.asFloat(RecordFieldId("maxpsfsidelobe")));

                           state.PeakResidualNoMask = max( state.PeakResidualNoMask, initRecord.asFloat(RecordFieldId("peakresidualnomask")));
                           state.MadRMS = max(state.MadRMS, initRecord.asFloat(RecordFieldId("madrms")));
                           state.NsigmaThreshold = initRecord.asFloat(RecordFieldId("nsigmathreshold"));

                           /*
                             It has been reset to -1.0.
                             If no masks have changed, it should remain at -1.0
                             If any mask has changed, the sum will come in, and should be added to this.
                           */
                           float thismasksum = initRecord.asFloat(RecordFieldId("masksum"));
                           if( thismasksum != -1.0 ) {
                               if ( state.MaskSum == -1.0 ) state.MaskSum = thismasksum;
                               else state.MaskSum += thismasksum;
                           }

                           if ( state.PrevPeakResidual == -1.0 ) state.PrevPeakResidual = state.PeakResidual;
                           if ( state.PrevPeakResidualNoMask == -1.0 ) state.PrevPeakResidualNoMask = state.PeakResidualNoMask;
                           if( state.IsCycleThresholdAuto == true ) updateCycleThreshold(state);

                           return dummy; } ) );
	}


	void grpcInteractiveCleanManager::mergeMinorCycleSummary( const Array<Double> &summary, grpcInteractiveCleanState &state ){

		IPosition cShp = state.SummaryMinor.shape();
		IPosition nShp = summary.shape();

		if( cShp.nelements() != 2 || cShp[0] != state.NSummaryFields ||
			nShp.nelements() != 2 || nShp[0] != state.NSummaryFields )
			throw(AipsError("Internal error in shape of global minor-cycle summary record"));

		state.SummaryMinor.resize( IPosition( 2, state.NSummaryFields, cShp[1]+nShp[1] ) ,true );

		for (unsigned int row = 0; row < nShp[1]; row++) {
			// iterations done
			state.SummaryMinor( IPosition(2,0,cShp[1]+row) ) = state.IterDone + summary(IPosition(2,0,row));
			// peak residual
			state.SummaryMinor( IPosition(2,1,cShp[1]+row) ) = summary(IPosition(2,1,row));
			// model flux
			state.SummaryMinor( IPosition(2,2,cShp[1]+row) ) = summary(IPosition(2,2,row));
			// cycle threshold
			state.SummaryMinor( IPosition(2,3,cShp[1]+row) ) = summary(IPosition(2,3,row));
			// mapper id
			state.SummaryMinor( IPosition(2,4,cShp[1]+row) ) = summary(IPosition(2,4,row));
			// chunk id (channel/stokes)
			state.SummaryMinor( IPosition(2,5,cShp[1]+row) ) = summary(IPosition(2,5,row));
		}
	}

	void grpcInteractiveCleanManager::mergeCycleExecutionRecord( Record& execRecord ){
        LogIO os( LogOrigin("grpcInteractiveCleanManager",__FUNCTION__,WHERE) );

        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {
                           mergeMinorCycleSummary( execRecord.asArrayDouble( RecordFieldId("summaryminor")), state );

                           state.IterDone += execRecord.asInt(RecordFieldId("iterdone"));

                           state.MaxCycleIterDone = max( state.MaxCycleIterDone, execRecord.asInt(RecordFieldId("maxcycleiterdone")) );

                           state.MinorCyclePeakResidual = max( state.PeakResidual, execRecord.asFloat(RecordFieldId("peakresidual")) );

                           state.UpdatedModelFlag |=execRecord.asBool( RecordFieldId("updatedmodelflag") );

                           os << "Completed " << state.IterDone << " iterations." << LogIO::POST;
                           //with peak residual "<< state.PeakResidual << LogIO::POST;
                           return dummy; } ) );
	}

	void grpcInteractiveCleanManager::changeStopFlag( bool stopEnabled ) {
        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {
                           state.StopFlag = stopEnabled;
                           return dummy;
                       } ) );
	}

    //====================================================================================================

    static bool isdir( const char *path ) {
        struct stat statbuf;
        int err = stat(path, &statbuf);
        if ( err == -1 ) return false;
        if ( S_ISDIR(statbuf.st_mode) ) return true;
        return false;
    }

    static std::string trim_trailing_slash( const char *str ) {
        char *temp = strdup(str);
        for ( int off = strlen(str) - 1; off >= 0; --off ) {
            if ( temp[off] == '/' ) temp[off] = '\0';
            else break;
        }
        std::string result = temp;
        free(temp);
        return result;
    }

    grpcInteractiveCleanGui::grpcInteractiveCleanGui( ) : viewer_pid(0), viewer_started(false) { }
    grpcInteractiveCleanGui::~grpcInteractiveCleanGui( ) {
        static const auto debug = getenv("GRPC_DEBUG");

        if ( ! viewer_started ) {
            if ( debug ) {
                std::cerr << "viewer shutdown required (" << viewer_uri << ")" << 
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
        } else {
            if ( debug ) {
                std::cerr << "sending shutdown message to viewer (" << viewer_uri << ")" << 
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }

            bool stopped = stop_viewer( );

            if ( debug ) {
                if ( stopped ) {
                    std::cerr << "viewer shutdown successful (" << viewer_uri << ")" << 
                        " (process " << getpid( ) << ", thread " << 
                        std::this_thread::get_id() << ")" << std::endl;
                } else {
                    std::cerr << "viewer shutdown failed (" << viewer_uri << ")" << 
                        " (process " << getpid( ) << ", thread " << 
                        std::this_thread::get_id() << ")" << std::endl;
                }
                fflush(stderr);
            }
        }
    }


    bool grpcInteractiveCleanGui::alive( ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if ( debug ) {
            std::cerr << "pinging viewer (" << viewer_uri << ")" << 
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        grpc::ClientContext context;
        ::google::protobuf::Empty resp;
        ::google::protobuf::Empty msg;
        auto ping = casatools::rpc::Ping::NewStub( grpc::CreateChannel( viewer_uri, grpc::InsecureChannelCredentials( ) ) );
        ::grpc::Status status = ping->now( &context, msg, &resp );
        bool ping_result = status.ok( );
        if ( debug ) {
            std::cerr << "ping result: " << (ping_result ? "OK" : "FAIL")<< 
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        if ( ping_result == false ) {
            int proc_status;
            waitpid( viewer_pid, &proc_status, WUNTRACED | WCONTINUED | WNOHANG );
            viewer_pid = 0;
            viewer_proxy.release( );
            viewer_started = false;
            if ( debug ) {
                std::cerr << "ping failed resetting state" << 
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
        }
        return ping_result;
    }

    bool grpcInteractiveCleanGui::launch( ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if ( viewer_started == false ) {
            // start the viewer process if it is not already running...
            if ( debug ) {
                std::cerr << "spawning viewer process" <<
                    " (process " << getpid( ) << ", thread " <<
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            return spawn_viewer( );
        } else {
	    if ( alive( ) ) {
	      if ( debug ) {
                std::cerr << "viewer process available" <<
                    " (process " << getpid( ) << ", thread " <<
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
	      }
	      return true;
	    } else {
	      if ( debug ) {
                std::cerr << "re-spawning viewer process" <<
                    " (process " << getpid( ) << ", thread " <<
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
	      }
	      return launch( );
	    }
        }
        return false;
    }

    void grpcInteractiveCleanGui::close_panel( int id ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if ( debug ) {
            std::cerr << "close_panel(" << id << ")" <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        if ( id != -1 && alive( ) ) {
            if ( debug ) {
                std::cerr << "close_panel(" << id << ") -- closing panel" <<
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            rpc::img::Id panel;
            grpc::ClientContext context;
            ::google::protobuf::Empty resp;
            panel.set_id(id);
            viewer_proxy->close( &context, panel, &resp );
        }
    }

    int grpcInteractiveCleanGui::open_panel( std::list<std::tuple<std::string,bool,bool>> images ) {
        static const auto debug = getenv("GRPC_DEBUG");
        if ( viewer_started == false ) {
            if ( launch( ) == false ) return -1;
        }
        if ( debug ) {
            std::cerr << "opening viewer panel" <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        grpc::ClientContext context;
        ::rpc::img::NewPanel np;
        rpc::img::Id resp;
        np.set_type("clean2");
        np.set_hidden(false);
        viewer_proxy->panel( &context, np, &resp );
        int result = resp.id( );

        if ( debug ) {
            std::cerr << "opened viewer panel " << result <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        // state for interactive masking in the new viewer panel
        clean_state.insert( std::pair<int,CleanState>(result, CleanState( )) );

        if ( debug ) {
            std::cerr << "created panel " << result <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        return result;
    }

    void grpcInteractiveCleanGui::unload( int id ) {
        grpc::ClientContext context;
        ::rpc::img::Id data;
        ::google::protobuf::Empty resp;
        data.set_id(id);
        viewer_proxy->unload( &context, data, &resp );
    }

    bool grpcInteractiveCleanGui::clone( const std::string &imageName, const std::string &newImageName ) {
        LogIO os(LogOrigin("grpcInteractiveCleanGui", __FUNCTION__, WHERE));

        try {
            PagedImage<Float> oldImage( imageName );
            PagedImage<Float> newImage( TiledShape( oldImage.shape(), oldImage.niceCursorShape()),
                                        oldImage.coordinates(), newImageName );
            newImage.set(0.0);
            newImage.table().flush(true, true);
        } catch (AipsError x) {
            os << LogIO::SEVERE << "Exception: " << x.getMesg() << LogIO::POST;
            return false;
        } 
        return true;
    }

    float grpcInteractiveCleanGui::maskSum(const std::string &maskname) {

        PagedImage<Float> mask( maskname );

        LatticeExprNode msum( sum( mask ) );
        float maskSum = msum.getFloat( );

        mask.unlock();
        mask.tempClose();

        return maskSum;
    }

    int grpcInteractiveCleanGui::interactivemask( int panel, const std::string &image, const std::string &mask,
                                                  int &niter, int &cycleniter, std::string &thresh,
                                                  std::string &cyclethresh, const bool forceReload ) {

        static const auto debug = getenv("GRPC_DEBUG");
        LogIO os( LogOrigin("grpcInteractiveCleanGui",__FUNCTION__,WHERE) );

        if ( debug ) {
            std::cerr << "starting interactivemask( " <<
                                 panel << ", " << image << ", " << mask << ", " <<
                                 niter << ", " << cycleniter << ", " << thresh << ", " <<
                                 cyclethresh << ", " << (forceReload ? "true" : "false") << ")" << 
                                 " (process " << getpid( ) << ", thread " << 
                                 std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        if ( viewer_started == false ) {
            // viewer should be started before calling interactivemask(...)
            os << LogIO::WARN << "Viewer GUI Not Available" << LogIO::POST;
            return 0;
        }

        auto state = clean_state.find(panel);
        if ( state == clean_state.end( ) ) {
            os << LogIO::WARN << "Invalid clean panel id used for interactive masking" << LogIO::POST;
            return 0;
        }

        if( Table::isReadable(mask) ) {
            if ( ! Table::isWritable(mask) ) {
                os << LogIO::WARN << "Mask image is not modifiable " << LogIO::POST;
                return 0;
            }
            // we should regrid here if image and mask do not match
        } else {
            clone(image, mask);
        }
    
        double startmask = maskSum(mask);

        if ( state->second.image_id == 0 || state->second.mask_id == 0 || forceReload ) {

            //Make sure image left after a "no more" is pressed is cleared
            if ( forceReload && state->second.image_id !=0 )
                state->second.prev_image_id = state->second.image_id;
            if ( forceReload && state->second.mask_id !=0 )
                state->second.prev_mask_id = state->second.mask_id;

            if ( state->second.prev_image_id ){
                if ( debug ) {
                    std::cerr << "preparing to unload prev_image_id " << state->second.prev_image_id << " (panel " << panel << ")" <<
                        " (process " << getpid( ) << ", thread " << 
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }
                unload( state->second.prev_image_id );
            }
            if ( state->second.prev_mask_id ) {
                if ( debug ) {
                    std::cerr << "preparing to unload prev_mask_id " << state->second.prev_mask_id << " (panel " << panel << ")" <<
                        " (process " << getpid( ) << ", thread " << 
                        std::this_thread::get_id() << ")" << std::endl;
                    fflush(stderr);
                }
                unload( state->second.prev_mask_id );
            }

            state->second.prev_image_id = 0;
            state->second.prev_mask_id = 0;

            {
                grpc::ClientContext context;
                ::rpc::img::NewData nd;
                rpc::img::Id resp;
                nd.mutable_panel( )->set_id(panel);
                nd.set_path(image);
                nd.set_type("raster");
                nd.set_scale(0);
                viewer_proxy->load( &context, nd, &resp );
                state->second.image_id = resp.id( );
            }
            {
                grpc::ClientContext context;
                ::rpc::img::NewData nd;
                rpc::img::Id resp;
                nd.mutable_panel( )->set_id(panel);
                nd.set_path(mask);
                nd.set_type("contour");
                nd.set_scale(0);
                viewer_proxy->load( &context, nd, &resp );
                state->second.mask_id = resp.id( );
            }

        } else {
            grpc::ClientContext context;
            ::rpc::img::Id id;
            ::google::protobuf::Empty resp;
            id.set_id(state->second.image_id);
            viewer_proxy->reload( &context, id, &resp );
            id.set_id(state->second.mask_id);
            viewer_proxy->reload( &context, id, &resp );
        }

        grpc::ClientContext context;
        ::rpc::img::InteractiveMaskOptions options;
        options.mutable_panel( )->set_id(state->first);
        options.set_niter(niter);
        options.set_cycleniter(cycleniter);
        options.set_threshold(thresh);
        options.set_cyclethreshold(cyclethresh);
        ::rpc::img::InteractiveMaskResult imresult;
        ::grpc::Status s = viewer_proxy->interactivemask( &context, options, &imresult );

        if ( ! s.ok( ) ) {
            std::cerr << "interactive mask failed: " << s.error_details( ) << std::endl;
            fflush(stderr);
        }

        niter = imresult.state( ).niter( );
        cycleniter = imresult.state( ).cycleniter( );
        thresh = imresult.state( ).threshold( );
        cyclethresh = imresult.state( ).cyclethreshold( );

        int result = 1;
        std::string action = imresult.action( );
        if ( action == "stop" ) result = 3;
        else if ( action == "no more" ) result = 2;
        else if ( action == "continue" ) result = 1;
        else {
            os << "ill-formed action result (" << action << ")" << LogIO::WARN << LogIO::POST;
            return 0;
        }

        state->second.prev_image_id = state->second.image_id;
        state->second.prev_mask_id = state->second.mask_id;
    
        state->second.image_id = 0;
        state->second.mask_id = 0;

        if ( debug ) {
            std::cerr << "set prev_image_id to " << state->second.prev_image_id << " (panel " << panel << ")" <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            std::cerr << "set prev_mask_id to " << state->second.prev_mask_id << " (panel " << panel << ")" <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        double endmask = maskSum(mask);

        if( startmask != endmask ) {
            result = -1 * result;
            LogIO os( LogOrigin("grpcInteractiveCleanGui",__FUNCTION__,WHERE) );
            os << "[" << mask << "] Mask modified from " << startmask << " pixels to " << endmask << " pixels " << LogIO::POST;
        }

        return result;
    }

    bool grpcInteractiveCleanGui::stop_viewer( ) {
        // viewer is not running...
        if ( ! viewer_started ) return false;
        static const auto debug = getenv("GRPC_DEBUG");
        if ( debug ) {
            std::cerr << "sending shutdown message to viewer (" << viewer_uri << ")" << 
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        // send shutdown message to viewer...
        grpc::ClientContext context;
        ::google::protobuf::Empty req;
        ::google::protobuf::Empty resp;
        auto shutdown = casatools::rpc::Shutdown::NewStub( grpc::CreateChannel( viewer_uri,
                                                                                grpc::InsecureChannelCredentials( ) ) );
        shutdown->now( &context, req, &resp );

        // wait on viewer (appimage) to exit...
        int status;
        pid_t w = waitpid( viewer_pid, &status, WUNTRACED | WCONTINUED );
        if ( w == -1 ){
            if ( debug ) {
                std::cerr << "viewer process waitpid failed " <<
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            // waitpid failed
            return false;
        } else if ( w == 0 ) {
            if ( debug ) {
                std::cerr << "viewer process not found " <<
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            return false;
        } else {
            if ( debug ) {
                std::cerr << "viewer process exited, status fetched " <<
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            return true;
        }

        viewer_pid = 0;
        viewer_proxy.release( );
        viewer_started = false;
        return true;
    }

    bool grpcInteractiveCleanGui::spawn_viewer( ) {
        static const auto debug = getenv("GRPC_DEBUG");

        std::string viewer_path = get_viewer_path( );
        if ( viewer_path.size( ) == 0 ) return false;

        // sanity check on viewer path...
        struct stat statbuf;
        if ( stat( viewer_path.c_str( ), &statbuf ) < 0 ) {
            // file (or dir) does not exist... e.g.
            //   >>>>>>registry available at 0.0.0.0:40939
            //   stopping registry<<<<<<
            return false;
        }

        std::string fifo = get_fifo( );
        if ( fifo.size( ) == 0 ) return false;

        // here we start the viewer in a very basic manner... we do not bother
        // with all of the theatrics needed to daemonize the launched process
        // (see https://stackoverflow.com/questions/17954432/creating-a-daemon-in-linux)
        // it could be that this should be done in the future, but for now we
        // will adopt the simple...

        int argc = 3;
        int logarg = argc;    // if a log file is specfied it comes last...
        std::string log_path = casatools::get_state( ).logPath( );
        if ( log_path.size( ) > 0 ) ++argc;

	    char **arguments = (char**) malloc(sizeof(char*) * (argc + 1));
        arguments[argc] = 0;
        arguments[0] = strdup(viewer_path.c_str( ));
        arguments[1] = (char*) malloc(sizeof(char) * (fifo.size( ) + 12));
        sprintf( arguments[1], "--server=%s", fifo.c_str( ) );
        arguments[2] = strdup("--oldregions");
        if ( log_path.size( ) > 0 ) {
            arguments[logarg] = (char*) malloc(sizeof(char) * (log_path.size( ) + 17));
            sprintf( arguments[logarg], "--casalogfile=%s", log_path.c_str( ) );
        }

        if ( debug ) {
            std::cerr << "forking viewer process: ";
            for (int i=0; i < argc; ++i) std::cout << arguments[i] << " ";
            std::cerr << " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        pid_t pid = fork( );

	    if ( pid == 0 ) {
            if ( debug ) {
                std::cerr << "execing viewer process: ";
                for (int i=0; i < argc; ++i) std::cout << arguments[i] << " ";
                std::cerr << " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            execvp( arguments[0], (char* const*) arguments );
            perror( "grpcInteractiveCleanGui::launch(...) child process exec failed" );
            exit(1);
	    }

	    for ( int i=0; i < argc; ++i ) free(arguments[i]);
	    free(arguments);

        if ( pid == -1 ) {
            perror( "grpcInteractiveCleanGui::launch(...) child process fork failed" );
            return false;
        }

        // perform a health check, after a delay...
        int status;
        sleep(2);
        pid_t w = waitpid( pid, &status, WUNTRACED | WCONTINUED | WNOHANG );
        if ( w == -1 ){
            if ( debug ) {
                std::cerr << "viewer process failed " <<
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            // waitpid failed
            return false;
        } else if ( w != 0 ) {
            if ( debug ) {
                std::cerr << "viewer process died " <<
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            // process exited
            if ( WIFEXITED(status) ) {
                printf("exited, status=%d\n", WEXITSTATUS(status));
            } else if (WIFSIGNALED(status)) {
                printf("killed by signal %d\n", WTERMSIG(status));
            } else if (WIFSTOPPED(status)) {
                printf("stopped by signal %d\n", WSTOPSIG(status));
            }
            return false;
        }

        if ( debug ) {
            std::cerr << "fetching viewer uri from " << fifo <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }
        char buffer[512];
        std::string uri_buffer;
        FILE *fp = fopen(fifo.c_str( ), "r");
        while ( fgets( buffer, sizeof(buffer), fp ) ) { uri_buffer = uri_buffer + buffer; }
        fclose(fp);
        trim(uri_buffer);

        // validate viewer uri...
        if ( ! std::regex_match( uri_buffer, std::regex("^([0-9]+\\.){3}[0-9]+:[0-9]+$") ) ) {
            //rework of regex required for IPv6...
            if ( debug ) {
                std::cerr << "bad viewer uri " << uri_buffer <<
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            return false;
        }

        if ( debug ) {
            std::cerr << "received viewer uri: " << uri_buffer <<
                " (process " << getpid( ) << ", thread " << 
                std::this_thread::get_id() << ")" << std::endl;
            fflush(stderr);
        }

        viewer_uri = uri_buffer;
        viewer_pid = pid;
        viewer_proxy = rpc::img::view::NewStub( grpc::CreateChannel( viewer_uri,
                                                                     grpc::InsecureChannelCredentials( ) ) );
        viewer_started = true;

        return true;
    }

    std::string grpcInteractiveCleanGui::get_python_path( ) {
        return casatools::get_state( ).pythonPath( );
    }

    std::string grpcInteractiveCleanGui::get_viewer_path( ) {
        std::string python_path = get_python_path( );
        if ( python_path.size( ) == 0 ) return std::string( );

        //*** python3 -m casaviewer --app-path
        char python_cmd[python_path.size( ) + 30];
        sprintf( python_cmd, "%s -m casaviewer --app-path", python_path.c_str( ) );
        std::array<char, 512> buffer;
        std::string result;
        std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(python_cmd, "r"), pclose);
        if ( ! pipe ) return std::string( );        //*** failed to start python
        while ( fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr ) {
            result += buffer.data();
        }
        trim(result);
        if ( result.size( ) == 0 ) return std::string( );
        return result;
    }

    std::string grpcInteractiveCleanGui::get_fifo( ) {
        static const char *env_tmpdir = getenv("TMPDIR");
        static std::string fifo_template = trim_trailing_slash(env_tmpdir && isdir(env_tmpdir) ? env_tmpdir : P_tmpdir) + "/vwr-XXXXXXXXXX";
        static int fifo_template_size = fifo_template.size( );
        char fifo_path[fifo_template_size+1];
        strncpy( fifo_path, fifo_template.c_str( ), fifo_template_size );
        fifo_path[fifo_template_size] = '\0';
        int fd = mkstemp(fifo_path);
        if ( fd == -1 ) throw std::runtime_error("mkstemp failed...");
        close( fd );
        unlink(fifo_path);
        mkfifo( fifo_path, 0666 );
        return fifo_path;
    }

    casacore::Record grpcInteractiveCleanManager::pauseForUserInteraction( ) {
		LogIO os( LogOrigin("grpcInteractiveCleanManager",__FUNCTION__,WHERE) );
        static const auto debug = getenv("GRPC_DEBUG");

        if ( clean_images.size( ) == 0 ) {
            // cannot open clean panel in viewer if not images are available...
            if ( debug ) {
                std::cerr << "no clean images available" <<
                    " (process " << getpid( ) << ", thread " << 
                    std::this_thread::get_id() << ")" << std::endl;
                fflush(stderr);
            }
            return Record( );
        }

        if ( clean_panel_id == -1 || ! gui.alive( ) ) {
            // open panel if it is not already open...
            clean_panel_id = gui.open_panel( clean_images );
        }
        
        int niter=0,cycleniter=0,iterdone;
        float threshold=0.0, cyclethreshold=0.0;
        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {
                           niter = state.Niter;
                           cycleniter = state.CycleNiter;
                           threshold = state.Threshold;
                           cyclethreshold = state.CycleThreshold;
                           iterdone = state.IterDone;
                           return dummy;
                       } ) );
                           
        std::string strthresh = std::to_string(threshold)+"Jy";
        std::string strcycthresh = std::to_string(cyclethreshold)+"Jy";
		  
        int iterleft = niter - iterdone;
        if( iterleft<0 ) iterleft=0;

        casacore::Vector<int> itsActionCodes(clean_images.size( ));
        itsActionCodes = 1.0;

        unsigned ind = 0;
        for ( auto it = clean_images.begin( ); it != clean_images.end( ); ++it, ++ind ) {
            if ( std::get<2>(*it) ) continue;
            if ( fabs(itsActionCodes[ind]) == 1.0 ) {
                std::string imageName = std::get<0>(*it) + ".residual" + ( std::get<1>(*it) ? ".tt0" : "" );
                std::string maskName = std::get<0>(*it) + ".mask";
                itsActionCodes[ind] = gui.interactivemask( clean_panel_id, imageName, maskName, iterleft,
                                                           cycleniter, strthresh, strcycthresh );
                if( itsActionCodes[ind] < 0 ) os << "[" << std::get<0>(*it) <<"] Mask changed interactively." << LogIO::POST;
                std::get<2>(*it) = (fabs(itsActionCodes[ind])==3);
            }
        }


        Quantity qa;
        casacore::Quantity::read(qa,strthresh);
        threshold = qa.getValue(Unit("Jy"));


	float oldcyclethreshold = cyclethreshold;
	Quantity qb;
	casacore::Quantity::read(qb,strcycthresh);
	cyclethreshold = qb.getValue(Unit("Jy"));

        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {
                           state.Niter = iterdone+iterleft;
                           state.CycleNiter = cycleniter;
                           state.Threshold = threshold;

                           if ( fabs(cyclethreshold) > 0 && fabs(oldcyclethreshold) > 0 &&
				fabs( cyclethreshold - oldcyclethreshold ) > 1e-06 )
                               state.CycleThreshold = cyclethreshold;
                           return dummy;
                       } ) );

		  Bool alldone=true;
		  for ( ind = 0; ind < clean_images.size( ); ++ind ) {
		      alldone = alldone & ( fabs(itsActionCodes[ind])==3 );
		  }
		  if( alldone==true ) changeStopFlag( true );

          Record returnRec;
          for( ind = 0; ind < clean_images.size( ); ind++ ){
              returnRec.define( RecordFieldId( String::toString(ind)), itsActionCodes[ind] );
          }

          return returnRec;
	}

    void grpcInteractiveCleanManager::closePanel( ) {
        gui.close_panel(clean_panel_id);
        clean_panel_id = -1;
        clean_images.clear( );
        access( (void*) 0,
                std::function< void* ( void*, grpcInteractiveCleanState& )>(
                       [&]( void *dummy, grpcInteractiveCleanState &state ) -> void* {
                           state.reset( );
                           return dummy; } ) );
    }

} //# NAMESPACE CASA - END
