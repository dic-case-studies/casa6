
// --- //

template <typename T>
Bool calanalysis::writeInput( const CalAnalysis::OUTPUT<T>& oOutput,
    const casacore::uInt& row, const casacore::uInt& col, ::casac::record& oRecIter ) {

  // Write the field

  unsigned int uiField = oOutput.uiField;
  std::string oFieldString( casacore::String::toString(uiField).c_str() );
  oRecIter.insert( std::string("field"), oFieldString );


  // Write the antenna 1

  unsigned int uiAntenna1 = oOutput.uiAntenna1;
  std::string oAntenna1String( casacore::String::toString(uiAntenna1).c_str() );
  oRecIter.insert( std::string("antenna1"), oAntenna1String );


  // Write the antenna 2

  int iAntenna2 = oOutput.iAntenna2;
  std::string oAntenna2String( casacore::String::toString(iAntenna2).c_str() );
  oRecIter.insert( std::string("antenna2"), oAntenna2String );


  // Write the feed

  std::string oFeedKey( "feed" );
  std::string oFeedValue( oOutput.oOut->operator()(row,col).oAxes.sFeed.c_str() );

  oRecIter.insert( oFeedKey, oFeedValue );


  // Write the user-defined iteration axis

  std::string oAxisIterKey;
  CalStats::AXIS eAxis = oOutput.oOut->operator()(row,col).oAxes.eAxisIterUserID;

  if ( eAxis == CalStats::FREQUENCY ) {
    oAxisIterKey = std::string( "frequency" );
  } else {
    oAxisIterKey = std::string( "time" );
  }

  double dAxisIterValue = oOutput.oOut->operator()(row,col).oAxes.dAxisIterUser;
  oRecIter.insert( oAxisIterKey, dAxisIterValue );


  // Write the RAP parameter

  std::string oRAPString;

  if ( oOutput.eRAP == CalAnalysis::REAL ) {
    oRAPString = std::string( "REAL" );
  } else if ( oOutput.eRAP == CalAnalysis::AMPLITUDE ) {
    oRAPString = std::string( "AMPLITUDE" );
  } else {
    oRAPString = std::string( "PHASE" );
  }

  oRecIter.insert( std::string("rap"), oRAPString );


  // Write the amplitude normalization boolean

  if ( oOutput.eRAP == CalAnalysis::AMPLITUDE ) {
    oRecIter.insert( std::string("norm"), (bool) oOutput.bNorm );
  }


  // Write the phase unwrapping boolean and the maximum phase jump parameter

  if ( oOutput.eRAP == CalAnalysis::PHASE ) {
    oRecIter.insert( std::string("unwrap"), (bool) oOutput.bUnwrap );
    oRecIter.insert( std::string("jumpMax"), (double) oOutput.dJumpMax );
  }


  // Return true

  return( true );

}

// --- //

template <typename T>
Bool calanalysis::writeData( const CalAnalysis::OUTPUT<T>& oOutput,
    const casacore::uInt& row, const casacore::uInt& col, ::casac::record& oRecIter ) {

  // Write the non-iteration axis

  std::string oAxisNonIterKey;
  CalStats::AXIS eAxisNon = oOutput.oOut->operator()(row,col).oAxes.eAxisNonIterID;

  if ( eAxisNon == CalStats::FREQUENCY ) {
    oAxisNonIterKey = std::string( "frequency" );
  } else {
    oAxisNonIterKey = std::string( "time" );
  }

  oRecIter.insert( std::string("abscissa"), oAxisNonIterKey );


  // Get the number of abscissae

  casacore::uInt uiNumAbs = oOutput.oOut->operator()(row,col).oData.oAbs.nelements();


  // Write the abscissae from the non-iteration axis (either times or
  // frequencies)

  std::vector<double> oAbs( uiNumAbs );

  for ( casacore::uInt a=0; a<uiNumAbs; a++ ) {
    oAbs[a] = oOutput.oOut->operator()(row,col).oData.oAbs[a];
  }
  
  oRecIter.insert( oAxisNonIterKey, oAbs );


  // Write the values (either reals, amplitudes, or phases)

  std::vector<double> oValue( uiNumAbs );

  for ( casacore::uInt a=0; a<uiNumAbs; a++ ) {
    oValue[a] = oOutput.oOut->operator()(row,col).oData.oValue[a];
  }

  oRecIter.insert( string("value"), oValue );


  // Write the value errors (either real, amplitude, or phase errors)

  std::vector<double> oValueErr( uiNumAbs );

  for ( casacore::uInt a=0; a<uiNumAbs; a++ ) {
    oValueErr[a] = oOutput.oOut->operator()(row,col).oData.oValueErr[a];
  }

  oRecIter.insert( string("valueErr"), oValueErr );


  // Write the flags

  std::vector<bool> oFlag( uiNumAbs );

  for ( casacore::uInt a=0; a<uiNumAbs; a++ ) {
    oFlag[a] = oOutput.oOut->operator()(row,col).oData.oFlag[a];
  }

  oRecIter.insert( string("flag"), oFlag );


  // Return true

  return( true );

}

// --- //

template <typename T>
Bool calanalysis::writeFit( const CalStats::ARG<T>& oArg,
    const CalAnalysis::OUTPUT<T>& oOutput, const casacore::uInt& row, const casacore::uInt& col,
    ::casac::record& oRecIter ) {

  // Write the fit order, type, and weight parameters

  oRecIter.insert( std::string("order"),
      CalStatsFitter::orderName(oArg.eOrder) );

  oRecIter.insert( std::string("type"), CalStatsFitter::typeName(oArg.eType) );

  oRecIter.insert( std::string("weight"),
      CalStatsFitter::weightName(oArg.eWeight) );


  // Write the fit validity flag

  bool valid = oOutput.oOut->operator()(row,col).oT->bValid;

  oRecIter.insert( std::string("validFit"), valid );


  // Write the fit reduced chi2 value

  double redchi2 = oOutput.oOut->operator()(row,col).oT->dRedChi2;

  oRecIter.insert( std::string("redChi2"), redchi2 );


  // Write the fit parameters

  casacore::uInt uiNumPars = oOutput.oOut->operator()(row,col).oT->oPars.nelements();
  std::vector<double> oPars( uiNumPars );

  for ( casacore::uInt p=0; p<uiNumPars; p++ ) {
    oPars[p] = oOutput.oOut->operator()(row,col).oT->oPars[p];
  }

  oRecIter.insert( std::string("pars"), oPars );


  // Write the fit variances

  std::vector<double> oVars( uiNumPars );

  for ( casacore::uInt p=0; p<uiNumPars; p++ ) {
    oVars[p] = oOutput.oOut->operator()(row,col).oT->oCovars(p,p);
  }

  oRecIter.insert( std::string("vars"), oVars );


  // Write the fit covariances

  std::vector<double> oCovars( uiNumPars * (uiNumPars-1) / 2 );

  for ( casacore::uInt pr=0,p=0; pr<uiNumPars; pr++ ) {
    for ( casacore::uInt pc=pr+1; pc<uiNumPars; pc++,p++ ) {
      oCovars[p] = oOutput.oOut->operator()(row,col).oT->oCovars(pr,pc);
    }
  }

  oRecIter.insert( std::string("covars"), oCovars );


  // Write the fit model

  casacore::uInt uiNumData = oOutput.oOut->operator()(row,col).oT->oModel.nelements();
  std::vector<double> oModel( uiNumData );

  for ( casacore::uInt d=0; d<uiNumData; d++ ) {
    oModel[d] = oOutput.oOut->operator()(row,col).oT->oModel[d];
  }

  oRecIter.insert( std::string("model"), oModel );


  // Write the fit residuals

  std::vector<double> oRes( uiNumData );

  for ( casacore::uInt d=0; d<uiNumData; d++ ) {
    oRes[d] = oOutput.oOut->operator()(row,col).oT->oRes[d];
  }

  oRecIter.insert( std::string("res"), oRes );


  // Write the fit residual variance

  casacore::Double dResVar = oOutput.oOut->operator()(row,col).oT->dResVar;

  oRecIter.insert( std::string("resVar"), dResVar );


  // Write the fit residual mean

  casacore::Double dResMean = oOutput.oOut->operator()(row,col).oT->dResMean;

  oRecIter.insert( std::string("resMean"), dResMean );


  // Return true

  return( true );

}
