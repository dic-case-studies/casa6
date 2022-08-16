
#include <iostream>
#include <mstransformer_cmpt.h>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogOrigin.h>
#include <casacore/casa/Exceptions/Error.h>
#include <mstransform/MSTransform/MSTransform.h>
#include <mstransform/MSTransform/MSTransformDataHandler.h>
#include <casacore/casa/Containers/RecordInterface.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/sstream.h>
#include <unistd.h>
#include <stdcasa/StdCasa/CasacSupport.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/measures/Measures/MeasIERS.h>


using namespace std;
using namespace casacore;
using namespace casa;

using namespace casacore;
namespace casac {

mstransformer::mstransformer(): logger_p(new LogIO(LogOrigin("mstransformer","",WHERE)))
{
	try
	{
		mstransformer_p = new MSTransform();

	} catch (AipsError x) {
		*logger_p << LogIO::DEBUG1
				<< "Stack Trace: " << x.getStackTrace()
				<< LogIO::POST;
		*logger_p 	<< LogIO::DEBUG1
					<< "Exception Reported: " << x.getMesg()
					<< LogIO::POST;
		RETHROW(x);
	}

}

mstransformer::~mstransformer()
{
	try
	{
		done();

	} catch (AipsError x) {
	    Table::relinquishAutoLocks(true);
		*logger_p << LogIO::DEBUG1 		<< "Stack Trace: " << x.getStackTrace()
										<< LogIO::POST;
		*logger_p << LogIO::DEBUG1 		<< "Exception Reported: " << x.getMesg()
										<< LogIO::POST;
		RETHROW(x);
	}
    Table::relinquishAutoLocks(true);
    MeasIERS::closeTables();

}

bool
mstransformer::done()
{
	try
	{
		if (mstransformer_p) {
			delete mstransformer_p;
			mstransformer_p = nullptr;
		}

	} catch (AipsError x) {
	    Table::relinquishAutoLocks(true);
		*logger_p << LogIO::DEBUG1 		<< "Stack Trace: " << x.getStackTrace()
										<< LogIO::POST;
		*logger_p << LogIO::DEBUG1 		<< "Exception Reported: " << x.getMesg()
										<< LogIO::POST;
		RETHROW(x);
	}
    Table::relinquishAutoLocks(true);
	return true;
}

bool
mstransformer::config(const ::casac::record& pars)
{
	try
	{
		if (!mstransformer_p) {
			mstransformer_p = new MSTransform();
		}

		if (!pars.empty()){
			Record mstpars = *toRecord(pars);
			return mstransformer_p->configure(mstpars);
		}
		else {
			return false;
		}

		return true;

	} catch(AipsError x) {
		*logger_p << LogIO::DEBUG1 		<< "Stack Trace: " << x.getStackTrace()
										<< LogIO::POST;
		*logger_p << LogIO::DEBUG1 		<< "Exception Reported: " << x.getMesg()
										<< LogIO::POST;
		RETHROW(x);
	}

}

bool
mstransformer::open()
{
	try
	{
		if (!mstransformer_p) {
			mstransformer_p = new MSTransform();
		}
		if (mstransformer_p) {
			return mstransformer_p->open();
		}

		return false;
	} catch(AipsError x) {
		*logger_p << LogIO::DEBUG1 		<< "Stack Trace: " << x.getStackTrace()
										<< LogIO::POST;
		*logger_p << LogIO::DEBUG1 		<< "Exception Reported: " << x.getMesg()
										<< LogIO::POST;
		RETHROW(x);
	}

}


::casac::record*
 mstransformer::run()
{
	casac::record *rstat(0);
	try
	{
		if(mstransformer_p){
			rstat =  fromRecord(mstransformer_p->run());
		}
		else{
			rstat = fromRecord(Record());
		}

		return rstat;
	} catch(AipsError x){
		*logger_p << LogIO::DEBUG1 		<< "Stack Trace: " << x.getStackTrace()
										<< LogIO::POST;
		*logger_p << LogIO::DEBUG1 		<< "Exception Reported: " << x.getMesg()
										<< LogIO::POST;
		RETHROW(x);

	}
}

bool
mstransformer::mergespwtables(const std::vector<std::string> &filenames)
{
	try
	{
		if (!mstransformer_p) {
			mstransformer_p = new MSTransform();
		}

		// Call the static method directly
		MSTransformDataHandler::mergeSpwSubTables(toVectorString(filenames));

	} catch(AipsError x){
		*logger_p << LogIO::DEBUG1 		<< "Stack Trace: " << x.getStackTrace()
										<< LogIO::POST;
		*logger_p << LogIO::DEBUG1 		<< "Exception Reported: " << x.getMesg()
										<< LogIO::POST;
		RETHROW(x);
	}

	return true;

}



} // casac namespace
