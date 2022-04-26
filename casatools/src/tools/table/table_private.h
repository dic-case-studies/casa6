//Privates for table_cmpt
//
//  These are internal variables the connect us to the objects
//
casacore::LogIO      *itsLog;

// shared with casac::tablerow
std::shared_ptr<casacore::TableProxy> itsTable;

//
// Private constructor so we can make components on the fly
//
table(casacore::TableProxy *myTable);
//
casacore::TableLock *getLockOptions(casac::record &lockoptions);
