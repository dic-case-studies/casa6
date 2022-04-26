public:
tablerow( );
private:

//Privates for table_cmpt
//
//  These are internal variables the connect us to the objects
//
casacore::LogIO                       *itsLog;
casacore::TableRowProxy               *itsRow;
std::shared_ptr<casacore::TableProxy> itsTable;

//
// Private constructor so we can make components on the fly
//
friend class casac::table;
tablerow(std::shared_ptr<casacore::TableProxy> myTable, const std::vector<std::string> &columnnames, bool exclude );
