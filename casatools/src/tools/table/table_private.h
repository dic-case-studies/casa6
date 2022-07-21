#include <utility>
friend class tablerow;

//Privates for table_cmpt
//
//  These are internal variables the connect us to the objects
//
casacore::LogIO      *itsLog;

// shared with casac::tablerow
std::shared_ptr<casacore::TableProxy> itsTable;

// rows created by this table
std::list<tablerow*> created_rows;
void remove_tablerow( tablerow *tr ) { created_rows.remove( tr ); }
void remove_all_tablerows( ) {
    // clear all outstanding tablerow objects
    for ( auto *tr : created_rows ) tr->reset( );
    created_rows.clear( );
}

//
// Private constructor so we can make components on the fly
//
table(casacore::TableProxy *myTable);
//
casacore::TableLock *getLockOptions(casac::record &lockoptions);

void _checkCorner(
    const std::vector<long>& corner, const casacore::String& name,
    const casacore::IPosition& shape,
    const std::pair<std::vector<long>, std::vector<long>> *const &blctrc=nullptr
);

