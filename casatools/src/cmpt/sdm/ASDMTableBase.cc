#include "ASDMTableBase.h"

using namespace casacore;

namespace asdm {

    ASDM_TABLE_BASE::ASDM_TABLE_BASE() {table_p_ = 0;}

    ASDM_TABLE_BASE::~ASDM_TABLE_BASE() {
        close();
    }

    void ASDM_TABLE_BASE::close() {
        if (table_p_ != 0) delete table_p_;
        table_p_ = 0;
    }

    casacore::Table* ASDM_TABLE_BASE::table_p() {return table_p_;}

    const string& ASDM_TABLE_BASE::name() const { return name_; }

    void ASDM_TABLE_BASE::buildAndAttachTable(casacore::MS* attachMS) {
        casacore::SetupNewTable tableSetup(attachMS->tableName() + "/" + String(name_), tableDesc(), casacore::Table::New);
        table_p_ = new casacore::Table(tableSetup, casacore::TableLock(casacore::TableLock::PermanentLockingWait));
        AlwaysAssert(table_p_, casacore::AipsError);
        attachMS->rwKeywordSet().defineTable(name_, *table_p_);
    }

} // end namespace asdm


