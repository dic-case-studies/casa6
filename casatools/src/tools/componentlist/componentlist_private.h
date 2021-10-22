// Private members and functions
//
casacore::LogIO         *itsLog;
casa::ComponentList *itsList;
casa::ComponentList *itsBin;

vector<int> checkIndices(int which,
                                     const casacore::String& function,
                                     const casacore::String& message) const;

vector<int> checkIndices(const vector<int>& which,
                                     const casacore::String& function,
                                     const casacore::String& message) const;

void _checkIndex(int which) const;

int checkFluxPol(const casacore::String &polString);
