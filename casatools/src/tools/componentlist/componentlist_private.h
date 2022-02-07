// Private members and functions
//
casacore::LogIO         *itsLog;
casa::ComponentList *itsList;
casa::ComponentList *itsBin;

vector<int> _checkIndices(int which,
                                     const casacore::String& function,
                                     const casacore::String& message) const;

vector<int> _checkIndices(const vector<int>& which,
                                     const casacore::String& function,
                                     const casacore::String& message) const;

void _checkIndex(int which) const;

int _checkFluxPol(const casacore::String &polString);
