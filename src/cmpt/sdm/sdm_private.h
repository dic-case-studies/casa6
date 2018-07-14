std::ostringstream infostream;
std::ostringstream errstream;
std::string sdm_path;

typedef struct SpectralWindowSummaryStruct {
  int		numChan;
  string	measFreqRef;
  asdm::Frequency     firstChan;
  asdm::Frequency	chanWidth;
  asdm::Frequency	refFreq;
} SpectralWindowSummary;

void antennaSummary(const asdm::ExecBlockRow* eb_p);
SpectralWindowSummary spectralWindowSummary(asdm::SpectralWindowRow * spw_p);
void mainSummary(asdm::ExecBlockRow* eb_p, int scanNumber, int subscanNumber);
void subscanSummary(asdm::ExecBlockRow* eb_p, int scanNumber);
void scanSummary(asdm::ExecBlockRow* eb_p);
void execBlockSummary(const asdm::ASDM& ds);
