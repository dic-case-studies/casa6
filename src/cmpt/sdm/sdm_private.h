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

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
// asdm2MS.cc
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

//
// A number of EnumSet to encode the different selection criteria.
//

EnumSet<CorrelationMode>         es_cm;
EnumSet<SpectralResolutionType>  es_srt;
EnumSet<TimeSampling>            es_ts;
Enum<CorrelationMode>            e_query_cm; 
EnumSet<AtmPhaseCorrection>      es_query_apc;    
EnumSet<AtmPhaseCorrection>      es_apc;

EnumSet<AtmPhaseCorrection> apcLiterals(const asdm::ASDM& ds);

//
// By default the resulting MS will not contain compressed columns
// unless the 'compress' option has been given.
// 
bool                             withCompression = false;

std::vector<int> dataDescriptionIdx2Idx;
std::map<asdm::MainRow*, int> stateIdx2Idx;
std::map<AtmPhaseCorrection, ASDM2MSFiller*> msFillers; // There will be one filler per value of the axis APC.
ASDM2MSFiller* msFiller;

bool hasCorrectedData(const EnumSet<AtmPhaseCorrection>& es);
bool hasUncorrectedData(const EnumSet<AtmPhaseCorrection>& es);

map<int, int> swIdx2Idx ;                       // A map which associates old and new index of Spectral Windows before/after reordering.
void fillSpectralWindow( asdm::ASDM* ds_p, map<unsigned int, double>& effectiveBwPerSpwId_m );
void fillEphemeris( asdm::ASDM* ds_p, uint64_t timeStepInNanoSecond, bool interpolate_ephemeris, string telescopeName );
void fillField( asdm::ASDM* ds_p, bool considerEphemeris );
void fillSysPower( const string asdmDirectory, asdm::ASDM* ds_p, bool ignoreTime, const vector<asdm::ScanRow *>& selectedScanRow_v,
                   map<AtmPhaseCorrection, ASDM2MSFiller*>& msFillers_m );
void fillMainLazily( const string& dsName, asdm::ASDM* ds_p, std::map<int, std::set<int> >& selected_eb_scan_m,
                     std::map<unsigned int , double>& effectiveBwPerDD_m, Enum<CorrelationMode> e_query_cm, bool checkdupints );
void fillState( asdm::MainRow* r_p );
void fillMain( asdm::MainRow* r_p, sdmbin::SDMBinData &sdmBinData, const sdmbin::VMSData *vmsData_p, UvwCoords &uvwCoords,
               std::map<unsigned int, double>& effectiveBwPerDD_m, bool complexData, bool mute,
               bool ac_xc_per_timestamp, bool skipFirstTime=false );

        
template<class T> void checkVectorSize( const string& vectorAttrName, const vector<T>& vectorAttr, const string& sizeAttrName,
                                        unsigned int sizeAttr, const string& tableName, unsigned int rowNumber );

bool gen_ms( const std::string &vis, bool createmms, const std::string &separationaxis, const variant &numsubms,    
             const std::string &corr_mode, const std::string &srt, const std::string &time_sampling,
             const std::string &ocorr_mode, bool compression, bool lazy, const std::string &asis,
             const std::string &wvr_corrected_data, const std::string &scans, bool ignore_time,
             bool process_syspower, bool process_caldevice, bool process_pointing, bool process_flags,
             double tbuff, bool applyflags, bool savecmds, const ::casac::variant& outfile, bool flagbackup,
             bool verbose, bool overwrite, bool showversion, const std::string &useversion, bool bdfflags,
             bool with_pointing_correction, bool remove_ref_undef, bool convert_ephem2geo, double polyephem_tabtimestep );
bool bd_flagger( const std::string &ms, const std::string &ocm, const std::string &scansOptionValue, bool lazy, bool processUncorrectedData );
void loadBDFlags( map<string, unsigned int>& abbrev2bitpos );
void processCorrelatorFlagsPerSlices( asdm::MainRow *mR_p, unsigned int iASDMIndex, const vector<int32_t> &mainRowIndex,
                                      const string &bdfPath, uint64_t bdfSliceSizeInMb, const vector<string> &antennas,
                                      const vector<pair<string, string> > &dataDescriptions, MSFlagEval &flagEval,
                                      unsigned &iMSRow, casacore::ArrayColumn<bool> &flag, casacore::ScalarColumn<bool> &flagRow,
                                      CorrelationModeMod::CorrelationMode ocorrelationMode );

public:
const bool verbose = true;

int sysPowerAntennaId(const asdm::SysPowerRow* row) const;
int sysPowerSpectralWindowId(const asdm::SysPowerRow* row) const;
int sysPowerFeedId(const asdm::SysPowerRow* row) const;
double sysPowerMidTimeInSeconds(const asdm::SysPowerRow* row) const;
double sysPowerIntervalInSeconds(const asdm::SysPowerRow* row) const;

void info(const string& message) const;
void error(const string& message, int status=1) const;
void warning (const string& message) const;
