#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <alma/ASDM/ASDM.h>
#include <alma/ASDM/AntennaRow.h>
#include <alma/ASDM/AntennaTable.h>
#include <alma/ASDM/ConfigDescriptionRow.h>
#include <alma/ASDM/ConfigDescriptionTable.h>
#include <alma/ASDM/DataDescriptionRow.h>
#include <alma/ASDM/DataDescriptionTable.h>
#include <alma/ASDM/ExecBlockRow.h>
#include <alma/ASDM/ExecBlockTable.h>
#include <alma/ASDM/MainRow.h>
#include <alma/ASDM/MainTable.h>
#include <alma/ASDM/PolarizationRow.h>
#include <alma/ASDM/PolarizationTable.h>
#include <alma/ASDM/ScanRow.h>
#include <alma/ASDM/ScanTable.h>
#include <alma/ASDM/SpectralWindowRow.h>
#include <alma/ASDM/SpectralWindowTable.h>
#include <alma/ASDM/StationRow.h>
#include <alma/ASDM/StationTable.h>
#include <alma/ASDM/SubscanRow.h>
#include <alma/ASDM/SubscanTable.h>
#include <alma/ASDM/Misc.h>

#include <alma/Enumerations/CAntennaMake.h>
#include <alma/Enumerations/CAtmPhaseCorrection.h>
#include <alma/Enumerations/CCorrelationMode.h>
#include <alma/Enumerations/CStokesParameter.h>
#include <alma/Enumerations/CFrequencyReferenceCode.h>
#include <alma/Enumerations/CScanIntent.h>
#include <alma/Enumerations/CSpectralResolutionType.h>
#include <alma/Enumerations/CSubscanIntent.h>
#include <alma/Enumerations/CTimeSampling.h>
#include <casa/Logging/StreamLogSink.h>
#include <casa/Logging/LogSink.h>

#include <alma/MS2ASDM/MS2ASDM.h>

#include <sdm_cmpt.h>

// A facility to get rid of blanks at start and end of a string.
//
string lrtrim(std::string& s,const std::string& drop = " ")
{
  std::string r=s.erase(s.find_last_not_of(drop)+1);
  return r.erase(0,r.find_first_not_of(drop));
}

using namespace std;

template<typename Enum, typename EnumHelper>
void output1 (typename  std::vector<Enum>::iterator begin, typename std::vector<Enum>::iterator end, ostringstream & oss) {
  if (begin == end) return;
  oss << ',' << EnumHelper::name(*begin);
  output1<Enum, EnumHelper>(begin+1, end, oss);
}

template<typename Enum, typename EnumHelper>
void output (typename std::vector<Enum>::iterator begin, typename std::vector<Enum>::iterator end, std::ostringstream & oss) {
  if (begin == end) return;
  oss << EnumHelper::name(*begin);
  output1<Enum, EnumHelper>(begin+1, end, oss);
}


template<typename T>
void output1 (typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end, ostringstream& oss) {
  if (begin == end) return;
  oss << "," << *begin;
  output1<T> (begin+1, end, oss);
}

template<typename T>
void output (typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end, ostringstream& oss) {
  if (begin == end) return;
  oss << *begin;
  output1<T> (begin+1, end, oss);
}

static bool notNull(int n) { return n != 0 ; }


namespace casac {
    sdm::sdm( const std::string &path ) {
        sdm_path = path;
    }

    bool sdm::fromms( const std::string &mspath, const std::string &datacolumn, const std::string &archiveid,
                      const std::string &rangeid, double subscanduration, double sbduration,
                      bool apcorrected, bool verbose) {
        struct stat path_stat;
        if ( stat( sdm_path.c_str( ), &path_stat ) != -1 ) {
            if ( S_ISREG(path_stat.st_mode) )
                throw runtime_error("SDM path exists and is a file");
            else if ( S_ISDIR(path_stat.st_mode) ) {
                DIR *dir = opendir(sdm_path.c_str( ));
                for ( int i=0; i <= 2 && readdir(dir); ++i ) {
                    if ( i == 2 ) throw runtime_error("SDM directory exists and is not empty");
                }
            }
        }
        casacore::MeasurementSet ms(mspath.c_str( ));
        casa::MS2ASDM m2a(ms);
        return m2a.writeASDM( sdm_path, datacolumn, archiveid, rangeid, verbose, subscanduration, sbduration, apcorrected );
    }

    std::string sdm::summarystr( ) {
        infostream.str("");
        infostream.clear( );
        infostream << "========================================================================================" << endl;
        infostream << "ASDM dataset :" << sdm_path << endl;
        infostream << "========================================================================================" << endl;

        asdm::ASDM ds;
        try {
            ds.setFromFile(sdm_path, asdm::ASDMParseOptions().loadTablesOnDemand(true).checkRowUniqueness(false));
            execBlockSummary(ds);
        } catch ( asdm::ConversionException ce ) {
            std::string result = std::string("ERROR[conversion-exception]: ") + ce.getMessage( );
            return result;
        }

        return infostream.str( );
    }

    void sdm::antennaSummary(const asdm::ExecBlockRow* eb_p) {
        asdm::ASDM& ds = eb_p->getTable().getContainer();
        asdm::AntennaTable& aT = ds.getAntenna();
        asdm::StationTable& sT = ds.getStation();

        // infostream.str("");
        const std::vector<asdm::Tag> antennaIds = eb_p->getAntennaId();
        asdm::AntennaRow * antenna_p = NULL ;
        asdm::StationRow * station_p = NULL ;
        infostream << endl;
        infostream << antennaIds.size() << " antennas have been used in this exec block." << endl;
        infostream << "        Id     Name         Make Station    Diameter         X              Y             Z" << endl;
        for (unsigned int i = 0; i < antennaIds.size(); i++) {
            antenna_p = aT.getRowByKey(antennaIds[i]);
            station_p = sT.getRowByKey(antenna_p->getStationId());
            std::vector<asdm::Length> position = station_p->getPosition();
            //infostream.fill('');
            infostream.width(12);infostream << antenna_p->getAntennaId() ;
            infostream.width(6); infostream.setf(ios::right); infostream   << antenna_p->getName() ;
            infostream.width(13); infostream  << CAntennaMake::name(antenna_p->getAntennaMake()) ;
            infostream.width(6); infostream   << station_p->getName() ;
            infostream.width(10); infostream.precision(10); infostream << antenna_p->getDishDiameter() ;
            infostream.width(15); infostream.precision(10); infostream << position[0] ;
            infostream.width(15); infostream.precision(10); infostream << position[1] ;
            infostream.width(15); infostream.precision(10); infostream << position[2] << endl;
        }
        // info(infostream.str());
    }

    sdm::SpectralWindowSummary sdm::spectralWindowSummary(asdm::SpectralWindowRow * spw_p) {
        SpectralWindowSummary result;

        result.numChan = spw_p->getNumChan();

        if (spw_p->isChanFreqStartExists())
            result.firstChan = spw_p->getChanFreqStart();
        else
            if (spw_p->isChanFreqArrayExists())
                result.firstChan = spw_p->getChanFreqArray()[0];
            else
                result.firstChan = asdm::Frequency(0.0);

        if (spw_p->isChanWidthArrayExists())
            result.chanWidth = spw_p->getChanWidthArray()[0];
        else
            if (spw_p->isChanWidthExists())
                result.chanWidth = spw_p->getChanWidth();
            else
                result.chanWidth = asdm::Frequency(0.0);

        if (spw_p->isMeasFreqRefExists())
            result.measFreqRef = CFrequencyReferenceCode::name(spw_p->getMeasFreqRef());
        else
            result.measFreqRef = "TOPO";

        result.refFreq = spw_p->getRefFreq();

        return result;
    }

    void sdm::mainSummary(asdm::ExecBlockRow* eb_p, int scanNumber, int subscanNumber) {

        asdm::ASDM& ds = eb_p->getTable().getContainer();

        asdm::Tag ebId = eb_p->getExecBlockId();

        const std::vector<asdm::MainRow *>& mains = ds.getMain().get();
        std::vector<asdm::MainRow *> eb_mains;

        for(asdm::MainRow* main_p: mains) {
            if ( main_p->getExecBlockId() == ebId && main_p->getScanNumber() == scanNumber &&
                 main_p->getSubscanNumber() == subscanNumber )
                eb_mains.push_back(main_p);
        }

        asdm::DataDescriptionTable& ddT = ds.getDataDescription();
        asdm::PolarizationTable& polT = ds.getPolarization();
        asdm::SpectralWindowTable& spwT = ds.getSpectralWindow();
        asdm::ConfigDescriptionTable& cfgDescT = ds.getConfigDescription();

        for ( asdm::MainRow* main_p: eb_mains ) {
            // infostream.str("");
            infostream << endl;
            infostream << "\t\t Binary data in " << main_p->getDataUID().getEntityId() << endl;
            infostream << "\t\t Number of integrations : " << main_p->getNumIntegration() << endl;
            infostream << "\t\t Time sampling : " << CTimeSampling::name(main_p->getTimeSampling()) << endl;
            asdm::ConfigDescriptionRow* cfgDesc_p = cfgDescT.getRowByKey(main_p->getConfigDescriptionId());
            infostream << "\t\t Correlation Mode : " << CCorrelationMode::name(cfgDesc_p->getCorrelationMode()) << endl;
            infostream << "\t\t Spectral resolution type : " << CSpectralResolutionType::name(cfgDesc_p->getSpectralType()) << endl;
            infostream << "\t\t Atmospheric phase correction : " ;
            std::vector<AtmPhaseCorrectionMod::AtmPhaseCorrection> apcs = cfgDesc_p->getAtmPhaseCorrection();
            output<AtmPhaseCorrectionMod::AtmPhaseCorrection, CAtmPhaseCorrection>(apcs.begin(), apcs.end(), infostream);
            infostream << endl;
            // info(infostream.str());

            std::vector<asdm::Tag> ddIds = cfgDesc_p->getDataDescriptionId();
            for ( asdm::Tag ddId: ddIds ) {
                asdm::DataDescriptionRow * dd_p = ddT.getRowByKey(ddId);
                asdm::SpectralWindowRow * spw_p = spwT.getRowByKey(dd_p->getSpectralWindowId());
                asdm::PolarizationRow * p_p = polT.getRowByKey(dd_p->getPolOrHoloId());
                // infostream.str("");
                SpectralWindowSummary spwSummary = spectralWindowSummary(spw_p);
                infostream << "\t\t " << spw_p->getSpectralWindowId() << " : numChan = " << spwSummary.numChan
                           << ", frame = " << spwSummary.measFreqRef
                           << ", firstChan = " << spwSummary.firstChan
                           << ", chandWidth = " << spwSummary.chanWidth
                           << " x "
                           << p_p->getPolarizationId() << " : corr = " ;
                std::vector<StokesParameterMod::StokesParameter> corrType = p_p->getCorrType();
                output<StokesParameterMod::StokesParameter, CStokesParameter>(corrType.begin(), corrType.end(), infostream);
                infostream << endl;
                // info(infostream.str());
            }
        }
    }

    void sdm::subscanSummary(asdm::ExecBlockRow* eb_p, int scanNumber) {

        asdm::ASDM& ds = eb_p->getTable().getContainer();
        asdm::Tag ebId = eb_p->getExecBlockId();

        const std::vector<asdm::SubscanRow *>& subscans = ds.getSubscan().get();
        std::vector<asdm::SubscanRow *> eb_subscans;
        for (asdm::SubscanRow * sscan_p: subscans) {
            if (sscan_p->getExecBlockId() == ebId && sscan_p->getScanNumber() == scanNumber)
                eb_subscans.push_back(sscan_p);
        }

        for (asdm::SubscanRow* sscan_p: eb_subscans) {
            // infostream.str("");
            infostream << "\tSubscan #" << sscan_p->getSubscanNumber()
                       << " from " << sscan_p->getStartTime().toFITS()
                       << " to " << sscan_p->getEndTime().toFITS()
                       << endl;
            infostream << "\t\tIntent : " << CSubscanIntent::name(sscan_p->getSubscanIntent()) << endl;
            infostream << "\t\tNumber of integrations : " << sscan_p->getNumIntegration() << endl;
            std::vector<int> numSubintegration = sscan_p->getNumSubintegration();
            if (find_if(numSubintegration.begin(), numSubintegration.end(), notNull) != numSubintegration.end()) {
                infostream << "\t\tNumber of subintegrations per integration : ";
                output<int>(numSubintegration.begin(), numSubintegration.end(), infostream);
                infostream << endl;
            }
            // info(infostream.str());

            mainSummary(eb_p, scanNumber, sscan_p->getSubscanNumber());
        }

    }

    void sdm::scanSummary(asdm::ExecBlockRow* eb_p) {

        asdm::ASDM& ds = eb_p->getTable().getContainer();
        asdm::Tag ebId = eb_p->getExecBlockId();

        const std::vector<asdm::MainRow *>& mains = ds.getMain().get();
        std::vector<asdm::MainRow *> eb_mains;

        for(asdm::MainRow* main: mains) {
            if ( main->getExecBlockId() == ebId) eb_mains.push_back(main);
        }

        const std::vector<asdm::ScanRow*>& scans = ds.getScan().get();
        std::vector<asdm::ScanRow *> eb_scans;
        for(asdm::ScanRow* scan: scans) {
            if ( scan->getExecBlockId() == ebId) eb_scans.push_back(scan);
        }

        // infostream.str("");
        infostream << endl;
        infostream << "Number of scans in this exec Block : " << eb_scans.size() << endl;
        // info(infostream.str());
        if (eb_scans.size() > 0) {
            for (asdm::ScanRow* scan_p: eb_scans) {
                // infostream.str("");
                infostream << endl;
                infostream << "scan #" << scan_p->getScanNumber()
                           << " from " << scan_p->getStartTime().toFITS()
                           << " to " <<  scan_p->getEndTime().toFITS()
                           << endl;

                std::vector<ScanIntentMod::ScanIntent> scis = scan_p->getScanIntent();
                if (scis.size() > 0) {
                    infostream << "\tIntents : ";
                    output<ScanIntentMod::ScanIntent, CScanIntent>(scis.begin(), scis.end(), infostream);
                    infostream << endl;
                }

                if ( scan_p->isFieldNameExists() ) {
                    std::vector<string> fields = scan_p->getFieldName();
                    if (fields.size() > 0) {
                        infostream << "\tFields : ";
                        output<string>(fields.begin(), fields.end(), infostream);
                        infostream << endl;
                    }
                }

                if ( scan_p->isSourceNameExists() ) {
                    infostream << "\tSources : " << scan_p->getSourceName() << endl;
                }
                // info(infostream.str());
                subscanSummary(eb_p, scan_p->getScanNumber());
            }
        }
    }

    void sdm::execBlockSummary(const asdm::ASDM& ds) {
        // infostream.str("");

        const std::vector<asdm::ExecBlockRow*>& ebs = ds.getExecBlock().get();
        for (unsigned int i = 0; i < ebs.size(); i++) {
            asdm::ExecBlockRow* eb_p = ebs[i];
            infostream << "\n";
            infostream << "Exec Block : " << eb_p->getExecBlockId() << endl;
            infostream << "Telescope : " << eb_p->getTelescopeName() << endl;
            infostream << "Configuration name : " << eb_p->getConfigName() << endl;
            infostream << "Observer name : " << eb_p->getObserverName() << endl;
            infostream << "The exec block started on " << eb_p->getStartTime().toFITS() << " and ended on " << eb_p->getEndTime().toFITS() << endl;
            if (eb_p->getAborted())
                infostream << "It was aborted." << endl;

            antennaSummary(eb_p);
            scanSummary(eb_p);
        }
    }

    sdm::~sdm( ) { }

}
