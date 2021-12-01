#if ! defined (CASA_VI_TEST_MsFactory_H_130114_1357)
#define CASA_VI_TEST_MsFactory_H_130114_1357

#include <casa/string.h>
#include <ms/MeasurementSets.h>
#include <msvis/MSVis/UtilJ.h>
#include <casacore/ms/MSOper/NewMSSimulator.h>
#include <casacore/measures/Measures/MeasTable.h>

#include <utility>

namespace casa {

    class ROVisibilityIterator;

namespace vi {

    class VisibilityIterator2;
    class VisBuffer2;

namespace test {

class FillState {

public:

    casacore::Int antenna1_p;
    casacore::Int antenna2_p;
    casacore::Int nAntennas_p;
    casacore::Int nFlagCategories_p;
    casacore::Int nChannels_p;
    casacore::Int nCorrelations_p;
    casacore::Int rowNumber_p;
    casacore::Int spectralWindow_p;
    casacore::Double time_p;
    casacore::Double timeDelta_p;
    casacore::Double uvw_p [3];

};


class GeneratorBase {

public:

    virtual ~GeneratorBase () {}
};

template <typename T>
class Generator : public GeneratorBase {
public:

    virtual ~Generator () {}

    virtual T operator() (const FillState &, casacore::Int channel, casacore::Int correlation) const = 0;

};

template <typename T>
class GenerateConstant : public Generator<T>{

public:

    GenerateConstant (const T & c) : c_p (c) {}

    T
    operator () (const FillState &, casacore::Int, casacore::Int ) const {
        return c_p;
    }

private:

    T c_p;
};

class GenerateAntenna1 : public Generator<casacore::Int> {

public:

    casacore::Int
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {
        return fillState.antenna1_p;
    }
};

class GenerateAntenna2 : public Generator<casacore::Int> {


public:

    casacore::Int
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {
        return fillState.antenna2_p;
    }
};

class GenerateDdi: public Generator<casacore::Int> {

public:

    casacore::Int
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {
        return fillState.spectralWindow_p;
    }
};

class GenerateExposure : public Generator<casacore::Double> {

public:

    casacore::Double
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {
        return fillState.timeDelta_p * .98;
    }
};

class GenerateFlag : public Generator<casacore::Bool> {

public:

    casacore::Bool
    operator () (const FillState & /*fillState*/, casacore::Int channel, casacore::Int correlation) const {

        // Generate a of T F T F ... for even rows and F T F T ... for odd rows

        casacore::Bool value = (channel % 2 == 0) == (correlation % 2 == 0);
        return value;
    }
};

class GenerateFlagRow : public Generator<casacore::Bool> {

public:

    casacore::Bool
    operator () (const FillState & fillState, casacore::Int /*channel*/, casacore::Int /*correlation*/) const {

        // Generate a of T F T F ... for even rows and F T F T ... for odd rows

        casacore::Bool value = (fillState.rowNumber_p % 3 == 0) || (fillState.rowNumber_p % 5 == 0);
        return value;
    }
};

class GenerateFlagCategory : public GeneratorBase {

public:

    virtual casacore::Bool
    operator () (const FillState & fillState, casacore::Int channel, casacore::Int correlation, casacore::Int /*category*/) const {

        casacore::Bool result;

        switch (correlation){

        case 0:

            result = fillState.rowNumber_p & 0x1;
            break;

        case 1:

            result = fillState.rowNumber_p & 0x2;
            break;

        case 2:

            result = channel & 0x1;
            break;

        case 3:

            result = channel & 0x2;
            break;

        default:

            Assert (false);
            result = false;
            break;

        }

        return result;
    }
};


template <typename T>
class GenerateUsingRow : public Generator<T> {


public:

    GenerateUsingRow (T c) : c_p (c) {}

    T
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {
        return fillState.rowNumber_p * 100 + c_p;
    }

private:

    const T c_p;
};


class GenerateTime : public Generator<casacore::Double> {


public:

    casacore::Double
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {
        return fillState.time_p;
    }
};

class GenerateTimeCentroid : public Generator<casacore::Double> {

public:

    casacore::Double
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {
        return fillState.time_p; // same as time for this model
    }
};

class GenerateTimeInterval : public Generator<casacore::Double> {

public:

    casacore::Double
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {
        return fillState.timeDelta_p;
    }
};

class GenerateUvw : public Generator <casacore::Vector <casacore::Double> > {

    casacore::Vector<casacore::Double>
    operator () (const FillState & fillState, casacore::Int, casacore::Int ) const {

        casacore::Vector <casacore::Double> result (3);
        result [0] = fillState.rowNumber_p * 10;
        result [1] = result [0] + 1;
        result [2] = result [0] + 2;
        return result;
    }
};

class GenerateVisibility : public Generator<casacore::Complex> {

public:

    GenerateVisibility (casacore::Int c) : c_p (c) {}

    casacore::Complex
    operator () (const FillState & fillState, casacore::Int channel, casacore::Int correlation) const {

        casacore::Float r = fillState.rowNumber_p * 10 + fillState.spectralWindow_p;
        casacore::Float i = channel * 100 + correlation * 10 + c_p;
        return casacore::Complex (r, i);
    }

private:

    const casacore::Int c_p;
};

class GenerateFloatData : public Generator<casacore::Float> {

public:

    GenerateFloatData (casacore::Int c) : c_p (c) {}

    casacore::Float
    operator () (const FillState & fillState, casacore::Int channel, casacore::Int correlation) const {

        casacore::Float r = fillState.rowNumber_p * 10 + 
            fillState.spectralWindow_p + channel * 100 + correlation * 1000;
        return r;
    }

private:

    const casacore::Int c_p;
};


class GenerateWeightSpectrum : public Generator<casacore::Float> {

public:

    casacore::Float
    operator () (const FillState & fillState, casacore::Int channel, casacore::Int correlation) const {

        casacore::Float r = fillState.rowNumber_p * 1000 + fillState.spectralWindow_p * 100 +
                  channel * 10 + correlation + 1;
        return r;
    }
};

class GenerateWeightSpectrumCorrected : public Generator<casacore::Float> {

public:

    casacore::Float
    operator () (const FillState & fillState, casacore::Int channel, casacore::Int correlation) const {

        casacore::Float r = fillState.rowNumber_p * 1000 + fillState.spectralWindow_p * 100 +
                  channel * 10 + correlation + 2;
        return r;
    }
};


class MsFactory {

public:

    MsFactory (const casacore::String & name);

    ~MsFactory ();

    std::pair<casacore::MeasurementSet *, casacore::rownr_t> createMs ();

    void addAntennas (casacore::Int nAntennas);

    //// void addColumn (casacore::MSMainEnums::PredefinedColumns columnId);
    void addWeightSpectrum (casacore::Bool addIt);
    void addFloatData (casacore::Bool addIt);

    void addCubeColumn (casacore::MSMainEnums::PredefinedColumns columnId,
                        const casacore::String & dataStorageMangerName);

    void removeColumn (casacore::MSMainEnums::PredefinedColumns columnId);

    void addFeeds (casacore::Int nFeeds);
    void addField (const casacore::String & name,
                   const casacore::MDirection & direction);
    void addSpectralWindow (const casacore::String & name,
                            casacore::Int nChannels,
                            casacore::Double frequency,
                            casacore::Double frequencyDelta,
                            const casacore::String & stokes,
                            casacore::Double endingTime = -1);
    void addSpectralWindows (int nSpectralWindows);

    void setDataGenerator (casacore::MSMainEnums::PredefinedColumns, GeneratorBase * generator);

    void setIncludeAutocorrelations (casacore::Bool includeThem);

    void setTimeInfo (casacore::Double startingTime,
                      casacore::Double endingTime,
                      casacore::Double interval);

protected:

    void addColumns ();
    void addDefaults ();
    void attachColumns ();
    void fillCollections (FillState & fillState);
    template <typename T>
    void
    fillCube (casacore::ArrayColumn<T> & column, const FillState & fillState,
              const GeneratorBase * generator);
    void fillCubes (FillState & fillState);
    void fillData ();
    void fillRows (FillState & fillState);

    template <typename T>
    void
    fillScalar (casacore::ScalarColumn<T> & column, const FillState & fillState,
                const GeneratorBase * generator);
    template <typename T>
    void
    fillCollection (casacore::ArrayColumn<T> & column, const FillState & fillState,
                    const GeneratorBase * generator);
    void fillScalars (FillState & fillState);
    void fillVisCubeCorrected (FillState & fillState);
    void fillVisCubeModel (FillState & fillState);
    void fillVisCubeObserved (FillState & fillState);
    void fillVisFloatData (FillState & fillState);
    void fillWeight (FillState & fillState);
    void fillSigma (FillState & fillState);
    void fillUvw (FillState & fillState);
    void fillWeightSpectrumCube (FillState & fillState);
    void fillWeightSpectrumCorrectedCube (FillState & fillState);
    void fillFlagCube (FillState & fillState);
    void fillFlagCategories (const FillState & fillState);

private:

    class Columns {

    public:

        Columns & operator= (const Columns & other);

        casacore::ScalarColumn<casacore::Int>    antenna1_p;
        casacore::ScalarColumn<casacore::Int>    antenna2_p;
        casacore::ScalarColumn<casacore::Int>    array_p;
        casacore::ArrayColumn<casacore::Complex> corrVis_p;
        casacore::ScalarColumn<casacore::Int>    dataDescriptionId_p;
        casacore::ScalarColumn<casacore::Double> exposure_p;
        casacore::ScalarColumn<casacore::Int>    feed1_p;
        casacore::ScalarColumn<casacore::Int>    feed2_p;
        casacore::ScalarColumn<casacore::Int>    field_p;
        casacore::ArrayColumn<casacore::Bool>    flagCategory_p;
        casacore::ScalarColumn<casacore::Bool>   flagRow_p;
        casacore::ArrayColumn<casacore::Bool>    flag_p;
        casacore::ArrayColumn<casacore::Float>   floatVis_p;
        casacore::ArrayColumn<casacore::Complex> modelVis_p;
        casacore::ScalarColumn<casacore::Int>    observation_p;
        casacore::ScalarColumn<casacore::Int>    processor_p;
        casacore::ScalarColumn<casacore::Int>    scan_p;
        casacore::ArrayColumn<casacore::Float>   sigma_p;
        casacore::ScalarColumn<casacore::Int>    state_p;
        casacore::ScalarColumn<casacore::Double> timeCentroid_p;
        casacore::ScalarColumn<casacore::Double> timeInterval_p;
        casacore::ScalarColumn<casacore::Double> time_p;
        casacore::ArrayColumn<casacore::Double>  uvw_p;
        casacore::ArrayColumn<casacore::Complex> vis_p;
        casacore::ArrayColumn<casacore::Float>   weightSpectrum_p;
        casacore::ArrayColumn<casacore::Float>   weightSpectrumCorrected_p;
        casacore::ArrayColumn<casacore::Float>   weight_p;

    };

    class Generators {

    public:

        typedef std::map<casacore::MSMainEnums::PredefinedColumns, const GeneratorBase *> GeneratorMap;

        Generators ()
        {
            generatorMap_p [casacore::MSMainEnums::ANTENNA1] = new GenerateAntenna1 ();
            generatorMap_p [casacore::MSMainEnums::ANTENNA2] = new GenerateAntenna2 ();
            generatorMap_p [casacore::MSMainEnums::ARRAY_ID] = new GenerateConstant<casacore::Int> (17);
            generatorMap_p [casacore::MSMainEnums::DATA_DESC_ID] = new GenerateDdi ();
            generatorMap_p [casacore::MSMainEnums::EXPOSURE] = new GenerateUsingRow<casacore::Double> (8);
            generatorMap_p [casacore::MSMainEnums::FEED1] = new GenerateConstant<casacore::Int> (0);
            generatorMap_p [casacore::MSMainEnums::FEED2] = new GenerateConstant<casacore::Int> (0);
            generatorMap_p [casacore::MSMainEnums::FIELD_ID] = new GenerateConstant<casacore::Int> (0);
            generatorMap_p [casacore::MSMainEnums::FLAG_CATEGORY] = new GenerateFlagCategory ();
            generatorMap_p [casacore::MSMainEnums::FLAG_ROW] = new GenerateFlagRow ();
            generatorMap_p [casacore::MSMainEnums::OBSERVATION_ID] = new GenerateUsingRow<casacore::Int> (7);
            generatorMap_p [casacore::MSMainEnums::PROCESSOR_ID] = new GenerateConstant<casacore::Int> (0);
            generatorMap_p [casacore::MSMainEnums::SCAN_NUMBER] = new GenerateUsingRow<casacore::Int> (5);
            generatorMap_p [casacore::MSMainEnums::STATE_ID] = new GenerateConstant<casacore::Int> (0);
            generatorMap_p [casacore::MSMainEnums::TIME_CENTROID] = new GenerateUsingRow<casacore::Double> (1);
            generatorMap_p [casacore::MSMainEnums::INTERVAL] = new GenerateUsingRow<casacore::Double> (2);
            generatorMap_p [casacore::MSMainEnums::TIME] = new GenerateTime ();
            generatorMap_p [casacore::MSMainEnums::FLAG] = new GenerateFlag ();
            generatorMap_p [casacore::MSMainEnums::SIGMA] = new GenerateUsingRow<casacore::Float> (3);
            generatorMap_p [casacore::MSMainEnums::SIGMA_SPECTRUM] = new GenerateConstant<casacore::Float> (0);
            generatorMap_p [casacore::MSMainEnums::UVW] = new GenerateUvw ();
            generatorMap_p [casacore::MSMainEnums::CORRECTED_DATA] = new GenerateVisibility (1);
            generatorMap_p [casacore::MSMainEnums::MODEL_DATA] = new GenerateVisibility (2);
            generatorMap_p [casacore::MSMainEnums::DATA] = new GenerateVisibility (0);
            generatorMap_p [casacore::MSMainEnums::FLOAT_DATA] = new GenerateFloatData (3);
            generatorMap_p [casacore::MSMainEnums::WEIGHT] = new GenerateUsingRow<casacore::Float> (4);
            generatorMap_p [casacore::MSMainEnums::WEIGHT_SPECTRUM] = new GenerateWeightSpectrum ();
            generatorMap_p [casacore::MSMainEnums::CORRECTED_WEIGHT_SPECTRUM] = new GenerateWeightSpectrumCorrected ();

        }

        ~Generators ()
        {
            for (GeneratorMap::const_iterator i = generatorMap_p.begin();
                 i != generatorMap_p.end();
                 i ++){

                delete i->second;
            }
        }

        const GeneratorBase *
        get (casacore::MSMainEnums::PredefinedColumns key) const
        {
            GeneratorMap::const_iterator i = generatorMap_p.find (key);

            ThrowIf (i == generatorMap_p.end(),
                     casacore::String::format ("No such key: %d", key));

            return i->second;
        }

        void
        set (casacore::MSMainEnums::PredefinedColumns key, const GeneratorBase * newValue)
        {
            GeneratorMap::iterator i = generatorMap_p.find (key);

            if (i != generatorMap_p.end()){
                delete i->second;
                i->second =  newValue;
            }
            else{
                generatorMap_p [key] = newValue;
            }

        }

    private:

        GeneratorMap generatorMap_p;

    };

    casacore::Bool addWeightSpectrum_p;
    casacore::Bool addFloatData_p;
    Columns columns_p;
    Generators generators_p;
    casacore::Bool includeAutocorrelations_p;
    casacore::MeasurementSet * ms_p;
    casacore::Int nAntennas_p;
    casacore::rownr_t nRows_p;
    std::unique_ptr<casacore::NewMSSimulator> simulator_p;
    casacore::Double timeEnd_p;
    casacore::Double timeInterval_p;
    casacore::Double timeStart_p;
    std::map<casacore::Int, casacore::Double> endingTimePerSpw_p;
};

inline MsFactory::MsFactory (const casacore::String & msName)
 : addWeightSpectrum_p (true),
   addFloatData_p (false),
   includeAutocorrelations_p (false),
   simulator_p (new casacore::NewMSSimulator (msName)),
   timeStart_p (-1)
{
    ms_p = new casacore::MeasurementSet (* simulator_p->getMs ()); //
}

inline MsFactory::~MsFactory ()
{
}


inline void MsFactory::addAntennas (casacore::Int nAntennas)
{
    casacore::Vector<casacore::Double> x (nAntennas), y (nAntennas), z (nAntennas),
                   diameter (nAntennas), offset (nAntennas);

    casacore::Vector<casacore::String> mount (nAntennas), name (nAntennas), pad (nAntennas);

    for (casacore::Int i = 0; i < nAntennas; i++){

        casacore::Double angle = ((i - 1) % 3) * (2 * 3.14159 / 3.0);
        casacore::Double radius = (i - 1) / 3.0 * 100;

        x [i] = radius * cos (angle);
        y [i] = radius * sin (angle);
        z [i] = 0;

        name [i] = casacore::String::format ("a%02d", i);
        pad [i] = casacore::String::format ("p%02d", i);
    }

    diameter = 10;
    offset = 0;
    mount = "ALT-AZ";

    casacore::MPosition vlaPosition;
    casacore::MeasTable::Observatory(vlaPosition, "VLA");

    simulator_p->initAnt ("Simulated", x, y, z, diameter, offset,
                          mount, name, pad, "local", vlaPosition);

    nAntennas_p = nAntennas;

}


inline void MsFactory::addDefaults ()
{
    // Configure fields if not present.

    casacore::Int nFields;
    casacore::Vector<casacore::String> sourceName;
    casacore::Vector<casacore::MDirection> sourceDirection;
    casacore::Vector<casacore::String> calCode;

    simulator_p->getFields (nFields, sourceName, sourceDirection, calCode);

    if (nFields == 0){

        casacore::Quantity ra (85.25, "deg");   // 05 41.0  -02 25 (horsehead nebula)
        casacore::Quantity dec (-2.417, "deg");
        casacore::MDirection direction (ra, dec);

        addField ("HorseHeadNebula", direction);
    }

    // Configure antennas if not present

    casacore::String telescope;
    casacore::Int nAntennas;
    casacore::Matrix<casacore::Double> xyz;
    casacore::Vector<casacore::Double> diameter;
    casacore::Vector<casacore::Double> offset;
    casacore::Vector<casacore::String> mount;
    casacore::Vector<casacore::String> name;
    casacore::Vector<casacore::String> pad;
    casacore::String coordinateSystem;
    casacore::MPosition referenceLocation;

    casacore::Bool ok = simulator_p->getAnt (telescope, nAntennas, & xyz, diameter, offset,
                                   mount, name, pad, coordinateSystem, referenceLocation);

    if (! ok){
        addAntennas (4);
    }

    // Configure feeds

    casacore::Vector<casacore::Double> x (2) , y (2);
    casacore::Vector<casacore::String> polarization (2);

    x[0] = 0;
    y[0] = .005;
    polarization [0] = "R";

    x[1] = 0;
    y[1] = .005;
    polarization [1] = "L";

    simulator_p->initFeeds ("", x, y, polarization);

    casacore::Int nSpectralWindows;
    casacore::Vector<casacore::String> spWindowNames;
    casacore::Vector<casacore::Int> nChannels;
    casacore::Vector<casacore::Quantity> startFrequencies;
    casacore::Vector<casacore::Quantity> frequencyDeltas;
    casacore::Vector<casacore::String> stokesString;

    ok = simulator_p->getSpWindows(nSpectralWindows,
                                   spWindowNames,
                                   nChannels,
                                   startFrequencies,
                                   frequencyDeltas,
                                   stokesString);

    if (! ok || nSpectralWindows == 0){

        addSpectralWindows (4);
    }

    if (timeStart_p < 0){
        timeStart_p = 0;
        timeInterval_p = 1;
        timeEnd_p = 15;
    }

}

inline void MsFactory::addField (const casacore::String & name,
                     const casacore::MDirection & direction)
{
    simulator_p->initFields (name, direction, "");
}

inline void MsFactory::addFeeds (casacore::Int nFeeds)
{
    casacore::Vector<casacore::Double> x (nFeeds, 0);
    casacore::Vector<casacore::String> polarization (nFeeds, "RR");
    simulator_p->initFeeds ("", x, x, polarization);
}

inline void MsFactory::addSpectralWindows (int nSpectralWindows)
{
    for (casacore::Int i = 0; i < nSpectralWindows; i++){

        casacore::String name = casacore::String::format ("sp%d", i);

        addSpectralWindow (name,
                           10 + i,
                           1e9 * (i + 1),
                           1e6 * (i + 1),
                           "LL RR RL LR");
    }
}

inline void MsFactory::addSpectralWindow (const casacore::String & name,
                              casacore::Int nChannels,
                              casacore::Double frequency,
                              casacore::Double frequencyDelta,
                              const casacore::String & stokes,
                              casacore::Double endingTime)
{
    simulator_p->initSpWindows (name,
                                nChannels,
                                casacore::Quantity (frequency, "Hz"),
                                casacore::Quantity (frequencyDelta, "Hz"),

                                casacore::MFrequency::TOPO,
                                stokes);

    // If so specified, mark that the ending time for this SPW is different
    if(endingTime != -1 )
    {
      // Get how many SPWs are present
      casacore::Int currentNSpectralWindows;
      casacore::Vector<casacore::String> currentSpWindowNames;
      casacore::Vector<casacore::Int> currentNChannels;
      casacore::Vector<casacore::Quantity> currentStartFrequencies;
      casacore::Vector<casacore::Quantity> currentFrequencyDeltas;
      casacore::Vector<casacore::String> currentStokesString;

      simulator_p->getSpWindows(currentNSpectralWindows,
                                currentSpWindowNames,
                                currentNChannels,
                                currentStartFrequencies,
                                currentFrequencyDeltas,
                                currentStokesString);

      // There are currentNSpectralWindows SPWs added so far.
      endingTimePerSpw_p[currentNSpectralWindows - 1] = endingTime;
    }
}

inline void MsFactory::addWeightSpectrum (casacore::Bool addIt)
{
    addWeightSpectrum_p = addIt;
}

inline void MsFactory::addFloatData (casacore::Bool addIt)
{
    addFloatData_p = addIt;
}

inline void MsFactory::addColumns ()
{
    if (addWeightSpectrum_p){
        addCubeColumn (casacore::MS::WEIGHT_SPECTRUM, "WeightSpectrumTiled");
    }
    if (addFloatData_p){
        addCubeColumn (casacore::MS::FLOAT_DATA, "FloatDataTiled");
    }
}

inline void MsFactory::addCubeColumn (casacore::MSMainEnums::PredefinedColumns columnId,
                          const casacore::String & dataStorageManagerName)
{
    casacore::TableDesc tableDescription;// = ms_p->actualTableDesc ();

    casacore::MS::addColumnToDesc (tableDescription, columnId, 2);
    casacore::IPosition tileShape (3, 4, 100, 100);

    TiledShapeStMan storageManager (dataStorageManagerName, tileShape);

    ms_p->addColumn (tableDescription, storageManager);

    ms_p->flush();
}

inline void MsFactory::removeColumn (casacore::MSMainEnums::PredefinedColumns columnId)
{
    casacore::String columnName = casacore::MS::columnName(columnId);

    ms_p->removeColumn(columnName);

    ms_p->flush();
}

inline void MsFactory::attachColumns ()
{
    const ColumnDescSet & cds = ms_p->tableDesc ().columnDescSet ();

    columns_p.antenna1_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::ANTENNA1));
    columns_p.antenna2_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::ANTENNA2));
    columns_p.array_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::ARRAY_ID));

    if (cds.isDefined ("CORRECTED_DATA")) {
        columns_p.corrVis_p.attach (* ms_p, "CORRECTED_DATA");
    }

    columns_p.dataDescriptionId_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::DATA_DESC_ID));
    columns_p.exposure_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::EXPOSURE));
    columns_p.field_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::FIELD_ID));
    columns_p.feed1_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::FEED1));
    columns_p.feed2_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::FEED2));
    columns_p.flag_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::FLAG));
    columns_p.flagCategory_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::FLAG_CATEGORY));
    columns_p.flagRow_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::FLAG_ROW));

    if (cds.isDefined (casacore::MS::columnName (casacore::MS::FLOAT_DATA))) {
        columns_p.floatVis_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::FLOAT_DATA));
        //floatDataFound_p = true;
    } else {
        //floatDataFound_p = false;
    }

    if (cds.isDefined ("MODEL_DATA")) {
        columns_p.modelVis_p.attach (* ms_p, "MODEL_DATA");
    }

    columns_p.observation_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::OBSERVATION_ID));
    columns_p.processor_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::PROCESSOR_ID));
    columns_p.scan_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::SCAN_NUMBER));
    columns_p.sigma_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::SIGMA));
    columns_p.state_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::STATE_ID));
    columns_p.time_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::TIME));
    columns_p.timeCentroid_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::TIME_CENTROID));
    columns_p.timeInterval_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::INTERVAL));
    columns_p.uvw_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::UVW));

    if (cds.isDefined (casacore::MS::columnName (casacore::MS::DATA))) {
        columns_p.vis_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::DATA));
    }

    columns_p.weight_p.attach (* ms_p, casacore::MS::columnName (casacore::MS::WEIGHT));

    if (cds.isDefined ("WEIGHT_SPECTRUM")) {
        columns_p.weightSpectrum_p.attach (* ms_p, "WEIGHT_SPECTRUM");
    }

    if (cds.isDefined ("CORRECTED_WEIGHT_SPECTRUM")) {
        columns_p.weightSpectrumCorrected_p.attach (* ms_p, "CORRECTED_WEIGHT_SPECTRUM");
    }

}

inline pair<casacore::MeasurementSet *, rownr_t> MsFactory::createMs ()
{
    addColumns ();

    fillData ();

    ms_p->flush();

    pair<casacore::MeasurementSet *, rownr_t> result = make_pair (ms_p, nRows_p);

    ms_p = 0; // give up all ownership and access

    return result;
}

inline void MsFactory::fillData ()
{
    addDefaults ();

    attachColumns ();

    FillState fillState;

    fillState.rowNumber_p = 0;
    fillState.nAntennas_p = nAntennas_p;
    fillState.nFlagCategories_p = 3;
    fillState.timeDelta_p = timeInterval_p;

    casacore::Double time = timeStart_p;

    casacore::Int nSpectralWindows;
    casacore::Vector<casacore::String> spWindowNames;
    casacore::Vector<casacore::Int> nChannels;
    casacore::Vector<casacore::Quantity> startFrequencies;
    casacore::Vector<casacore::Quantity> frequencyDeltas;
    casacore::Vector<casacore::String> stokesString;

    simulator_p->getSpWindows(nSpectralWindows,
                              spWindowNames,
                              nChannels,
                              startFrequencies,
                              frequencyDeltas,
                              stokesString);

    while (time < timeEnd_p){
        fillState.time_p = time;

        for (casacore::Int j = 0; j < nSpectralWindows; j++){
            // Only add this time to this SPW if either this SPW does not have
            // a custom end time or if the timestamp is less than this
            // SPW custom end time 
            bool addThisTime = (endingTimePerSpw_p.find(j) == endingTimePerSpw_p.end()) ? 
                true : time < endingTimePerSpw_p[j];
            if(addThisTime)
            {
                fillState.spectralWindow_p = j;
                fillState.nChannels_p = nChannels [j];
                vector<casacore::String> stokesComponents = utilj::split (stokesString [j], " ", true);
                fillState.nCorrelations_p = stokesComponents.size();
                fillRows (fillState);
            }
        }

        time += timeInterval_p;
    }

    printf ("\n---Total of %d rows filled\n", fillState.rowNumber_p);
    printf ("---Time range %f to %f\n", timeStart_p, timeEnd_p);

    nRows_p = fillState.rowNumber_p;
}

inline void MsFactory::fillRows (FillState & fillState)
{
    // Extend the MS to have one row per every unique pairing
    // of antennas.

    casacore::Int n = fillState.nAntennas_p + (includeAutocorrelations_p ? 1 : 0);
    casacore::Int nNewRows = (n * (n - 1)) / 2;
    ms_p->addRow(nNewRows);

    // Fill in a row for every unique antenna pairing.


    for (casacore::Int a1 = 0; a1 < fillState.nAntennas_p; a1 ++){

        fillState.antenna1_p = a1;

        casacore::Int firstA2 = includeAutocorrelations_p ? a1 : a1 + 1;

        for (casacore::Int a2 = firstA2; a2 < fillState.nAntennas_p; a2 ++){

            fillState.antenna2_p = a2;

            fillCubes (fillState);

            fillFlagCategories (fillState);

            fillCollections (fillState);

            fillScalars (fillState);

            fillState.rowNumber_p ++;
        }
    }
}

template <typename T>
void MsFactory::fillCube (ArrayColumn<T> & column, const FillState & fillState,
                     const GeneratorBase * generatorBase)
{

    const Generator<T> * generator = dynamic_cast <const Generator<T> *> (generatorBase);

    ThrowIf (generator == 0, "Bad return type on generator");

    casacore::Matrix <T> cell (casacore::IPosition (2, fillState.nCorrelations_p, fillState.nChannels_p));

    for (casacore::Int channel = 0; channel < fillState.nChannels_p; channel ++){

        for (casacore::Int correlation = 0;
             correlation < fillState.nCorrelations_p;
             correlation ++){

            cell (correlation, channel) = (* generator) (fillState, channel, correlation);
        }
    }

    column.put (fillState.rowNumber_p, cell);
}

template <typename T>
void MsFactory::fillScalar (ScalarColumn<T> & column, const FillState & fillState,
                       const GeneratorBase * generatorBase)
{
    const Generator<T> * generator = dynamic_cast <const Generator<T> *> (generatorBase);

    ThrowIf (generator == 0, "Bad return type on generator");

    ThrowIf (column.isNull(),
             casacore::String::format ("Column is not attached (%s)", typeid (column).name()));

    T value = (* generator) (fillState, -1, -1);

    column.put (fillState.rowNumber_p, value);
}


inline void MsFactory::fillCubes (FillState & fillState)
{
    fillVisCubeCorrected (fillState);
    fillVisCubeModel (fillState);
    fillVisCubeObserved (fillState);
    fillVisFloatData (fillState);

    fillWeightSpectrumCube (fillState);
    fillWeightSpectrumCorrectedCube (fillState);
    fillFlagCube (fillState);
}

inline void MsFactory::fillFlagCategories (const FillState & fillState)
{

    const GenerateFlagCategory * generator =
            dynamic_cast <const GenerateFlagCategory *> (generators_p.get(casacore::MSMainEnums::FLAG_CATEGORY));

    ThrowIf (generator == 0, "Bad return type on generator");

    Cube <casacore::Bool> cell (casacore::IPosition (3, fillState.nCorrelations_p, fillState.nChannels_p,
                                  fillState.nFlagCategories_p));

    for (casacore::Int channel = 0; channel < fillState.nChannels_p; channel ++){

        for (casacore::Int correlation = 0;
             correlation < fillState.nCorrelations_p;
             correlation ++){

            for (casacore::Int category = 0;
                 category < fillState.nFlagCategories_p;
                 category ++){

                casacore::Bool value = (* generator) (fillState, channel, correlation, category);

                cell (correlation, channel, category) = value;
            }
        }
    }

    columns_p.flagCategory_p.put (fillState.rowNumber_p, cell);
}


inline void MsFactory::fillUvw (FillState & fillState)
{
    const GeneratorBase * generatorBase = generators_p.get (casacore::MSMainEnums::UVW);
    const Generator<casacore::Vector <casacore::Double> > * generator =
            dynamic_cast<const Generator<casacore::Vector <casacore::Double> > *> (generatorBase);

    casacore::Vector<casacore::Double> uvw (3);
    uvw = (* generator) (fillState, -1, -1);

    columns_p.uvw_p.put (fillState.rowNumber_p, uvw);
}


inline void MsFactory::fillVisCubeCorrected (FillState & fillState)
{
    if (! columns_p.corrVis_p.isNull ()){

        fillCube (columns_p.corrVis_p, fillState, generators_p.get(casacore::MSMainEnums::CORRECTED_DATA));
    }
}

inline void MsFactory::fillWeightSpectrumCube (FillState & fillState)
{
    if (! columns_p.weightSpectrum_p.isNull()){

        fillCube (columns_p.weightSpectrum_p, fillState, generators_p.get(casacore::MSMainEnums::WEIGHT_SPECTRUM));
    }
}

inline void MsFactory::fillWeightSpectrumCorrectedCube (FillState & fillState)
{
    if (! columns_p.weightSpectrumCorrected_p.isNull()){

        fillCube (columns_p.weightSpectrumCorrected_p, fillState,
                  generators_p.get(casacore::MSMainEnums::CORRECTED_WEIGHT_SPECTRUM));
    }
}

inline void MsFactory::fillFlagCube (FillState & fillState)
{
    fillCube (columns_p.flag_p, fillState, generators_p.get(casacore::MSMainEnums::FLAG));
}


inline void MsFactory::fillVisCubeModel (FillState & fillState)
{
    if (! columns_p.modelVis_p.isNull ()){

        fillCube (columns_p.modelVis_p, fillState, generators_p.get(casacore::MSMainEnums::MODEL_DATA));
    }
}

inline void MsFactory::fillVisFloatData (FillState & fillState)
{
    if (! columns_p.floatVis_p.isNull ()){
        fillCube (columns_p.floatVis_p, fillState, generators_p.get(casacore::MSMainEnums::FLOAT_DATA));
    }
}

inline void MsFactory::fillVisCubeObserved (FillState & fillState)
{
  if (! columns_p.vis_p.isNull ()){
    fillCube (columns_p.vis_p, fillState, generators_p.get (casacore::MSMainEnums::DATA));
  }
}

inline void MsFactory::fillCollections (FillState & fillState)
{
    fillWeight (fillState);
    fillSigma (fillState);
    fillUvw (fillState);
}

inline void MsFactory::fillScalars (FillState & fillState)
{
    fillScalar (columns_p.antenna1_p, fillState, generators_p.get (casacore::MSMainEnums::ANTENNA1));
    fillScalar (columns_p.antenna2_p, fillState, generators_p.get (casacore::MSMainEnums::ANTENNA2));
    fillScalar (columns_p.array_p, fillState, generators_p.get (casacore::MSMainEnums::ARRAY_ID));
    fillScalar (columns_p.dataDescriptionId_p, fillState, generators_p.get (casacore::MSMainEnums::DATA_DESC_ID));
    fillScalar (columns_p.exposure_p, fillState, generators_p.get (casacore::MSMainEnums::EXPOSURE));
    fillScalar (columns_p.feed1_p, fillState, generators_p.get (casacore::MSMainEnums::FEED1));
    fillScalar (columns_p.feed2_p, fillState, generators_p.get (casacore::MSMainEnums::FEED2));
    fillScalar (columns_p.field_p, fillState, generators_p.get (casacore::MSMainEnums::FIELD_ID));
    fillScalar (columns_p.flagRow_p, fillState, generators_p.get (casacore::MSMainEnums::FLAG_ROW));
    fillScalar (columns_p.observation_p, fillState, generators_p.get (casacore::MSMainEnums::OBSERVATION_ID));
    fillScalar (columns_p.processor_p, fillState, generators_p.get (casacore::MSMainEnums::PROCESSOR_ID));
    fillScalar (columns_p.scan_p, fillState, generators_p.get (casacore::MSMainEnums::SCAN_NUMBER));
    fillScalar (columns_p.state_p, fillState, generators_p.get (casacore::MSMainEnums::STATE_ID));
    fillScalar (columns_p.timeCentroid_p, fillState, generators_p.get (casacore::MSMainEnums::TIME_CENTROID));
    fillScalar (columns_p.timeInterval_p, fillState, generators_p.get (casacore::MSMainEnums::INTERVAL));
    fillScalar (columns_p.time_p, fillState, generators_p.get (casacore::MSMainEnums::TIME));
}

inline void MsFactory::fillSigma (FillState & fillState)
{
    casacore::Vector<Float> sigmas (fillState.nCorrelations_p);
    const Generator<Float> * generator =
            dynamic_cast <const Generator<Float> *> (generators_p.get (casacore::MSMainEnums::SIGMA));

    for (casacore::Int i = 0; i < fillState.nCorrelations_p; i ++){
        sigmas (i) = (* generator) (fillState, -1, i);
    }


    columns_p.sigma_p.put (fillState.rowNumber_p, sigmas);

}

inline void MsFactory::fillWeight (FillState & fillState)
{
    casacore::Vector<Float> weights (fillState.nCorrelations_p);
    const GeneratorBase * generatorBase = generators_p.get (casacore::MSMainEnums::WEIGHT);
    const Generator<Float> * generator =
        dynamic_cast <const Generator<Float> *> (generatorBase);

    for (casacore::Int i = 0; i < fillState.nCorrelations_p; i ++){
        weights (i) = (* generator) (fillState, -1, i);
    }

    columns_p.weight_p.put (fillState.rowNumber_p, weights);
}

inline void MsFactory::setDataGenerator (casacore::MSMainEnums::PredefinedColumns column, GeneratorBase * generator)
{
    generators_p.set (column, generator);
}

inline void MsFactory::setIncludeAutocorrelations (casacore::Bool includeThem)
{
    includeAutocorrelations_p = includeThem;
}

inline void MsFactory::setTimeInfo (casacore::Double startingTime, casacore::Double endingTime, casacore::Double interval)
{
    timeStart_p = startingTime;
    timeEnd_p = endingTime;
    timeInterval_p = interval;
}


} // end namespace test

} // end namespace vi

} // end namespace casa

#endif // ! defined (CASA_VI_TEST_MsFactory_H_130114_1357)
