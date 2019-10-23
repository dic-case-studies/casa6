#include <iostream>
#include <iomanip>
#include <sstream>
#include <bitset>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>
#include <iterator>
#include <utility>
#include <stdlib.h>

using namespace std;

#include <alma/ASDM/ASDMAll.h>

using namespace asdm;

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayUtil.h>
#include <casa/Arrays/ArrayLogical.h>

#include <alma/Enumerations/CAxisName.h>
using namespace AxisNameMod;
#include <alma/Enumerations/CProcessorType.h>
using namespace ProcessorTypeMod;

#include <alma/ASDMBinaries/SDMDataObject.h>
#include <alma/ASDMBinaries/SDMDataObjectReader.h>
#include <alma/ASDMBinaries/SDMDataObjectStreamReader.h>

using namespace asdmbinaries;
using namespace casacore;

#include "asdm2MSGeneric.h"

#include <ScansParser.h>

#include <alma/Enumerations/CScanIntent.h>
using namespace ScanIntentMod;
#include <alma/Enumerations/CSubscanIntent.h>
using namespace SubscanIntentMod;

#include <casacore/casa/System/AppState.h>

static bool ddebug = false;

/*
** A simplistic tracing toolbox.
*/
static bool debug = (getenv("ASDM_DEBUG") != NULL);

static vector<char> logIndent;
#define LOGENTER(name) if (debug) { std::for_each(logIndent.begin(), logIndent.end(), [](char v) { cout << v; }); logIndent.push_back('\t'); cout << #name ": entering" << endl; }
#define LOGEXIT(name)  if (debug) { logIndent.pop_back(); std::for_each(logIndent.begin(), logIndent.end(), [](char v) { cout << v; } ); cout << #name ": exiting" << endl; }
#define LOG(msg) if (debug) { std::for_each(logIndent.begin(), logIndent.end(), [](char v) { cout << v; } ); cout << msg << endl; }

/*
** Two string streams used all over the applications in general to prepare messages sent to the logging facilities.
*/
#include <casa/Logging/StreamLogSink.h>
#include <casa/Logging/LogSink.h>
#include <sdm_cmpt.h>

namespace casac {

    /*
    ** A class with static variables and methods to inform about the presence of axes in the description of attachment (attribute axes)
    ** It's meant to be used as kind of global datatructure which can be refered to from any point of the code (maybe a singleton would have been more elegant).
    **
    ** 1) Initialize the dataStructure with a call to the set method
    ** 2) Then query the dataStructure with hasXXX() calls from any location in the code.
    */
    class CorrelatorFlagsAxes {
    public:
        static void set(const std::vector<AxisName>& axes);
        static bool hasBAL();
        static void hasBAL (bool b);
        static bool hasANT();
        static void hasANT (bool b);
        static bool hasBAB();
        static void hasBAB(bool b);
        static bool hasSPW();
        static void hasSPW(bool b);
        static bool hasPOL();
        static void hasPOL(bool b);

        static void clear();

    private:
        static bool hasBAL_;
        static bool hasANT_;
        static bool hasBAB_;
        static bool hasSPW_;
        static bool hasPOL_;

        CorrelatorFlagsAxes() {;}
    };

    bool CorrelatorFlagsAxes::hasBAL_ = false;
    bool CorrelatorFlagsAxes::hasANT_ = false;
    bool CorrelatorFlagsAxes::hasBAB_ = false;
    bool CorrelatorFlagsAxes::hasSPW_ = false;
    bool CorrelatorFlagsAxes::hasPOL_ = false;

    void CorrelatorFlagsAxes::set(const vector<AxisName>& axes) {
        CorrelatorFlagsAxes::clear();
        for (unsigned int i = 0; i < axes.size(); i++) {
            switch (axes[i]) {
                case AxisName::BAL : hasBAB_ = true; break;
                case AxisName::ANT : hasANT_ = true; break;
                case AxisName::BAB : hasBAB_ = true; break;
                case AxisName::SPW : hasSPW_ = true; break;
                case AxisName::POL : hasPOL_ = true; break;
                default: break;
            }
        }
    }

    bool CorrelatorFlagsAxes::hasBAL() {return hasBAL_;}
    void CorrelatorFlagsAxes::hasBAL(bool b) {hasBAL_=b;}

    bool CorrelatorFlagsAxes::hasANT() {return hasANT_;}
    void CorrelatorFlagsAxes::hasANT(bool b) {hasANT_=b;}

    bool CorrelatorFlagsAxes::hasBAB() {return hasBAB_;}
    void CorrelatorFlagsAxes::hasBAB(bool b) {hasBAB_=b;}

    bool CorrelatorFlagsAxes::hasSPW() {return hasSPW_;}
    void CorrelatorFlagsAxes::hasSPW(bool b) {hasSPW_=b;}

    bool CorrelatorFlagsAxes::hasPOL() {return hasPOL_;}
    void CorrelatorFlagsAxes::hasPOL(bool b) {hasPOL_=b;}

    void CorrelatorFlagsAxes::clear() {
        hasBAL_ = false;
        hasANT_ = false;
        hasBAB_ = false;
        hasSPW_ = false;
        hasPOL_ = false;
    }
    /*
    ** A class to encode the exceptions occuring during the usage of one instance of class ProcessFlags.
    */
    class ProcessFlagsException {
    private :
        string message;

    public :
        ProcessFlagsException ();
        ProcessFlagsException (const string & message ) ;
        virtual ~ProcessFlagsException ();
        const string &getMessage();
    };

    ProcessFlagsException::ProcessFlagsException () : message("ProcessFlagsException") {;}

    ProcessFlagsException::ProcessFlagsException (const string & message) : message("ProcessFlagsException: " + message) {;}

    ProcessFlagsException::~ProcessFlagsException () {;}

    const string& ProcessFlagsException::getMessage() {
        return message;
    }
    // end of class ProcessFlagsException

    /*
    ** A class to perform linearized index computations between a matrix lower left and upper right corners. This class was designed during the work on http://jira.alma.cl/browse/PRTSIR-9678
    **
    */
    class Transposer {
    private:
        unsigned int n_;
        vector<unsigned int> ur2llTransposed_v_;

        unsigned int nnm1o2_() const;
        Transposer() {;};

    public :
        Transposer(unsigned int n);
        unsigned int transposed(unsigned int n) const;
    };

    unsigned int Transposer::nnm1o2_() const {
        return n_*(n_-1)/2;
    }

    Transposer::Transposer(unsigned int n) {
        if (n == 0 )
            throw ProcessFlagsException("Internal error. The Tranposer constructor was called with n == 0.");
        n_ = n;
        vector<unsigned int> v_v(n_);
        vector<vector<unsigned int> > m_vv(n_, v_v);

        ur2llTransposed_v_.resize(nnm1o2_(), 0);

        unsigned int k = 0;

        // Lower left traversal and index linearization
        k = 0;
        for (unsigned int i = 0; i < n_-1; i++)
            for (unsigned int j = i+1; j < n_; j++)
                m_vv[i][j] = k++;

        // Upper right traversal and linearized index transposition
        k = 0;
        for (unsigned int i = 1; i < n_; i++)
            for (unsigned int j = 0; j < i ; j++)
                ur2llTransposed_v_[k++]=m_vv[j][i];
    }

    unsigned int Transposer::transposed(unsigned int n) const {
        if ( n >= nnm1o2_() ) {
            ostringstream oss;
            oss << "Internal error.The method Transposer::transposed was called with an invalid value '"
                << n << "'. The maximum allowed is '" << nnm1o2_() << "==" << n_ << "*" << n_-1 << "/2'";
            throw ProcessFlagsException(oss.str());
        }
        else
            return  ur2llTransposed_v_[n];
    }


    /*
    ** Some typedefs for quite complex data structures used in many places
    */
    typedef pair <unsigned int, unsigned int> FLAG_SHAPE;
    typedef vector<char> FLAG_V;
    typedef pair<FLAG_SHAPE, FLAG_V> FLAG_CELL;

    /*
    ** A generic class to implement a machine "consuming" elements of the "flags" attachments in the BDF parts.
    ** The parameter T is meant to handle the base type of the flags.
    */
    template <typename T >
    class BDFFlagConsumer {
    private:
        const T* items_p;
        unsigned int numItems;
        unsigned int i0, i1;
        BDFFlagConsumer(){;}

    public:
        BDFFlagConsumer(const T* items_p, unsigned int numItems) : items_p(items_p), numItems(numItems) {
            LOGENTER("BDFFlagConsumer::BDFFlagConsumer(const T* items_p, unsigned int numItems)");
            i0 = 0; i1 = 1;
            LOGEXIT("BDFFlagConsumer::BDFFlagConsumer(const T* items_p, unsigned int numItems)");
        }

        pair<const T* , const T*>  consume(unsigned int n) {
            LOGENTER("BDFFlagConsumer::pair<const T* , const T*>  consume(unsigned int n)");
            if (i0 >= numItems) return make_pair<const T*, const T*>(NULL, NULL);
            i1 = i0 + n;
            if (i1 > numItems)
                return make_pair<const T*, const T*>(NULL, NULL);

            const T* p0 = items_p + i0;
            const T* p1 = items_p + i1;
            i0 = i1;
            LOGEXIT("BDFFlagsConsumer::pair<const T* , const T*>  consume(unsigned int n)");
            return make_pair(p0, p1);
        }
    }; // end of class BDFFlagConsumer.

    /*
    ** A class whose instances are meant to be used as functors returning the evaluation of
    ** flags to be written into the MS given the mask flag provided on the command line.
    **
    */
    class MSFlagEval {
    private:
        MSFlagEval() ;
        FLAGSTYPE mask;
        FLAGSTYPE ignore;

    public:
        MSFlagEval(FLAGSTYPE mask, FLAGSTYPE ignore = 0 ); // 0 plays the role of a kind of neutral element as a matter of flagging.
        ~MSFlagEval();
        char operator()(FLAGSTYPE);
        MSFlagEval& operator=(const MSFlagEval& other);
    };

    MSFlagEval::MSFlagEval(FLAGSTYPE mask, FLAGSTYPE ignore):mask(mask),ignore(ignore) {;}
    MSFlagEval::~MSFlagEval() {;}
    //char MSFlagEval::operator()(FLAGSTYPE f) { return ((f == ignore) || (f & mask & 0x0000FFFF) == 0) ? 0 : 1; }
    char MSFlagEval::operator()(FLAGSTYPE f) { return ((f == ignore) || ((f & mask) == 0)) ? 0 : 1; } // As of 2015.8 correlator
    // uses the definitions for binary
    // contained in BinaryFlags
    // see https://bugs.nrao.edu/browse/CAS-8491
    // There is no reason to continue to
    // ignore the 16 upper bits.
    MSFlagEval & MSFlagEval::operator=(const MSFlagEval& other) {
        if ( this != &other ) {
            mask = other.mask;
        }
        return *this;
    }

    // end of class MSFlagEval

    /*
    ** A class to define a machine accumulating MS Flags cells and their descriptors as the BDF Flags are traversed
    ** the BDFlag consumer.
    */
    template <typename T>
    class MSFlagAccumulator {
    private:
        unsigned int numIntegration_;
        unsigned int numBAL_;
        unsigned int numDD_;

        unsigned int currIntegration_;
        unsigned int currBAL_;
        unsigned int currDD_;

        vector < vector < vector < FLAG_CELL > > > flagCell_vvv;  // The 1st axis ( left ) is TIME
        // The 2nd axis ( middle ) is BAL or ANT
        // The 3rd axis ( right ) is DD

        vector<T> MSFlags_v;

        vector< FLAG_SHAPE * > flagShapes_p_v;
        vector< FLAG_V * > flagValues_p_v;

        unsigned int numFlaggedRows;
        MSFlagAccumulator();

    public:
        MSFlagAccumulator(unsigned int numIntegration,
                          unsigned int numBAL,
                          unsigned int numDD) {
            LOGENTER("MSFlagAccumulator::MSFlagAccumulator()");
            this->numBAL_ = numBAL;
            this->numIntegration_ = numIntegration;
            this->numDD_ = numDD;
            vector<FLAG_CELL> v(numDD_);
            vector <vector<FLAG_CELL> > vv(numBAL_, v);
            flagCell_vvv.assign(numIntegration_, vv);

            currIntegration_ = 0;
            currBAL_ = 0;
            currDD_ = 0;

            numFlaggedRows = 0;
            LOGEXIT("MSFlagAccumulator::MSFlagAccumulator()");
        }

        virtual ~MSFlagAccumulator() {;}
        void accumulate(unsigned int numChan, unsigned int numPol, const vector<T>& values) {
            LOGENTER("MSFlagAccumulator::accumulate(unsigned int numChan, unsigned int numPol, T* values)");
            vector<T> values_ = values;
            flagCell_vvv[currIntegration_][currBAL_][currDD_] = make_pair(make_pair(numChan, numPol), values_);
            LOGEXIT("MSFlagAccumulator::accumulate(unsigned int numChan, unsigned int numPol, T* values)");
        }

        const vector < vector < vector < FLAG_CELL > > > &flagCell() {  // The 1st axis ( left ) is TIME
            // The 2nd axis ( middle ) is BAL or ANT
            // The 3rd axis ( right ) is DD
            return flagCell_vvv;
        }

        unsigned int numIntegration() const {
            return numIntegration_;
        }

        unsigned int numBAL() const {
            return numBAL_;
        }

        unsigned int numDD() const {
            return numDD_;
        }

        void nextBAL() {
            LOGENTER("MSFlagAccumulator::nextBAL()");
            currBAL_++;
            LOGEXIT("MSFlagAccumulator::nextBAL()");
        }

        void resetBAL() {
            LOGENTER("MSFlagAccumulator::resetBAL()");
            currBAL_=0;
            LOGEXIT("MSFlagAccumulator::resetBAL()");
        }

        void nextDD() {
            LOGENTER("MSFlagAccumulator::nextDD()");
            currDD_++;
            LOGEXIT("MSFlagAccumulator::nextDD()");
        }

        void resetDD() {
            LOGENTER("MSFlagAccumulator::resetDD()");
            currDD_=0;
            LOGEXIT("MSFlagAccumulator::resetDD()");
        }

        void nextIntegration() {
            LOGENTER("MSFlagAccumulator::nextIntegration()");
            currIntegration_++;
            LOGEXIT("MSFlagAccumulator::nextIntegration()");
        }

        void resetIntegration() {
            LOGENTER("MSFlagAccumulator::resetIntegration()");
            currIntegration_=0;
            LOGEXIT("MSFlagAccumulator::resetIntegration()");
        }

        pair< vector<FLAG_SHAPE *>*,
              vector<FLAG_V *>* >  orderedByDDTIMBAL(bool skipFirstIntegration=false) {
            LOGENTER("MSFlagAccumulator::orderedByDDTIMBAL(bool MSORDER=true)");

            unsigned int numIntUsed = skipFirstIntegration ? (numIntegration_-1) : numIntegration_;
            flagShapes_p_v.clear();
            flagShapes_p_v.resize(numIntUsed*numBAL_*numDD_);
            flagValues_p_v.clear();
            flagValues_p_v.resize(numIntUsed*numBAL_*numDD_);
    
            unsigned int k = 0;
            for (unsigned int iDD = 0; iDD < numDD_; iDD++) {
                for (unsigned int iIntegration = 0; iIntegration < numIntegration_; iIntegration++) {
                    if (skipFirstIntegration && (iIntegration==0)) continue;
	
                    for (unsigned int iBAL = 0; iBAL < numBAL_; iBAL++) {
                        flagShapes_p_v[k] = &(flagCell_vvv[iIntegration][iBAL][iDD].first);
                        flagValues_p_v[k] = &(flagCell_vvv[iIntegration][iBAL][iDD].second);
                        k++;
                    }
                }
            }
    
            LOGEXIT("MSFlagAccumulator::orderedByDDTIMBAL(bool MSORDER=true)");
            return make_pair< vector<FLAG_SHAPE *> *, vector<FLAG_V * > *> (&flagShapes_p_v, &flagValues_p_v);
        }

        void info (ostream& os) {
            LOGENTER("MSFlagAccumulator::info");
            os << "I have " << numIntegration_ * numBAL_ * numDD_  << " MS flag cells." << endl;
            LOGEXIT("MSFlagAccumulator::info");
        }

        void dump (ostream& os, bool MSORDER=true) {
            LOGENTER("MSFlagAccumulator::dump");

            pair< vector<FLAG_SHAPE *>*, vector<FLAG_V *>* > cds = orderedByDDTIMBAL(MSORDER);

            vector<FLAG_SHAPE *>* shapes_p_v_p = cds.first;
//          vector<FLAG_V*>* values_p_v_p = cds.second;
            unsigned int numCell = shapes_p_v_p->size();
            unsigned rn = 0;
            for (unsigned int i = 0; i < numCell; i++) {
                os << rn++ << " - shape=[" << shapes_p_v_p->at(i)->first << "," << shapes_p_v_p->at(i)->second << "]" << endl;
            }
            LOGEXIT("MSFlagAccumulator::dump");
        }
    }; // end of class MSFlagAccumulator


    /*
    ** The following functions drive the processing of the BDF flags
    */
    static void traverseBAB( bool sameAntenna, const vector<SDMDataObject::Baseband> &basebands,
                             const pair<unsigned int, const FLAGSTYPE*> &flagsPair,
                             BDFFlagConsumer<FLAGSTYPE> &consumer, MSFlagEval &flagEval,
                             MSFlagAccumulator<char> &accumulator ) {

        LOGENTER("traverseBAB");


        unsigned int		 numFlags = flagsPair.first;

        accumulator.resetDD();
        for (SDMDataObject::Baseband bab: basebands) {
            const vector<SDMDataObject::SpectralWindow>& spws = bab.spectralWindows();
            bool firstSPW = true;
            vector<StokesParameterMod::StokesParameter>::const_iterator ppIter, ppBegin, ppEnd;
            unsigned int numPolProducts;
            pair<const FLAGSTYPE*, const FLAGSTYPE*> range;

            for(SDMDataObject::SpectralWindow spw: spws) {
                if ((firstSPW && !CorrelatorFlagsAxes::hasSPW()) || CorrelatorFlagsAxes::hasSPW()) {
                    /*
                     * Consider the spectral window if 
                     **  (it's the first spw in the vector and SPW does not belong to the list of axes for flags) or (if SPW belongs to the list of axes for flags).
                     */
                    ppIter = ppBegin = sameAntenna ? spw.sdPolProducts().begin() : spw.crossPolProducts().begin();
                    ppEnd = sameAntenna ? spw.sdPolProducts().end() : spw.crossPolProducts().end();
                    numPolProducts = ppEnd - ppBegin;
                }

                unsigned int flagsCellNumPolProducts =  (numPolProducts == 3 && sameAntenna) ? 4 : numPolProducts;

                vector<char>  MSFlagsCell(spw.numSpectralPoint()*flagsCellNumPolProducts, (char) 0);
                if (numFlags) {
                    if ((firstSPW && !CorrelatorFlagsAxes::hasSPW()) || CorrelatorFlagsAxes::hasSPW()) {
                        /*
                        ** Consume flags if and only if (it's the first spw in the vector and SPW does not belong to the list of axes for flags) or (if SPW belongs to the list of axes for flags).
                        */
                        range  = consumer.consume(numPolProducts);
                    }
                    unsigned int k = 0;
                    for (unsigned int i = 0; i < spw.numSpectralPoint(); i++){
                        unsigned int kk = 0;
                        for (const FLAGSTYPE* p = range.first; p != range.second; p++){
                            MSFlagsCell[k] = flagEval(*p) ; k++;
                            if (ddebug && *p) cout << *p << ":" << (short) flagEval(*p) << " ";
                            kk++;
                            if ((flagsCellNumPolProducts == 4) && sameAntenna && ( kk == 1 )) { // If we are in a case like XX XY YX, don't forget to repeat the value of index 1 .
                                MSFlagsCell[k] =  MSFlagsCell[k-1]; k++;
                            }
                        }
                    }
                    if (ddebug) cout << endl;
                }
                firstSPW = false;
                accumulator.accumulate(spw.numSpectralPoint(), flagsCellNumPolProducts, MSFlagsCell);
                accumulator.nextDD();
            }
        }
        LOGEXIT("traverseBAB");
    }

    static void traverseANT( const vector<SDMDataObject::Baseband> &basebands,
                             const vector<string> &antennas,
                             const pair<unsigned int, const FLAGSTYPE*> &flagsPair,
                             BDFFlagConsumer<FLAGSTYPE> &consumer, MSFlagEval &flagEval,
                             MSFlagAccumulator<char> &accumulator) {

        LOGENTER("traverseANT");
        accumulator.resetBAL();
        for (unsigned int i = 0; i < antennas.size() ; i++) {
            traverseBAB(true, basebands, flagsPair, consumer, flagEval, accumulator);
            accumulator.nextBAL();
        }
        LOGEXIT("traverseANT");
    }

    static void traverseBAL( const vector<SDMDataObject::Baseband>& basebands,
                             const vector<string> &antennas,
                             const pair<unsigned int, const FLAGSTYPE*> &flagsPair,
                             BDFFlagConsumer<FLAGSTYPE> &consumer,
                             MSFlagEval &flagEval, MSFlagAccumulator<char> &accumulator ) {

        LOGENTER("traverseBAL");
        accumulator.resetBAL();
        for (unsigned int i = 1; i < antennas.size(); i++)
            //for (unsigned int j = i+1; j < antennas.size(); j++) {
            for (unsigned int j = 0; j < i; j++) {
                traverseBAB(false, basebands, flagsPair, consumer, flagEval, accumulator);
                accumulator.nextBAL();
            }
        LOGEXIT("traverseBAL");
    }

    static void traverseALMACorrelatorFlagsAxes( const vector<SDMDataObject::Baseband> &basebands,
                                                 const vector<string> &antennas,
                                                 const pair<unsigned int, const FLAGSTYPE*> &flagsPair,
                                                 MSFlagEval &flagEval, 
                                                 CorrelationModeMod::CorrelationMode correlationMode,
                                                 MSFlagAccumulator<char> &autoAccumulator,
                                                 MSFlagAccumulator<char> &crossAccumulator ) {

        LOGENTER("traverseALMACorrelatorFlagsAxes");

        const FLAGSTYPE*	flags_p	 = flagsPair.second;
        unsigned int		numFlags = flagsPair.first;
        BDFFlagConsumer<FLAGSTYPE> consumer(flags_p, numFlags);

        // Attention the next two calls must be done in *this* order (BAL then ANT) !!!!
        if (correlationMode != CorrelationModeMod::AUTO_ONLY)
            traverseBAL(basebands, antennas,  flagsPair, consumer, flagEval, crossAccumulator);
        if (correlationMode != CorrelationModeMod::CROSS_ONLY)
            traverseANT(basebands, antennas, flagsPair, consumer, flagEval, autoAccumulator);

        LOGEXIT("traverseALMACorrelatorFlagsAxes");
    }

    static void traverseALMARadiometerFlagsAxes( unsigned int numTime,
                                                 const vector<SDMDataObject::Baseband> &basebands,
                                                 const vector<string> &antennas,
                                                 const pair<unsigned int, const FLAGSTYPE*> &flagsPair,
                                                 MSFlagEval &flagEval, MSFlagAccumulator<char> &accumulator) {

        LOGENTER("traverseALMARadiometerFlagsAxes");

        const FLAGSTYPE*	flags_p	 = flagsPair.second;
        unsigned int		numFlags = flagsPair.first;
        BDFFlagConsumer<FLAGSTYPE> consumer(flags_p, numFlags);

        // accumulator.resetIntegration(); This call has been transfered to the caller 01/02/2016 Michel Caillat
        for (unsigned int iTime = 0; iTime < numTime; iTime++) {
            traverseANT(basebands, antennas, flagsPair, consumer, flagEval, accumulator);
            accumulator.nextIntegration();
        }
        LOGEXIT("traverseALMARadiometerFlagsAxes");
    }

    static bool isNotNull(char f){
        return f != 0;
    }

    /**
     *
     * Populates one cell of the columns FLAG and FLAG_ROW from the content of what's been found in the BDF flags attachment.
     *
     * (function used for Correlator data)
     *
     * @parameter flagsShape_p a pointer to the shape of the flag in the MS Main row number iRow0.
     * @parameter flag_v_p a pointer to a vector of flags as a flattened version of the 2D natural represenation of flags.
     * @parameter iRow0 the row number in the MS Main table
     * @parameter flag the column FLAG in the MS Main table
     * @parameter flagRow the column FLAG_ROW in the MS Main table
     */
    static bool putCell( FLAG_SHAPE* flagShape_p, FLAG_V* flag_v_p, uInt iRow0,
                         ArrayColumn<Bool> &flag, ScalarColumn<Bool> flagRow ) {

        LOGENTER("putCell");
        uInt numChan = flagShape_p->first;
        uInt numCorr = flagShape_p->second;

        bool cellFlagged = false;
        bool flagged = false;

        Matrix<Bool> flagCell(IPosition(2, numCorr, numChan));
        if (debug) {
            cout << "irow0 = " << iRow0 << endl;
            cout << "expecting a cell of shape numCorr=" << numCorr << ", numChan=" << numChan << endl;
            cout << "actual shape is " << flag.shape(iRow0) << endl;
        }
        flag.get((uInt)iRow0, flagCell);
        bool allSet = true;
        //unsigned int notNullBDF =  count_if(flag_v_p->begin(), flag_v_p->end(), isNotNull);
        //if (notNullBDF) cout << "row " << iRow0 << " - putCell in front of " << flag_v_p->size() << " elements of which " << notNullBDF << " are not null";
        int notNull = 0;
        int k = 0;
        if (ddebug) cout << "Put into MS :" ;
        for (uInt i = 0;  i < numChan; i++) {
            for (uInt j = 0; j < numCorr; j++) {
                flagged = flag_v_p->at(k++) != (char) 0;
                if (ddebug) cout << flagged << " " ;
                if (flagged) notNull++;
                flagCell(j, i) = flagged; // flagCell(j, i) || flagged;    // Let's OR the content of flag with what's found in the BDF flags.
                cellFlagged = cellFlagged || flagged;
                allSet = allSet && flagged;
            }
        }
        if (ddebug) cout << endl;
        flag.put((uInt)iRow0, flagCell);
        //if (notNull)  cout << "row " << iRow0 <<  " - putCell counted actually" << notNull << " non null flags so that cellFlagged = " << cellFlagged << endl;

        flagRow.put((uInt)iRow0, flagRow.get((uInt)iRow0) || allSet);  // Let's OR the content of flagRow with what's found in the BDF flags.

        LOGEXIT("putCell");
        return cellFlagged;
    }

    /*
     * Used for correlator data.
     */
    static pair<uInt, uInt> mergeAndPut( CorrelationModeMod::CorrelationMode correlationMode,
                                         CorrelationModeMod::CorrelationMode ocorrelationMode,
                                         MSFlagAccumulator<char>& autoAccumulator,
                                         MSFlagAccumulator<char>& crossAccumulator,
                                         uInt iRow0,
                                         ArrayColumn<Bool>& flag,
                                         ScalarColumn<Bool> flagRow) {
        LOGENTER("mergeAnPut");

        unsigned int	numIntegration = autoAccumulator.numIntegration();
        unsigned int	numANT	       = autoAccumulator.numBAL();
        unsigned int  numBAL         = crossAccumulator.numBAL();
        unsigned int	numDD	       = autoAccumulator.numDD();

        if (numIntegration != crossAccumulator.numIntegration()	||
            (numANT * (numANT - 1 )) / 2 != numBAL			||
            numDD			   != crossAccumulator.numDD() )
            throw ProcessFlagsException("The accumulators of cross and auto data flags do not have compatible dimensions");

        pair< vector<FLAG_SHAPE * >*, vector<FLAG_V * >* > autoFlagPair  = autoAccumulator.orderedByDDTIMBAL();
        pair< vector<FLAG_SHAPE * >*, vector<FLAG_V * >* > crossFlagPair = crossAccumulator.orderedByDDTIMBAL();

        vector<pair <unsigned int, unsigned int> *>* autoFlagShapes_p_v_p = autoFlagPair.first;
        vector<vector<char> *>* autoFlagValues_p_v_p = autoFlagPair.second;

        vector<pair <unsigned int, unsigned int> *>* crossFlagShapes_p_v_p = crossFlagPair.first;
        vector<vector<char> *>* crossFlagValues_p_v_p = crossFlagPair.second;

        unsigned int		kAuto	       = 0;
        unsigned int		kCross	       = 0;

        Transposer t(numANT);  // Yes we need that since the filler has performed a transposition on the base line order !

        unsigned int		numFlaggedRows = 0; // Maintain a count of the number of flagged rows.

        //
        // For each Data Description...
        for (unsigned int iDD = 0; iDD < numDD; iDD++) {
            //
            // ... put the flags for auto correlation firstly
            if (correlationMode != CorrelationModeMod::CROSS_ONLY && ocorrelationMode != CorrelationModeMod::CROSS_ONLY) {
                if (debug) cout << "AUTO numInt =" << numIntegration << ", numDD = " << numDD << ", iDD = " << iDD << ", " << autoFlagShapes_p_v_p->at(0)->first << ", " << autoFlagShapes_p_v_p->at(0)->second << endl;
                for (unsigned int iTIMEANT = 0; iTIMEANT < numIntegration * numANT; iTIMEANT++) {
                    if (putCell(autoFlagShapes_p_v_p->at(kAuto),
                                autoFlagValues_p_v_p->at(kAuto),
                                iRow0,
                                flag,
                                flagRow))
                        numFlaggedRows++;
                    //if (count_if(autoFlagValues_p_v_p->at(kAuto)->begin(), autoFlagValues_p_v_p->at(kAuto)->end(), isNotNull) != 0)
                    //  numFlaggedRows++;

                    iRow0++;
                    kAuto++;
                }
            }

            //
            // ... put the flags for cross correlation then
            if (correlationMode != CorrelationModeMod::AUTO_ONLY && ocorrelationMode != CorrelationModeMod::AUTO_ONLY) {
                if (debug) cout << "CROSS " << numIntegration << ", " << iDD << ", " << crossFlagShapes_p_v_p->at(0)->first << ", " << crossFlagShapes_p_v_p->at(0)->second << endl;
                for (unsigned int iIntegration = 0; iIntegration < numIntegration; iIntegration++) {
                    for (unsigned int iBAL = 0; iBAL < numBAL; iBAL++) {
                        if (putCell(crossFlagShapes_p_v_p->at(kCross),
                                    crossFlagValues_p_v_p->at(kCross),
                                    iRow0 + t.transposed(iBAL), // This fixes the bug related in http://jira.alma.cl/browse/PRTSIR-9678
                                    flag,
                                    flagRow))
                            numFlaggedRows++;
                        kCross++;
                    }
                    iRow0 += numBAL;
                }
            }
        }
        LOGEXIT("mergeAndPut");
        return make_pair(iRow0, numFlaggedRows);
    }

    /*
     * Used for Radiometer data.
     */
    static pair<uInt, uInt> put( MSFlagAccumulator<char>& accumulator, uInt iRow0,
                                 ArrayColumn<Bool>& flag, ScalarColumn<Bool> flagRow,
                                 bool skipFirstIntegration=false) {
        LOGENTER("put");
        pair< vector<FLAG_SHAPE * >*,
              vector<FLAG_V * >* > cds = accumulator.orderedByDDTIMBAL(skipFirstIntegration); 

        vector<pair <unsigned int, unsigned int> *>* shapes_p_v_p = cds.first;
        vector<vector<char> *>* values_p_v_p = cds.second;
        unsigned int numCells = shapes_p_v_p->size();

        uInt numFlaggedRows = 0;
        for (unsigned int i = 0; i < numCells; i++) {
            FLAG_V * flag_v_p = values_p_v_p->at(i);
            int numChan = shapes_p_v_p->at(i)->first;
            int numCorr = shapes_p_v_p->at(i)->second;

            bool cellFlagged = false;
            bool flagged = false;

            Matrix<Bool> flagCell(IPosition(2, numCorr, numChan));
            flag.get((uInt)iRow0, flagCell);

            bool allSet = true;

            int k = 0;
            for (int i = 0;  i < numChan; i++)
                for (int j = 0; j < numCorr; j++) {
                    flagged = flag_v_p->at(k++) != (char) 0;

                    flagCell(j, i) = flagged ; // flagCell(j, i) || flagged;    // Let's OR the content of flag with what's found in the BDF flags.
                    cellFlagged = cellFlagged || flagged;
                    allSet = allSet && flagged;
                }
            if (cellFlagged) numFlaggedRows ++;

            flag.put((uInt)iRow0, flagCell);
            //    flagRow.put((uInt)iRow0, flagRow.get((uInt)iRow0) || allSet);  // Let's OR the content of flagRow with what's found in the BDF flags.
            flagRow.put((uInt)iRow0, allSet);  // Let's OR the content of flagRow with what's found in the BDF flags.
            iRow0++;
        }

        LOGEXIT("put");
        return make_pair(iRow0, numFlaggedRows);
    }

    void sdm::loadBDFlags(map<string, unsigned int>& abbrev2bitpos) {
        LOGENTER("void loadBDFlags(map<string, unsigned int>& abbrev2bitpos");

        string bdflagsFilename = "bdflags.abbrev.txt";
        char * rootDir_p;

        string bdflagsPath = "";
        const std::list<std::string> &state_path = AppStateSource::fetch( ).dataPath( );
        for ( std::list<std::string>::const_iterator it=state_path.begin(); bdflagsPath.size( ) == 0 && it != state_path.end(); ++it ) {
            std::string path = *it + "/alma/asdm/" + bdflagsFilename;
            if ( file_exists(path) ) bdflagsPath = path;
        }
        if ( bdflagsPath.size( ) == 0 ) {
            std::ostringstream errstream;
            errstream << "The file '" << bdflagsPath << "' containing the collection of BDF flags can't be found." << endl;
            error(errstream.str());
        }
        
        ifstream bdflagsAbbrev;
        bdflagsAbbrev.clear(ifstream::failbit | ifstream::badbit );
        string line;
        unsigned int bitposition = 0;
        try {
            bdflagsAbbrev.open(bdflagsPath.c_str(), fstream::in);
            while (! bdflagsAbbrev.eof()) {
                getline(bdflagsAbbrev, line);
                if (line.size() > 0)
                    abbrev2bitpos[line] = bitposition++;
            }
            bdflagsAbbrev.close();
        } catch (ifstream::failure e) {
            errstream.str("");
            errstream << "I/O error with the file containing the collection of BDF flags, the message was '" << e.what() << "'." << endl;
            error(errstream.str());
        }
        LOGEXIT("void loadBDFlags(map<string, unsigned int>& abbrev2bitpos");
    }

    typedef map<string, unsigned int> s2ui;
    typedef pair<const string, unsigned int> s_ui;

    /** 
     * This function tries to allocate the number of integrations into slices
     * based on the average size of one integration and the requested maximum size of one slice.
     * The sum of all of the values in the returned vector is always equal to nIntegrations.
     * @parameter nIntegrations the number of integrations in the BDF
     * @parameter bdfSize the size of the BDF, in bytes
     * @parameter approxSizeInMemory the approximate size that one wants for one slice, in bytes.
     */
    static vector<unsigned int> getIntegrationSlices(int nIntegrations, uint64_t bdfSize, uint64_t approxSizeInMemory) {
        if (debug) cout << "getIntegrationSlices: entering" << endl;
        vector<unsigned int> result;
        if (nIntegrations > 0) {
            if (bdfSize < approxSizeInMemory) {
                result.push_back(nIntegrations);
            } else {
                uint64_t avgIntSize = bdfSize/nIntegrations;
                unsigned int nIntPerSlice = approxSizeInMemory/avgIntSize;
                unsigned int nSlice = nIntegrations/nIntPerSlice;
                result.resize(nSlice,nIntPerSlice);
                unsigned int residualInts = nIntegrations % nIntPerSlice;
                if (residualInts > (nSlice*nIntPerSlice/5)) {
                    // if the number of left over integrations is more than 1/5 of what's in the other slices
                    // then this makes an acceptable last slice - no need to redistribute
                    result.push_back(residualInts);
                } else {
                    // final slice would too small, redistribute it to the other slices
                    while(residualInts > 0) {
                        for (unsigned int i = 0; residualInts > 0 && i < result.size(); i++) {
                            result[i]++; 
                            residualInts--;
                        }
                    }
                }
            }
        }
        if (debug) cout << "getItegrationSlices: exiting" << endl;
        return result;
    }

    static unsigned int flagsSizeIncludingSPWAndPOL( unsigned int numAntenna,
                                                     CorrelationModeMod::CorrelationMode correlationMode,
                                                     const SDMDataObject::DataStruct & dataStruct) {
        unsigned int autoresult = 0, crossresult = 0;

        const vector<SDMDataObject::Baseband>& basebands = dataStruct.basebands();
        for (SDMDataObject::Baseband baseband : basebands)
            for (SDMDataObject::SpectralWindow spw : baseband.spectralWindows()) {
                if (correlationMode != CorrelationModeMod::CROSS_ONLY) {
                    autoresult += spw.sdPolProducts().end() - spw.sdPolProducts().begin();
                }

                if (correlationMode != CorrelationModeMod::AUTO_ONLY) {
                    crossresult += spw.crossPolProducts().end() - spw.crossPolProducts().begin();
                }
            }
        return (autoresult * numAntenna + crossresult * numAntenna * (numAntenna - 1) / 2);
    }

    static void processCorrelatorFlags( unsigned int numIntegration, const string &bdfPath,
                                        const vector<string> &antennas,
                                        const vector<pair<string, string> > &dataDescriptions,
                                        MSFlagEval &flagEval, uInt &iMSRow,
                                        ArrayColumn<Bool> &flag, ScalarColumn<Bool> &flagRow,
                                        CorrelationModeMod::CorrelationMode ocorrelationMode ) {

        std::regex  ALMACorrelatorFlagsAxesRegex("(BAL )?ANT (BAB )?(SPW )?(POL )?");

        uInt numFlaggedRows = 0;
        SDMDataObjectStreamReader sdosr;

        sdosr.open(bdfPath);
        CorrelationModeMod::CorrelationMode	correlationMode = sdosr.correlationMode();
        const SDMDataObject::BinaryPart & flagsBP = sdosr.dataStruct().flags();
        const vector<AxisName> & flagsAxes = flagsBP.axes();
        CorrelatorFlagsAxes::set(flagsAxes);

        ostringstream oss;
        //hack06 hack06_instance(oss);
        //for_each(flagsAxes.begin(), flagsAxes.end(), hack06_instance);
        for_each(flagsAxes.begin(), flagsAxes.end(), [&](AxisName ax){oss << CAxisName::toString(ax) << " ";});
        string axes_s = oss.str();
        if (axes_s.find("SPW ") == std::string::npos) {
            /*
             * This is a workaround to take into account the possbility that the description of the axes of the flags attachment (attribute axes) *does not*
             * contain the 'SPW' (expectedly between 'BAB' and 'POL') while the flags attachment contain values for *all* the spw of each baseband.
             * The test is based on the number of flags values announced in the flags description (attribute size) which is compared with the calculation of
             * the number of flags *if 'SPW' had present in the list of axes*. If the equality is verified then presence of "SPW" is considered as true.
             *
             * This work was done in the context of http://jira.alma.cl/browse/PRTSPR-21531
             */
            unsigned int announcedFlagSize = flagsBP.size();
            unsigned int calculatedFlagSize = flagsSizeIncludingSPWAndPOL(antennas.size(), correlationMode, sdosr.dataStruct());
            if ( announcedFlagSize == calculatedFlagSize ) {
                CorrelatorFlagsAxes::hasSPW(true);
                if (debug) cout << announcedFlagSize << "==" << calculatedFlagSize << " -> forcing SPW to be present." << endl;
            }
            else {
                CorrelatorFlagsAxes::hasSPW(true);
                if (debug) cout << announcedFlagSize << "!=" << calculatedFlagSize << " -> considering SPW as absent." << endl;
            }
        }

        // Check the validity of the sequence of flags axes.
        if (!std::regex_match(axes_s, ALMACorrelatorFlagsAxesRegex)) {
            throw ProcessFlagsException("'" + oss.str() + "' is not a valid sequence of flags axes for an ALMA correlator.");
        }

        unsigned int				numBAL		= antennas.size() * (antennas.size() - 1) / 2;
        unsigned int				numDD		= dataDescriptions.size();

        MSFlagAccumulator<char> autoAccumulator(numIntegration, antennas.size(), numDD);
        MSFlagAccumulator<char> crossAccumulator(numIntegration, numBAL, numDD);
        if ( correlationMode != CorrelationModeMod::CROSS_ONLY )
            autoAccumulator.resetIntegration();
        if ( correlationMode != CorrelationModeMod::AUTO_ONLY )
            crossAccumulator.resetIntegration();

        while (sdosr.hasSubset()) {
            const FLAGSTYPE * flags_p;

            unsigned int numFlags = sdosr.getSubset().flags(flags_p);
            pair<unsigned int, const FLAGSTYPE *> flagsPair(numFlags, flags_p);
            traverseALMACorrelatorFlagsAxes(sdosr.dataStruct().basebands(),
                                            antennas,
                                            flagsPair,
                                            flagEval,
                                            correlationMode,
                                            autoAccumulator,
                                            crossAccumulator);
            if ( correlationMode != CorrelationModeMod::CROSS_ONLY )
                autoAccumulator.nextIntegration();
            if ( correlationMode != CorrelationModeMod::AUTO_ONLY )
                crossAccumulator.nextIntegration();
        }

        pair<uInt, uInt> putReturn = mergeAndPut(correlationMode,
                                                 ocorrelationMode,
                                                 autoAccumulator,
                                                 crossAccumulator,
                                                 iMSRow,
                                                 flag,
                                                 flagRow);
        iMSRow = putReturn.first;
        numFlaggedRows += putReturn.second;
        sdosr.close();

        return;
    }

    void sdm::processCorrelatorFlagsPerSlices( MainRow *mR_p, unsigned int iASDMIndex,
                                               const vector<int32_t> &mainRowIndex,
                                               const string &bdfPath, uint64_t bdfSliceSizeInMb,
                                               const vector<string> &antennas,
                                               const vector<pair<string, string> > &dataDescriptions,
                                               MSFlagEval &flagEval, uInt &iMSRow,
                                               ArrayColumn<Bool> &flag, ScalarColumn<Bool> &flagRow,
                                               CorrelationModeMod::CorrelationMode ocorrelationMode ) {

        // Regular expressions for the correct sequences of axes in the flags in the case of ALMA data.
        std::regex  ALMACorrelatorFlagsAxesRegex("(BAL )?ANT (BAB )?(SPW )?(POL )?");

        uInt numFlaggedRows = 0;
        SDMDataObjectStreamReader sdosr;

        sdosr.open(bdfPath);
        CorrelationModeMod::CorrelationMode	correlationMode = sdosr.correlationMode();
        const SDMDataObject::BinaryPart & flagsBP = sdosr.dataStruct().flags();
        const vector<AxisName> & flagsAxes = flagsBP.axes();
        CorrelatorFlagsAxes::set(flagsAxes);

        ostringstream oss;
        //hack06 hack06_instance(oss);
        //for_each(flagsAxes.begin(), flagsAxes.end(), hack06_instance);
        for_each(flagsAxes.begin(), flagsAxes.end(), [&](AxisName ax){oss << CAxisName::toString(ax) << " ";});
        string axes_s = oss.str();
        if (axes_s.find("SPW ") == std::string::npos) {
            /*
             * This is a workaround to take into account the possbility that the description of the axes of the flags attachment (attribute axes) *does not*
             * contain the 'SPW' (expectedly between 'BAB' and 'POL') while the flags attachment contain values for *all* the spw of each baseband.
             * The test is based on the number of flags values announced in the flags description (attribute size) which is compared with the calculation of
             * the number of flags *if 'SPW' had present in the list of axes*. If the equality is verified then presence of "SPW" is considered as true.
             *
             * This work was done in the context of http://jira.alma.cl/browse/PRTSPR-21531
             */
            unsigned int announcedFlagSize = flagsBP.size();
            unsigned int calculatedFlagSize = flagsSizeIncludingSPWAndPOL(antennas.size(), correlationMode, sdosr.dataStruct());
            if ( announcedFlagSize == calculatedFlagSize ) {
                CorrelatorFlagsAxes::hasSPW(true);
                if (debug) cout << announcedFlagSize << "==" << calculatedFlagSize << " -> forcing SPW to be present." << endl;
            }
            else {
                CorrelatorFlagsAxes::hasSPW(false);
                if (debug) cout << announcedFlagSize << "!=" << calculatedFlagSize << " -> considering SPW as absent." << endl;
            }
        }

        // Check the validity of the sequence of flags axes.
        if (!std::regex_match(axes_s, ALMACorrelatorFlagsAxesRegex)) {
            throw ProcessFlagsException("'" + oss.str() + "' is not a valid sequence of flags axes for an ALMA correlator.");
        }

        int N = mR_p->getNumIntegration();
        uint64_t bdfSize = mR_p->getDataSize();
        vector<unsigned int> integrationSlices(getIntegrationSlices(N, bdfSize, bdfSliceSizeInMb*1024*1024));

        int32_t		numberOfIntegrations	 = 0;
        int32_t		numberOfReadIntegrations = 0;

        unsigned int				numBAL		= antennas.size() * (antennas.size() - 1) / 2;
        unsigned int				numDD		= dataDescriptions.size();

        for (unsigned int j = 0; j < integrationSlices.size(); j++) {
            if (debug) cout << "PLATOON" << endl;
            numberOfIntegrations = integrationSlices[j];
            if (debug) cout << "------------> " << numberOfIntegrations << endl;
            if (numberOfIntegrations) {
                MSFlagAccumulator<char> autoAccumulator(numberOfIntegrations, antennas.size(), numDD);
                MSFlagAccumulator<char> crossAccumulator(numberOfIntegrations, numBAL, numDD);
                if ( correlationMode != CorrelationModeMod::CROSS_ONLY )
                    autoAccumulator.resetIntegration();
                if ( correlationMode != CorrelationModeMod::AUTO_ONLY )
                    crossAccumulator.resetIntegration();

                for (int iIntegration = 0; iIntegration < numberOfIntegrations; iIntegration++) {
                    const FLAGSTYPE * flags_p;
                    unsigned int numFlags = sdosr.getSubset().flags(flags_p);
                    pair<unsigned int, const FLAGSTYPE *> flagsPair(numFlags, flags_p);
                    traverseALMACorrelatorFlagsAxes(sdosr.dataStruct().basebands(),
                                                    antennas,
                                                    flagsPair,
                                                    flagEval,
                                                    correlationMode,
                                                    autoAccumulator,
                                                    crossAccumulator);
                    if ( correlationMode != CorrelationModeMod::CROSS_ONLY )
                        autoAccumulator.nextIntegration();
                    if ( correlationMode != CorrelationModeMod::AUTO_ONLY )
                        crossAccumulator.nextIntegration();
                }
                pair<uInt, uInt> putReturn = mergeAndPut(correlationMode,
                                                         ocorrelationMode,
                                                         autoAccumulator,
                                                         crossAccumulator,
                                                         iMSRow,
                                                         flag,
                                                         flagRow);
                iMSRow = putReturn.first;
                numFlaggedRows += putReturn.second;
                numberOfReadIntegrations += numberOfIntegrations;

                infostream.str("");

                infostream << "ASDM Main row #" << mainRowIndex[iASDMIndex] << " - " << numberOfReadIntegrations   << "/" << N << " integrations done so far.";
                info(infostream.str());
            }
        }

        // this should no longer be necessary, but keep this block in place just in case
        if (debug) cout << "REMAINING" << endl;
        uint32_t numberOfRemainingIntegrations = N - numberOfReadIntegrations;
        if (numberOfRemainingIntegrations) {
            MSFlagAccumulator<char> autoAccumulator(numberOfRemainingIntegrations, antennas.size(), numDD);
            MSFlagAccumulator<char> crossAccumulator(numberOfRemainingIntegrations, numBAL, numDD);
            if ( correlationMode != CorrelationModeMod::CROSS_ONLY )
                autoAccumulator.resetIntegration();
            if ( correlationMode != CorrelationModeMod::AUTO_ONLY )
                crossAccumulator.resetIntegration();
            for (unsigned int iIntegration = 0; iIntegration < numberOfRemainingIntegrations; iIntegration++) {
                const FLAGSTYPE * flags_p;
                unsigned int numFlags = sdosr.getSubset().flags(flags_p);
                pair<unsigned int, const FLAGSTYPE *> flagsPair(numFlags, flags_p);
                traverseALMACorrelatorFlagsAxes(sdosr.dataStruct().basebands(),
                                                antennas,
                                                flagsPair,
                                                flagEval,
                                                correlationMode,
                                                autoAccumulator,
                                                crossAccumulator);
                if ( correlationMode != CorrelationModeMod::CROSS_ONLY )
                    autoAccumulator.nextIntegration();
                if ( correlationMode != CorrelationModeMod::AUTO_ONLY )
                    crossAccumulator.nextIntegration();
            }
            pair<uInt, uInt> putReturn = mergeAndPut(correlationMode,
                                                     ocorrelationMode,
                                                     autoAccumulator,
                                                     crossAccumulator,
                                                     iMSRow,
                                                     flag,
                                                     flagRow);
            iMSRow = putReturn.first;
            numFlaggedRows += putReturn.second;
            numberOfReadIntegrations += numberOfIntegrations;

            infostream.str("");

            infostream << "ASDM Main row #" << mainRowIndex[iASDMIndex]   << " - " << numberOfReadIntegrations << "/" << N << " integrations done so far.";
            info(infostream.str());
        }
        sdosr.close();
    }


    // 'lazy' parameter must be set to true with the option if the MS has been produced by asdm2MS with the option --lazy !!!
    bool sdm::bd_flagger( const std::string &ms, const std::string &ocm,
                          const std::string &scansOptionValue, bool lazy,
                          bool processUncorrectedData ) {

        const string flagcond_opt = "ALL";

        LOGENTER("sdm::bdf_flagger( " + ms + " , ... )");

        string dsName = sdm_path;
        string msName = ms;
        bitset<32> flagmask;
        map<string, unsigned int> abbrev2bitpos;

        uint64_t bdfSliceSizeInMb = 500; // The default size of the BDF slice hold in memory.
        bool checkdupints = true; // hidden option, used to turn off duplicate integration checks for RADIOMETER data during unit tests.
        CorrelationModeMod::CorrelationMode ocorrelationMode = CorrelationModeMod::CROSS_AND_AUTO;

        // Load the BDF flags abbreviations.
        loadBDFlags(abbrev2bitpos);

        //string abbrevList;
        //hack02 hack02_instance(abbrevList);
        //for_each (abbrev2bitpos.begin(), abbrev2bitpos.end(), hack02_instance);
        string abbrevList = accumulate( abbrev2bitpos.begin(), abbrev2bitpos.end(), string( ),
                                        []( const std::string &acc, const map<string, unsigned int>::value_type &v ) {
                                            return acc + "\n\t* " + v.first;
                                        });

        // bdflags2MS processing "-f ALL" (FLAGCOND) flag
        for ( s_ui p: abbrev2bitpos ) {
            flagmask.set(p.second, true);
            infostream << " " << p.first;
        }

        if (ocm == "ca")
            ocorrelationMode = CorrelationModeMod::CROSS_AND_AUTO;
        else if (ocm == "ao")
            ocorrelationMode =  CorrelationModeMod::AUTO_ONLY;
        else if (ocm == "co")
            ocorrelationMode = CorrelationModeMod::CROSS_ONLY;
        else {
            errstream.str("");
            errstream << "The value given to the option ocm ('" << ocm << "') is invalid. Expecting one of {'ca', 'ao'}" << endl;
            error(errstream.str());
        }

        // The tables will be loaded in memory only on demand.
        ASDMParseOptions parseOpt; parseOpt.loadTablesOnDemand(true);
        ASDM ds;

        // Open the dataset.
        try {
            ds.setFromFile(sdm_path, parseOpt);
        } catch (ConversionException e) {
            errstream.str("");
            errstream << e.getMessage() << endl;
            error(errstream.str());
        } catch (...) {
            errstream.str("");
            errstream << "Unexpected error while trying to access the dataset." << endl;
            error(errstream.str());
        }

        //
        // Reject the dataset if it does not contain ALMA data.
        //
        // In order to decide if it contains ALMA data, let's consider simply
        // the first row of the ExecBlock table.
        //
        const vector<ExecBlockRow *>& ebRs = ds.getExecBlock().get();
        if (ebRs.size() == 0) {
            errstream.str("");
            errstream << "The execblock table has no row; I can't determine the origin of the data." << endl;
            error(errstream.str());
        }
        ExecBlockRow * ebR = ebRs[0];
        string telescopeName = str_toupper(trim_copy(ebR->getTelescopeName()));
        vector<string> telescopeNames = {"ALMA", "OSF", "AOS"};
        if (find(telescopeNames.begin(), telescopeNames.end(), telescopeName) == telescopeNames.end()) {
            errstream.str("");
            errstream << "This dataset announces telescopeName == '" << telescopeName << "', which is not ALMA. Flags can't be processed." << endl;
            error(errstream.str());
        }

        if (debug) {
            cout << "checkdupints : " << checkdupints << endl;
        }

        //
        // Selection of the kind of data - uncorrected or corrected - to consider.
        //
        infostream.str("");
        infostream << "only " << (processUncorrectedData ? "uncorrected" : "corrected") << " data will be considered." << endl;
        info(infostream.str());

        //
        // Selection of the scans to consider.
        //
        vector<ScanRow *>	scanRow_v	       = ds.getScan().get();
        map<int, set<int> > all_eb_scan_m;
        for (vector<ScanRow *>::size_type i = 0; i < scanRow_v.size(); i++)
            all_eb_scan_m[scanRow_v[i]->getExecBlockId().getTagValue()].insert(scanRow_v[i]->getScanNumber());

        vector<ScanRow *>	selectedScanRow_v;
        map<int, set<int> >   selected_eb_scan_m;

        string scansOptionInfo;
        if (scansOptionValue.size()>0) {
            map<int, set<int> > eb_scan_m;
            int status = scansParser(scansOptionValue, eb_scan_m);

            if (status == 0) {
                errstream.str("");
                errstream << "'" << scansOptionValue << "' is an invalid scans selection." << endl;
                error(errstream.str());
            }

            vector<ScanRow *> scanRow_v = ds.getScan().get();
            map<int, set<int> >::iterator iter_m = eb_scan_m.find(-1);

            if (iter_m != eb_scan_m.end())
                for (map<int, set<int> >::iterator iterr_m = all_eb_scan_m.begin(); iterr_m != all_eb_scan_m.end(); iterr_m++) {
                    if ((iter_m->second).empty())
                        selected_eb_scan_m[iterr_m->first] = iterr_m->second;
                    else
                        selected_eb_scan_m[iterr_m->first] = SetAndSet<int>(iter_m->second, iterr_m->second);
                }
            for (map<int, set<int> >::iterator iterr_m = all_eb_scan_m.begin(); iterr_m != all_eb_scan_m.end(); iterr_m++)
                if ((iter_m=eb_scan_m.find(iterr_m->first)) != eb_scan_m.end()) {
                    if ((iter_m->second).empty())
                        selected_eb_scan_m[iterr_m->first].insert((iterr_m->second).begin(), (iterr_m->second).end());
                    else {
                        set<int> s = SetAndSet<int>(iter_m->second, iterr_m->second);
                        selected_eb_scan_m[iterr_m->first].insert(s.begin(), s.end());
                    }
                }

            ostringstream	oss;
            oss << "The following scans will be processed : " << endl;
            for (map<int, set<int> >::const_iterator iter_m = selected_eb_scan_m.begin(); iter_m != selected_eb_scan_m.end(); iter_m++) {
                oss << "eb#" << iter_m->first << " -> " << displaySet<int>(iter_m->second) << endl;

                Tag		execBlockTag  = Tag(iter_m->first, TagType::ExecBlock);
                for (set<int>::const_iterator iter_s = iter_m->second.begin();
                     iter_s		     != iter_m->second.end();
                     iter_s++)
                    selectedScanRow_v.push_back(ds.getScan().getRowByKey(execBlockTag, *iter_s));
            }

            scansOptionInfo = oss.str();
        } else {
            selectedScanRow_v = ds.getScan().get();
            selected_eb_scan_m = all_eb_scan_m;
            scansOptionInfo = "All scans of all exec blocks will be processed \n";
        }
        infostream.str("");
        infostream << scansOptionInfo;
        info(infostream.str());

        infostream.str("");
        infostream << "The MS FLAG column will be calculated for the correlation mode '" << CCorrelationMode::toString(ocorrelationMode) << "'." << endl;
        info(infostream.str());

        // Open the MS MAIN table of the measurement set in Update mode.

        Table mainTable;
        // ArrayColumn<Complex>	data;
        ArrayColumn<Bool>	flag;
        ScalarColumn<Bool>	flagRow;

        try {
            mainTable = Table (ms, Table::Update);
            // data.attach(mainTable, "DATA");
            flag.attach(mainTable, "FLAG");
            flagRow.attach(mainTable, "FLAG_ROW");
        } catch (AipsError e) {
            errstream.str("");
            errstream << e.getMesg() << endl;
            error(errstream.str());
        }

        // Regular expressions for the correct sequences of axes in the flags in the case of ALMA data.
        std::regex ALMARadiometerPackedFlagsAxesRegex("(TIM ANT )?(BAB BIN POL )?");
        std::regex ALMARadiometerFlagsAxesRegex("(ANT )?(BAB BIN POL )?");

        ConfigDescriptionTable &	cfgT = ds.getConfigDescription();
        DataDescriptionTable &	ddT  = ds.getDataDescription();

        SDMDataObjectStreamReader sdosr;
        SDMDataObjectReader sdor;

        unsigned long fm2Ulong = 0;
        if (processUncorrectedData) {                  // A bit of manipulation on the mask in the case of uncorrected data...
            bitset<32> UDflagmask(flagmask.to_ulong());  // ... because ...
            UDflagmask.reset(abbrev2bitpos["WVR_APC"]);  // the mask must NOT have the WVR_APC bit set to 1 !
            fm2Ulong = UDflagmask.to_ulong();            // (see https://bugs.nrao.edu/browse/CAS-8491)
        } else {
            fm2Ulong = flagmask.to_ulong();
        }
        MSFlagEval flagEval(fm2Ulong, processUncorrectedData ? 3 : 0); // If we process uncorrected data we ignore the combination
        // 3 = WVR_APC | INTEGRATION_FULLY_BLANKED

        //
        //
        // Consider only the Main rows whose execBlockId and scanNumber attributes correspond to the selection and which
        // contain the data with appropriate atmospheric phase correction characteristics.
        //
        AtmPhaseCorrectionMod::AtmPhaseCorrection queriedAPC = processUncorrectedData ? AtmPhaseCorrectionMod::AP_UNCORRECTED : AtmPhaseCorrectionMod::AP_CORRECTED;
        vector<MainRow*> v;
        vector<int32_t> mainRowIndex;
        const vector<MainRow *>& temp = ds.getMain().get();

        for ( vector<MainRow *>::const_iterator iter_v = temp.begin(); iter_v != temp.end(); iter_v++) {
            map<int, set<int> >::iterator iter_m = selected_eb_scan_m.find((*iter_v)->getExecBlockId().getTagValue());
            if ( iter_m != selected_eb_scan_m.end() && iter_m->second.find((*iter_v)->getScanNumber()) != iter_m->second.end() ) {
                bool toBeProcessed  = cfgT.getRowByKey((*iter_v)->getConfigDescriptionId())->getProcessorType() == RADIOMETER ; // RADIOMETER data are always put in both (UN/CORRECTED) MS.
                if (!toBeProcessed) {
                    vector<AtmPhaseCorrectionMod::AtmPhaseCorrection > apc_v = cfgT.getRowByKey((*iter_v)->getConfigDescriptionId())->getAtmPhaseCorrection();
                    toBeProcessed =   find(apc_v.begin(), apc_v.end(), queriedAPC) != apc_v.end();
                }
                if (toBeProcessed) {
                    mainRowIndex.push_back(iter_v - temp.begin());
                    v.push_back(*iter_v);
                }
            }
        }

        infostream.str("");
        infostream << "The dataset has " << temp.size() << " main(s)...";
        infostream << v.size() << " of them in the selected exec blocks / scans." << endl;
        info(infostream.str());

        uInt	iMSRow	    = 0;	// Row index in the MS Main table.
        uInt	iMSRowBegin = 0;	// Index of the first row in the MS Main table of the slice corresponding to one row in the ASDM Main table.

        unsigned int numFlaggedRowsTotal = 0;

        // used in checking for duplicate integrations in the WVR (Radiometer) case
        // This holds the most recent last integration time for each configDescriptionId - but only for Radiometer data.
        map<Tag, double> lastTimeMap;

        iMSRowBegin = iMSRow;
        unsigned int iASDMIndex = 0;
        for (MainRow * mR: v) {
            //
            // Retrieve metadata informations.
            //
            unsigned int numIntegration = mR->getNumIntegration();

            ConfigDescriptionRow * cfgR = cfgT.getRowByKey(mR->getConfigDescriptionId());

            vector<Tag> antennaTags = cfgR->getAntennaId();
            //vector<string> antennas;
            //hack04 hack04_instance(antennas);
            //for_each(antennaTags.begin(), antennaTags.end(), hack04_instance);
            vector<string> antennas = accumulate( antennaTags.begin(), antennaTags.end(), vector<string>( ),
                                                  []( vector<string> &acc, const Tag &tag ) {
                                                      acc.push_back(tag.toString( ));
                                                      return acc;
                                                  });

            vector<Tag> dataDescriptionTags = cfgR->getDataDescriptionId();
            //vector<DataDescriptionRow *> ddRs;
            //hack03 hack03_instance(ddT,ddRs);
            //for_each(dataDescriptionTags.begin(), dataDescriptionTags.end(), hack03_instance);
            vector<DataDescriptionRow*> ddRs = accumulate( dataDescriptionTags.begin(), dataDescriptionTags.end(), vector<DataDescriptionRow *>( ),
                                                           [&]( vector<DataDescriptionRow *> &acc, const Tag &tag ) {
                                                               acc.push_back( ddT.getRowByKey(tag) );
                                                               return acc;
                                                           });

            unsigned int numDD = ddRs.size();

            //vector<pair<string, string> > dataDescriptions;
            //hack05 hack05_instance(dataDescriptions);
            //for_each(ddRs.begin(),
            //     ddRs.end(),
            //     hack05_instance);
            vector<pair<string, string> > dataDescriptions = accumulate( ddRs.begin(), ddRs.end(), vector<pair<string, string> >( ),
                                                                         []( vector<pair<string,string> > &acc, const DataDescriptionRow *row ) {
                                                                             acc.push_back(make_pair(row->getSpectralWindowId( ).toString( ),row->getPolOrHoloId( ).toString( )));
                                                                             return acc;
                                                                         });

            string dataUID = mR->getDataUID().getEntityId().toString();
            replace(dataUID.begin(),dataUID.end(),':','_');
            replace(dataUID.begin(),dataUID.end(),'/','_');
            string bdfPath = dsName+"/ASDMBinary/"+dataUID;
            if (debug) cout << "BDF " << bdfPath << endl;
            ProcessorType pt = cfgR->getProcessorType();
            uInt numFlaggedRows = 0;
            try {
                infostream.str("");
                infostream << "ASDM Main row #" << mainRowIndex[iASDMIndex]
                           << " (scan #" << mR->getScanNumber()
                           <<", subscan #" <<  mR->getSubscanNumber()
                           <<", " << CProcessorType::toString(pt)
                           <<", " << mR->getConfigDescriptionId().toString() << ")"
                           << " - BDF '" << bdfPath << "' - Size " << mR->getDataSize() << " bytes." <<  endl;
                info(infostream.str());

                //
                // Check that the triple (execBlockId,scanNumber, subscanNumber) refers to an existing entry in the subscan table
                // if it's not the case ignore this row of the Main table to have a behaviour similar to the filler's one.
                //
                ASDM&			ds	   = mR -> getTable() . getContainer();
                ScanRow*		scanR_p	   = ds.getScan().getRowByKey(mR -> getExecBlockId(),	mR -> getScanNumber());
                vector<ScanIntent>	scanIntent = scanR_p -> getScanIntent();
                SubscanRow*		sscanR_p   = ds.getSubscan().getRowByKey(mR -> getExecBlockId(),
                                                                         mR -> getScanNumber(),
                                                                         mR -> getSubscanNumber());
                if (sscanR_p == 0) {
                    infostream << "Could not find a row in the Subscan table for the following key value (execBlockId="
                               << mR->getExecBlockId().toString()
                               <<", scanNumber="<< mR->getScanNumber()
                               <<", subscanNum=" << mR->getSubscanNumber() << "). Aborting. "
                               << endl;
                    info(infostream.str());
                    continue;  // goto the next main row.
                }
                // End of check

                switch (pt) {
                    case ProcessorTypeMod::CORRELATOR :
                    {
                        //ddebug =  (mR -> getScanNumber() == 2) ;
                        if (lazy)
                            processCorrelatorFlags( numIntegration, bdfPath, antennas, dataDescriptions,
                                                    flagEval, iMSRow, flag, flagRow, ocorrelationMode );
                        else
                            processCorrelatorFlagsPerSlices( mR, iASDMIndex, mainRowIndex, bdfPath, bdfSliceSizeInMb,
                                                             antennas, dataDescriptions, flagEval, iMSRow, flag,
                                                             flagRow, ocorrelationMode);
                    }
                    break;

                    case ProcessorTypeMod::RADIOMETER :
                    {
                        const SDMDataObject& sdo = sdor.read(bdfPath);

                        // should this correlationMode be skipped? This is probably always AUTO_ONLY, but this is a general test just in case
                        if ((sdo.correlationMode()==CorrelationModeMod::AUTO_ONLY and ocorrelationMode==CorrelationModeMod::CROSS_ONLY) || 
                            (sdo.correlationMode()==CorrelationModeMod::CROSS_ONLY and ocorrelationMode==CorrelationModeMod::AUTO_ONLY)) {
                            if (debug) {
                                cout << "Skipped file " << bdfPath << " due to output correlation mode selection" << endl;
                            }
                        } else {
                            const SDMDataObject::BinaryPart & flagsBP = sdo.dataStruct().flags();
                            const vector<AxisName> & flagsAxes = flagsBP.axes();
                            ostringstream oss;
                            //hack06 hack06_instance(oss);
                            //for_each(flagsAxes.begin(), flagsAxes.end(), hack06_instance);
                            for_each(flagsAxes.begin(), flagsAxes.end(), [&](AxisName ax){oss << CAxisName::toString(ax) << " ";});

                            // Check the validity of the sequence of flags axes (depending on the fact that data are packed or not).
                            // also check tos ee if the first integration should be skipped - only done for packed data
                            // must be known before accumulator is instantiated.
                            bool skipFirstIntegration = false;
                            if (sdo.hasPackedData()) {
                                if (!std::regex_match(oss.str(), ALMARadiometerPackedFlagsAxesRegex))
                                    throw ProcessFlagsException("'" + oss.str() + "' is not a valid sequence of flags axes for an ALMA radiometer.");

                                // determine if the first integration should be skipped due a duplicate time from a subscan.
                                const SDMDataSubset& sdmDataSubset = sdo.sdmDataSubsets()[0];
                                int64_t  deltaTime = sdmDataSubset.interval() / sdo.numTime();
                                int64_t startTime = (int64_t)sdmDataSubset.time() -  (int64_t)sdmDataSubset.interval()/2LL + deltaTime/2LL;

                                // should the first integration be skipped? Any actual skipping happens later.
                                skipFirstIntegration = checkdupints && lastTimeMap[mR->getConfigDescriptionId()] == ArrayTime(startTime).getMJD();
                                if (debug && skipFirstIntegration) {
                                    cout << "Duplicate time seen in Row : " << mainRowIndex[iASDMIndex]
                                         << " cdId : " << mR->getConfigDescriptionId()
                                         << " " << mR->getDataUID().getEntityId().toString()
                                         << " numTime : " << sdo.numTime()
                                         << endl;
                                }
                                lastTimeMap[mR->getConfigDescriptionId()] = ArrayTime(startTime+(sdo.numTime()-1)*deltaTime).getMJD();
                            } else {
                                if (!std::regex_match(oss.str(), ALMARadiometerFlagsAxesRegex))
                                    throw ProcessFlagsException("'" + oss.str() + "' is not a valid sequence of flags axes for an ALMA radiometer.");
                            }
                            unsigned int numIntegrations = sdo.hasPackedData() ? sdo.numTime() : sdo.sdmDataSubsets().size();
                            if (numIntegrations != numIntegration) {
                                infostream << "(the number of integrations actually read in the BDF (numIntegrations = " << numIntegrations << ") is different from the value announced in the Main table (numIntegration = " << numIntegration
                                           << "). Using " << numIntegrations << ")";
                            }
                            MSFlagAccumulator<char> accumulator(numIntegrations, antennas.size(), numDD);

                            const FLAGSTYPE *	flags_p	 = NULL;
                            unsigned int	numFlags = 0;

                            if (sdo.hasPackedData()) {
                                numFlags = sdo.tpDataSubset().flags(flags_p);
                                pair<unsigned int, const FLAGSTYPE *> flagsPair(numFlags, flags_p);
                                accumulator.resetIntegration();
                                traverseALMARadiometerFlagsAxes(numIntegrations, sdo.dataStruct().basebands(), antennas, flagsPair, flagEval, accumulator);
                            } else {
                                // duplicate integrations are not checked or skipped here
                                const vector<SDMDataSubset>& sdmDataSubsets = sdo.sdmDataSubsets();
                                accumulator.resetIntegration();
                                for (unsigned int iIntegration = 0; iIntegration < numIntegrations; iIntegration++) {
                                    numFlags = sdmDataSubsets[iIntegration].flags(flags_p);
                                    pair<unsigned int, const FLAGSTYPE *> flagsPair(numFlags, flags_p);
                                    traverseALMARadiometerFlagsAxes(1, sdo.dataStruct().basebands(), antennas, flagsPair, flagEval, accumulator);
                                }
                                // MAFlagAccumulator::flagCell() just returns a data member, it doesn't change anything
                                // in the object, and here, the returned value (x) is not used.  So this can be
                                // safely commented out.
                                // const vector < vector < vector <FLAG_CELL> > >& x = accumulator.flagCell();
                            }
                            infostream.str("");

                            infostream << "ASDM Main row #" << mainRowIndex[iASDMIndex] << " - " << numIntegrations  << "/"
                                       << numIntegrations << " integrations done so far.";
                            info(infostream.str());

                            pair<uInt, uInt> putReturn = put(accumulator, iMSRow, flag, flagRow,skipFirstIntegration);
                            //	  accumulator.dump(cout, true);
                            iMSRow = putReturn.first;
                            numFlaggedRows = putReturn.second;
                        }
                        sdor.done();
                    }
                    break;

                    default:
                        throw ProcessFlagsException("Unrecognized processor type.");
                }
                infostream.str("");
                infostream << "ASDM Main row #" << mainRowIndex[iASDMIndex]
                           << " (scan #" << mR->getScanNumber()
                           <<", subscan #" <<  mR->getSubscanNumber()
                           <<", " << CProcessorType::toString(pt)
                           <<", " << mR->getConfigDescriptionId().toString() << ")"
                           << " - BDF '" << bdfPath << "' - Size " << mR->getDataSize()
                           << " bytes,  produced " << numFlaggedRows
                           << " flagged rows in the MS Main table rows "
                           << iMSRowBegin << " to " << (iMSRow - 1) << endl << endl << endl;
                info(infostream.str());
                iMSRowBegin = iMSRow;
                numFlaggedRowsTotal += numFlaggedRows;
            } catch (AipsError e) {
                info(infostream.str());
                infostream.str("");
                infostream << e.getMesg() << endl;
                info(infostream.str());
                return false;
            } catch (ProcessFlagsException e) {
                info(infostream.str());
                infostream.str("");
                infostream << e.getMessage() << " , bdf path = " << bdfPath << ", processor type = " << CProcessorType::toString(pt) << endl;
                info(infostream.str());
            } catch (SDMDataObjectParserException e) {
                info(infostream.str());
                infostream.str("");
                infostream << e.getMessage() << endl;
                info(infostream.str());
            } catch (SDMDataObjectException e) {
                info(infostream.str());
                infostream.str("");
                infostream << e.getMessage() << endl;
                info(infostream.str());
            } catch (SDMDataObjectReaderException e) {
                info(infostream.str());
                infostream.str("");
                infostream << e.getMessage() << endl;
                info(infostream.str());
            } catch (SDMDataObjectStreamReaderException e) {
                info(infostream.str());
                infostream.str("");
                infostream << e.getMessage() << endl;
                info(infostream.str());
            }
            iASDMIndex++;
        }
        infostream.str("");
        if (mainTable.nrow() > 0) { // this is a paranoid test...
            infostream << numFlaggedRowsTotal << " rows have been flagged in the " << mainTable.nrow()
                       << " of the MS Main table. "
                       << setprecision(4)
                       << ((float) numFlaggedRowsTotal) / mainTable.nrow() * 100.0 << "%." << endl;
            info(infostream.str());
        }
        mainTable.flush();
        LOGEXIT("int main (int argc, char * argv[])");
        return true;
    }
}
