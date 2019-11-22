#ifndef _ASDMTABLEBASE_H_
#define _ASDMTABLEBASE_H_

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/ArrColDesc.h>
#include <tables/Tables/TableRow.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Containers/RecordDesc.h>
#include <casa/Containers/RecordField.h>
#include <casa/Containers/RecordInterface.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/TableLock.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <casa/Exceptions/Error.h>
#include <casa/Utilities/Assert.h>

#include <alma/ASDM/ASDM.h>
#include <alma/ASDM/ArrayTime.h>
#include <alma/ASDM/ArrayTimeInterval.h>

namespace asdm {

    template<typename T>
    class ASDM_TABLE_SINGLETON {
    protected:
        ASDM_TABLE_SINGLETON() {;}
        ~ASDM_TABLE_SINGLETON() {;}
    public:
        static T *instance() {
            if ( NULL == instance_)
                instance_ = new T;
            return (static_cast<T*>(instance_));
        }
    private:
        static T *instance_;
    };

    template<typename T> T *ASDM_TABLE_SINGLETON<T>::instance_ = NULL;
    class ASDM_TABLE_BASE {
    protected:
        ASDM_TABLE_BASE();
        virtual ~ASDM_TABLE_BASE();
        std::string name_;
        casacore::Table* table_p_;
    public:
        casacore::Table* table_p();
        const std::string& name() const;
        virtual const casacore::TableDesc& tableDesc() const = 0;
        void buildAndAttachTable(casacore::MS* attachMS);
        virtual void fill(const ASDM& asdm) = 0;
        void close();

        template<typename T, typename U>  casacore::Vector<U> basic2CASA1D(const std::vector<T>& v) {
            casacore::Vector<U> result;
            if (v.size() == 0)
                return result;
            
            result.resize(v.size());
            for (unsigned int i = 0; i < v.size(); i++)
                result(i) = v.at(i);
            return result;
        }

        template<typename T, typename U>  casacore::Matrix<U> basic2CASA2D(const std::vector<std::vector <T> >& v) {
            casacore::Matrix<U> result;
            if (v.size() == 0 || v.at(0).size() == 0)
                return result;

            result.resize(v.size(), v.at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    result(i,j) = v.at(i).at(j);
            return result;
        }

        template<typename T, typename U>  casacore::Cube<U> basic2CASA3D(const std::vector<std::vector <std::vector <T> > >& v) {
            casacore::Cube<U> result;
            if (v.size() == 0 || v.at(0).size() == 0 || v.at(0).at(0).size() == 0)
                return result;

            result.resize(v.size(), v.at(0).size(), v.at(0).at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    for (unsigned int k = 0; k < v.at(0).at(0).size(); k++)
                        result(i,j,k) = v.at(i).at(j).at(k);
            return result;
        }

        template<typename T, typename U>  casacore::Vector<U> ext2CASA1D(const std::vector<T>& v) {
            casacore::Vector<U> result;
            if (v.size() == 0)
                return result;

            result.resize(v.size());
            for (unsigned int i = 0; i < v.size(); i++)
                result(i) = v.at(i).get();
            return result;
        }

        template<typename T, typename U>  casacore::Matrix<U> ext2CASA2D(const std::vector<std::vector <T> >& v) {
            casacore::Matrix<U> result;
            if (v.size() == 0 || v.at(0).size() == 0)
                return result;

            result.resize(v.size(), v.at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    result(i,j) = v.at(i).at(j).get();
            return result;
        }

        template<typename T, typename U>  casacore::Cube<U> ext2CASA3D(const std::vector<std::vector <std::vector <T> > >& v) {
            casacore::Cube<U> result;
            if (v.size() == 0 || v.at(0).size() == 0 || v.at(0).at(0).size() == 0)
                return result;

            result.resize(v.size(), v.at(0).size(), v.at(0).at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    for (unsigned int k = 0; k < v.at(0).at(0).size(); k++)
                        result(i,j,k) = v.at(i).at(j).at(k).get();
            return result;
        }

        template<typename T, typename U>  casacore::Vector<U> _2CASAString1D(const std::vector<T>& v) {
            casacore::Vector<U> result;
            if (v.size() == 0)
                return result;

            result.resize(v.size());
            for (unsigned int i = 0; i < v.size(); i++)
                result(i) = v.at(i).toString();
            return result;
        }

        template<typename T, typename U>  casacore::Matrix<U> _2CASAString2D(const std::vector<std::vector <T> >& v) {
            casacore::Matrix<U> result;
            if (v.size() == 0 || v.at(0).size() == 0)
                return result;

            result.resize(v.size(), v.at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    result(i,j) = v.at(i).at(j).toString();
            return result;
        }

        template<typename T, typename U>  casacore::Cube<U> _2CASAString3D(const std::vector<std::vector <std::vector <T> > >& v) {
            casacore::Cube<U> result;
            if (v.size() == 0 || v.at(0).size() == 0 || v.at(0).at(0).size() == 0)
                return result;

            result.resize(v.size(), v.at(0).size(), v.at(0).at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    for (unsigned int k = 0; k < v.at(0).at(0).size(); k++)
                        result(i,j,k) = v.at(i).at(j).at(k).toString();
            return result;
        }

        template<typename enumT, typename CenumT>  casacore::Vector<casacore::String> enum2CASA1D (const std::vector<enumT>& v) {
            casacore::Vector<casacore::String> result;
            if (v.size() == 0) return result;

            result.resize(v.size());
            for (unsigned int i = 0; i < v.size(); i++)
                result(i) = CenumT::name(v.at(i));
            return result;
        }

        template<typename enumT, typename CenumT>  casacore::Matrix<casacore::String> enum2CASA2D (const std::vector<std::vector<enumT> >& v) {
            casacore::Matrix<casacore::String> result;
            if (v.size() == 0 || v.at(0).size() == 0) return result;

            result.resize(v.size(), v.at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0 ; j < v.at(0).size(); j++)
                    result(i,j) = CenumT::name(v.at(i).at(j));
            return result;
        }

        template<typename enumT, typename CenumT>  casacore::Cube<casacore::String> enum2CASA3D (const std::vector<std::vector<std::vector<enumT> > >& v) {
            casacore::Cube<casacore::String> result;
            if (v.size() == 0 || v.at(0).size() == 0 || v.at(0).at(0) == 0) return result;

            result.resize(v.size(), v.at(0).size(), v.at(0).at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0 ; j < v.at(0).size(); j++)
                    for (unsigned int k = 0 ; k < v.at(0).at(0).size(); k++)
                        result(i,j,k) = CenumT::name(v.at(i).at(j).at(k));
            return result;
        }

        template<typename T, typename U> casacore::Vector<U> interval2CASA1D(const std::vector<T>& v) {
            casacore::Vector<U> result;
            if (v.size()==0) return result;

            result.resize(v.size());
            for (unsigned int i = 0; i < v.size(); i++)
                result(i) = v.at(i).get()/1e09;
            return result;
        }

        template<typename T, typename U> casacore::Matrix<U> interval2CASA2D(const std::vector<std::vector<T> >& v) {
            casacore::Matrix<U> result;
            if (v.size()==0 || v.at(0).size()) return result;

            result.resize(v.size(), v.at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    result(i,j) = v.at(i).at(j).get()/1.e09;
            return result;
        }

        template<typename T, typename U> casacore::Cube<U> interval2CASA3D(const std::vector<std::vector< std::vector<T> > >& v) {
            casacore::Cube<U> result;
            if (v.size() == 0 || v.at(0).size() == 0 || v.at(0).at(0).size() == 0) return result;

            result.resize(v.size(), v.at(0).size(), v.at(0).at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    for (unsigned int k = 0; k < v.at(0).at(0).size(); k++)
                        result(i,j,k) = v.at(i).at(j).at(k).get()/1.e09;
            return result;
        }

        template<typename U>  casacore::Vector<U> at2CASA1D(const std::vector<ArrayTime>& v) {
            casacore::Vector<U> result;
            if (v.size()==0) return result;

            result.resize(v.size());
            for (unsigned int i = 0; i < v.size(); i++)
                result(i) = v.at(i).get()/1e09;
            return result;
        }

        template<typename U> casacore::Matrix<U> at2CASA2D(const std::vector<std::vector<ArrayTime> >& v) {
            casacore::Matrix<U> result;
            if (v.size()==0 || v.at(0).size()) return result;

            result.resize(v.size(), v.at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    result(i,j) = v.at(i).at(j).get()/1.e09;
            return result;
        }

        template<typename U>  casacore::Cube<U> at2CASA3D(const std::vector<std::vector< std::vector<ArrayTime> > >& v) {
            casacore::Cube<U> result;
            if (v.size() == 0 || v.at(0).size() == 0 || v.at(0).at(0).size() == 0) return result;

            result.resize(v.size(), v.at(0).size(), v.at(0).at(0).size());
            for (unsigned int i = 0; i < v.size(); i++)
                for (unsigned int j = 0; j < v.at(0).size(); j++)
                    for (unsigned int k = 0; k < v.at(0).at(0).size(); k++)
                        result(i,j,k) = v.at(i).at(j).at(k).get()/1.e09;
            return result;
        }

        template<typename U> casacore::Vector<U> ati2CASA1D(const ArrayTimeInterval& ati) {
            casacore::Vector<U> result(2);
            result(0) = ((double) ati.getStart().get()) / ArrayTime::unitsInASecond;
            result(1) = ((double) ati.getDuration().get()) / ArrayTime::unitsInASecond;
            return result;
        }

        template<typename U> casacore::Matrix<U> ati2CASA2D(const std::vector<ArrayTimeInterval>& v) {
            casacore::Matrix<U> result;
            if (v.size() == 0) return result;

            result.resize(v.size(), 2);
            for (std::vector<ArrayTimeInterval>::size_type i = 0; i < v.size(); i++) {
                result(i, 0) = ((double) v[i].getStart().get()) / ArrayTime::unitsInASecond;
                result(i, 1) = ((double) v[i].getDuration().get()) / ArrayTime::unitsInASecond;
            }
            return result;
        }

        template<typename U> casacore::Cube<U> ati2CASA3D(const std::vector<std::vector<ArrayTimeInterval> >& v) {
            casacore::Cube<U> result;
            if (v.size() == 0 || v.at(0).size() == 0) return result;

            result.resize(v.size(), v.at(0).size(), 2);
            for (std::vector<std::vector<ArrayTimeInterval> >::size_type i = 0; i < v.size(); i++) {
                for (std::vector<ArrayTimeInterval>::size_type j = 0; j < v.at(0).size(); j++) {
                    result(i, j, 0) = ((double) v[i][j].getStart().get()) / ArrayTime::unitsInASecond;
                    result(i, j, 1) = ((double) v[i][j].getDuration().get()) / ArrayTime::unitsInASecond;
                }
            }
            return result;
        }

    };

} // end namespace asdm

#endif // _ASDMTABLEBASE_H_
