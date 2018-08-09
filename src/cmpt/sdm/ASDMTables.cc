
#include "ASDMTables.h"

#include <stdint.h>
#include <alma/ASDM/MainTable.h>
#include <alma/ASDM/MainRow.h>
#include <alma/ASDM/AlmaRadiometerTable.h>
#include <alma/ASDM/AlmaRadiometerRow.h>
#include <alma/ASDM/AnnotationTable.h>
#include <alma/ASDM/AnnotationRow.h>
#include <alma/ASDM/AntennaTable.h>
#include <alma/ASDM/AntennaRow.h>
#include <alma/ASDM/CalAmpliTable.h>
#include <alma/ASDM/CalAmpliRow.h>
#include <alma/ASDM/CalAntennaSolutionsTable.h>
#include <alma/ASDM/CalAntennaSolutionsRow.h>
#include <alma/ASDM/CalAppPhaseTable.h>
#include <alma/ASDM/CalAppPhaseRow.h>
#include <alma/ASDM/CalAtmosphereTable.h>
#include <alma/ASDM/CalAtmosphereRow.h>
#include <alma/ASDM/CalBandpassTable.h>
#include <alma/ASDM/CalBandpassRow.h>
#include <alma/ASDM/CalCurveTable.h>
#include <alma/ASDM/CalCurveRow.h>
#include <alma/ASDM/CalDataTable.h>
#include <alma/ASDM/CalDataRow.h>
#include <alma/ASDM/CalDelayTable.h>
#include <alma/ASDM/CalDelayRow.h>
#include <alma/ASDM/CalDeviceTable.h>
#include <alma/ASDM/CalDeviceRow.h>
#include <alma/ASDM/CalFluxTable.h>
#include <alma/ASDM/CalFluxRow.h>
#include <alma/ASDM/CalFocusTable.h>
#include <alma/ASDM/CalFocusRow.h>
#include <alma/ASDM/CalFocusModelTable.h>
#include <alma/ASDM/CalFocusModelRow.h>
#include <alma/ASDM/CalGainTable.h>
#include <alma/ASDM/CalGainRow.h>
#include <alma/ASDM/CalHolographyTable.h>
#include <alma/ASDM/CalHolographyRow.h>
#include <alma/ASDM/CalPhaseTable.h>
#include <alma/ASDM/CalPhaseRow.h>
#include <alma/ASDM/CalPointingTable.h>
#include <alma/ASDM/CalPointingRow.h>
#include <alma/ASDM/CalPointingModelTable.h>
#include <alma/ASDM/CalPointingModelRow.h>
#include <alma/ASDM/CalPositionTable.h>
#include <alma/ASDM/CalPositionRow.h>
#include <alma/ASDM/CalPrimaryBeamTable.h>
#include <alma/ASDM/CalPrimaryBeamRow.h>
#include <alma/ASDM/CalReductionTable.h>
#include <alma/ASDM/CalReductionRow.h>
#include <alma/ASDM/CalSeeingTable.h>
#include <alma/ASDM/CalSeeingRow.h>
#include <alma/ASDM/CalWVRTable.h>
#include <alma/ASDM/CalWVRRow.h>
#include <alma/ASDM/ConfigDescriptionTable.h>
#include <alma/ASDM/ConfigDescriptionRow.h>
#include <alma/ASDM/CorrelatorModeTable.h>
#include <alma/ASDM/CorrelatorModeRow.h>
#include <alma/ASDM/DataDescriptionTable.h>
#include <alma/ASDM/DataDescriptionRow.h>
#include <alma/ASDM/DelayModelTable.h>
#include <alma/ASDM/DelayModelRow.h>
#include <alma/ASDM/DelayModelFixedParametersTable.h>
#include <alma/ASDM/DelayModelFixedParametersRow.h>
#include <alma/ASDM/DelayModelVariableParametersTable.h>
#include <alma/ASDM/DelayModelVariableParametersRow.h>
#include <alma/ASDM/DopplerTable.h>
#include <alma/ASDM/DopplerRow.h>
#include <alma/ASDM/EphemerisTable.h>
#include <alma/ASDM/EphemerisRow.h>
#include <alma/ASDM/ExecBlockTable.h>
#include <alma/ASDM/ExecBlockRow.h>
#include <alma/ASDM/FeedTable.h>
#include <alma/ASDM/FeedRow.h>
#include <alma/ASDM/FieldTable.h>
#include <alma/ASDM/FieldRow.h>
#include <alma/ASDM/FlagTable.h>
#include <alma/ASDM/FlagRow.h>
#include <alma/ASDM/FlagCmdTable.h>
#include <alma/ASDM/FlagCmdRow.h>
#include <alma/ASDM/FocusTable.h>
#include <alma/ASDM/FocusRow.h>
#include <alma/ASDM/FocusModelTable.h>
#include <alma/ASDM/FocusModelRow.h>
#include <alma/ASDM/FreqOffsetTable.h>
#include <alma/ASDM/FreqOffsetRow.h>
#include <alma/ASDM/GainTrackingTable.h>
#include <alma/ASDM/GainTrackingRow.h>
#include <alma/ASDM/HistoryTable.h>
#include <alma/ASDM/HistoryRow.h>
#include <alma/ASDM/HolographyTable.h>
#include <alma/ASDM/HolographyRow.h>
#include <alma/ASDM/ObservationTable.h>
#include <alma/ASDM/ObservationRow.h>
#include <alma/ASDM/PointingTable.h>
#include <alma/ASDM/PointingRow.h>
#include <alma/ASDM/PointingModelTable.h>
#include <alma/ASDM/PointingModelRow.h>
#include <alma/ASDM/PolarizationTable.h>
#include <alma/ASDM/PolarizationRow.h>
#include <alma/ASDM/ProcessorTable.h>
#include <alma/ASDM/ProcessorRow.h>
#include <alma/ASDM/PulsarTable.h>
#include <alma/ASDM/PulsarRow.h>
#include <alma/ASDM/ReceiverTable.h>
#include <alma/ASDM/ReceiverRow.h>
#include <alma/ASDM/SBSummaryTable.h>
#include <alma/ASDM/SBSummaryRow.h>
#include <alma/ASDM/ScaleTable.h>
#include <alma/ASDM/ScaleRow.h>
#include <alma/ASDM/ScanTable.h>
#include <alma/ASDM/ScanRow.h>
#include <alma/ASDM/SeeingTable.h>
#include <alma/ASDM/SeeingRow.h>
#include <alma/ASDM/SourceTable.h>
#include <alma/ASDM/SourceRow.h>
#include <alma/ASDM/SpectralWindowTable.h>
#include <alma/ASDM/SpectralWindowRow.h>
#include <alma/ASDM/SquareLawDetectorTable.h>
#include <alma/ASDM/SquareLawDetectorRow.h>
#include <alma/ASDM/StateTable.h>
#include <alma/ASDM/StateRow.h>
#include <alma/ASDM/StationTable.h>
#include <alma/ASDM/StationRow.h>
#include <alma/ASDM/SubscanTable.h>
#include <alma/ASDM/SubscanRow.h>
#include <alma/ASDM/SwitchCycleTable.h>
#include <alma/ASDM/SwitchCycleRow.h>
#include <alma/ASDM/SysCalTable.h>
#include <alma/ASDM/SysCalRow.h>
#include <alma/ASDM/SysPowerTable.h>
#include <alma/ASDM/SysPowerRow.h>
#include <alma/ASDM/TotalPowerTable.h>
#include <alma/ASDM/TotalPowerRow.h>
#include <alma/ASDM/WVMCalTable.h>
#include <alma/ASDM/WVMCalRow.h>
#include <alma/ASDM/WeatherTable.h>
#include <alma/ASDM/WeatherRow.h>
using namespace casacore;
using namespace TimeSamplingMod;
using namespace BasebandNameMod;
using namespace AntennaMakeMod;
using namespace AntennaTypeMod;
using namespace AtmPhaseCorrectionMod;
using namespace ReceiverBandMod;
using namespace BasebandNameMod;
using namespace PolarizationTypeMod;
using namespace AtmPhaseCorrectionMod;
using namespace BasebandNameMod;
using namespace ReceiverBandMod;
using namespace PolarizationTypeMod;
using namespace BasebandNameMod;
using namespace ReceiverBandMod;
using namespace BasebandNameMod;
using namespace PolarizationTypeMod;
using namespace SyscalMethodMod;
using namespace BasebandNameMod;
using namespace NetSidebandMod;
using namespace AtmPhaseCorrectionMod;
using namespace CalCurveTypeMod;
using namespace ReceiverBandMod;
using namespace PolarizationTypeMod;
using namespace AtmPhaseCorrectionMod;
using namespace CalCurveTypeMod;
using namespace ReceiverBandMod;
using namespace PolarizationTypeMod;
using namespace CalDataOriginMod;
using namespace CalTypeMod;
using namespace AssociatedCalNatureMod;
using namespace ScanIntentMod;
using namespace AtmPhaseCorrectionMod;
using namespace BasebandNameMod;
using namespace ReceiverBandMod;
using namespace PolarizationTypeMod;
using namespace ReceiverSidebandMod;
using namespace CalibrationDeviceMod;
using namespace FluxCalibrationMethodMod;
using namespace StokesParameterMod;
using namespace DirectionReferenceCodeMod;
using namespace SourceModelMod;
using namespace ReceiverBandMod;
using namespace AtmPhaseCorrectionMod;
using namespace FocusMethodMod;
using namespace PolarizationTypeMod;
using namespace ReceiverBandMod;
using namespace PolarizationTypeMod;
using namespace AntennaMakeMod;
using namespace AntennaMakeMod;
using namespace PolarizationTypeMod;
using namespace ReceiverBandMod;
using namespace BasebandNameMod;
using namespace ReceiverBandMod;
using namespace AtmPhaseCorrectionMod;
using namespace PolarizationTypeMod;
using namespace ReceiverBandMod;
using namespace AntennaMakeMod;
using namespace AtmPhaseCorrectionMod;
using namespace PointingModelModeMod;
using namespace PointingMethodMod;
using namespace PolarizationTypeMod;
using namespace ReceiverBandMod;
using namespace AntennaMakeMod;
using namespace PointingModelModeMod;
using namespace PolarizationTypeMod;
using namespace AtmPhaseCorrectionMod;
using namespace PositionMethodMod;
using namespace ReceiverBandMod;
using namespace ReceiverBandMod;
using namespace AntennaMakeMod;
using namespace PolarizationTypeMod;
using namespace PrimaryBeamDescriptionMod;
using namespace InvalidatingConditionMod;
using namespace AtmPhaseCorrectionMod;
using namespace WVRMethodMod;
using namespace CorrelationModeMod;
using namespace AtmPhaseCorrectionMod;
using namespace ProcessorTypeMod;
using namespace SpectralResolutionTypeMod;
using namespace SpectralResolutionTypeMod;
using namespace BasebandNameMod;
using namespace AccumModeMod;
using namespace AxisNameMod;
using namespace FilterModeMod;
using namespace CorrelatorNameMod;
using namespace PolarizationTypeMod;
using namespace DifferenceTypeMod;
using namespace DifferenceTypeMod;
using namespace DopplerReferenceCodeMod;
using namespace PolarizationTypeMod;
using namespace DirectionReferenceCodeMod;
using namespace PolarizationTypeMod;
using namespace PolarizationTypeMod;
using namespace ReceiverBandMod;
using namespace PolarizationTypeMod;
using namespace HolographyChannelTypeMod;
using namespace DirectionReferenceCodeMod;
using namespace PolarizationTypeMod;
using namespace ReceiverBandMod;
using namespace StokesParameterMod;
using namespace PolarizationTypeMod;
using namespace ProcessorTypeMod;
using namespace ProcessorSubTypeMod;
using namespace ReceiverBandMod;
using namespace ReceiverSidebandMod;
using namespace NetSidebandMod;
using namespace ReceiverBandMod;
using namespace SBTypeMod;
using namespace DirectionReferenceCodeMod;
using namespace TimeScaleMod;
using namespace DataScaleMod;
using namespace DataScaleMod;
using namespace WeightTypeMod;
using namespace ScanIntentMod;
using namespace CalDataOriginMod;
using namespace CalibrationFunctionMod;
using namespace CalibrationSetMod;
using namespace AntennaMotionPatternMod;
using namespace DirectionReferenceCodeMod;
using namespace SourceModelMod;
using namespace FrequencyReferenceCodeMod;
using namespace StokesParameterMod;
using namespace RadialVelocityReferenceCodeMod;
using namespace BasebandNameMod;
using namespace NetSidebandMod;
using namespace SidebandProcessingModeMod;
using namespace WindowFunctionMod;
using namespace CorrelationBitMod;
using namespace FrequencyReferenceCodeMod;
using namespace SpectralResolutionTypeMod;
using namespace DetectorBandTypeMod;
using namespace CalibrationDeviceMod;
using namespace StationTypeMod;
using namespace SubscanIntentMod;
using namespace SwitchingModeMod;
using namespace CorrelatorCalibrationMod;
using namespace DirectionReferenceCodeMod;
using namespace WVRMethodMod;

/*
 * ALMA - Atacama Large Millimeter Array
 * (c) European Southern Observatory, 2002
 * (c) Associated Universities Inc., 2002
 * Copyright by ESO (in the framework of the ALMA collaboration),
 * Copyright by AUI (in the framework of the ALMA collaboration),
 * All rights reserved.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY, without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307  USA
 *
 * Warning!
 *  -------------------------------------------------------------------- 
 * | This is generated code!  Do not modify this file.                  |
 * | If you do, all changes will be lost when the file is re-generated. |
 *  --------------------------------------------------------------------
 *
 * File ASDMTables.cpp
 */


namespace casac {

    ASDM_MAIN::ASDM_MAIN() {
        name_ = "ASDM_MAIN";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Main table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("time", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("configDescriptionId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("fieldId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("timeSampling", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("interval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numIntegration", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("scanNumber", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("subscanNumber", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("dataSize", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("dataUID", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("stateId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockId", "blabla"));
  		
  		  		
    }

    ASDM_MAIN::~ASDM_MAIN() {
    }

    const TableDesc& ASDM_MAIN::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_MAIN::fill(const ASDM& asdm) {
        vector<MainRow*> rows = asdm.getMain().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<double> time(*table_p_, "time");             
  		
        ScalarColumn<String> configDescriptionId(*table_p_, "configDescriptionId");             
  		
        ScalarColumn<String> fieldId(*table_p_, "fieldId");             
  		
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ScalarColumn<String> timeSampling(*table_p_, "timeSampling");             
  		
        ScalarColumn<double> interval(*table_p_, "interval");             
  		
        ScalarColumn<int> numIntegration(*table_p_, "numIntegration");             
  		
        ScalarColumn<int> scanNumber(*table_p_, "scanNumber");             
  		
        ScalarColumn<int> subscanNumber(*table_p_, "subscanNumber");             
  		
        ScalarColumn<int> dataSize(*table_p_, "dataSize");             
  		
        ScalarColumn<String> dataUID(*table_p_, "dataUID");             
  		
        ArrayColumn<String> stateId(*table_p_, "stateId");             
  		
        ScalarColumn<String> execBlockId(*table_p_, "execBlockId");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            time.put(rowIndex, rows.at(i)->getTime().get()/(1.0e9));
	

	
            configDescriptionId.put(rowIndex, rows.at(i)->getConfigDescriptionId().toString());
	

	
            fieldId.put(rowIndex, rows.at(i)->getFieldId().toString());
	

		
	
            numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            timeSampling.put(rowIndex, CTimeSampling::name(rows.at(i)->getTimeSampling()));
	

	
            interval.put(rowIndex, rows.at(i)->getInterval().get()/(1.0e9));
	

	
            numIntegration.put(rowIndex, rows.at(i)->getNumIntegration());
	

	
            scanNumber.put(rowIndex, rows.at(i)->getScanNumber());
	

	
            subscanNumber.put(rowIndex, rows.at(i)->getSubscanNumber());
	

	
            dataSize.put(rowIndex, rows.at(i)->getDataSize());
	

	
            dataUID.put(rowIndex, rows.at(i)->getDataUID().toString());
	

	
            stateId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getStateId()));
	

	
            execBlockId.put(rowIndex, rows.at(i)->getExecBlockId().toString());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_ALMARADIOMETER::ASDM_ALMARADIOMETER() {
        name_ = "ASDM_ALMARADIOMETER";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset AlmaRadiometer table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("almaRadiometerId", "blabla"));
  		
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("spectralWindowId", "blabla"));
  		  		
    }

    ASDM_ALMARADIOMETER::~ASDM_ALMARADIOMETER() {
    }

    const TableDesc& ASDM_ALMARADIOMETER::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
	
    void ASDM_ALMARADIOMETER::fill(const ASDM& asdm) {
        vector<AlmaRadiometerRow*> rows = asdm.getAlmaRadiometer().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> almaRadiometerId(*table_p_, "almaRadiometerId");             
  		
  		
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ArrayColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            almaRadiometerId.put(rowIndex, rows.at(i)->getAlmaRadiometerId().toString());
	

		
		
	
            if (rows.at(i)->isNumAntennaExists())
                numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            if (rows.at(i)->isSpectralWindowIdExists())
                spectralWindowId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getSpectralWindowId()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_ANNOTATION::ASDM_ANNOTATION() {
        name_ = "ASDM_ANNOTATION";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Annotation table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("annotationId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("time", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("issue", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("details", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBaseband", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("interval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("dValue", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("vdValue", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("vvdValues", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("llValue", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("vllValue", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("vvllValue", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("antennaId", "blabla"));
  		  		
    }

    ASDM_ANNOTATION::~ASDM_ANNOTATION() {
    }

    const TableDesc& ASDM_ANNOTATION::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_ANNOTATION::fill(const ASDM& asdm) {
        vector<AnnotationRow*> rows = asdm.getAnnotation().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> annotationId(*table_p_, "annotationId");             
  		
  		
        ScalarColumn<double> time(*table_p_, "time");             
  		
        ScalarColumn<String> issue(*table_p_, "issue");             
  		
        ScalarColumn<String> details(*table_p_, "details");             
  		
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ArrayColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<int> numBaseband(*table_p_, "numBaseband");             
  		
        ScalarColumn<double> interval(*table_p_, "interval");             
  		
        ScalarColumn<double> dValue(*table_p_, "dValue");             
  		
        ArrayColumn<double> vdValue(*table_p_, "vdValue");             
  		
        ArrayColumn<double> vvdValues(*table_p_, "vvdValues");             
  		
        ScalarColumn<int> llValue(*table_p_, "llValue");             
  		
        ArrayColumn<int> vllValue(*table_p_, "vllValue");             
  		
        ArrayColumn<int> vvllValue(*table_p_, "vvllValue");             
  		
        ArrayColumn<String> antennaId(*table_p_, "antennaId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            annotationId.put(rowIndex, rows.at(i)->getAnnotationId().toString());
	

		
	
            time.put(rowIndex, rows.at(i)->getTime().get()/(1.0e9));
	

	
            issue.put(rowIndex, rows.at(i)->getIssue());
	

	
            details.put(rowIndex, rows.at(i)->getDetails());
	

		
	
            if (rows.at(i)->isNumAntennaExists())
                numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            if (rows.at(i)->isBasebandNameExists())
                basebandName.put(rowIndex, enum2CASA1D<BasebandName,CBasebandName>(rows.at(i)->getBasebandName()));
	

	
            if (rows.at(i)->isNumBasebandExists())
                numBaseband.put(rowIndex, rows.at(i)->getNumBaseband());
	

	
            if (rows.at(i)->isIntervalExists())
                interval.put(rowIndex, rows.at(i)->getInterval().get()/(1.0e9));
	

	
            if (rows.at(i)->isDValueExists())
                dValue.put(rowIndex, rows.at(i)->getDValue());
	

	
            if (rows.at(i)->isVdValueExists())
                vdValue.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getVdValue()));
	

	
            if (rows.at(i)->isVvdValuesExists())
                vvdValues.put(rowIndex, basic2CASA2D<double,double>(rows.at(i)->getVvdValues()));
	

	
            if (rows.at(i)->isLlValueExists())
                llValue.put(rowIndex, rows.at(i)->getLlValue());
	

	
            if (rows.at(i)->isVllValueExists())
                vllValue.put(rowIndex, basic2CASA1D<int64_t,int>(rows.at(i)->getVllValue()));
	

	
            if (rows.at(i)->isVvllValueExists())
                vvllValue.put(rowIndex, basic2CASA2D<int64_t,int>(rows.at(i)->getVvllValue()));
	

	
            if (rows.at(i)->isAntennaIdExists())
                antennaId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getAntennaId()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_ANTENNA::ASDM_ANTENNA() {
        name_ = "ASDM_ANTENNA";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Antenna table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("name", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaMake", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("dishDiameter", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("position", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("offset", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("time", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("stationId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("assocAntennaId", "blabla"));
  		  		
    }

    ASDM_ANTENNA::~ASDM_ANTENNA() {
    }

    const TableDesc& ASDM_ANTENNA::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_ANTENNA::fill(const ASDM& asdm) {
        vector<AntennaRow*> rows = asdm.getAntenna().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
  		
        ScalarColumn<String> name(*table_p_, "name");             
  		
        ScalarColumn<String> antennaMake(*table_p_, "antennaMake");             
  		
        ScalarColumn<String> antennaType(*table_p_, "antennaType");             
  		
        ScalarColumn<double> dishDiameter(*table_p_, "dishDiameter");             
  		
        ArrayColumn<double> position(*table_p_, "position");             
  		
        ArrayColumn<double> offset(*table_p_, "offset");             
  		
        ScalarColumn<double> time(*table_p_, "time");             
  		
        ScalarColumn<String> stationId(*table_p_, "stationId");             
  		
  		
        ScalarColumn<String> assocAntennaId(*table_p_, "assocAntennaId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

		
	
            name.put(rowIndex, rows.at(i)->getName());
	

	
            antennaMake.put(rowIndex, CAntennaMake::name(rows.at(i)->getAntennaMake()));
	

	
            antennaType.put(rowIndex, CAntennaType::name(rows.at(i)->getAntennaType()));
	

	
            dishDiameter.put(rowIndex, rows.at(i)->getDishDiameter().get());
	

	
            position.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getPosition()));
	

	
            offset.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getOffset()));
	

	
            time.put(rowIndex, rows.at(i)->getTime().get()/(1.0e9));
	

	
            stationId.put(rowIndex, rows.at(i)->getStationId().toString());
	

		
	
            if (rows.at(i)->isAssocAntennaIdExists())
                assocAntennaId.put(rowIndex, rows.at(i)->getAssocAntennaId().toString());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALAMPLI::ASDM_CALAMPLI() {
        name_ = "ASDM_CALAMPLI";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalAmpli table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("apertureEfficiency", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("apertureEfficiencyError", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("correctionValidity", "blabla"));
  		  		
    }

    ASDM_CALAMPLI::~ASDM_CALAMPLI() {
    }

    const TableDesc& ASDM_CALAMPLI::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
			
		
			
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALAMPLI::fill(const ASDM& asdm) {
        vector<CalAmpliRow*> rows = asdm.getCalAmpli().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ArrayColumn<float> apertureEfficiency(*table_p_, "apertureEfficiency");             
  		
        ArrayColumn<float> apertureEfficiencyError(*table_p_, "apertureEfficiencyError");             
  		
  		
        ScalarColumn<bool> correctionValidity(*table_p_, "correctionValidity");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            basebandName.put(rowIndex, CBasebandName::name(rows.at(i)->getBasebandName()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            apertureEfficiency.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getApertureEfficiency()));
	

	
            apertureEfficiencyError.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getApertureEfficiencyError()));
	

		
	
            if (rows.at(i)->isCorrectionValidityExists())
                correctionValidity.put(rowIndex, rows.at(i)->getCorrectionValidity());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALANTENNASOLUTIONS::ASDM_CALANTENNASOLUTIONS() {
        name_ = "ASDM_CALANTENNASOLUTIONS";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalAntennaSolutions table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("refAntennaName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("direction", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("integrationTime", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("correctionValidity", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phaseAnt", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phaseAntRMS", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("amplitudeAnt", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("amplitudeAntRMS", "blabla"));
  		
  		  		
    }

    ASDM_CALANTENNASOLUTIONS::~ASDM_CALANTENNASOLUTIONS() {
    }

    const TableDesc& ASDM_CALANTENNASOLUTIONS::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALANTENNASOLUTIONS::fill(const ASDM& asdm) {
        vector<CalAntennaSolutionsRow*> rows = asdm.getCalAntennaSolutions().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ScalarColumn<String> refAntennaName(*table_p_, "refAntennaName");             
  		
        ArrayColumn<double> direction(*table_p_, "direction");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<double> integrationTime(*table_p_, "integrationTime");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ScalarColumn<bool> correctionValidity(*table_p_, "correctionValidity");             
  		
        ArrayColumn<float> phaseAnt(*table_p_, "phaseAnt");             
  		
        ArrayColumn<float> phaseAntRMS(*table_p_, "phaseAntRMS");             
  		
        ArrayColumn<float> amplitudeAnt(*table_p_, "amplitudeAnt");             
  		
        ArrayColumn<float> amplitudeAntRMS(*table_p_, "amplitudeAntRMS");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            basebandName.put(rowIndex, CBasebandName::name(rows.at(i)->getBasebandName()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            refAntennaName.put(rowIndex, rows.at(i)->getRefAntennaName());
	

	
            direction.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getDirection()));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            integrationTime.put(rowIndex, rows.at(i)->getIntegrationTime().get()/(1.0e9));
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            correctionValidity.put(rowIndex, rows.at(i)->getCorrectionValidity());
	

	
            phaseAnt.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getPhaseAnt()));
	

	
            phaseAntRMS.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getPhaseAntRMS()));
	

	
            amplitudeAnt.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getAmplitudeAnt()));
	

	
            amplitudeAntRMS.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getAmplitudeAntRMS()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALAPPPHASE::ASDM_CALAPPPHASE() {
        name_ = "ASDM_CALAPPPHASE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalAppPhase table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("scanNumber", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("adjustTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("adjustToken", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("phasingMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPhasedAntennas", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("phasedAntennas", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("refAntennaIndex", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("candRefAntennaIndex", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("phasePacking", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptors", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numChannels", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPhaseValues", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phaseValues", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numCompare", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numEfficiencies", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("compareArray", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("efficiencyIndices", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("efficiencies", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("quality", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("phasedSumAntenna", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("typeSupports", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numSupports", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phaseSupports", "blabla"));
  		  		
    }

    ASDM_CALAPPPHASE::~ASDM_CALAPPPHASE() {
    }

    const TableDesc& ASDM_CALAPPPHASE::tableDesc() const {
        return tableDesc_;
    }


		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALAPPPHASE::fill(const ASDM& asdm) {
        vector<CalAppPhaseRow*> rows = asdm.getCalAppPhase().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<int> scanNumber(*table_p_, "scanNumber");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<double> adjustTime(*table_p_, "adjustTime");             
  		
        ScalarColumn<String> adjustToken(*table_p_, "adjustToken");             
  		
        ScalarColumn<String> phasingMode(*table_p_, "phasingMode");             
  		
        ScalarColumn<int> numPhasedAntennas(*table_p_, "numPhasedAntennas");             
  		
        ArrayColumn<String> phasedAntennas(*table_p_, "phasedAntennas");             
  		
        ScalarColumn<int> refAntennaIndex(*table_p_, "refAntennaIndex");             
  		
        ScalarColumn<int> candRefAntennaIndex(*table_p_, "candRefAntennaIndex");             
  		
        ScalarColumn<String> phasePacking(*table_p_, "phasePacking");             
  		
        ScalarColumn<int> numReceptors(*table_p_, "numReceptors");             
  		
        ScalarColumn<int> numChannels(*table_p_, "numChannels");             
  		
        ScalarColumn<int> numPhaseValues(*table_p_, "numPhaseValues");             
  		
        ArrayColumn<float> phaseValues(*table_p_, "phaseValues");             
  		
        ScalarColumn<int> numCompare(*table_p_, "numCompare");             
  		
        ScalarColumn<int> numEfficiencies(*table_p_, "numEfficiencies");             
  		
        ArrayColumn<String> compareArray(*table_p_, "compareArray");             
  		
        ArrayColumn<int> efficiencyIndices(*table_p_, "efficiencyIndices");             
  		
        ArrayColumn<float> efficiencies(*table_p_, "efficiencies");             
  		
        ArrayColumn<float> quality(*table_p_, "quality");             
  		
        ScalarColumn<String> phasedSumAntenna(*table_p_, "phasedSumAntenna");             
  		
  		
        ScalarColumn<String> typeSupports(*table_p_, "typeSupports");             
  		
        ScalarColumn<int> numSupports(*table_p_, "numSupports");             
  		
        ArrayColumn<float> phaseSupports(*table_p_, "phaseSupports");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            basebandName.put(rowIndex, CBasebandName::name(rows.at(i)->getBasebandName()));
	

	
            scanNumber.put(rowIndex, rows.at(i)->getScanNumber());
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            adjustTime.put(rowIndex, rows.at(i)->getAdjustTime().get()/(1.0e9));
	

	
            adjustToken.put(rowIndex, rows.at(i)->getAdjustToken());
	

	
            phasingMode.put(rowIndex, rows.at(i)->getPhasingMode());
	

	
            numPhasedAntennas.put(rowIndex, rows.at(i)->getNumPhasedAntennas());
	

	
            phasedAntennas.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getPhasedAntennas()));
	

	
            refAntennaIndex.put(rowIndex, rows.at(i)->getRefAntennaIndex());
	

	
            candRefAntennaIndex.put(rowIndex, rows.at(i)->getCandRefAntennaIndex());
	

	
            phasePacking.put(rowIndex, rows.at(i)->getPhasePacking());
	

	
            numReceptors.put(rowIndex, rows.at(i)->getNumReceptors());
	

	
            numChannels.put(rowIndex, rows.at(i)->getNumChannels());
	

	
            numPhaseValues.put(rowIndex, rows.at(i)->getNumPhaseValues());
	

	
            phaseValues.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getPhaseValues()));
	

	
            numCompare.put(rowIndex, rows.at(i)->getNumCompare());
	

	
            numEfficiencies.put(rowIndex, rows.at(i)->getNumEfficiencies());
	

	
            compareArray.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCompareArray()));
	

	
            efficiencyIndices.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getEfficiencyIndices()));
	

	
            efficiencies.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getEfficiencies()));
	

	
            quality.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getQuality()));
	

	
            phasedSumAntenna.put(rowIndex, rows.at(i)->getPhasedSumAntenna());
	

		
	
            if (rows.at(i)->isTypeSupportsExists())
                typeSupports.put(rowIndex, rows.at(i)->getTypeSupports());
	

	
            if (rows.at(i)->isNumSupportsExists())
                numSupports.put(rowIndex, rows.at(i)->getNumSupports());
	

	
            if (rows.at(i)->isPhaseSupportsExists())
                phaseSupports.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getPhaseSupports()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALATMOSPHERE::ASDM_CALATMOSPHERE() {
        name_ = "ASDM_CALATMOSPHERE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalAtmosphere table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numFreq", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numLoad", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("forwardEffSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("groundPressure", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("groundRelHumidity", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencySpectrum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("groundTemperature", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("powerSkySpectrum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("powerLoadSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("syscalType", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tAtmSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tRecSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tSysSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("tauSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tAtm", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tRec", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tSys", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("tau", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("water", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("waterError", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("alphaSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("forwardEfficiency", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("forwardEfficiencyError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("sbGain", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("sbGainError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("sbGainSpectrum", "blabla"));
  		  		
    }

    ASDM_CALATMOSPHERE::~ASDM_CALATMOSPHERE() {
    }

    const TableDesc& ASDM_CALATMOSPHERE::tableDesc() const {
        return tableDesc_;
    }


		
			
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALATMOSPHERE::fill(const ASDM& asdm) {
        vector<CalAtmosphereRow*> rows = asdm.getCalAtmosphere().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<int> numFreq(*table_p_, "numFreq");             
  		
        ScalarColumn<int> numLoad(*table_p_, "numLoad");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<float> forwardEffSpectrum(*table_p_, "forwardEffSpectrum");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<double> groundPressure(*table_p_, "groundPressure");             
  		
        ScalarColumn<double> groundRelHumidity(*table_p_, "groundRelHumidity");             
  		
        ArrayColumn<double> frequencySpectrum(*table_p_, "frequencySpectrum");             
  		
        ScalarColumn<double> groundTemperature(*table_p_, "groundTemperature");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<float> powerSkySpectrum(*table_p_, "powerSkySpectrum");             
  		
        ArrayColumn<float> powerLoadSpectrum(*table_p_, "powerLoadSpectrum");             
  		
        ScalarColumn<String> syscalType(*table_p_, "syscalType");             
  		
        ArrayColumn<double> tAtmSpectrum(*table_p_, "tAtmSpectrum");             
  		
        ArrayColumn<double> tRecSpectrum(*table_p_, "tRecSpectrum");             
  		
        ArrayColumn<double> tSysSpectrum(*table_p_, "tSysSpectrum");             
  		
        ArrayColumn<float> tauSpectrum(*table_p_, "tauSpectrum");             
  		
        ArrayColumn<double> tAtm(*table_p_, "tAtm");             
  		
        ArrayColumn<double> tRec(*table_p_, "tRec");             
  		
        ArrayColumn<double> tSys(*table_p_, "tSys");             
  		
        ArrayColumn<float> tau(*table_p_, "tau");             
  		
        ArrayColumn<double> water(*table_p_, "water");             
  		
        ArrayColumn<double> waterError(*table_p_, "waterError");             
  		
  		
        ArrayColumn<float> alphaSpectrum(*table_p_, "alphaSpectrum");             
  		
        ArrayColumn<float> forwardEfficiency(*table_p_, "forwardEfficiency");             
  		
        ArrayColumn<double> forwardEfficiencyError(*table_p_, "forwardEfficiencyError");             
  		
        ArrayColumn<float> sbGain(*table_p_, "sbGain");             
  		
        ArrayColumn<float> sbGainError(*table_p_, "sbGainError");             
  		
        ArrayColumn<float> sbGainSpectrum(*table_p_, "sbGainSpectrum");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            basebandName.put(rowIndex, CBasebandName::name(rows.at(i)->getBasebandName()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            numFreq.put(rowIndex, rows.at(i)->getNumFreq());
	

	
            numLoad.put(rowIndex, rows.at(i)->getNumLoad());
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            forwardEffSpectrum.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getForwardEffSpectrum()));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            groundPressure.put(rowIndex, rows.at(i)->getGroundPressure().get());
	

	
            groundRelHumidity.put(rowIndex, rows.at(i)->getGroundRelHumidity().get());
	

	
            frequencySpectrum.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencySpectrum()));
	

	
            groundTemperature.put(rowIndex, rows.at(i)->getGroundTemperature().get());
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            powerSkySpectrum.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getPowerSkySpectrum()));
	

	
            powerLoadSpectrum.put(rowIndex, basic2CASA3D<float,float>(rows.at(i)->getPowerLoadSpectrum()));
	

	
            syscalType.put(rowIndex, CSyscalMethod::name(rows.at(i)->getSyscalType()));
	

	
            tAtmSpectrum.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getTAtmSpectrum()));
	

	
            tRecSpectrum.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getTRecSpectrum()));
	

	
            tSysSpectrum.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getTSysSpectrum()));
	

	
            tauSpectrum.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getTauSpectrum()));
	

	
            tAtm.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getTAtm()));
	

	
            tRec.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getTRec()));
	

	
            tSys.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getTSys()));
	

	
            tau.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getTau()));
	

	
            water.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getWater()));
	

	
            waterError.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getWaterError()));
	

		
	
            if (rows.at(i)->isAlphaSpectrumExists())
                alphaSpectrum.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getAlphaSpectrum()));
	

	
            if (rows.at(i)->isForwardEfficiencyExists())
                forwardEfficiency.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getForwardEfficiency()));
	

	
            if (rows.at(i)->isForwardEfficiencyErrorExists())
                forwardEfficiencyError.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getForwardEfficiencyError()));
	

	
            if (rows.at(i)->isSbGainExists())
                sbGain.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getSbGain()));
	

	
            if (rows.at(i)->isSbGainErrorExists())
                sbGainError.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getSbGainError()));
	

	
            if (rows.at(i)->isSbGainSpectrumExists())
                sbGainSpectrum.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getSbGainSpectrum()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALBANDPASS::ASDM_CALBANDPASS() {
        name_ = "ASDM_CALBANDPASS";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalBandpass table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sideband", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("typeCurve", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPoly", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("antennaNames", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("refAntennaName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("freqLimits", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("curve", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("reducedChiSquared", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBaseline", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numFreq", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("rms", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numSpectralWindow", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("chanFreqStart", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("chanFreqStep", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("numSpectralWindowChan", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("spectrum", "blabla"));
  		  		
    }

    ASDM_CALBANDPASS::~ASDM_CALBANDPASS() {
    }

    const TableDesc& ASDM_CALBANDPASS::tableDesc() const {
        return tableDesc_;
    }


		
			
			
		
			
			
		
			
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALBANDPASS::fill(const ASDM& asdm) {
        vector<CalBandpassRow*> rows = asdm.getCalBandpass().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<String> sideband(*table_p_, "sideband");             
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> typeCurve(*table_p_, "typeCurve");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ScalarColumn<int> numPoly(*table_p_, "numPoly");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<String> antennaNames(*table_p_, "antennaNames");             
  		
        ScalarColumn<String> refAntennaName(*table_p_, "refAntennaName");             
  		
        ArrayColumn<double> freqLimits(*table_p_, "freqLimits");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<float> curve(*table_p_, "curve");             
  		
        ArrayColumn<double> reducedChiSquared(*table_p_, "reducedChiSquared");             
  		
  		
        ScalarColumn<int> numBaseline(*table_p_, "numBaseline");             
  		
        ScalarColumn<int> numFreq(*table_p_, "numFreq");             
  		
        ArrayColumn<float> rms(*table_p_, "rms");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<int> numSpectralWindow(*table_p_, "numSpectralWindow");             
  		
        ArrayColumn<double> chanFreqStart(*table_p_, "chanFreqStart");             
  		
        ArrayColumn<double> chanFreqStep(*table_p_, "chanFreqStep");             
  		
        ArrayColumn<int> numSpectralWindowChan(*table_p_, "numSpectralWindowChan");             
  		
        ArrayColumn<float> spectrum(*table_p_, "spectrum");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            basebandName.put(rowIndex, CBasebandName::name(rows.at(i)->getBasebandName()));
	

	
            sideband.put(rowIndex, CNetSideband::name(rows.at(i)->getSideband()));
	

	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            typeCurve.put(rowIndex, CCalCurveType::name(rows.at(i)->getTypeCurve()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            numPoly.put(rowIndex, rows.at(i)->getNumPoly());
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            antennaNames.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getAntennaNames()));
	

	
            refAntennaName.put(rowIndex, rows.at(i)->getRefAntennaName());
	

	
            freqLimits.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFreqLimits()));
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            curve.put(rowIndex, basic2CASA3D<float,float>(rows.at(i)->getCurve()));
	

	
            reducedChiSquared.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getReducedChiSquared()));
	

		
	
            if (rows.at(i)->isNumBaselineExists())
                numBaseline.put(rowIndex, rows.at(i)->getNumBaseline());
	

	
            if (rows.at(i)->isNumFreqExists())
                numFreq.put(rowIndex, rows.at(i)->getNumFreq());
	

	
            if (rows.at(i)->isRmsExists())
                rms.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getRms()));
	

	
            if (rows.at(i)->isFrequencyRangeExists())
                frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            if (rows.at(i)->isNumSpectralWindowExists())
                numSpectralWindow.put(rowIndex, rows.at(i)->getNumSpectralWindow());
	

	
            if (rows.at(i)->isChanFreqStartExists())
                chanFreqStart.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getChanFreqStart()));
	

	
            if (rows.at(i)->isChanFreqStepExists())
                chanFreqStep.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getChanFreqStep()));
	

	
            if (rows.at(i)->isNumSpectralWindowChanExists())
                numSpectralWindowChan.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getNumSpectralWindowChan()));
	

	
            if (rows.at(i)->isSpectrumExists())
                spectrum.put(rowIndex, basic2CASA3D<float,float>(rows.at(i)->getSpectrum()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALCURVE::ASDM_CALCURVE() {
        name_ = "ASDM_CALCURVE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalCurve table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("typeCurve", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPoly", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("antennaNames", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("refAntennaName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("curve", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("reducedChiSquared", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBaseline", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("rms", "blabla"));
  		  		
    }

    ASDM_CALCURVE::~ASDM_CALCURVE() {
    }

    const TableDesc& ASDM_CALCURVE::tableDesc() const {
        return tableDesc_;
    }


		
			
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALCURVE::fill(const ASDM& asdm) {
        vector<CalCurveRow*> rows = asdm.getCalCurve().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> typeCurve(*table_p_, "typeCurve");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ScalarColumn<int> numPoly(*table_p_, "numPoly");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<String> antennaNames(*table_p_, "antennaNames");             
  		
        ScalarColumn<String> refAntennaName(*table_p_, "refAntennaName");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<float> curve(*table_p_, "curve");             
  		
        ArrayColumn<double> reducedChiSquared(*table_p_, "reducedChiSquared");             
  		
  		
        ScalarColumn<int> numBaseline(*table_p_, "numBaseline");             
  		
        ArrayColumn<float> rms(*table_p_, "rms");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            typeCurve.put(rowIndex, CCalCurveType::name(rows.at(i)->getTypeCurve()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            numPoly.put(rowIndex, rows.at(i)->getNumPoly());
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            antennaNames.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getAntennaNames()));
	

	
            refAntennaName.put(rowIndex, rows.at(i)->getRefAntennaName());
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            curve.put(rowIndex, basic2CASA3D<float,float>(rows.at(i)->getCurve()));
	

	
            reducedChiSquared.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getReducedChiSquared()));
	

		
	
            if (rows.at(i)->isNumBaselineExists())
                numBaseline.put(rowIndex, rows.at(i)->getNumBaseline());
	

	
            if (rows.at(i)->isRmsExists())
                rms.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getRms()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALDATA::ASDM_CALDATA() {
        name_ = "ASDM_CALDATA";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalData table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startTimeObserved", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endTimeObserved", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numScan", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("scanSet", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("assocCalDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("assocCalNature", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("fieldName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("sourceName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("sourceCode", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("scanIntent", "blabla"));
  		  		
    }

    ASDM_CALDATA::~ASDM_CALDATA() {
    }

    const TableDesc& ASDM_CALDATA::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
			
		
	
    void ASDM_CALDATA::fill(const ASDM& asdm) {
        vector<CalDataRow*> rows = asdm.getCalData().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
  		
        ScalarColumn<double> startTimeObserved(*table_p_, "startTimeObserved");             
  		
        ScalarColumn<double> endTimeObserved(*table_p_, "endTimeObserved");             
  		
        ScalarColumn<String> execBlockUID(*table_p_, "execBlockUID");             
  		
        ScalarColumn<String> calDataType(*table_p_, "calDataType");             
  		
        ScalarColumn<String> calType(*table_p_, "calType");             
  		
        ScalarColumn<int> numScan(*table_p_, "numScan");             
  		
        ArrayColumn<int> scanSet(*table_p_, "scanSet");             
  		
  		
        ScalarColumn<String> assocCalDataId(*table_p_, "assocCalDataId");             
  		
        ScalarColumn<String> assocCalNature(*table_p_, "assocCalNature");             
  		
        ArrayColumn<String> fieldName(*table_p_, "fieldName");             
  		
        ArrayColumn<String> sourceName(*table_p_, "sourceName");             
  		
        ArrayColumn<String> sourceCode(*table_p_, "sourceCode");             
  		
        ArrayColumn<String> scanIntent(*table_p_, "scanIntent");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

		
	
            startTimeObserved.put(rowIndex, rows.at(i)->getStartTimeObserved().get()/(1.0e9));
	

	
            endTimeObserved.put(rowIndex, rows.at(i)->getEndTimeObserved().get()/(1.0e9));
	

	
            execBlockUID.put(rowIndex, rows.at(i)->getExecBlockUID().toString());
	

	
            calDataType.put(rowIndex, CCalDataOrigin::name(rows.at(i)->getCalDataType()));
	

	
            calType.put(rowIndex, CCalType::name(rows.at(i)->getCalType()));
	

	
            numScan.put(rowIndex, rows.at(i)->getNumScan());
	

	
            scanSet.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getScanSet()));
	

		
	
            if (rows.at(i)->isAssocCalDataIdExists())
                assocCalDataId.put(rowIndex, rows.at(i)->getAssocCalDataId().toString());
	

	
            if (rows.at(i)->isAssocCalNatureExists())
                assocCalNature.put(rowIndex, CAssociatedCalNature::name(rows.at(i)->getAssocCalNature()));
	

	
            if (rows.at(i)->isFieldNameExists())
                fieldName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getFieldName()));
	

	
            if (rows.at(i)->isSourceNameExists())
                sourceName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getSourceName()));
	

	
            if (rows.at(i)->isSourceCodeExists())
                sourceCode.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getSourceCode()));
	

	
            if (rows.at(i)->isScanIntentExists())
                scanIntent.put(rowIndex, enum2CASA1D<ScanIntent,CScanIntent>(rows.at(i)->getScanIntent()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALDELAY::ASDM_CALDELAY() {
        name_ = "ASDM_CALDELAY";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalDelay table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("refAntennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("delayError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("delayOffset", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("reducedChiSquared", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("appliedDelay", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("crossDelayOffset", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("crossDelayOffsetError", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numSideband", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("refFreq", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("refFreqPhase", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("sidebands", "blabla"));
  		  		
    }

    ASDM_CALDELAY::~ASDM_CALDELAY() {
    }

    const TableDesc& ASDM_CALDELAY::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
	
    void ASDM_CALDELAY::fill(const ASDM& asdm) {
        vector<CalDelayRow*> rows = asdm.getCalDelay().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<String> refAntennaName(*table_p_, "refAntennaName");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<double> delayError(*table_p_, "delayError");             
  		
        ArrayColumn<double> delayOffset(*table_p_, "delayOffset");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<double> reducedChiSquared(*table_p_, "reducedChiSquared");             
  		
        ArrayColumn<double> appliedDelay(*table_p_, "appliedDelay");             
  		
  		
        ScalarColumn<double> crossDelayOffset(*table_p_, "crossDelayOffset");             
  		
        ScalarColumn<double> crossDelayOffsetError(*table_p_, "crossDelayOffsetError");             
  		
        ScalarColumn<int> numSideband(*table_p_, "numSideband");             
  		
        ArrayColumn<double> refFreq(*table_p_, "refFreq");             
  		
        ArrayColumn<double> refFreqPhase(*table_p_, "refFreqPhase");             
  		
        ArrayColumn<String> sidebands(*table_p_, "sidebands");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            basebandName.put(rowIndex, CBasebandName::name(rows.at(i)->getBasebandName()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            refAntennaName.put(rowIndex, rows.at(i)->getRefAntennaName());
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            delayError.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getDelayError()));
	

	
            delayOffset.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getDelayOffset()));
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            reducedChiSquared.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getReducedChiSquared()));
	

	
            appliedDelay.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getAppliedDelay()));
	

		
	
            if (rows.at(i)->isCrossDelayOffsetExists())
                crossDelayOffset.put(rowIndex, rows.at(i)->getCrossDelayOffset());
	

	
            if (rows.at(i)->isCrossDelayOffsetErrorExists())
                crossDelayOffsetError.put(rowIndex, rows.at(i)->getCrossDelayOffsetError());
	

	
            if (rows.at(i)->isNumSidebandExists())
                numSideband.put(rowIndex, rows.at(i)->getNumSideband());
	

	
            if (rows.at(i)->isRefFreqExists())
                refFreq.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getRefFreq()));
	

	
            if (rows.at(i)->isRefFreqPhaseExists())
                refFreqPhase.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getRefFreqPhase()));
	

	
            if (rows.at(i)->isSidebandsExists())
                sidebands.put(rowIndex, enum2CASA1D<ReceiverSideband,CReceiverSideband>(rows.at(i)->getSidebands()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALDEVICE::ASDM_CALDEVICE() {
        name_ = "ASDM_CALDEVICE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalDevice table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("feedId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numCalload", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("calLoadNames", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("calEff", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("noiseCal", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("coupledNoiseCal", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("temperatureLoad", "blabla"));
  		  		
    }

    ASDM_CALDEVICE::~ASDM_CALDEVICE() {
    }

    const TableDesc& ASDM_CALDEVICE::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALDEVICE::fill(const ASDM& asdm) {
        vector<CalDeviceRow*> rows = asdm.getCalDevice().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
        ScalarColumn<int> feedId(*table_p_, "feedId");             
  		
  		
        ScalarColumn<int> numCalload(*table_p_, "numCalload");             
  		
        ArrayColumn<String> calLoadNames(*table_p_, "calLoadNames");             
  		
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<float> calEff(*table_p_, "calEff");             
  		
        ArrayColumn<double> noiseCal(*table_p_, "noiseCal");             
  		
        ArrayColumn<float> coupledNoiseCal(*table_p_, "coupledNoiseCal");             
  		
        ArrayColumn<double> temperatureLoad(*table_p_, "temperatureLoad");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

	
            feedId.put(rowIndex, rows.at(i)->getFeedId());
	

		
	
            numCalload.put(rowIndex, rows.at(i)->getNumCalload());
	

	
            calLoadNames.put(rowIndex, enum2CASA1D<CalibrationDevice,CCalibrationDevice>(rows.at(i)->getCalLoadNames()));
	

		
	
            if (rows.at(i)->isNumReceptorExists())
                numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            if (rows.at(i)->isCalEffExists())
                calEff.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getCalEff()));
	

	
            if (rows.at(i)->isNoiseCalExists())
                noiseCal.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getNoiseCal()));
	

	
            if (rows.at(i)->isCoupledNoiseCalExists())
                coupledNoiseCal.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getCoupledNoiseCal()));
	

	
            if (rows.at(i)->isTemperatureLoadExists())
                temperatureLoad.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getTemperatureLoad()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALFLUX::ASDM_CALFLUX() {
        name_ = "ASDM_CALFLUX";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalFlux table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sourceName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numFrequencyRanges", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numStokes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRanges", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("fluxMethod", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("flux", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("fluxError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("stokes", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("direction", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("directionCode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("directionEquinox", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("PA", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("PAError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("size", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("sizeError", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sourceModel", "blabla"));
  		  		
    }

    ASDM_CALFLUX::~ASDM_CALFLUX() {
    }

    const TableDesc& ASDM_CALFLUX::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
	
    void ASDM_CALFLUX::fill(const ASDM& asdm) {
        vector<CalFluxRow*> rows = asdm.getCalFlux().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> sourceName(*table_p_, "sourceName");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<int> numFrequencyRanges(*table_p_, "numFrequencyRanges");             
  		
        ScalarColumn<int> numStokes(*table_p_, "numStokes");             
  		
        ArrayColumn<double> frequencyRanges(*table_p_, "frequencyRanges");             
  		
        ScalarColumn<String> fluxMethod(*table_p_, "fluxMethod");             
  		
        ArrayColumn<double> flux(*table_p_, "flux");             
  		
        ArrayColumn<double> fluxError(*table_p_, "fluxError");             
  		
        ArrayColumn<String> stokes(*table_p_, "stokes");             
  		
  		
        ArrayColumn<double> direction(*table_p_, "direction");             
  		
        ScalarColumn<String> directionCode(*table_p_, "directionCode");             
  		
        ScalarColumn<double> directionEquinox(*table_p_, "directionEquinox");             
  		
        ArrayColumn<double> PA(*table_p_, "PA");             
  		
        ArrayColumn<double> PAError(*table_p_, "PAError");             
  		
        ArrayColumn<double> size(*table_p_, "size");             
  		
        ArrayColumn<double> sizeError(*table_p_, "sizeError");             
  		
        ScalarColumn<String> sourceModel(*table_p_, "sourceModel");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            sourceName.put(rowIndex, rows.at(i)->getSourceName());
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            numFrequencyRanges.put(rowIndex, rows.at(i)->getNumFrequencyRanges());
	

	
            numStokes.put(rowIndex, rows.at(i)->getNumStokes());
	

	
            frequencyRanges.put(rowIndex, ext2CASA2D<Frequency,double>(rows.at(i)->getFrequencyRanges()));
	

	
            fluxMethod.put(rowIndex, CFluxCalibrationMethod::name(rows.at(i)->getFluxMethod()));
	

	
            flux.put(rowIndex, basic2CASA2D<double,double>(rows.at(i)->getFlux()));
	

	
            fluxError.put(rowIndex, basic2CASA2D<double,double>(rows.at(i)->getFluxError()));
	

	
            stokes.put(rowIndex, enum2CASA1D<StokesParameter,CStokesParameter>(rows.at(i)->getStokes()));
	

		
	
            if (rows.at(i)->isDirectionExists())
                direction.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getDirection()));
	

	
            if (rows.at(i)->isDirectionCodeExists())
                directionCode.put(rowIndex, CDirectionReferenceCode::name(rows.at(i)->getDirectionCode()));
	

	
            if (rows.at(i)->isDirectionEquinoxExists())
                directionEquinox.put(rowIndex, rows.at(i)->getDirectionEquinox().get());
	

	
            if (rows.at(i)->isPAExists())
                PA.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getPA()));
	

	
            if (rows.at(i)->isPAErrorExists())
                PAError.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getPAError()));
	

	
            if (rows.at(i)->isSizeExists())
                size.put(rowIndex, ext2CASA3D<Angle,double>(rows.at(i)->getSize()));
	

	
            if (rows.at(i)->isSizeErrorExists())
                sizeError.put(rowIndex, ext2CASA3D<Angle,double>(rows.at(i)->getSizeError()));
	

	
            if (rows.at(i)->isSourceModelExists())
                sourceModel.put(rowIndex, CSourceModel::name(rows.at(i)->getSourceModel()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALFOCUS::ASDM_CALFOCUS() {
        name_ = "ASDM_CALFOCUS";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalFocus table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("ambientTemperature", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("focusMethod", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("pointingDirection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("wereFixed", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("offset", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("offsetError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("offsetWasTied", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("reducedChiSquared", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("position", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("polarizationsAveraged", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("focusCurveWidth", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("focusCurveWidthError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("focusCurveWasFixed", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("offIntensity", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("offIntensityError", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("offIntensityWasFixed", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("peakIntensity", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("peakIntensityError", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("peakIntensityWasFixed", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("astigmPlus", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("astigmPlusError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("astigmMult", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("astigmMultError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("illumOffset", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("illumOffsetError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("fitRMS", "blabla"));
  		  		
    }

    ASDM_CALFOCUS::~ASDM_CALFOCUS() {
    }

    const TableDesc& ASDM_CALFOCUS::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
			
		
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALFOCUS::fill(const ASDM& asdm) {
        vector<CalFocusRow*> rows = asdm.getCalFocus().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<double> ambientTemperature(*table_p_, "ambientTemperature");             
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> focusMethod(*table_p_, "focusMethod");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ArrayColumn<double> pointingDirection(*table_p_, "pointingDirection");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<bool> wereFixed(*table_p_, "wereFixed");             
  		
        ArrayColumn<double> offset(*table_p_, "offset");             
  		
        ArrayColumn<double> offsetError(*table_p_, "offsetError");             
  		
        ArrayColumn<bool> offsetWasTied(*table_p_, "offsetWasTied");             
  		
        ArrayColumn<double> reducedChiSquared(*table_p_, "reducedChiSquared");             
  		
        ArrayColumn<double> position(*table_p_, "position");             
  		
  		
        ScalarColumn<bool> polarizationsAveraged(*table_p_, "polarizationsAveraged");             
  		
        ArrayColumn<double> focusCurveWidth(*table_p_, "focusCurveWidth");             
  		
        ArrayColumn<double> focusCurveWidthError(*table_p_, "focusCurveWidthError");             
  		
        ArrayColumn<bool> focusCurveWasFixed(*table_p_, "focusCurveWasFixed");             
  		
        ArrayColumn<double> offIntensity(*table_p_, "offIntensity");             
  		
        ArrayColumn<double> offIntensityError(*table_p_, "offIntensityError");             
  		
        ScalarColumn<bool> offIntensityWasFixed(*table_p_, "offIntensityWasFixed");             
  		
        ArrayColumn<double> peakIntensity(*table_p_, "peakIntensity");             
  		
        ArrayColumn<double> peakIntensityError(*table_p_, "peakIntensityError");             
  		
        ScalarColumn<bool> peakIntensityWasFixed(*table_p_, "peakIntensityWasFixed");             
  		
        ArrayColumn<double> astigmPlus(*table_p_, "astigmPlus");             
  		
        ArrayColumn<double> astigmPlusError(*table_p_, "astigmPlusError");             
  		
        ArrayColumn<double> astigmMult(*table_p_, "astigmMult");             
  		
        ArrayColumn<double> astigmMultError(*table_p_, "astigmMultError");             
  		
        ArrayColumn<double> illumOffset(*table_p_, "illumOffset");             
  		
        ArrayColumn<double> illumOffsetError(*table_p_, "illumOffsetError");             
  		
        ArrayColumn<double> fitRMS(*table_p_, "fitRMS");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            ambientTemperature.put(rowIndex, rows.at(i)->getAmbientTemperature().get());
	

	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            focusMethod.put(rowIndex, CFocusMethod::name(rows.at(i)->getFocusMethod()));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            pointingDirection.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getPointingDirection()));
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            wereFixed.put(rowIndex, basic2CASA1D<bool,bool>(rows.at(i)->getWereFixed()));
	

	
            offset.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getOffset()));
	

	
            offsetError.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getOffsetError()));
	

	
            offsetWasTied.put(rowIndex, basic2CASA2D<bool,bool>(rows.at(i)->getOffsetWasTied()));
	

	
            reducedChiSquared.put(rowIndex, basic2CASA2D<double,double>(rows.at(i)->getReducedChiSquared()));
	

	
            position.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getPosition()));
	

		
	
            if (rows.at(i)->isPolarizationsAveragedExists())
                polarizationsAveraged.put(rowIndex, rows.at(i)->getPolarizationsAveraged());
	

	
            if (rows.at(i)->isFocusCurveWidthExists())
                focusCurveWidth.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getFocusCurveWidth()));
	

	
            if (rows.at(i)->isFocusCurveWidthErrorExists())
                focusCurveWidthError.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getFocusCurveWidthError()));
	

	
            if (rows.at(i)->isFocusCurveWasFixedExists())
                focusCurveWasFixed.put(rowIndex, basic2CASA1D<bool,bool>(rows.at(i)->getFocusCurveWasFixed()));
	

	
            if (rows.at(i)->isOffIntensityExists())
                offIntensity.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getOffIntensity()));
	

	
            if (rows.at(i)->isOffIntensityErrorExists())
                offIntensityError.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getOffIntensityError()));
	

	
            if (rows.at(i)->isOffIntensityWasFixedExists())
                offIntensityWasFixed.put(rowIndex, rows.at(i)->getOffIntensityWasFixed());
	

	
            if (rows.at(i)->isPeakIntensityExists())
                peakIntensity.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getPeakIntensity()));
	

	
            if (rows.at(i)->isPeakIntensityErrorExists())
                peakIntensityError.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getPeakIntensityError()));
	

	
            if (rows.at(i)->isPeakIntensityWasFixedExists())
                peakIntensityWasFixed.put(rowIndex, rows.at(i)->getPeakIntensityWasFixed());
	

	
            if (rows.at(i)->isAstigmPlusExists())
                astigmPlus.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getAstigmPlus()));
	

	
            if (rows.at(i)->isAstigmPlusErrorExists())
                astigmPlusError.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getAstigmPlusError()));
	

	
            if (rows.at(i)->isAstigmMultExists())
                astigmMult.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getAstigmMult()));
	

	
            if (rows.at(i)->isAstigmMultErrorExists())
                astigmMultError.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getAstigmMultError()));
	

	
            if (rows.at(i)->isIllumOffsetExists())
                illumOffset.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getIllumOffset()));
	

	
            if (rows.at(i)->isIllumOffsetErrorExists())
                illumOffsetError.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getIllumOffsetError()));
	

	
            if (rows.at(i)->isFitRMSExists())
                fitRMS.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getFitRMS()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALFOCUSMODEL::ASDM_CALFOCUSMODEL() {
        name_ = "ASDM_CALFOCUSMODEL";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalFocusModel table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("polarizationType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaMake", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numCoeff", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numSourceObs", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("coeffName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("coeffFormula", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("coeffValue", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("coeffError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("coeffFixed", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("focusModel", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("focusRMS", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("reducedChiSquared", "blabla"));
  		
  		  		
    }

    ASDM_CALFOCUSMODEL::~ASDM_CALFOCUSMODEL() {
    }

    const TableDesc& ASDM_CALFOCUSMODEL::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALFOCUSMODEL::fill(const ASDM& asdm) {
        vector<CalFocusModelRow*> rows = asdm.getCalFocusModel().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> polarizationType(*table_p_, "polarizationType");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<String> antennaMake(*table_p_, "antennaMake");             
  		
        ScalarColumn<int> numCoeff(*table_p_, "numCoeff");             
  		
        ScalarColumn<int> numSourceObs(*table_p_, "numSourceObs");             
  		
        ArrayColumn<String> coeffName(*table_p_, "coeffName");             
  		
        ArrayColumn<String> coeffFormula(*table_p_, "coeffFormula");             
  		
        ArrayColumn<float> coeffValue(*table_p_, "coeffValue");             
  		
        ArrayColumn<float> coeffError(*table_p_, "coeffError");             
  		
        ArrayColumn<bool> coeffFixed(*table_p_, "coeffFixed");             
  		
        ScalarColumn<String> focusModel(*table_p_, "focusModel");             
  		
        ArrayColumn<double> focusRMS(*table_p_, "focusRMS");             
  		
        ScalarColumn<double> reducedChiSquared(*table_p_, "reducedChiSquared");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            polarizationType.put(rowIndex, CPolarizationType::name(rows.at(i)->getPolarizationType()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            antennaMake.put(rowIndex, CAntennaMake::name(rows.at(i)->getAntennaMake()));
	

	
            numCoeff.put(rowIndex, rows.at(i)->getNumCoeff());
	

	
            numSourceObs.put(rowIndex, rows.at(i)->getNumSourceObs());
	

	
            coeffName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCoeffName()));
	

	
            coeffFormula.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCoeffFormula()));
	

	
            coeffValue.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getCoeffValue()));
	

	
            coeffError.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getCoeffError()));
	

	
            coeffFixed.put(rowIndex, basic2CASA1D<bool,bool>(rows.at(i)->getCoeffFixed()));
	

	
            focusModel.put(rowIndex, rows.at(i)->getFocusModel());
	

	
            focusRMS.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getFocusRMS()));
	

	
            reducedChiSquared.put(rowIndex, rows.at(i)->getReducedChiSquared());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALGAIN::ASDM_CALGAIN() {
        name_ = "ASDM_CALGAIN";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalGain table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("gain", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("gainValid", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("fit", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("fitWeight", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("totalGainValid", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("totalFit", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("totalFitWeight", "blabla"));
  		
  		  		
    }

    ASDM_CALGAIN::~ASDM_CALGAIN() {
    }

    const TableDesc& ASDM_CALGAIN::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALGAIN::fill(const ASDM& asdm) {
        vector<CalGainRow*> rows = asdm.getCalGain().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<float> gain(*table_p_, "gain");             
  		
        ScalarColumn<bool> gainValid(*table_p_, "gainValid");             
  		
        ScalarColumn<float> fit(*table_p_, "fit");             
  		
        ScalarColumn<float> fitWeight(*table_p_, "fitWeight");             
  		
        ScalarColumn<bool> totalGainValid(*table_p_, "totalGainValid");             
  		
        ScalarColumn<float> totalFit(*table_p_, "totalFit");             
  		
        ScalarColumn<float> totalFitWeight(*table_p_, "totalFitWeight");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            gain.put(rowIndex, rows.at(i)->getGain());
	

	
            gainValid.put(rowIndex, rows.at(i)->getGainValid());
	

	
            fit.put(rowIndex, rows.at(i)->getFit());
	

	
            fitWeight.put(rowIndex, rows.at(i)->getFitWeight());
	

	
            totalGainValid.put(rowIndex, rows.at(i)->getTotalGainValid());
	

	
            totalFit.put(rowIndex, rows.at(i)->getTotalFit());
	

	
            totalFitWeight.put(rowIndex, rows.at(i)->getTotalFitWeight());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALHOLOGRAPHY::ASDM_CALHOLOGRAPHY() {
        name_ = "ASDM_CALHOLOGRAPHY";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalHolography table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaMake", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("ambientTemperature", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("focusPosition", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("illuminationTaper", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPanelModes", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("beamMapUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("rawRMS", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("weightedRMS", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("surfaceMapUID", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("direction", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numScrew", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("screwName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("screwMotion", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("screwMotionError", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("gravCorrection", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("gravOptRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("tempCorrection", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tempOptRange", "blabla"));
  		  		
    }

    ASDM_CALHOLOGRAPHY::~ASDM_CALHOLOGRAPHY() {
    }

    const TableDesc& ASDM_CALHOLOGRAPHY::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALHOLOGRAPHY::fill(const ASDM& asdm) {
        vector<CalHolographyRow*> rows = asdm.getCalHolography().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<String> antennaMake(*table_p_, "antennaMake");             
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<double> ambientTemperature(*table_p_, "ambientTemperature");             
  		
        ArrayColumn<double> focusPosition(*table_p_, "focusPosition");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<double> illuminationTaper(*table_p_, "illuminationTaper");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ScalarColumn<int> numPanelModes(*table_p_, "numPanelModes");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> beamMapUID(*table_p_, "beamMapUID");             
  		
        ScalarColumn<double> rawRMS(*table_p_, "rawRMS");             
  		
        ScalarColumn<double> weightedRMS(*table_p_, "weightedRMS");             
  		
        ScalarColumn<String> surfaceMapUID(*table_p_, "surfaceMapUID");             
  		
        ArrayColumn<double> direction(*table_p_, "direction");             
  		
  		
        ScalarColumn<int> numScrew(*table_p_, "numScrew");             
  		
        ArrayColumn<String> screwName(*table_p_, "screwName");             
  		
        ArrayColumn<double> screwMotion(*table_p_, "screwMotion");             
  		
        ArrayColumn<double> screwMotionError(*table_p_, "screwMotionError");             
  		
        ScalarColumn<bool> gravCorrection(*table_p_, "gravCorrection");             
  		
        ArrayColumn<double> gravOptRange(*table_p_, "gravOptRange");             
  		
        ScalarColumn<bool> tempCorrection(*table_p_, "tempCorrection");             
  		
        ArrayColumn<double> tempOptRange(*table_p_, "tempOptRange");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            antennaMake.put(rowIndex, CAntennaMake::name(rows.at(i)->getAntennaMake()));
	

	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            ambientTemperature.put(rowIndex, rows.at(i)->getAmbientTemperature().get());
	

	
            focusPosition.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getFocusPosition()));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            illuminationTaper.put(rowIndex, rows.at(i)->getIlluminationTaper());
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            numPanelModes.put(rowIndex, rows.at(i)->getNumPanelModes());
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            beamMapUID.put(rowIndex, rows.at(i)->getBeamMapUID().toString());
	

	
            rawRMS.put(rowIndex, rows.at(i)->getRawRMS().get());
	

	
            weightedRMS.put(rowIndex, rows.at(i)->getWeightedRMS().get());
	

	
            surfaceMapUID.put(rowIndex, rows.at(i)->getSurfaceMapUID().toString());
	

	
            direction.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getDirection()));
	

		
	
            if (rows.at(i)->isNumScrewExists())
                numScrew.put(rowIndex, rows.at(i)->getNumScrew());
	

	
            if (rows.at(i)->isScrewNameExists())
                screwName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getScrewName()));
	

	
            if (rows.at(i)->isScrewMotionExists())
                screwMotion.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getScrewMotion()));
	

	
            if (rows.at(i)->isScrewMotionErrorExists())
                screwMotionError.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getScrewMotionError()));
	

	
            if (rows.at(i)->isGravCorrectionExists())
                gravCorrection.put(rowIndex, rows.at(i)->getGravCorrection());
	

	
            if (rows.at(i)->isGravOptRangeExists())
                gravOptRange.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getGravOptRange()));
	

	
            if (rows.at(i)->isTempCorrectionExists())
                tempCorrection.put(rowIndex, rows.at(i)->getTempCorrection());
	

	
            if (rows.at(i)->isTempOptRangeExists())
                tempOptRange.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getTempOptRange()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALPHASE::ASDM_CALPHASE() {
        name_ = "ASDM_CALPHASE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalPhase table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBaseline", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("ampli", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("antennaNames", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("baselineLengths", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("decorrelationFactor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("direction", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("integrationTime", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phase", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phaseRMS", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("statPhaseRMS", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("correctionValidity", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("singleAntennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("refAntennaName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phaseAnt", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phaseAntRMS", "blabla"));
  		  		
    }

    ASDM_CALPHASE::~ASDM_CALPHASE() {
    }

    const TableDesc& ASDM_CALPHASE::tableDesc() const {
        return tableDesc_;
    }


		
			
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALPHASE::fill(const ASDM& asdm) {
        vector<CalPhaseRow*> rows = asdm.getCalPhase().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<int> numBaseline(*table_p_, "numBaseline");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<float> ampli(*table_p_, "ampli");             
  		
        ArrayColumn<String> antennaNames(*table_p_, "antennaNames");             
  		
        ArrayColumn<double> baselineLengths(*table_p_, "baselineLengths");             
  		
        ArrayColumn<float> decorrelationFactor(*table_p_, "decorrelationFactor");             
  		
        ArrayColumn<double> direction(*table_p_, "direction");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<double> integrationTime(*table_p_, "integrationTime");             
  		
        ArrayColumn<float> phase(*table_p_, "phase");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<float> phaseRMS(*table_p_, "phaseRMS");             
  		
        ArrayColumn<float> statPhaseRMS(*table_p_, "statPhaseRMS");             
  		
  		
        ArrayColumn<bool> correctionValidity(*table_p_, "correctionValidity");             
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ArrayColumn<String> singleAntennaName(*table_p_, "singleAntennaName");             
  		
        ScalarColumn<String> refAntennaName(*table_p_, "refAntennaName");             
  		
        ArrayColumn<float> phaseAnt(*table_p_, "phaseAnt");             
  		
        ArrayColumn<float> phaseAntRMS(*table_p_, "phaseAntRMS");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            basebandName.put(rowIndex, CBasebandName::name(rows.at(i)->getBasebandName()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            numBaseline.put(rowIndex, rows.at(i)->getNumBaseline());
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            ampli.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getAmpli()));
	

	
            antennaNames.put(rowIndex, basic2CASA2D<string,String>(rows.at(i)->getAntennaNames()));
	

	
            baselineLengths.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getBaselineLengths()));
	

	
            decorrelationFactor.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getDecorrelationFactor()));
	

	
            direction.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getDirection()));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            integrationTime.put(rowIndex, rows.at(i)->getIntegrationTime().get()/(1.0e9));
	

	
            phase.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getPhase()));
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            phaseRMS.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getPhaseRMS()));
	

	
            statPhaseRMS.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getStatPhaseRMS()));
	

		
	
            if (rows.at(i)->isCorrectionValidityExists())
                correctionValidity.put(rowIndex, basic2CASA1D<bool,bool>(rows.at(i)->getCorrectionValidity()));
	

	
            if (rows.at(i)->isNumAntennaExists())
                numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            if (rows.at(i)->isSingleAntennaNameExists())
                singleAntennaName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getSingleAntennaName()));
	

	
            if (rows.at(i)->isRefAntennaNameExists())
                refAntennaName.put(rowIndex, rows.at(i)->getRefAntennaName());
	

	
            if (rows.at(i)->isPhaseAntExists())
                phaseAnt.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getPhaseAnt()));
	

	
            if (rows.at(i)->isPhaseAntRMSExists())
                phaseAntRMS.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getPhaseAntRMS()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALPOINTING::ASDM_CALPOINTING() {
        name_ = "ASDM_CALPOINTING";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalPointing table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("ambientTemperature", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaMake", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("direction", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("pointingModelMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("pointingMethod", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("collOffsetRelative", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("collOffsetAbsolute", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("collError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("collOffsetTied", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("reducedChiSquared", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("averagedPolarizations", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("beamPA", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("beamPAError", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("beamPAWasFixed", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("beamWidth", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("beamWidthError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("beamWidthWasFixed", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("offIntensity", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("offIntensityError", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("offIntensityWasFixed", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("peakIntensity", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("peakIntensityError", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("peakIntensityWasFixed", "blabla"));
  		  		
    }

    ASDM_CALPOINTING::~ASDM_CALPOINTING() {
    }

    const TableDesc& ASDM_CALPOINTING::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALPOINTING::fill(const ASDM& asdm) {
        vector<CalPointingRow*> rows = asdm.getCalPointing().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<double> ambientTemperature(*table_p_, "ambientTemperature");             
  		
        ScalarColumn<String> antennaMake(*table_p_, "antennaMake");             
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ArrayColumn<double> direction(*table_p_, "direction");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<String> pointingModelMode(*table_p_, "pointingModelMode");             
  		
        ScalarColumn<String> pointingMethod(*table_p_, "pointingMethod");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<double> collOffsetRelative(*table_p_, "collOffsetRelative");             
  		
        ArrayColumn<double> collOffsetAbsolute(*table_p_, "collOffsetAbsolute");             
  		
        ArrayColumn<double> collError(*table_p_, "collError");             
  		
        ArrayColumn<bool> collOffsetTied(*table_p_, "collOffsetTied");             
  		
        ArrayColumn<double> reducedChiSquared(*table_p_, "reducedChiSquared");             
  		
  		
        ScalarColumn<bool> averagedPolarizations(*table_p_, "averagedPolarizations");             
  		
        ArrayColumn<double> beamPA(*table_p_, "beamPA");             
  		
        ArrayColumn<double> beamPAError(*table_p_, "beamPAError");             
  		
        ScalarColumn<bool> beamPAWasFixed(*table_p_, "beamPAWasFixed");             
  		
        ArrayColumn<double> beamWidth(*table_p_, "beamWidth");             
  		
        ArrayColumn<double> beamWidthError(*table_p_, "beamWidthError");             
  		
        ArrayColumn<bool> beamWidthWasFixed(*table_p_, "beamWidthWasFixed");             
  		
        ArrayColumn<double> offIntensity(*table_p_, "offIntensity");             
  		
        ArrayColumn<double> offIntensityError(*table_p_, "offIntensityError");             
  		
        ScalarColumn<bool> offIntensityWasFixed(*table_p_, "offIntensityWasFixed");             
  		
        ArrayColumn<double> peakIntensity(*table_p_, "peakIntensity");             
  		
        ArrayColumn<double> peakIntensityError(*table_p_, "peakIntensityError");             
  		
        ScalarColumn<bool> peakIntensityWasFixed(*table_p_, "peakIntensityWasFixed");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            ambientTemperature.put(rowIndex, rows.at(i)->getAmbientTemperature().get());
	

	
            antennaMake.put(rowIndex, CAntennaMake::name(rows.at(i)->getAntennaMake()));
	

	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            direction.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getDirection()));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            pointingModelMode.put(rowIndex, CPointingModelMode::name(rows.at(i)->getPointingModelMode()));
	

	
            pointingMethod.put(rowIndex, CPointingMethod::name(rows.at(i)->getPointingMethod()));
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            collOffsetRelative.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getCollOffsetRelative()));
	

	
            collOffsetAbsolute.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getCollOffsetAbsolute()));
	

	
            collError.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getCollError()));
	

	
            collOffsetTied.put(rowIndex, basic2CASA2D<bool,bool>(rows.at(i)->getCollOffsetTied()));
	

	
            reducedChiSquared.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getReducedChiSquared()));
	

		
	
            if (rows.at(i)->isAveragedPolarizationsExists())
                averagedPolarizations.put(rowIndex, rows.at(i)->getAveragedPolarizations());
	

	
            if (rows.at(i)->isBeamPAExists())
                beamPA.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getBeamPA()));
	

	
            if (rows.at(i)->isBeamPAErrorExists())
                beamPAError.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getBeamPAError()));
	

	
            if (rows.at(i)->isBeamPAWasFixedExists())
                beamPAWasFixed.put(rowIndex, rows.at(i)->getBeamPAWasFixed());
	

	
            if (rows.at(i)->isBeamWidthExists())
                beamWidth.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getBeamWidth()));
	

	
            if (rows.at(i)->isBeamWidthErrorExists())
                beamWidthError.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getBeamWidthError()));
	

	
            if (rows.at(i)->isBeamWidthWasFixedExists())
                beamWidthWasFixed.put(rowIndex, basic2CASA1D<bool,bool>(rows.at(i)->getBeamWidthWasFixed()));
	

	
            if (rows.at(i)->isOffIntensityExists())
                offIntensity.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getOffIntensity()));
	

	
            if (rows.at(i)->isOffIntensityErrorExists())
                offIntensityError.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getOffIntensityError()));
	

	
            if (rows.at(i)->isOffIntensityWasFixedExists())
                offIntensityWasFixed.put(rowIndex, rows.at(i)->getOffIntensityWasFixed());
	

	
            if (rows.at(i)->isPeakIntensityExists())
                peakIntensity.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getPeakIntensity()));
	

	
            if (rows.at(i)->isPeakIntensityErrorExists())
                peakIntensityError.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getPeakIntensityError()));
	

	
            if (rows.at(i)->isPeakIntensityWasFixedExists())
                peakIntensityWasFixed.put(rowIndex, rows.at(i)->getPeakIntensityWasFixed());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALPOINTINGMODEL::ASDM_CALPOINTINGMODEL() {
        name_ = "ASDM_CALPOINTINGMODEL";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalPointingModel table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaMake", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("pointingModelMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("polarizationType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numCoeff", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("coeffName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("coeffVal", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("coeffError", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("coeffFixed", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("azimuthRMS", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("elevationRms", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("skyRMS", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("reducedChiSquared", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numObs", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("coeffFormula", "blabla"));
  		  		
    }

    ASDM_CALPOINTINGMODEL::~ASDM_CALPOINTINGMODEL() {
    }

    const TableDesc& ASDM_CALPOINTINGMODEL::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
		
			
		
			
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALPOINTINGMODEL::fill(const ASDM& asdm) {
        vector<CalPointingModelRow*> rows = asdm.getCalPointingModel().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<String> antennaMake(*table_p_, "antennaMake");             
  		
        ScalarColumn<String> pointingModelMode(*table_p_, "pointingModelMode");             
  		
        ScalarColumn<String> polarizationType(*table_p_, "polarizationType");             
  		
        ScalarColumn<int> numCoeff(*table_p_, "numCoeff");             
  		
        ArrayColumn<String> coeffName(*table_p_, "coeffName");             
  		
        ArrayColumn<float> coeffVal(*table_p_, "coeffVal");             
  		
        ArrayColumn<float> coeffError(*table_p_, "coeffError");             
  		
        ArrayColumn<bool> coeffFixed(*table_p_, "coeffFixed");             
  		
        ScalarColumn<double> azimuthRMS(*table_p_, "azimuthRMS");             
  		
        ScalarColumn<double> elevationRms(*table_p_, "elevationRms");             
  		
        ScalarColumn<double> skyRMS(*table_p_, "skyRMS");             
  		
        ScalarColumn<double> reducedChiSquared(*table_p_, "reducedChiSquared");             
  		
  		
        ScalarColumn<int> numObs(*table_p_, "numObs");             
  		
        ArrayColumn<String> coeffFormula(*table_p_, "coeffFormula");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            antennaMake.put(rowIndex, CAntennaMake::name(rows.at(i)->getAntennaMake()));
	

	
            pointingModelMode.put(rowIndex, CPointingModelMode::name(rows.at(i)->getPointingModelMode()));
	

	
            polarizationType.put(rowIndex, CPolarizationType::name(rows.at(i)->getPolarizationType()));
	

	
            numCoeff.put(rowIndex, rows.at(i)->getNumCoeff());
	

	
            coeffName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCoeffName()));
	

	
            coeffVal.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getCoeffVal()));
	

	
            coeffError.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getCoeffError()));
	

	
            coeffFixed.put(rowIndex, basic2CASA1D<bool,bool>(rows.at(i)->getCoeffFixed()));
	

	
            azimuthRMS.put(rowIndex, rows.at(i)->getAzimuthRMS().get());
	

	
            elevationRms.put(rowIndex, rows.at(i)->getElevationRms().get());
	

	
            skyRMS.put(rowIndex, rows.at(i)->getSkyRMS().get());
	

	
            reducedChiSquared.put(rowIndex, rows.at(i)->getReducedChiSquared());
	

		
	
            if (rows.at(i)->isNumObsExists())
                numObs.put(rowIndex, rows.at(i)->getNumObs());
	

	
            if (rows.at(i)->isCoeffFormulaExists())
                coeffFormula.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCoeffFormula()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALPOSITION::ASDM_CALPOSITION() {
        name_ = "ASDM_CALPOSITION";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalPosition table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("antennaPosition", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("stationName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("stationPosition", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("positionMethod", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("refAntennaNames", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("axesOffset", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("axesOffsetErr", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("axesOffsetFixed", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("positionOffset", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("positionErr", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("reducedChiSquared", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("delayRms", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("phaseRms", "blabla"));
  		  		
    }

    ASDM_CALPOSITION::~ASDM_CALPOSITION() {
    }

    const TableDesc& ASDM_CALPOSITION::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALPOSITION::fill(const ASDM& asdm) {
        vector<CalPositionRow*> rows = asdm.getCalPosition().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ArrayColumn<double> antennaPosition(*table_p_, "antennaPosition");             
  		
        ScalarColumn<String> stationName(*table_p_, "stationName");             
  		
        ArrayColumn<double> stationPosition(*table_p_, "stationPosition");             
  		
        ScalarColumn<String> positionMethod(*table_p_, "positionMethod");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ArrayColumn<String> refAntennaNames(*table_p_, "refAntennaNames");             
  		
        ScalarColumn<double> axesOffset(*table_p_, "axesOffset");             
  		
        ScalarColumn<double> axesOffsetErr(*table_p_, "axesOffsetErr");             
  		
        ScalarColumn<bool> axesOffsetFixed(*table_p_, "axesOffsetFixed");             
  		
        ArrayColumn<double> positionOffset(*table_p_, "positionOffset");             
  		
        ArrayColumn<double> positionErr(*table_p_, "positionErr");             
  		
        ScalarColumn<double> reducedChiSquared(*table_p_, "reducedChiSquared");             
  		
  		
        ScalarColumn<double> delayRms(*table_p_, "delayRms");             
  		
        ScalarColumn<double> phaseRms(*table_p_, "phaseRms");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            antennaPosition.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getAntennaPosition()));
	

	
            stationName.put(rowIndex, rows.at(i)->getStationName());
	

	
            stationPosition.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getStationPosition()));
	

	
            positionMethod.put(rowIndex, CPositionMethod::name(rows.at(i)->getPositionMethod()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            refAntennaNames.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getRefAntennaNames()));
	

	
            axesOffset.put(rowIndex, rows.at(i)->getAxesOffset().get());
	

	
            axesOffsetErr.put(rowIndex, rows.at(i)->getAxesOffsetErr().get());
	

	
            axesOffsetFixed.put(rowIndex, rows.at(i)->getAxesOffsetFixed());
	

	
            positionOffset.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getPositionOffset()));
	

	
            positionErr.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getPositionErr()));
	

	
            reducedChiSquared.put(rowIndex, rows.at(i)->getReducedChiSquared());
	

		
	
            if (rows.at(i)->isDelayRmsExists())
                delayRms.put(rowIndex, rows.at(i)->getDelayRms());
	

	
            if (rows.at(i)->isPhaseRmsExists())
                phaseRms.put(rowIndex, rows.at(i)->getPhaseRms().get());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALPRIMARYBEAM::ASDM_CALPRIMARYBEAM() {
        name_ = "ASDM_CALPRIMARYBEAM";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalPrimaryBeam table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaMake", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numSubband", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("mainBeamEfficiency", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("beamDescriptionUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("relativeAmplitudeRms", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("direction", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("minValidDirection", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("maxValidDirection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("descriptionType", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("imageChannelNumber", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("imageNominalFrequency", "blabla"));
  		
  		  		
    }

    ASDM_CALPRIMARYBEAM::~ASDM_CALPRIMARYBEAM() {
    }

    const TableDesc& ASDM_CALPRIMARYBEAM::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
	
    void ASDM_CALPRIMARYBEAM::fill(const ASDM& asdm) {
        vector<CalPrimaryBeamRow*> rows = asdm.getCalPrimaryBeam().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<String> antennaMake(*table_p_, "antennaMake");             
  		
        ScalarColumn<int> numSubband(*table_p_, "numSubband");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<double> mainBeamEfficiency(*table_p_, "mainBeamEfficiency");             
  		
        ScalarColumn<String> beamDescriptionUID(*table_p_, "beamDescriptionUID");             
  		
        ScalarColumn<float> relativeAmplitudeRms(*table_p_, "relativeAmplitudeRms");             
  		
        ArrayColumn<double> direction(*table_p_, "direction");             
  		
        ArrayColumn<double> minValidDirection(*table_p_, "minValidDirection");             
  		
        ArrayColumn<double> maxValidDirection(*table_p_, "maxValidDirection");             
  		
        ScalarColumn<String> descriptionType(*table_p_, "descriptionType");             
  		
        ArrayColumn<int> imageChannelNumber(*table_p_, "imageChannelNumber");             
  		
        ArrayColumn<double> imageNominalFrequency(*table_p_, "imageNominalFrequency");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            antennaMake.put(rowIndex, CAntennaMake::name(rows.at(i)->getAntennaMake()));
	

	
            numSubband.put(rowIndex, rows.at(i)->getNumSubband());
	

	
            frequencyRange.put(rowIndex, ext2CASA2D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            mainBeamEfficiency.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getMainBeamEfficiency()));
	

	
            beamDescriptionUID.put(rowIndex, rows.at(i)->getBeamDescriptionUID().toString());
	

	
            relativeAmplitudeRms.put(rowIndex, rows.at(i)->getRelativeAmplitudeRms());
	

	
            direction.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getDirection()));
	

	
            minValidDirection.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getMinValidDirection()));
	

	
            maxValidDirection.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getMaxValidDirection()));
	

	
            descriptionType.put(rowIndex, CPrimaryBeamDescription::name(rows.at(i)->getDescriptionType()));
	

	
            imageChannelNumber.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getImageChannelNumber()));
	

	
            imageNominalFrequency.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getImageNominalFrequency()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALREDUCTION::ASDM_CALREDUCTION() {
        name_ = "ASDM_CALREDUCTION";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalReduction table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numApplied", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("appliedCalibrations", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numParam", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("paramSet", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numInvalidConditions", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("invalidConditions", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("timeReduced", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("messages", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("software", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("softwareVersion", "blabla"));
  		
  		  		
    }

    ASDM_CALREDUCTION::~ASDM_CALREDUCTION() {
    }

    const TableDesc& ASDM_CALREDUCTION::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALREDUCTION::fill(const ASDM& asdm) {
        vector<CalReductionRow*> rows = asdm.getCalReduction().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<int> numApplied(*table_p_, "numApplied");             
  		
        ArrayColumn<String> appliedCalibrations(*table_p_, "appliedCalibrations");             
  		
        ScalarColumn<int> numParam(*table_p_, "numParam");             
  		
        ArrayColumn<String> paramSet(*table_p_, "paramSet");             
  		
        ScalarColumn<int> numInvalidConditions(*table_p_, "numInvalidConditions");             
  		
        ArrayColumn<String> invalidConditions(*table_p_, "invalidConditions");             
  		
        ScalarColumn<double> timeReduced(*table_p_, "timeReduced");             
  		
        ScalarColumn<String> messages(*table_p_, "messages");             
  		
        ScalarColumn<String> software(*table_p_, "software");             
  		
        ScalarColumn<String> softwareVersion(*table_p_, "softwareVersion");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            numApplied.put(rowIndex, rows.at(i)->getNumApplied());
	

	
            appliedCalibrations.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getAppliedCalibrations()));
	

	
            numParam.put(rowIndex, rows.at(i)->getNumParam());
	

	
            paramSet.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getParamSet()));
	

	
            numInvalidConditions.put(rowIndex, rows.at(i)->getNumInvalidConditions());
	

	
            invalidConditions.put(rowIndex, enum2CASA1D<InvalidatingCondition,CInvalidatingCondition>(rows.at(i)->getInvalidConditions()));
	

	
            timeReduced.put(rowIndex, rows.at(i)->getTimeReduced().get()/(1.0e9));
	

	
            messages.put(rowIndex, rows.at(i)->getMessages());
	

	
            software.put(rowIndex, rows.at(i)->getSoftware());
	

	
            softwareVersion.put(rowIndex, rows.at(i)->getSoftwareVersion());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALSEEING::ASDM_CALSEEING() {
        name_ = "ASDM_CALSEEING";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalSeeing table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyRange", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("integrationTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBaseLengths", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("baselineLengths", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("phaseRMS", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("seeing", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("seeingError", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("exponent", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("outerScale", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("outerScaleRMS", "blabla"));
  		  		
    }

    ASDM_CALSEEING::~ASDM_CALSEEING() {
    }

    const TableDesc& ASDM_CALSEEING::tableDesc() const {
        return tableDesc_;
    }


		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALSEEING::fill(const ASDM& asdm) {
        vector<CalSeeingRow*> rows = asdm.getCalSeeing().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ArrayColumn<double> frequencyRange(*table_p_, "frequencyRange");             
  		
        ScalarColumn<double> integrationTime(*table_p_, "integrationTime");             
  		
        ScalarColumn<int> numBaseLengths(*table_p_, "numBaseLengths");             
  		
        ArrayColumn<double> baselineLengths(*table_p_, "baselineLengths");             
  		
        ArrayColumn<double> phaseRMS(*table_p_, "phaseRMS");             
  		
        ScalarColumn<double> seeing(*table_p_, "seeing");             
  		
        ScalarColumn<double> seeingError(*table_p_, "seeingError");             
  		
  		
        ScalarColumn<float> exponent(*table_p_, "exponent");             
  		
        ScalarColumn<double> outerScale(*table_p_, "outerScale");             
  		
        ScalarColumn<double> outerScaleRMS(*table_p_, "outerScaleRMS");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            atmPhaseCorrection.put(rowIndex, CAtmPhaseCorrection::name(rows.at(i)->getAtmPhaseCorrection()));
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            frequencyRange.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyRange()));
	

	
            integrationTime.put(rowIndex, rows.at(i)->getIntegrationTime().get()/(1.0e9));
	

	
            numBaseLengths.put(rowIndex, rows.at(i)->getNumBaseLengths());
	

	
            baselineLengths.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getBaselineLengths()));
	

	
            phaseRMS.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getPhaseRMS()));
	

	
            seeing.put(rowIndex, rows.at(i)->getSeeing().get());
	

	
            seeingError.put(rowIndex, rows.at(i)->getSeeingError().get());
	

		
	
            if (rows.at(i)->isExponentExists())
                exponent.put(rowIndex, rows.at(i)->getExponent());
	

	
            if (rows.at(i)->isOuterScaleExists())
                outerScale.put(rowIndex, rows.at(i)->getOuterScale().get());
	

	
            if (rows.at(i)->isOuterScaleRMSExists())
                outerScaleRMS.put(rowIndex, rows.at(i)->getOuterScaleRMS().get());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CALWVR::ASDM_CALWVR() {
        name_ = "ASDM_CALWVR";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CalWVR table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDataId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calReductionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endValidTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("wvrMethod", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numInputAntennas", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("inputAntennaNames", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numChan", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("chanFreq", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("chanWidth", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("refTemp", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPoly", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("pathCoeff", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("polyFreqLimits", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("wetPath", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("dryPath", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("water", "blabla"));
  		
  		  		
    }

    ASDM_CALWVR::~ASDM_CALWVR() {
    }

    const TableDesc& ASDM_CALWVR::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_CALWVR::fill(const ASDM& asdm) {
        vector<CalWVRRow*> rows = asdm.getCalWVR().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaName(*table_p_, "antennaName");             
  		
        ScalarColumn<String> calDataId(*table_p_, "calDataId");             
  		
        ScalarColumn<String> calReductionId(*table_p_, "calReductionId");             
  		
  		
        ScalarColumn<double> startValidTime(*table_p_, "startValidTime");             
  		
        ScalarColumn<double> endValidTime(*table_p_, "endValidTime");             
  		
        ScalarColumn<String> wvrMethod(*table_p_, "wvrMethod");             
  		
        ScalarColumn<int> numInputAntennas(*table_p_, "numInputAntennas");             
  		
        ArrayColumn<String> inputAntennaNames(*table_p_, "inputAntennaNames");             
  		
        ScalarColumn<int> numChan(*table_p_, "numChan");             
  		
        ArrayColumn<double> chanFreq(*table_p_, "chanFreq");             
  		
        ArrayColumn<double> chanWidth(*table_p_, "chanWidth");             
  		
        ArrayColumn<double> refTemp(*table_p_, "refTemp");             
  		
        ScalarColumn<int> numPoly(*table_p_, "numPoly");             
  		
        ArrayColumn<float> pathCoeff(*table_p_, "pathCoeff");             
  		
        ArrayColumn<double> polyFreqLimits(*table_p_, "polyFreqLimits");             
  		
        ArrayColumn<float> wetPath(*table_p_, "wetPath");             
  		
        ArrayColumn<float> dryPath(*table_p_, "dryPath");             
  		
        ScalarColumn<double> water(*table_p_, "water");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaName.put(rowIndex, rows.at(i)->getAntennaName());
	

	
            calDataId.put(rowIndex, rows.at(i)->getCalDataId().toString());
	

	
            calReductionId.put(rowIndex, rows.at(i)->getCalReductionId().toString());
	

		
	
            startValidTime.put(rowIndex, rows.at(i)->getStartValidTime().get()/(1.0e9));
	

	
            endValidTime.put(rowIndex, rows.at(i)->getEndValidTime().get()/(1.0e9));
	

	
            wvrMethod.put(rowIndex, CWVRMethod::name(rows.at(i)->getWvrMethod()));
	

	
            numInputAntennas.put(rowIndex, rows.at(i)->getNumInputAntennas());
	

	
            inputAntennaNames.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getInputAntennaNames()));
	

	
            numChan.put(rowIndex, rows.at(i)->getNumChan());
	

	
            chanFreq.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getChanFreq()));
	

	
            chanWidth.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getChanWidth()));
	

	
            refTemp.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getRefTemp()));
	

	
            numPoly.put(rowIndex, rows.at(i)->getNumPoly());
	

	
            pathCoeff.put(rowIndex, basic2CASA3D<float,float>(rows.at(i)->getPathCoeff()));
	

	
            polyFreqLimits.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getPolyFreqLimits()));
	

	
            wetPath.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getWetPath()));
	

	
            dryPath.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getDryPath()));
	

	
            water.put(rowIndex, rows.at(i)->getWater().get());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CONFIGDESCRIPTION::ASDM_CONFIGDESCRIPTION() {
        name_ = "ASDM_CONFIGDESCRIPTION";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset ConfigDescription table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("configDescriptionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numDataDescription", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numFeed", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("correlationMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAtmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("atmPhaseCorrection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("processorType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralType", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("feedId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("switchCycleId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("dataDescriptionId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("processorId", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("phasedArrayList", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAssocValues", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("assocNature", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("assocConfigDescriptionId", "blabla"));
  		  		
    }

    ASDM_CONFIGDESCRIPTION::~ASDM_CONFIGDESCRIPTION() {
    }

    const TableDesc& ASDM_CONFIGDESCRIPTION::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
			
		
			
		
			
			
		
	
    void ASDM_CONFIGDESCRIPTION::fill(const ASDM& asdm) {
        vector<ConfigDescriptionRow*> rows = asdm.getConfigDescription().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> configDescriptionId(*table_p_, "configDescriptionId");             
  		
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ScalarColumn<int> numDataDescription(*table_p_, "numDataDescription");             
  		
        ScalarColumn<int> numFeed(*table_p_, "numFeed");             
  		
        ScalarColumn<String> correlationMode(*table_p_, "correlationMode");             
  		
        ScalarColumn<int> numAtmPhaseCorrection(*table_p_, "numAtmPhaseCorrection");             
  		
        ArrayColumn<String> atmPhaseCorrection(*table_p_, "atmPhaseCorrection");             
  		
        ScalarColumn<String> processorType(*table_p_, "processorType");             
  		
        ScalarColumn<String> spectralType(*table_p_, "spectralType");             
  		
        ArrayColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ArrayColumn<int> feedId(*table_p_, "feedId");             
  		
        ArrayColumn<String> switchCycleId(*table_p_, "switchCycleId");             
  		
        ArrayColumn<String> dataDescriptionId(*table_p_, "dataDescriptionId");             
  		
        ScalarColumn<String> processorId(*table_p_, "processorId");             
  		
  		
        ArrayColumn<int> phasedArrayList(*table_p_, "phasedArrayList");             
  		
        ScalarColumn<int> numAssocValues(*table_p_, "numAssocValues");             
  		
        ArrayColumn<String> assocNature(*table_p_, "assocNature");             
  		
        ArrayColumn<String> assocConfigDescriptionId(*table_p_, "assocConfigDescriptionId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            configDescriptionId.put(rowIndex, rows.at(i)->getConfigDescriptionId().toString());
	

		
	
            numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            numDataDescription.put(rowIndex, rows.at(i)->getNumDataDescription());
	

	
            numFeed.put(rowIndex, rows.at(i)->getNumFeed());
	

	
            correlationMode.put(rowIndex, CCorrelationMode::name(rows.at(i)->getCorrelationMode()));
	

	
            numAtmPhaseCorrection.put(rowIndex, rows.at(i)->getNumAtmPhaseCorrection());
	

	
            atmPhaseCorrection.put(rowIndex, enum2CASA1D<AtmPhaseCorrection,CAtmPhaseCorrection>(rows.at(i)->getAtmPhaseCorrection()));
	

	
            processorType.put(rowIndex, CProcessorType::name(rows.at(i)->getProcessorType()));
	

	
            spectralType.put(rowIndex, CSpectralResolutionType::name(rows.at(i)->getSpectralType()));
	

	
            antennaId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getAntennaId()));
	

	
            feedId.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getFeedId()));
	

	
            switchCycleId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getSwitchCycleId()));
	

	
            dataDescriptionId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getDataDescriptionId()));
	

	
            processorId.put(rowIndex, rows.at(i)->getProcessorId().toString());
	

		
	
            if (rows.at(i)->isPhasedArrayListExists())
                phasedArrayList.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getPhasedArrayList()));
	

	
            if (rows.at(i)->isNumAssocValuesExists())
                numAssocValues.put(rowIndex, rows.at(i)->getNumAssocValues());
	

	
            if (rows.at(i)->isAssocNatureExists())
                assocNature.put(rowIndex, enum2CASA1D<SpectralResolutionType,CSpectralResolutionType>(rows.at(i)->getAssocNature()));
	

	
            if (rows.at(i)->isAssocConfigDescriptionIdExists())
                assocConfigDescriptionId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getAssocConfigDescriptionId()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_CORRELATORMODE::ASDM_CORRELATORMODE() {
        name_ = "ASDM_CORRELATORMODE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset CorrelatorMode table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("correlatorModeId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBaseband", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("basebandNames", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("basebandConfig", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("accumMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("binMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAxes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("axesOrderArray", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("filterMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("correlatorName", "blabla"));
  		
  		  		
    }

    ASDM_CORRELATORMODE::~ASDM_CORRELATORMODE() {
    }

    const TableDesc& ASDM_CORRELATORMODE::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
			
		
			
			
		
			
		
			
		
			
			
		
			
			
		
			
			
		
	
    void ASDM_CORRELATORMODE::fill(const ASDM& asdm) {
        vector<CorrelatorModeRow*> rows = asdm.getCorrelatorMode().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> correlatorModeId(*table_p_, "correlatorModeId");             
  		
  		
        ScalarColumn<int> numBaseband(*table_p_, "numBaseband");             
  		
        ArrayColumn<String> basebandNames(*table_p_, "basebandNames");             
  		
        ArrayColumn<int> basebandConfig(*table_p_, "basebandConfig");             
  		
        ScalarColumn<String> accumMode(*table_p_, "accumMode");             
  		
        ScalarColumn<int> binMode(*table_p_, "binMode");             
  		
        ScalarColumn<int> numAxes(*table_p_, "numAxes");             
  		
        ArrayColumn<String> axesOrderArray(*table_p_, "axesOrderArray");             
  		
        ArrayColumn<String> filterMode(*table_p_, "filterMode");             
  		
        ScalarColumn<String> correlatorName(*table_p_, "correlatorName");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            correlatorModeId.put(rowIndex, rows.at(i)->getCorrelatorModeId().toString());
	

		
	
            numBaseband.put(rowIndex, rows.at(i)->getNumBaseband());
	

	
            basebandNames.put(rowIndex, enum2CASA1D<BasebandName,CBasebandName>(rows.at(i)->getBasebandNames()));
	

	
            basebandConfig.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getBasebandConfig()));
	

	
            accumMode.put(rowIndex, CAccumMode::name(rows.at(i)->getAccumMode()));
	

	
            binMode.put(rowIndex, rows.at(i)->getBinMode());
	

	
            numAxes.put(rowIndex, rows.at(i)->getNumAxes());
	

	
            axesOrderArray.put(rowIndex, enum2CASA1D<AxisName,CAxisName>(rows.at(i)->getAxesOrderArray()));
	

	
            filterMode.put(rowIndex, enum2CASA1D<FilterMode,CFilterMode>(rows.at(i)->getFilterMode()));
	

	
            correlatorName.put(rowIndex, CCorrelatorName::name(rows.at(i)->getCorrelatorName()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_DATADESCRIPTION::ASDM_DATADESCRIPTION() {
        name_ = "ASDM_DATADESCRIPTION";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset DataDescription table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("dataDescriptionId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("polOrHoloId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("pulsarId", "blabla"));
  		  		
    }

    ASDM_DATADESCRIPTION::~ASDM_DATADESCRIPTION() {
    }

    const TableDesc& ASDM_DATADESCRIPTION::tableDesc() const {
        return tableDesc_;
    }


		
			
		
	
    void ASDM_DATADESCRIPTION::fill(const ASDM& asdm) {
        vector<DataDescriptionRow*> rows = asdm.getDataDescription().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> dataDescriptionId(*table_p_, "dataDescriptionId");             
  		
  		
        ScalarColumn<String> polOrHoloId(*table_p_, "polOrHoloId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
  		
        ScalarColumn<String> pulsarId(*table_p_, "pulsarId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            dataDescriptionId.put(rowIndex, rows.at(i)->getDataDescriptionId().toString());
	

		
	
            polOrHoloId.put(rowIndex, rows.at(i)->getPolOrHoloId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

		
	
            if (rows.at(i)->isPulsarIdExists())
                pulsarId.put(rowIndex, rows.at(i)->getPulsarId().toString());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_DELAYMODEL::ASDM_DELAYMODEL() {
        name_ = "ASDM_DELAYMODEL";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset DelayModel table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPoly", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("phaseDelay", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("phaseDelayRate", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("groupDelay", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("groupDelayRate", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("fieldId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("timeOrigin", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("atmosphericGroupDelay", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("atmosphericGroupDelayRate", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("geometricDelay", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("geometricDelayRate", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numLO", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("LOOffset", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("LOOffsetRate", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("dispersiveDelay", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("dispersiveDelayRate", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("atmosphericDryDelay", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("atmosphericWetDelay", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("padDelay", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("antennaDelay", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationType", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("electronicDelay", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("electronicDelayRate", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("receiverDelay", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("IFDelay", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("LODelay", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("crossPolarizationDelay", "blabla"));
  		  		
    }

    ASDM_DELAYMODEL::~ASDM_DELAYMODEL() {
    }

    const TableDesc& ASDM_DELAYMODEL::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_DELAYMODEL::fill(const ASDM& asdm) {
        vector<DelayModelRow*> rows = asdm.getDelayModel().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
        ScalarColumn<int> numPoly(*table_p_, "numPoly");             
  		
        ArrayColumn<double> phaseDelay(*table_p_, "phaseDelay");             
  		
        ArrayColumn<double> phaseDelayRate(*table_p_, "phaseDelayRate");             
  		
        ArrayColumn<double> groupDelay(*table_p_, "groupDelay");             
  		
        ArrayColumn<double> groupDelayRate(*table_p_, "groupDelayRate");             
  		
        ScalarColumn<String> fieldId(*table_p_, "fieldId");             
  		
  		
        ScalarColumn<double> timeOrigin(*table_p_, "timeOrigin");             
  		
        ScalarColumn<double> atmosphericGroupDelay(*table_p_, "atmosphericGroupDelay");             
  		
        ScalarColumn<double> atmosphericGroupDelayRate(*table_p_, "atmosphericGroupDelayRate");             
  		
        ScalarColumn<double> geometricDelay(*table_p_, "geometricDelay");             
  		
        ScalarColumn<double> geometricDelayRate(*table_p_, "geometricDelayRate");             
  		
        ScalarColumn<int> numLO(*table_p_, "numLO");             
  		
        ArrayColumn<double> LOOffset(*table_p_, "LOOffset");             
  		
        ArrayColumn<double> LOOffsetRate(*table_p_, "LOOffsetRate");             
  		
        ScalarColumn<double> dispersiveDelay(*table_p_, "dispersiveDelay");             
  		
        ScalarColumn<double> dispersiveDelayRate(*table_p_, "dispersiveDelayRate");             
  		
        ScalarColumn<double> atmosphericDryDelay(*table_p_, "atmosphericDryDelay");             
  		
        ScalarColumn<double> atmosphericWetDelay(*table_p_, "atmosphericWetDelay");             
  		
        ScalarColumn<double> padDelay(*table_p_, "padDelay");             
  		
        ScalarColumn<double> antennaDelay(*table_p_, "antennaDelay");             
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<String> polarizationType(*table_p_, "polarizationType");             
  		
        ArrayColumn<double> electronicDelay(*table_p_, "electronicDelay");             
  		
        ArrayColumn<double> electronicDelayRate(*table_p_, "electronicDelayRate");             
  		
        ArrayColumn<double> receiverDelay(*table_p_, "receiverDelay");             
  		
        ArrayColumn<double> IFDelay(*table_p_, "IFDelay");             
  		
        ArrayColumn<double> LODelay(*table_p_, "LODelay");             
  		
        ScalarColumn<double> crossPolarizationDelay(*table_p_, "crossPolarizationDelay");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
	
            numPoly.put(rowIndex, rows.at(i)->getNumPoly());
	

	
            phaseDelay.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getPhaseDelay()));
	

	
            phaseDelayRate.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getPhaseDelayRate()));
	

	
            groupDelay.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getGroupDelay()));
	

	
            groupDelayRate.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getGroupDelayRate()));
	

	
            fieldId.put(rowIndex, rows.at(i)->getFieldId().toString());
	

		
	
            if (rows.at(i)->isTimeOriginExists())
                timeOrigin.put(rowIndex, rows.at(i)->getTimeOrigin().get()/(1.0e9));
	

	
            if (rows.at(i)->isAtmosphericGroupDelayExists())
                atmosphericGroupDelay.put(rowIndex, rows.at(i)->getAtmosphericGroupDelay());
	

	
            if (rows.at(i)->isAtmosphericGroupDelayRateExists())
                atmosphericGroupDelayRate.put(rowIndex, rows.at(i)->getAtmosphericGroupDelayRate());
	

	
            if (rows.at(i)->isGeometricDelayExists())
                geometricDelay.put(rowIndex, rows.at(i)->getGeometricDelay());
	

	
            if (rows.at(i)->isGeometricDelayRateExists())
                geometricDelayRate.put(rowIndex, rows.at(i)->getGeometricDelayRate());
	

	
            if (rows.at(i)->isNumLOExists())
                numLO.put(rowIndex, rows.at(i)->getNumLO());
	

	
            if (rows.at(i)->isLOOffsetExists())
                LOOffset.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getLOOffset()));
	

	
            if (rows.at(i)->isLOOffsetRateExists())
                LOOffsetRate.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getLOOffsetRate()));
	

	
            if (rows.at(i)->isDispersiveDelayExists())
                dispersiveDelay.put(rowIndex, rows.at(i)->getDispersiveDelay());
	

	
            if (rows.at(i)->isDispersiveDelayRateExists())
                dispersiveDelayRate.put(rowIndex, rows.at(i)->getDispersiveDelayRate());
	

	
            if (rows.at(i)->isAtmosphericDryDelayExists())
                atmosphericDryDelay.put(rowIndex, rows.at(i)->getAtmosphericDryDelay());
	

	
            if (rows.at(i)->isAtmosphericWetDelayExists())
                atmosphericWetDelay.put(rowIndex, rows.at(i)->getAtmosphericWetDelay());
	

	
            if (rows.at(i)->isPadDelayExists())
                padDelay.put(rowIndex, rows.at(i)->getPadDelay());
	

	
            if (rows.at(i)->isAntennaDelayExists())
                antennaDelay.put(rowIndex, rows.at(i)->getAntennaDelay());
	

	
            if (rows.at(i)->isNumReceptorExists())
                numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            if (rows.at(i)->isPolarizationTypeExists())
                polarizationType.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationType()));
	

	
            if (rows.at(i)->isElectronicDelayExists())
                electronicDelay.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getElectronicDelay()));
	

	
            if (rows.at(i)->isElectronicDelayRateExists())
                electronicDelayRate.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getElectronicDelayRate()));
	

	
            if (rows.at(i)->isReceiverDelayExists())
                receiverDelay.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getReceiverDelay()));
	

	
            if (rows.at(i)->isIFDelayExists())
                IFDelay.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getIFDelay()));
	

	
            if (rows.at(i)->isLODelayExists())
                LODelay.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getLODelay()));
	

	
            if (rows.at(i)->isCrossPolarizationDelayExists())
                crossPolarizationDelay.put(rowIndex, rows.at(i)->getCrossPolarizationDelay());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_DELAYMODELFIXEDPARAMETERS::ASDM_DELAYMODELFIXEDPARAMETERS() {
        name_ = "ASDM_DELAYMODELFIXEDPARAMETERS";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset DelayModelFixedParameters table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("delayModelFixedParametersId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("delayModelVersion", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("gaussConstant", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("newtonianConstant", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("gravity", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("earthFlattening", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("earthRadius", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("moonEarthMassRatio", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("ephemerisEpoch", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("earthTideLag", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("earthGM", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("moonGM", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("sunGM", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("loveNumberH", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("loveNumberL", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("precessionConstant", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("lightTime1AU", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("speedOfLight", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("delayModelFlags", "blabla"));
  		  		
    }

    ASDM_DELAYMODELFIXEDPARAMETERS::~ASDM_DELAYMODELFIXEDPARAMETERS() {
    }

    const TableDesc& ASDM_DELAYMODELFIXEDPARAMETERS::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_DELAYMODELFIXEDPARAMETERS::fill(const ASDM& asdm) {
        vector<DelayModelFixedParametersRow*> rows = asdm.getDelayModelFixedParameters().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> delayModelFixedParametersId(*table_p_, "delayModelFixedParametersId");             
  		
  		
        ScalarColumn<String> delayModelVersion(*table_p_, "delayModelVersion");             
  		
        ScalarColumn<String> execBlockId(*table_p_, "execBlockId");             
  		
  		
        ScalarColumn<double> gaussConstant(*table_p_, "gaussConstant");             
  		
        ScalarColumn<double> newtonianConstant(*table_p_, "newtonianConstant");             
  		
        ScalarColumn<double> gravity(*table_p_, "gravity");             
  		
        ScalarColumn<double> earthFlattening(*table_p_, "earthFlattening");             
  		
        ScalarColumn<double> earthRadius(*table_p_, "earthRadius");             
  		
        ScalarColumn<double> moonEarthMassRatio(*table_p_, "moonEarthMassRatio");             
  		
        ScalarColumn<String> ephemerisEpoch(*table_p_, "ephemerisEpoch");             
  		
        ScalarColumn<double> earthTideLag(*table_p_, "earthTideLag");             
  		
        ScalarColumn<double> earthGM(*table_p_, "earthGM");             
  		
        ScalarColumn<double> moonGM(*table_p_, "moonGM");             
  		
        ScalarColumn<double> sunGM(*table_p_, "sunGM");             
  		
        ScalarColumn<double> loveNumberH(*table_p_, "loveNumberH");             
  		
        ScalarColumn<double> loveNumberL(*table_p_, "loveNumberL");             
  		
        ScalarColumn<double> precessionConstant(*table_p_, "precessionConstant");             
  		
        ScalarColumn<double> lightTime1AU(*table_p_, "lightTime1AU");             
  		
        ScalarColumn<double> speedOfLight(*table_p_, "speedOfLight");             
  		
        ScalarColumn<String> delayModelFlags(*table_p_, "delayModelFlags");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            delayModelFixedParametersId.put(rowIndex, rows.at(i)->getDelayModelFixedParametersId().toString());
	

		
	
            delayModelVersion.put(rowIndex, rows.at(i)->getDelayModelVersion());
	

	
            execBlockId.put(rowIndex, rows.at(i)->getExecBlockId().toString());
	

		
	
            if (rows.at(i)->isGaussConstantExists())
                gaussConstant.put(rowIndex, rows.at(i)->getGaussConstant().get());
	

	
            if (rows.at(i)->isNewtonianConstantExists())
                newtonianConstant.put(rowIndex, rows.at(i)->getNewtonianConstant());
	

	
            if (rows.at(i)->isGravityExists())
                gravity.put(rowIndex, rows.at(i)->getGravity());
	

	
            if (rows.at(i)->isEarthFlatteningExists())
                earthFlattening.put(rowIndex, rows.at(i)->getEarthFlattening());
	

	
            if (rows.at(i)->isEarthRadiusExists())
                earthRadius.put(rowIndex, rows.at(i)->getEarthRadius().get());
	

	
            if (rows.at(i)->isMoonEarthMassRatioExists())
                moonEarthMassRatio.put(rowIndex, rows.at(i)->getMoonEarthMassRatio());
	

	
            if (rows.at(i)->isEphemerisEpochExists())
                ephemerisEpoch.put(rowIndex, rows.at(i)->getEphemerisEpoch());
	

	
            if (rows.at(i)->isEarthTideLagExists())
                earthTideLag.put(rowIndex, rows.at(i)->getEarthTideLag());
	

	
            if (rows.at(i)->isEarthGMExists())
                earthGM.put(rowIndex, rows.at(i)->getEarthGM());
	

	
            if (rows.at(i)->isMoonGMExists())
                moonGM.put(rowIndex, rows.at(i)->getMoonGM());
	

	
            if (rows.at(i)->isSunGMExists())
                sunGM.put(rowIndex, rows.at(i)->getSunGM());
	

	
            if (rows.at(i)->isLoveNumberHExists())
                loveNumberH.put(rowIndex, rows.at(i)->getLoveNumberH());
	

	
            if (rows.at(i)->isLoveNumberLExists())
                loveNumberL.put(rowIndex, rows.at(i)->getLoveNumberL());
	

	
            if (rows.at(i)->isPrecessionConstantExists())
                precessionConstant.put(rowIndex, rows.at(i)->getPrecessionConstant().get());
	

	
            if (rows.at(i)->isLightTime1AUExists())
                lightTime1AU.put(rowIndex, rows.at(i)->getLightTime1AU());
	

	
            if (rows.at(i)->isSpeedOfLightExists())
                speedOfLight.put(rowIndex, rows.at(i)->getSpeedOfLight().get());
	

	
            if (rows.at(i)->isDelayModelFlagsExists())
                delayModelFlags.put(rowIndex, rows.at(i)->getDelayModelFlags());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_DELAYMODELVARIABLEPARAMETERS::ASDM_DELAYMODELVARIABLEPARAMETERS() {
        name_ = "ASDM_DELAYMODELVARIABLEPARAMETERS";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset DelayModelVariableParameters table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("delayModelVariableParametersId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("time", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("ut1_utc", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("iat_utc", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("timeType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("gstAtUt0", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("earthRotationRate", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("polarOffsets", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("polarOffsetsType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("delayModelFixedParametersId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("nutationInLongitude", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("nutationInLongitudeRate", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("nutationInObliquity", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("nutationInObliquityRate", "blabla"));
  		  		
    }

    ASDM_DELAYMODELVARIABLEPARAMETERS::~ASDM_DELAYMODELVARIABLEPARAMETERS() {
    }

    const TableDesc& ASDM_DELAYMODELVARIABLEPARAMETERS::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_DELAYMODELVARIABLEPARAMETERS::fill(const ASDM& asdm) {
        vector<DelayModelVariableParametersRow*> rows = asdm.getDelayModelVariableParameters().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> delayModelVariableParametersId(*table_p_, "delayModelVariableParametersId");             
  		
  		
        ScalarColumn<double> time(*table_p_, "time");             
  		
        ScalarColumn<double> ut1_utc(*table_p_, "ut1_utc");             
  		
        ScalarColumn<double> iat_utc(*table_p_, "iat_utc");             
  		
        ScalarColumn<String> timeType(*table_p_, "timeType");             
  		
        ScalarColumn<double> gstAtUt0(*table_p_, "gstAtUt0");             
  		
        ScalarColumn<double> earthRotationRate(*table_p_, "earthRotationRate");             
  		
        ArrayColumn<double> polarOffsets(*table_p_, "polarOffsets");             
  		
        ScalarColumn<String> polarOffsetsType(*table_p_, "polarOffsetsType");             
  		
        ScalarColumn<String> delayModelFixedParametersId(*table_p_, "delayModelFixedParametersId");             
  		
  		
        ScalarColumn<double> nutationInLongitude(*table_p_, "nutationInLongitude");             
  		
        ScalarColumn<double> nutationInLongitudeRate(*table_p_, "nutationInLongitudeRate");             
  		
        ScalarColumn<double> nutationInObliquity(*table_p_, "nutationInObliquity");             
  		
        ScalarColumn<double> nutationInObliquityRate(*table_p_, "nutationInObliquityRate");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            delayModelVariableParametersId.put(rowIndex, rows.at(i)->getDelayModelVariableParametersId().toString());
	

		
	
            time.put(rowIndex, rows.at(i)->getTime().get()/(1.0e9));
	

	
            ut1_utc.put(rowIndex, rows.at(i)->getUt1_utc());
	

	
            iat_utc.put(rowIndex, rows.at(i)->getIat_utc());
	

	
            timeType.put(rowIndex, CDifferenceType::name(rows.at(i)->getTimeType()));
	

	
            gstAtUt0.put(rowIndex, rows.at(i)->getGstAtUt0().get());
	

	
            earthRotationRate.put(rowIndex, rows.at(i)->getEarthRotationRate().get());
	

	
            polarOffsets.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getPolarOffsets()));
	

	
            polarOffsetsType.put(rowIndex, CDifferenceType::name(rows.at(i)->getPolarOffsetsType()));
	

	
            delayModelFixedParametersId.put(rowIndex, rows.at(i)->getDelayModelFixedParametersId().toString());
	

		
	
            if (rows.at(i)->isNutationInLongitudeExists())
                nutationInLongitude.put(rowIndex, rows.at(i)->getNutationInLongitude().get());
	

	
            if (rows.at(i)->isNutationInLongitudeRateExists())
                nutationInLongitudeRate.put(rowIndex, rows.at(i)->getNutationInLongitudeRate().get());
	

	
            if (rows.at(i)->isNutationInObliquityExists())
                nutationInObliquity.put(rowIndex, rows.at(i)->getNutationInObliquity().get());
	

	
            if (rows.at(i)->isNutationInObliquityRateExists())
                nutationInObliquityRate.put(rowIndex, rows.at(i)->getNutationInObliquityRate().get());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_DOPPLER::ASDM_DOPPLER() {
        name_ = "ASDM_DOPPLER";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Doppler table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("dopplerId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("sourceId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("transitionIndex", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("velDef", "blabla"));
  		
  		  		
    }

    ASDM_DOPPLER::~ASDM_DOPPLER() {
    }

    const TableDesc& ASDM_DOPPLER::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
	
    void ASDM_DOPPLER::fill(const ASDM& asdm) {
        vector<DopplerRow*> rows = asdm.getDoppler().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<int> dopplerId(*table_p_, "dopplerId");             
  		
        ScalarColumn<int> sourceId(*table_p_, "sourceId");             
  		
  		
        ScalarColumn<int> transitionIndex(*table_p_, "transitionIndex");             
  		
        ScalarColumn<String> velDef(*table_p_, "velDef");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            dopplerId.put(rowIndex, rows.at(i)->getDopplerId());
	

	
            sourceId.put(rowIndex, rows.at(i)->getSourceId());
	

		
	
            transitionIndex.put(rowIndex, rows.at(i)->getTransitionIndex());
	

	
            velDef.put(rowIndex, CDopplerReferenceCode::name(rows.at(i)->getVelDef()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_EPHEMERIS::ASDM_EPHEMERIS() {
        name_ = "ASDM_EPHEMERIS";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Ephemeris table";
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("ephemerisId", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("observerLocation", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("equinoxEquator", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPolyDir", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("dir", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPolyDist", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("distance", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("timeOrigin", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("origin", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPolyRadVel", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("radVel", "blabla"));
  		  		
    }

    ASDM_EPHEMERIS::~ASDM_EPHEMERIS() {
    }

    const TableDesc& ASDM_EPHEMERIS::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_EPHEMERIS::fill(const ASDM& asdm) {
        vector<EphemerisRow*> rows = asdm.getEphemeris().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
        ScalarColumn<int> ephemerisId(*table_p_, "ephemerisId");             
  		
  		
        ArrayColumn<double> observerLocation(*table_p_, "observerLocation");             
  		
        ScalarColumn<double> equinoxEquator(*table_p_, "equinoxEquator");             
  		
        ScalarColumn<int> numPolyDir(*table_p_, "numPolyDir");             
  		
        ArrayColumn<double> dir(*table_p_, "dir");             
  		
        ScalarColumn<int> numPolyDist(*table_p_, "numPolyDist");             
  		
        ArrayColumn<double> distance(*table_p_, "distance");             
  		
        ScalarColumn<double> timeOrigin(*table_p_, "timeOrigin");             
  		
        ScalarColumn<String> origin(*table_p_, "origin");             
  		
  		
        ScalarColumn<int> numPolyRadVel(*table_p_, "numPolyRadVel");             
  		
        ArrayColumn<double> radVel(*table_p_, "radVel");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

	
            ephemerisId.put(rowIndex, rows.at(i)->getEphemerisId());
	

		
	
            observerLocation.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getObserverLocation()));
	

	
            equinoxEquator.put(rowIndex, rows.at(i)->getEquinoxEquator());
	

	
            numPolyDir.put(rowIndex, rows.at(i)->getNumPolyDir());
	

	
            dir.put(rowIndex, basic2CASA2D<double,double>(rows.at(i)->getDir()));
	

	
            numPolyDist.put(rowIndex, rows.at(i)->getNumPolyDist());
	

	
            distance.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getDistance()));
	

	
            timeOrigin.put(rowIndex, rows.at(i)->getTimeOrigin().get()/(1.0e9));
	

	
            origin.put(rowIndex, rows.at(i)->getOrigin());
	

		
	
            if (rows.at(i)->isNumPolyRadVelExists())
                numPolyRadVel.put(rowIndex, rows.at(i)->getNumPolyRadVel());
	

	
            if (rows.at(i)->isRadVelExists())
                radVel.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getRadVel()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_EXECBLOCK::ASDM_EXECBLOCK() {
        name_ = "ASDM_EXECBLOCK";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset ExecBlock table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("execBlockNum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("projectUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("configName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("telescopeName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("observerName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numObservingLog", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("observingLog", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sessionReference", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("baseRangeMin", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("baseRangeMax", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("baseRmsMinor", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("baseRmsMajor", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("basePa", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("aborted", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sBSummaryId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("releaseDate", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("schedulerMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("siteAltitude", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("siteLongitude", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("siteLatitude", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("observingScript", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("observingScriptUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("scaleId", "blabla"));
  		  		
    }

    ASDM_EXECBLOCK::~ASDM_EXECBLOCK() {
    }

    const TableDesc& ASDM_EXECBLOCK::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_EXECBLOCK::fill(const ASDM& asdm) {
        vector<ExecBlockRow*> rows = asdm.getExecBlock().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> execBlockId(*table_p_, "execBlockId");             
  		
  		
        ScalarColumn<double> startTime(*table_p_, "startTime");             
  		
        ScalarColumn<double> endTime(*table_p_, "endTime");             
  		
        ScalarColumn<int> execBlockNum(*table_p_, "execBlockNum");             
  		
        ScalarColumn<String> execBlockUID(*table_p_, "execBlockUID");             
  		
        ScalarColumn<String> projectUID(*table_p_, "projectUID");             
  		
        ScalarColumn<String> configName(*table_p_, "configName");             
  		
        ScalarColumn<String> telescopeName(*table_p_, "telescopeName");             
  		
        ScalarColumn<String> observerName(*table_p_, "observerName");             
  		
        ScalarColumn<int> numObservingLog(*table_p_, "numObservingLog");             
  		
        ArrayColumn<String> observingLog(*table_p_, "observingLog");             
  		
        ScalarColumn<String> sessionReference(*table_p_, "sessionReference");             
  		
        ScalarColumn<double> baseRangeMin(*table_p_, "baseRangeMin");             
  		
        ScalarColumn<double> baseRangeMax(*table_p_, "baseRangeMax");             
  		
        ScalarColumn<double> baseRmsMinor(*table_p_, "baseRmsMinor");             
  		
        ScalarColumn<double> baseRmsMajor(*table_p_, "baseRmsMajor");             
  		
        ScalarColumn<double> basePa(*table_p_, "basePa");             
  		
        ScalarColumn<bool> aborted(*table_p_, "aborted");             
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ArrayColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> sBSummaryId(*table_p_, "sBSummaryId");             
  		
  		
        ScalarColumn<double> releaseDate(*table_p_, "releaseDate");             
  		
        ScalarColumn<String> schedulerMode(*table_p_, "schedulerMode");             
  		
        ScalarColumn<double> siteAltitude(*table_p_, "siteAltitude");             
  		
        ScalarColumn<double> siteLongitude(*table_p_, "siteLongitude");             
  		
        ScalarColumn<double> siteLatitude(*table_p_, "siteLatitude");             
  		
        ScalarColumn<String> observingScript(*table_p_, "observingScript");             
  		
        ScalarColumn<String> observingScriptUID(*table_p_, "observingScriptUID");             
  		
        ScalarColumn<String> scaleId(*table_p_, "scaleId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            execBlockId.put(rowIndex, rows.at(i)->getExecBlockId().toString());
	

		
	
            startTime.put(rowIndex, rows.at(i)->getStartTime().get()/(1.0e9));
	

	
            endTime.put(rowIndex, rows.at(i)->getEndTime().get()/(1.0e9));
	

	
            execBlockNum.put(rowIndex, rows.at(i)->getExecBlockNum());
	

	
            execBlockUID.put(rowIndex, rows.at(i)->getExecBlockUID().toString());
	

	
            projectUID.put(rowIndex, rows.at(i)->getProjectUID().toString());
	

	
            configName.put(rowIndex, rows.at(i)->getConfigName());
	

	
            telescopeName.put(rowIndex, rows.at(i)->getTelescopeName());
	

	
            observerName.put(rowIndex, rows.at(i)->getObserverName());
	

	
            numObservingLog.put(rowIndex, rows.at(i)->getNumObservingLog());
	

	
            observingLog.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getObservingLog()));
	

	
            sessionReference.put(rowIndex, rows.at(i)->getSessionReference().toString());
	

	
            baseRangeMin.put(rowIndex, rows.at(i)->getBaseRangeMin().get());
	

	
            baseRangeMax.put(rowIndex, rows.at(i)->getBaseRangeMax().get());
	

	
            baseRmsMinor.put(rowIndex, rows.at(i)->getBaseRmsMinor().get());
	

	
            baseRmsMajor.put(rowIndex, rows.at(i)->getBaseRmsMajor().get());
	

	
            basePa.put(rowIndex, rows.at(i)->getBasePa().get());
	

	
            aborted.put(rowIndex, rows.at(i)->getAborted());
	

	
            numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            antennaId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getAntennaId()));
	

	
            sBSummaryId.put(rowIndex, rows.at(i)->getSBSummaryId().toString());
	

		
	
            if (rows.at(i)->isReleaseDateExists())
                releaseDate.put(rowIndex, rows.at(i)->getReleaseDate().get()/(1.0e9));
	

	
            if (rows.at(i)->isSchedulerModeExists())
                schedulerMode.put(rowIndex, rows.at(i)->getSchedulerMode());
	

	
            if (rows.at(i)->isSiteAltitudeExists())
                siteAltitude.put(rowIndex, rows.at(i)->getSiteAltitude().get());
	

	
            if (rows.at(i)->isSiteLongitudeExists())
                siteLongitude.put(rowIndex, rows.at(i)->getSiteLongitude().get());
	

	
            if (rows.at(i)->isSiteLatitudeExists())
                siteLatitude.put(rowIndex, rows.at(i)->getSiteLatitude().get());
	

	
            if (rows.at(i)->isObservingScriptExists())
                observingScript.put(rowIndex, rows.at(i)->getObservingScript());
	

	
            if (rows.at(i)->isObservingScriptUIDExists())
                observingScriptUID.put(rowIndex, rows.at(i)->getObservingScriptUID().toString());
	

	
            if (rows.at(i)->isScaleIdExists())
                scaleId.put(rowIndex, rows.at(i)->getScaleId().toString());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_FEED::ASDM_FEED() {
        name_ = "ASDM_FEED";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Feed table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("feedId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("beamOffset", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("focusReference", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationTypes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polResponse", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("receptorAngle", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("receiverId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("feedNum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("illumOffset", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("position", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("skyCoupling", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numChan", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("skyCouplingSpectrum", "blabla"));
  		  		
    }

    ASDM_FEED::~ASDM_FEED() {
    }

    const TableDesc& ASDM_FEED::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_FEED::fill(const ASDM& asdm) {
        vector<FeedRow*> rows = asdm.getFeed().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
        ScalarColumn<int> feedId(*table_p_, "feedId");             
  		
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<double> beamOffset(*table_p_, "beamOffset");             
  		
        ArrayColumn<double> focusReference(*table_p_, "focusReference");             
  		
        ArrayColumn<String> polarizationTypes(*table_p_, "polarizationTypes");             
  		
        ArrayColumn<String> polResponse(*table_p_, "polResponse");             
  		
        ArrayColumn<double> receptorAngle(*table_p_, "receptorAngle");             
  		
        ArrayColumn<int> receiverId(*table_p_, "receiverId");             
  		
  		
        ScalarColumn<int> feedNum(*table_p_, "feedNum");             
  		
        ArrayColumn<double> illumOffset(*table_p_, "illumOffset");             
  		
        ArrayColumn<double> position(*table_p_, "position");             
  		
        ScalarColumn<float> skyCoupling(*table_p_, "skyCoupling");             
  		
        ScalarColumn<int> numChan(*table_p_, "numChan");             
  		
        ArrayColumn<float> skyCouplingSpectrum(*table_p_, "skyCouplingSpectrum");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

	
            feedId.put(rowIndex, rows.at(i)->getFeedId());
	

		
	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            beamOffset.put(rowIndex, basic2CASA2D<double,double>(rows.at(i)->getBeamOffset()));
	

	
            focusReference.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getFocusReference()));
	

	
            polarizationTypes.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationTypes()));
	

	
            polResponse.put(rowIndex, _2CASAString2D<asdm::Complex,String>(rows.at(i)->getPolResponse()));
	

	
            receptorAngle.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getReceptorAngle()));
	

	
            receiverId.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getReceiverId()));
	

		
	
            if (rows.at(i)->isFeedNumExists())
                feedNum.put(rowIndex, rows.at(i)->getFeedNum());
	

	
            if (rows.at(i)->isIllumOffsetExists())
                illumOffset.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getIllumOffset()));
	

	
            if (rows.at(i)->isPositionExists())
                position.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getPosition()));
	

	
            if (rows.at(i)->isSkyCouplingExists())
                skyCoupling.put(rowIndex, rows.at(i)->getSkyCoupling());
	

	
            if (rows.at(i)->isNumChanExists())
                numChan.put(rowIndex, rows.at(i)->getNumChan());
	

	
            if (rows.at(i)->isSkyCouplingSpectrumExists())
                skyCouplingSpectrum.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getSkyCouplingSpectrum()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_FIELD::ASDM_FIELD() {
        name_ = "ASDM_FIELD";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Field table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("fieldId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("fieldName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPoly", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("delayDir", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("phaseDir", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("referenceDir", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("time", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("code", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("directionCode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("directionEquinox", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("assocNature", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("ephemerisId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("sourceId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("assocFieldId", "blabla"));
  		  		
    }

    ASDM_FIELD::~ASDM_FIELD() {
    }

    const TableDesc& ASDM_FIELD::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
	
    void ASDM_FIELD::fill(const ASDM& asdm) {
        vector<FieldRow*> rows = asdm.getField().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> fieldId(*table_p_, "fieldId");             
  		
  		
        ScalarColumn<String> fieldName(*table_p_, "fieldName");             
  		
        ScalarColumn<int> numPoly(*table_p_, "numPoly");             
  		
        ArrayColumn<double> delayDir(*table_p_, "delayDir");             
  		
        ArrayColumn<double> phaseDir(*table_p_, "phaseDir");             
  		
        ArrayColumn<double> referenceDir(*table_p_, "referenceDir");             
  		
  		
        ScalarColumn<double> time(*table_p_, "time");             
  		
        ScalarColumn<String> code(*table_p_, "code");             
  		
        ScalarColumn<String> directionCode(*table_p_, "directionCode");             
  		
        ScalarColumn<double> directionEquinox(*table_p_, "directionEquinox");             
  		
        ScalarColumn<String> assocNature(*table_p_, "assocNature");             
  		
        ScalarColumn<int> ephemerisId(*table_p_, "ephemerisId");             
  		
        ScalarColumn<int> sourceId(*table_p_, "sourceId");             
  		
        ScalarColumn<String> assocFieldId(*table_p_, "assocFieldId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            fieldId.put(rowIndex, rows.at(i)->getFieldId().toString());
	

		
	
            fieldName.put(rowIndex, rows.at(i)->getFieldName());
	

	
            numPoly.put(rowIndex, rows.at(i)->getNumPoly());
	

	
            delayDir.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getDelayDir()));
	

	
            phaseDir.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getPhaseDir()));
	

	
            referenceDir.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getReferenceDir()));
	

		
	
            if (rows.at(i)->isTimeExists())
                time.put(rowIndex, rows.at(i)->getTime().get()/(1.0e9));
	

	
            if (rows.at(i)->isCodeExists())
                code.put(rowIndex, rows.at(i)->getCode());
	

	
            if (rows.at(i)->isDirectionCodeExists())
                directionCode.put(rowIndex, CDirectionReferenceCode::name(rows.at(i)->getDirectionCode()));
	

	
            if (rows.at(i)->isDirectionEquinoxExists())
                directionEquinox.put(rowIndex, rows.at(i)->getDirectionEquinox().get()/(1.0e9));
	

	
            if (rows.at(i)->isAssocNatureExists())
                assocNature.put(rowIndex, rows.at(i)->getAssocNature());
	

	
            if (rows.at(i)->isEphemerisIdExists())
                ephemerisId.put(rowIndex, rows.at(i)->getEphemerisId());
	

	
            if (rows.at(i)->isSourceIdExists())
                sourceId.put(rowIndex, rows.at(i)->getSourceId());
	

	
            if (rows.at(i)->isAssocFieldIdExists())
                assocFieldId.put(rowIndex, rows.at(i)->getAssocFieldId().toString());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_FLAG::ASDM_FLAG() {
        name_ = "ASDM_FLAG";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Flag table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("flagId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("reason", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAntenna", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("antennaId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPolarizationType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numSpectralWindow", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPairedAntenna", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numChan", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationType", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("channel", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("pairedAntennaId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("spectralWindowId", "blabla"));
  		  		
    }

    ASDM_FLAG::~ASDM_FLAG() {
    }

    const TableDesc& ASDM_FLAG::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
	
    void ASDM_FLAG::fill(const ASDM& asdm) {
        vector<FlagRow*> rows = asdm.getFlag().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> flagId(*table_p_, "flagId");             
  		
  		
        ScalarColumn<double> startTime(*table_p_, "startTime");             
  		
        ScalarColumn<double> endTime(*table_p_, "endTime");             
  		
        ScalarColumn<String> reason(*table_p_, "reason");             
  		
        ScalarColumn<int> numAntenna(*table_p_, "numAntenna");             
  		
        ArrayColumn<String> antennaId(*table_p_, "antennaId");             
  		
  		
        ScalarColumn<int> numPolarizationType(*table_p_, "numPolarizationType");             
  		
        ScalarColumn<int> numSpectralWindow(*table_p_, "numSpectralWindow");             
  		
        ScalarColumn<int> numPairedAntenna(*table_p_, "numPairedAntenna");             
  		
        ScalarColumn<int> numChan(*table_p_, "numChan");             
  		
        ArrayColumn<String> polarizationType(*table_p_, "polarizationType");             
  		
        ArrayColumn<int> channel(*table_p_, "channel");             
  		
        ArrayColumn<String> pairedAntennaId(*table_p_, "pairedAntennaId");             
  		
        ArrayColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            flagId.put(rowIndex, rows.at(i)->getFlagId().toString());
	

		
	
            startTime.put(rowIndex, rows.at(i)->getStartTime().get()/(1.0e9));
	

	
            endTime.put(rowIndex, rows.at(i)->getEndTime().get()/(1.0e9));
	

	
            reason.put(rowIndex, rows.at(i)->getReason());
	

	
            numAntenna.put(rowIndex, rows.at(i)->getNumAntenna());
	

	
            antennaId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getAntennaId()));
	

		
	
            if (rows.at(i)->isNumPolarizationTypeExists())
                numPolarizationType.put(rowIndex, rows.at(i)->getNumPolarizationType());
	

	
            if (rows.at(i)->isNumSpectralWindowExists())
                numSpectralWindow.put(rowIndex, rows.at(i)->getNumSpectralWindow());
	

	
            if (rows.at(i)->isNumPairedAntennaExists())
                numPairedAntenna.put(rowIndex, rows.at(i)->getNumPairedAntenna());
	

	
            if (rows.at(i)->isNumChanExists())
                numChan.put(rowIndex, rows.at(i)->getNumChan());
	

	
            if (rows.at(i)->isPolarizationTypeExists())
                polarizationType.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationType()));
	

	
            if (rows.at(i)->isChannelExists())
                channel.put(rowIndex, basic2CASA2D<int,int>(rows.at(i)->getChannel()));
	

	
            if (rows.at(i)->isPairedAntennaIdExists())
                pairedAntennaId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getPairedAntennaId()));
	

	
            if (rows.at(i)->isSpectralWindowIdExists())
                spectralWindowId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getSpectralWindowId()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_FLAGCMD::ASDM_FLAGCMD() {
        name_ = "ASDM_FLAGCMD";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset FlagCmd table";
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("type", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("reason", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("level", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("severity", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("applied", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("command", "blabla"));
  		
  		  		
    }

    ASDM_FLAGCMD::~ASDM_FLAGCMD() {
    }

    const TableDesc& ASDM_FLAGCMD::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_FLAGCMD::fill(const ASDM& asdm) {
        vector<FlagCmdRow*> rows = asdm.getFlagCmd().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
        ScalarColumn<String> type(*table_p_, "type");             
  		
        ScalarColumn<String> reason(*table_p_, "reason");             
  		
        ScalarColumn<int> level(*table_p_, "level");             
  		
        ScalarColumn<int> severity(*table_p_, "severity");             
  		
        ScalarColumn<bool> applied(*table_p_, "applied");             
  		
        ScalarColumn<String> command(*table_p_, "command");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
	
            type.put(rowIndex, rows.at(i)->getType());
	

	
            reason.put(rowIndex, rows.at(i)->getReason());
	

	
            level.put(rowIndex, rows.at(i)->getLevel());
	

	
            severity.put(rowIndex, rows.at(i)->getSeverity());
	

	
            applied.put(rowIndex, rows.at(i)->getApplied());
	

	
            command.put(rowIndex, rows.at(i)->getCommand());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_FOCUS::ASDM_FOCUS() {
        name_ = "ASDM_FOCUS";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Focus table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("focusTracking", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("focusOffset", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("focusRotationOffset", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("focusModelId", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("measuredFocusPosition", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("measuredFocusRotation", "blabla"));
  		  		
    }

    ASDM_FOCUS::~ASDM_FOCUS() {
    }

    const TableDesc& ASDM_FOCUS::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_FOCUS::fill(const ASDM& asdm) {
        vector<FocusRow*> rows = asdm.getFocus().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
        ScalarColumn<bool> focusTracking(*table_p_, "focusTracking");             
  		
        ArrayColumn<double> focusOffset(*table_p_, "focusOffset");             
  		
        ArrayColumn<double> focusRotationOffset(*table_p_, "focusRotationOffset");             
  		
        ScalarColumn<int> focusModelId(*table_p_, "focusModelId");             
  		
  		
        ArrayColumn<double> measuredFocusPosition(*table_p_, "measuredFocusPosition");             
  		
        ArrayColumn<double> measuredFocusRotation(*table_p_, "measuredFocusRotation");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
	
            focusTracking.put(rowIndex, rows.at(i)->getFocusTracking());
	

	
            focusOffset.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getFocusOffset()));
	

	
            focusRotationOffset.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getFocusRotationOffset()));
	

	
            focusModelId.put(rowIndex, rows.at(i)->getFocusModelId());
	

		
	
            if (rows.at(i)->isMeasuredFocusPositionExists())
                measuredFocusPosition.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getMeasuredFocusPosition()));
	

	
            if (rows.at(i)->isMeasuredFocusRotationExists())
                measuredFocusRotation.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getMeasuredFocusRotation()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_FOCUSMODEL::ASDM_FOCUSMODEL() {
        name_ = "ASDM_FOCUSMODEL";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset FocusModel table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("focusModelId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("polarizationType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numCoeff", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("coeffName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("coeffFormula", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("coeffVal", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("assocNature", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("assocFocusModelId", "blabla"));
  		
  		  		
    }

    ASDM_FOCUSMODEL::~ASDM_FOCUSMODEL() {
    }

    const TableDesc& ASDM_FOCUSMODEL::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_FOCUSMODEL::fill(const ASDM& asdm) {
        vector<FocusModelRow*> rows = asdm.getFocusModel().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<int> focusModelId(*table_p_, "focusModelId");             
  		
  		
        ScalarColumn<String> polarizationType(*table_p_, "polarizationType");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<int> numCoeff(*table_p_, "numCoeff");             
  		
        ArrayColumn<String> coeffName(*table_p_, "coeffName");             
  		
        ArrayColumn<String> coeffFormula(*table_p_, "coeffFormula");             
  		
        ArrayColumn<float> coeffVal(*table_p_, "coeffVal");             
  		
        ScalarColumn<String> assocNature(*table_p_, "assocNature");             
  		
        ScalarColumn<int> assocFocusModelId(*table_p_, "assocFocusModelId");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            focusModelId.put(rowIndex, rows.at(i)->getFocusModelId());
	

		
	
            polarizationType.put(rowIndex, CPolarizationType::name(rows.at(i)->getPolarizationType()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            numCoeff.put(rowIndex, rows.at(i)->getNumCoeff());
	

	
            coeffName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCoeffName()));
	

	
            coeffFormula.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCoeffFormula()));
	

	
            coeffVal.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getCoeffVal()));
	

	
            assocNature.put(rowIndex, rows.at(i)->getAssocNature());
	

	
            assocFocusModelId.put(rowIndex, rows.at(i)->getAssocFocusModelId());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_FREQOFFSET::ASDM_FREQOFFSET() {
        name_ = "ASDM_FREQOFFSET";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset FreqOffset table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("feedId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("offset", "blabla"));
  		
  		  		
    }

    ASDM_FREQOFFSET::~ASDM_FREQOFFSET() {
    }

    const TableDesc& ASDM_FREQOFFSET::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
	
    void ASDM_FREQOFFSET::fill(const ASDM& asdm) {
        vector<FreqOffsetRow*> rows = asdm.getFreqOffset().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
        ScalarColumn<int> feedId(*table_p_, "feedId");             
  		
  		
        ScalarColumn<double> offset(*table_p_, "offset");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

	
            feedId.put(rowIndex, rows.at(i)->getFeedId());
	

		
	
            offset.put(rowIndex, rows.at(i)->getOffset().get());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_GAINTRACKING::ASDM_GAINTRACKING() {
        name_ = "ASDM_GAINTRACKING";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset GainTracking table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("feedId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("attenuator", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("polarizationType", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("samplingLevel", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAttFreq", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("attFreq", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("attSpectrum", "blabla"));
  		  		
    }

    ASDM_GAINTRACKING::~ASDM_GAINTRACKING() {
    }

    const TableDesc& ASDM_GAINTRACKING::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_GAINTRACKING::fill(const ASDM& asdm) {
        vector<GainTrackingRow*> rows = asdm.getGainTracking().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
        ScalarColumn<int> feedId(*table_p_, "feedId");             
  		
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ArrayColumn<float> attenuator(*table_p_, "attenuator");             
  		
        ArrayColumn<String> polarizationType(*table_p_, "polarizationType");             
  		
  		
        ScalarColumn<float> samplingLevel(*table_p_, "samplingLevel");             
  		
        ScalarColumn<int> numAttFreq(*table_p_, "numAttFreq");             
  		
        ArrayColumn<double> attFreq(*table_p_, "attFreq");             
  		
        ArrayColumn<String> attSpectrum(*table_p_, "attSpectrum");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

	
            feedId.put(rowIndex, rows.at(i)->getFeedId());
	

		
	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            attenuator.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getAttenuator()));
	

	
            polarizationType.put(rowIndex, enum2CASA1D<PolarizationType,CPolarizationType>(rows.at(i)->getPolarizationType()));
	

		
	
            if (rows.at(i)->isSamplingLevelExists())
                samplingLevel.put(rowIndex, rows.at(i)->getSamplingLevel());
	

	
            if (rows.at(i)->isNumAttFreqExists())
                numAttFreq.put(rowIndex, rows.at(i)->getNumAttFreq());
	

	
            if (rows.at(i)->isAttFreqExists())
                attFreq.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getAttFreq()));
	

	
            if (rows.at(i)->isAttSpectrumExists())
                attSpectrum.put(rowIndex, _2CASAString1D<asdm::Complex,String>(rows.at(i)->getAttSpectrum()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_HISTORY::ASDM_HISTORY() {
        name_ = "ASDM_HISTORY";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset History table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("time", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("message", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("priority", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("origin", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("objectId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("application", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("cliCommand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("appParms", "blabla"));
  		
  		  		
    }

    ASDM_HISTORY::~ASDM_HISTORY() {
    }

    const TableDesc& ASDM_HISTORY::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_HISTORY::fill(const ASDM& asdm) {
        vector<HistoryRow*> rows = asdm.getHistory().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> execBlockId(*table_p_, "execBlockId");             
  		
        ScalarColumn<double> time(*table_p_, "time");             
  		
  		
        ScalarColumn<String> message(*table_p_, "message");             
  		
        ScalarColumn<String> priority(*table_p_, "priority");             
  		
        ScalarColumn<String> origin(*table_p_, "origin");             
  		
        ScalarColumn<String> objectId(*table_p_, "objectId");             
  		
        ScalarColumn<String> application(*table_p_, "application");             
  		
        ScalarColumn<String> cliCommand(*table_p_, "cliCommand");             
  		
        ScalarColumn<String> appParms(*table_p_, "appParms");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            execBlockId.put(rowIndex, rows.at(i)->getExecBlockId().toString());
	

	
            time.put(rowIndex, rows.at(i)->getTime().get()/(1.0e9));
	

		
	
            message.put(rowIndex, rows.at(i)->getMessage());
	

	
            priority.put(rowIndex, rows.at(i)->getPriority());
	

	
            origin.put(rowIndex, rows.at(i)->getOrigin());
	

	
            objectId.put(rowIndex, rows.at(i)->getObjectId());
	

	
            application.put(rowIndex, rows.at(i)->getApplication());
	

	
            cliCommand.put(rowIndex, rows.at(i)->getCliCommand());
	

	
            appParms.put(rowIndex, rows.at(i)->getAppParms());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_HOLOGRAPHY::ASDM_HOLOGRAPHY() {
        name_ = "ASDM_HOLOGRAPHY";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Holography table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("holographyId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("distance", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("focus", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numCorr", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("type", "blabla"));
  		
  		  		
    }

    ASDM_HOLOGRAPHY::~ASDM_HOLOGRAPHY() {
    }

    const TableDesc& ASDM_HOLOGRAPHY::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
			
		
	
    void ASDM_HOLOGRAPHY::fill(const ASDM& asdm) {
        vector<HolographyRow*> rows = asdm.getHolography().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> holographyId(*table_p_, "holographyId");             
  		
  		
        ScalarColumn<double> distance(*table_p_, "distance");             
  		
        ScalarColumn<double> focus(*table_p_, "focus");             
  		
        ScalarColumn<int> numCorr(*table_p_, "numCorr");             
  		
        ArrayColumn<String> type(*table_p_, "type");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            holographyId.put(rowIndex, rows.at(i)->getHolographyId().toString());
	

		
	
            distance.put(rowIndex, rows.at(i)->getDistance().get());
	

	
            focus.put(rowIndex, rows.at(i)->getFocus().get());
	

	
            numCorr.put(rowIndex, rows.at(i)->getNumCorr());
	

	
            type.put(rowIndex, enum2CASA1D<HolographyChannelType,CHolographyChannelType>(rows.at(i)->getType()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_OBSERVATION::ASDM_OBSERVATION() {
        name_ = "ASDM_OBSERVATION";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Observation table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("observationId", "blabla"));
  		
  		
  		  		
    }

    ASDM_OBSERVATION::~ASDM_OBSERVATION() {
    }

    const TableDesc& ASDM_OBSERVATION::tableDesc() const {
        return tableDesc_;
    }


		
			
		
	
    void ASDM_OBSERVATION::fill(const ASDM& asdm) {
        vector<ObservationRow*> rows = asdm.getObservation().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> observationId(*table_p_, "observationId");             
  		
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            observationId.put(rowIndex, rows.at(i)->getObservationId().toString());
	

		
		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_POINTING::ASDM_POINTING() {
        name_ = "ASDM_POINTING";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Pointing table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numSample", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("encoder", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("pointingTracking", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("usePolynomials", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("timeOrigin", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numTerm", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("pointingDirection", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("target", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("offset", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("pointingModelId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("overTheTop", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("sourceOffset", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sourceOffsetReferenceCode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("sourceOffsetEquinox", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("sampledTimeInterval", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("atmosphericCorrection", "blabla"));
  		  		
    }

    ASDM_POINTING::~ASDM_POINTING() {
    }

    const TableDesc& ASDM_POINTING::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
	
    void ASDM_POINTING::fill(const ASDM& asdm) {
        vector<PointingRow*> rows = asdm.getPointing().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
        ScalarColumn<int> numSample(*table_p_, "numSample");             
  		
        ArrayColumn<double> encoder(*table_p_, "encoder");             
  		
        ScalarColumn<bool> pointingTracking(*table_p_, "pointingTracking");             
  		
        ScalarColumn<bool> usePolynomials(*table_p_, "usePolynomials");             
  		
        ScalarColumn<double> timeOrigin(*table_p_, "timeOrigin");             
  		
        ScalarColumn<int> numTerm(*table_p_, "numTerm");             
  		
        ArrayColumn<double> pointingDirection(*table_p_, "pointingDirection");             
  		
        ArrayColumn<double> target(*table_p_, "target");             
  		
        ArrayColumn<double> offset(*table_p_, "offset");             
  		
        ScalarColumn<int> pointingModelId(*table_p_, "pointingModelId");             
  		
  		
        ScalarColumn<bool> overTheTop(*table_p_, "overTheTop");             
  		
        ArrayColumn<double> sourceOffset(*table_p_, "sourceOffset");             
  		
        ScalarColumn<String> sourceOffsetReferenceCode(*table_p_, "sourceOffsetReferenceCode");             
  		
        ScalarColumn<double> sourceOffsetEquinox(*table_p_, "sourceOffsetEquinox");             
  		
        ArrayColumn<double> sampledTimeInterval(*table_p_, "sampledTimeInterval");             
  		
        ArrayColumn<double> atmosphericCorrection(*table_p_, "atmosphericCorrection");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
	
            numSample.put(rowIndex, rows.at(i)->getNumSample());
	

	
            encoder.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getEncoder()));
	

	
            pointingTracking.put(rowIndex, rows.at(i)->getPointingTracking());
	

	
            usePolynomials.put(rowIndex, rows.at(i)->getUsePolynomials());
	

	
            timeOrigin.put(rowIndex, rows.at(i)->getTimeOrigin().get()/(1.0e9));
	

	
            numTerm.put(rowIndex, rows.at(i)->getNumTerm());
	

	
            pointingDirection.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getPointingDirection()));
	

	
            target.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getTarget()));
	

	
            offset.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getOffset()));
	

	
            pointingModelId.put(rowIndex, rows.at(i)->getPointingModelId());
	

		
	
            if (rows.at(i)->isOverTheTopExists())
                overTheTop.put(rowIndex, rows.at(i)->getOverTheTop());
	

	
            if (rows.at(i)->isSourceOffsetExists())
                sourceOffset.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getSourceOffset()));
	

	
            if (rows.at(i)->isSourceOffsetReferenceCodeExists())
                sourceOffsetReferenceCode.put(rowIndex, CDirectionReferenceCode::name(rows.at(i)->getSourceOffsetReferenceCode()));
	

	
            if (rows.at(i)->isSourceOffsetEquinoxExists())
                sourceOffsetEquinox.put(rowIndex, rows.at(i)->getSourceOffsetEquinox().get()/(1.0e9));
	

	
            if (rows.at(i)->isSampledTimeIntervalExists())
                sampledTimeInterval.put(rowIndex, ati2CASA2D<double>(rows.at(i)->getSampledTimeInterval()));
	

	
            if (rows.at(i)->isAtmosphericCorrectionExists())
                atmosphericCorrection.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getAtmosphericCorrection()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_POINTINGMODEL::ASDM_POINTINGMODEL() {
        name_ = "ASDM_POINTINGMODEL";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset PointingModel table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("pointingModelId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numCoeff", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("coeffName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("coeffVal", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("polarizationType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("assocNature", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("assocPointingModelId", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("coeffFormula", "blabla"));
  		  		
    }

    ASDM_POINTINGMODEL::~ASDM_POINTINGMODEL() {
    }

    const TableDesc& ASDM_POINTINGMODEL::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
		
	
    void ASDM_POINTINGMODEL::fill(const ASDM& asdm) {
        vector<PointingModelRow*> rows = asdm.getPointingModel().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<int> pointingModelId(*table_p_, "pointingModelId");             
  		
  		
        ScalarColumn<int> numCoeff(*table_p_, "numCoeff");             
  		
        ArrayColumn<String> coeffName(*table_p_, "coeffName");             
  		
        ArrayColumn<float> coeffVal(*table_p_, "coeffVal");             
  		
        ScalarColumn<String> polarizationType(*table_p_, "polarizationType");             
  		
        ScalarColumn<String> receiverBand(*table_p_, "receiverBand");             
  		
        ScalarColumn<String> assocNature(*table_p_, "assocNature");             
  		
        ScalarColumn<int> assocPointingModelId(*table_p_, "assocPointingModelId");             
  		
  		
        ArrayColumn<String> coeffFormula(*table_p_, "coeffFormula");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            pointingModelId.put(rowIndex, rows.at(i)->getPointingModelId());
	

		
	
            numCoeff.put(rowIndex, rows.at(i)->getNumCoeff());
	

	
            coeffName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCoeffName()));
	

	
            coeffVal.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getCoeffVal()));
	

	
            polarizationType.put(rowIndex, CPolarizationType::name(rows.at(i)->getPolarizationType()));
	

	
            receiverBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getReceiverBand()));
	

	
            assocNature.put(rowIndex, rows.at(i)->getAssocNature());
	

	
            assocPointingModelId.put(rowIndex, rows.at(i)->getAssocPointingModelId());
	

		
	
            if (rows.at(i)->isCoeffFormulaExists())
                coeffFormula.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getCoeffFormula()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_POLARIZATION::ASDM_POLARIZATION() {
        name_ = "ASDM_POLARIZATION";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Polarization table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("polarizationId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numCorr", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("corrType", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("corrProduct", "blabla"));
  		
  		  		
    }

    ASDM_POLARIZATION::~ASDM_POLARIZATION() {
    }

    const TableDesc& ASDM_POLARIZATION::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
			
			
		
	
    void ASDM_POLARIZATION::fill(const ASDM& asdm) {
        vector<PolarizationRow*> rows = asdm.getPolarization().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> polarizationId(*table_p_, "polarizationId");             
  		
  		
        ScalarColumn<int> numCorr(*table_p_, "numCorr");             
  		
        ArrayColumn<String> corrType(*table_p_, "corrType");             
  		
        ArrayColumn<String> corrProduct(*table_p_, "corrProduct");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            polarizationId.put(rowIndex, rows.at(i)->getPolarizationId().toString());
	

		
	
            numCorr.put(rowIndex, rows.at(i)->getNumCorr());
	

	
            corrType.put(rowIndex, enum2CASA1D<StokesParameter,CStokesParameter>(rows.at(i)->getCorrType()));
	

	
            corrProduct.put(rowIndex, enum2CASA2D<PolarizationType,CPolarizationType>(rows.at(i)->getCorrProduct()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_PROCESSOR::ASDM_PROCESSOR() {
        name_ = "ASDM_PROCESSOR";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Processor table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("processorId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("modeId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("processorType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("processorSubType", "blabla"));
  		
  		  		
    }

    ASDM_PROCESSOR::~ASDM_PROCESSOR() {
    }

    const TableDesc& ASDM_PROCESSOR::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
			
			
		
	
    void ASDM_PROCESSOR::fill(const ASDM& asdm) {
        vector<ProcessorRow*> rows = asdm.getProcessor().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> processorId(*table_p_, "processorId");             
  		
  		
        ScalarColumn<String> modeId(*table_p_, "modeId");             
  		
        ScalarColumn<String> processorType(*table_p_, "processorType");             
  		
        ScalarColumn<String> processorSubType(*table_p_, "processorSubType");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            processorId.put(rowIndex, rows.at(i)->getProcessorId().toString());
	

		
	
            modeId.put(rowIndex, rows.at(i)->getModeId().toString());
	

	
            processorType.put(rowIndex, CProcessorType::name(rows.at(i)->getProcessorType()));
	

	
            processorSubType.put(rowIndex, CProcessorSubType::name(rows.at(i)->getProcessorSubType()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_PULSAR::ASDM_PULSAR() {
        name_ = "ASDM_PULSAR";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Pulsar table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("pulsarId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("refTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("refPulseFreq", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("refPhase", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBin", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPoly", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("phasePoly", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("timeSpan", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("startPhaseBin", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("endPhaseBin", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("dispersionMeasure", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("refFrequency", "blabla"));
  		  		
    }

    ASDM_PULSAR::~ASDM_PULSAR() {
    }

    const TableDesc& ASDM_PULSAR::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_PULSAR::fill(const ASDM& asdm) {
        vector<PulsarRow*> rows = asdm.getPulsar().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> pulsarId(*table_p_, "pulsarId");             
  		
  		
        ScalarColumn<double> refTime(*table_p_, "refTime");             
  		
        ScalarColumn<double> refPulseFreq(*table_p_, "refPulseFreq");             
  		
        ScalarColumn<double> refPhase(*table_p_, "refPhase");             
  		
        ScalarColumn<int> numBin(*table_p_, "numBin");             
  		
  		
        ScalarColumn<int> numPoly(*table_p_, "numPoly");             
  		
        ArrayColumn<double> phasePoly(*table_p_, "phasePoly");             
  		
        ScalarColumn<double> timeSpan(*table_p_, "timeSpan");             
  		
        ArrayColumn<float> startPhaseBin(*table_p_, "startPhaseBin");             
  		
        ArrayColumn<float> endPhaseBin(*table_p_, "endPhaseBin");             
  		
        ScalarColumn<double> dispersionMeasure(*table_p_, "dispersionMeasure");             
  		
        ScalarColumn<double> refFrequency(*table_p_, "refFrequency");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            pulsarId.put(rowIndex, rows.at(i)->getPulsarId().toString());
	

		
	
            refTime.put(rowIndex, rows.at(i)->getRefTime().get()/(1.0e9));
	

	
            refPulseFreq.put(rowIndex, rows.at(i)->getRefPulseFreq().get());
	

	
            refPhase.put(rowIndex, rows.at(i)->getRefPhase());
	

	
            numBin.put(rowIndex, rows.at(i)->getNumBin());
	

		
	
            if (rows.at(i)->isNumPolyExists())
                numPoly.put(rowIndex, rows.at(i)->getNumPoly());
	

	
            if (rows.at(i)->isPhasePolyExists())
                phasePoly.put(rowIndex, basic2CASA1D<double,double>(rows.at(i)->getPhasePoly()));
	

	
            if (rows.at(i)->isTimeSpanExists())
                timeSpan.put(rowIndex, rows.at(i)->getTimeSpan().get()/(1.0e9));
	

	
            if (rows.at(i)->isStartPhaseBinExists())
                startPhaseBin.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getStartPhaseBin()));
	

	
            if (rows.at(i)->isEndPhaseBinExists())
                endPhaseBin.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getEndPhaseBin()));
	

	
            if (rows.at(i)->isDispersionMeasureExists())
                dispersionMeasure.put(rowIndex, rows.at(i)->getDispersionMeasure());
	

	
            if (rows.at(i)->isRefFrequencyExists())
                refFrequency.put(rowIndex, rows.at(i)->getRefFrequency().get());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_RECEIVER::ASDM_RECEIVER() {
        name_ = "ASDM_RECEIVER";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Receiver table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("receiverId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("name", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numLO", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("frequencyBand", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("freqLO", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("receiverSideband", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("sidebandLO", "blabla"));
  		
  		  		
    }

    ASDM_RECEIVER::~ASDM_RECEIVER() {
    }

    const TableDesc& ASDM_RECEIVER::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
			
		
			
			
		
	
    void ASDM_RECEIVER::fill(const ASDM& asdm) {
        vector<ReceiverRow*> rows = asdm.getReceiver().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<int> receiverId(*table_p_, "receiverId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
        ScalarColumn<String> name(*table_p_, "name");             
  		
        ScalarColumn<int> numLO(*table_p_, "numLO");             
  		
        ScalarColumn<String> frequencyBand(*table_p_, "frequencyBand");             
  		
        ArrayColumn<double> freqLO(*table_p_, "freqLO");             
  		
        ScalarColumn<String> receiverSideband(*table_p_, "receiverSideband");             
  		
        ArrayColumn<String> sidebandLO(*table_p_, "sidebandLO");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            receiverId.put(rowIndex, rows.at(i)->getReceiverId());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
	
            name.put(rowIndex, rows.at(i)->getName());
	

	
            numLO.put(rowIndex, rows.at(i)->getNumLO());
	

	
            frequencyBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getFrequencyBand()));
	

	
            freqLO.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFreqLO()));
	

	
            receiverSideband.put(rowIndex, CReceiverSideband::name(rows.at(i)->getReceiverSideband()));
	

	
            sidebandLO.put(rowIndex, enum2CASA1D<NetSideband,CNetSideband>(rows.at(i)->getSidebandLO()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SBSUMMARY::ASDM_SBSUMMARY() {
        name_ = "ASDM_SBSUMMARY";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset SBSummary table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sBSummaryId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sbSummaryUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("projectUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("obsUnitSetUID", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("frequency", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("frequencyBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sbType", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("sbDuration", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numObservingMode", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("observingMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numberRepeats", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numScienceGoal", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("scienceGoal", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numWeatherConstraint", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("weatherConstraint", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("centerDirection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("centerDirectionCode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("centerDirectionEquinox", "blabla"));
  		  		
    }

    ASDM_SBSUMMARY::~ASDM_SBSUMMARY() {
    }

    const TableDesc& ASDM_SBSUMMARY::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
	
    void ASDM_SBSUMMARY::fill(const ASDM& asdm) {
        vector<SBSummaryRow*> rows = asdm.getSBSummary().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> sBSummaryId(*table_p_, "sBSummaryId");             
  		
  		
        ScalarColumn<String> sbSummaryUID(*table_p_, "sbSummaryUID");             
  		
        ScalarColumn<String> projectUID(*table_p_, "projectUID");             
  		
        ScalarColumn<String> obsUnitSetUID(*table_p_, "obsUnitSetUID");             
  		
        ScalarColumn<double> frequency(*table_p_, "frequency");             
  		
        ScalarColumn<String> frequencyBand(*table_p_, "frequencyBand");             
  		
        ScalarColumn<String> sbType(*table_p_, "sbType");             
  		
        ScalarColumn<double> sbDuration(*table_p_, "sbDuration");             
  		
        ScalarColumn<int> numObservingMode(*table_p_, "numObservingMode");             
  		
        ArrayColumn<String> observingMode(*table_p_, "observingMode");             
  		
        ScalarColumn<int> numberRepeats(*table_p_, "numberRepeats");             
  		
        ScalarColumn<int> numScienceGoal(*table_p_, "numScienceGoal");             
  		
        ArrayColumn<String> scienceGoal(*table_p_, "scienceGoal");             
  		
        ScalarColumn<int> numWeatherConstraint(*table_p_, "numWeatherConstraint");             
  		
        ArrayColumn<String> weatherConstraint(*table_p_, "weatherConstraint");             
  		
  		
        ArrayColumn<double> centerDirection(*table_p_, "centerDirection");             
  		
        ScalarColumn<String> centerDirectionCode(*table_p_, "centerDirectionCode");             
  		
        ScalarColumn<double> centerDirectionEquinox(*table_p_, "centerDirectionEquinox");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            sBSummaryId.put(rowIndex, rows.at(i)->getSBSummaryId().toString());
	

		
	
            sbSummaryUID.put(rowIndex, rows.at(i)->getSbSummaryUID().toString());
	

	
            projectUID.put(rowIndex, rows.at(i)->getProjectUID().toString());
	

	
            obsUnitSetUID.put(rowIndex, rows.at(i)->getObsUnitSetUID().toString());
	

	
            frequency.put(rowIndex, rows.at(i)->getFrequency());
	

	
            frequencyBand.put(rowIndex, CReceiverBand::name(rows.at(i)->getFrequencyBand()));
	

	
            sbType.put(rowIndex, CSBType::name(rows.at(i)->getSbType()));
	

	
            sbDuration.put(rowIndex, rows.at(i)->getSbDuration().get()/(1.0e9));
	

	
            numObservingMode.put(rowIndex, rows.at(i)->getNumObservingMode());
	

	
            observingMode.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getObservingMode()));
	

	
            numberRepeats.put(rowIndex, rows.at(i)->getNumberRepeats());
	

	
            numScienceGoal.put(rowIndex, rows.at(i)->getNumScienceGoal());
	

	
            scienceGoal.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getScienceGoal()));
	

	
            numWeatherConstraint.put(rowIndex, rows.at(i)->getNumWeatherConstraint());
	

	
            weatherConstraint.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getWeatherConstraint()));
	

		
	
            if (rows.at(i)->isCenterDirectionExists())
                centerDirection.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getCenterDirection()));
	

	
            if (rows.at(i)->isCenterDirectionCodeExists())
                centerDirectionCode.put(rowIndex, CDirectionReferenceCode::name(rows.at(i)->getCenterDirectionCode()));
	

	
            if (rows.at(i)->isCenterDirectionEquinoxExists())
                centerDirectionEquinox.put(rowIndex, rows.at(i)->getCenterDirectionEquinox().get()/(1.0e9));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SCALE::ASDM_SCALE() {
        name_ = "ASDM_SCALE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Scale table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("scaleId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("timeScale", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("crossDataScale", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("autoDataScale", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("weightType", "blabla"));
  		
  		  		
    }

    ASDM_SCALE::~ASDM_SCALE() {
    }

    const TableDesc& ASDM_SCALE::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
			
		
			
			
		
			
			
		
	
    void ASDM_SCALE::fill(const ASDM& asdm) {
        vector<ScaleRow*> rows = asdm.getScale().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> scaleId(*table_p_, "scaleId");             
  		
  		
        ScalarColumn<String> timeScale(*table_p_, "timeScale");             
  		
        ScalarColumn<String> crossDataScale(*table_p_, "crossDataScale");             
  		
        ScalarColumn<String> autoDataScale(*table_p_, "autoDataScale");             
  		
        ScalarColumn<String> weightType(*table_p_, "weightType");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            scaleId.put(rowIndex, rows.at(i)->getScaleId().toString());
	

		
	
            timeScale.put(rowIndex, CTimeScale::name(rows.at(i)->getTimeScale()));
	

	
            crossDataScale.put(rowIndex, CDataScale::name(rows.at(i)->getCrossDataScale()));
	

	
            autoDataScale.put(rowIndex, CDataScale::name(rows.at(i)->getAutoDataScale()));
	

	
            weightType.put(rowIndex, CWeightType::name(rows.at(i)->getWeightType()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SCAN::ASDM_SCAN() {
        name_ = "ASDM_SCAN";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Scan table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("scanNumber", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numIntent", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numSubscan", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("scanIntent", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("calDataType", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("calibrationOnLine", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("calibrationFunction", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("calibrationSet", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("calPattern", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numField", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("fieldName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sourceName", "blabla"));
  		  		
    }

    ASDM_SCAN::~ASDM_SCAN() {
    }

    const TableDesc& ASDM_SCAN::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
			
		
			
			
		
			
			
		
			
		
			
		
			
		
	
    void ASDM_SCAN::fill(const ASDM& asdm) {
        vector<ScanRow*> rows = asdm.getScan().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> execBlockId(*table_p_, "execBlockId");             
  		
        ScalarColumn<int> scanNumber(*table_p_, "scanNumber");             
  		
  		
        ScalarColumn<double> startTime(*table_p_, "startTime");             
  		
        ScalarColumn<double> endTime(*table_p_, "endTime");             
  		
        ScalarColumn<int> numIntent(*table_p_, "numIntent");             
  		
        ScalarColumn<int> numSubscan(*table_p_, "numSubscan");             
  		
        ArrayColumn<String> scanIntent(*table_p_, "scanIntent");             
  		
        ArrayColumn<String> calDataType(*table_p_, "calDataType");             
  		
        ArrayColumn<bool> calibrationOnLine(*table_p_, "calibrationOnLine");             
  		
  		
        ArrayColumn<String> calibrationFunction(*table_p_, "calibrationFunction");             
  		
        ArrayColumn<String> calibrationSet(*table_p_, "calibrationSet");             
  		
        ArrayColumn<String> calPattern(*table_p_, "calPattern");             
  		
        ScalarColumn<int> numField(*table_p_, "numField");             
  		
        ArrayColumn<String> fieldName(*table_p_, "fieldName");             
  		
        ScalarColumn<String> sourceName(*table_p_, "sourceName");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            execBlockId.put(rowIndex, rows.at(i)->getExecBlockId().toString());
	

	
            scanNumber.put(rowIndex, rows.at(i)->getScanNumber());
	

		
	
            startTime.put(rowIndex, rows.at(i)->getStartTime().get()/(1.0e9));
	

	
            endTime.put(rowIndex, rows.at(i)->getEndTime().get()/(1.0e9));
	

	
            numIntent.put(rowIndex, rows.at(i)->getNumIntent());
	

	
            numSubscan.put(rowIndex, rows.at(i)->getNumSubscan());
	

	
            scanIntent.put(rowIndex, enum2CASA1D<ScanIntent,CScanIntent>(rows.at(i)->getScanIntent()));
	

	
            calDataType.put(rowIndex, enum2CASA1D<CalDataOrigin,CCalDataOrigin>(rows.at(i)->getCalDataType()));
	

	
            calibrationOnLine.put(rowIndex, basic2CASA1D<bool,bool>(rows.at(i)->getCalibrationOnLine()));
	

		
	
            if (rows.at(i)->isCalibrationFunctionExists())
                calibrationFunction.put(rowIndex, enum2CASA1D<CalibrationFunction,CCalibrationFunction>(rows.at(i)->getCalibrationFunction()));
	

	
            if (rows.at(i)->isCalibrationSetExists())
                calibrationSet.put(rowIndex, enum2CASA1D<CalibrationSet,CCalibrationSet>(rows.at(i)->getCalibrationSet()));
	

	
            if (rows.at(i)->isCalPatternExists())
                calPattern.put(rowIndex, enum2CASA1D<AntennaMotionPattern,CAntennaMotionPattern>(rows.at(i)->getCalPattern()));
	

	
            if (rows.at(i)->isNumFieldExists())
                numField.put(rowIndex, rows.at(i)->getNumField());
	

	
            if (rows.at(i)->isFieldNameExists())
                fieldName.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getFieldName()));
	

	
            if (rows.at(i)->isSourceNameExists())
                sourceName.put(rowIndex, rows.at(i)->getSourceName());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SEEING::ASDM_SEEING() {
        name_ = "ASDM_SEEING";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Seeing table";
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBaseLength", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("baseLength", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("phaseRms", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("seeing", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("exponent", "blabla"));
  		
  		  		
    }

    ASDM_SEEING::~ASDM_SEEING() {
    }

    const TableDesc& ASDM_SEEING::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_SEEING::fill(const ASDM& asdm) {
        vector<SeeingRow*> rows = asdm.getSeeing().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
        ScalarColumn<int> numBaseLength(*table_p_, "numBaseLength");             
  		
        ArrayColumn<double> baseLength(*table_p_, "baseLength");             
  		
        ArrayColumn<double> phaseRms(*table_p_, "phaseRms");             
  		
        ScalarColumn<float> seeing(*table_p_, "seeing");             
  		
        ScalarColumn<float> exponent(*table_p_, "exponent");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
	
            numBaseLength.put(rowIndex, rows.at(i)->getNumBaseLength());
	

	
            baseLength.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getBaseLength()));
	

	
            phaseRms.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getPhaseRms()));
	

	
            seeing.put(rowIndex, rows.at(i)->getSeeing());
	

	
            exponent.put(rowIndex, rows.at(i)->getExponent());
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SOURCE::ASDM_SOURCE() {
        name_ = "ASDM_SOURCE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Source table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("sourceId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("code", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("direction", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("properMotion", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sourceName", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("directionCode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("directionEquinox", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("calibrationGroup", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("catalog", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("deltaVel", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("position", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numLines", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("transition", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("restFrequency", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("sysVel", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("rangeVel", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sourceModel", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("frequencyRefCode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numFreq", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numStokes", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequency", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("frequencyInterval", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("stokesParameter", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("flux", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("fluxErr", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("positionAngle", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("positionAngleErr", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("size", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("sizeErr", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("velRefCode", "blabla"));
  		  		
    }

    ASDM_SOURCE::~ASDM_SOURCE() {
    }

    const TableDesc& ASDM_SOURCE::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
	
    void ASDM_SOURCE::fill(const ASDM& asdm) {
        vector<SourceRow*> rows = asdm.getSource().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<int> sourceId(*table_p_, "sourceId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
  		
        ScalarColumn<String> code(*table_p_, "code");             
  		
        ArrayColumn<double> direction(*table_p_, "direction");             
  		
        ArrayColumn<double> properMotion(*table_p_, "properMotion");             
  		
        ScalarColumn<String> sourceName(*table_p_, "sourceName");             
  		
  		
        ScalarColumn<String> directionCode(*table_p_, "directionCode");             
  		
        ScalarColumn<double> directionEquinox(*table_p_, "directionEquinox");             
  		
        ScalarColumn<int> calibrationGroup(*table_p_, "calibrationGroup");             
  		
        ScalarColumn<String> catalog(*table_p_, "catalog");             
  		
        ScalarColumn<double> deltaVel(*table_p_, "deltaVel");             
  		
        ArrayColumn<double> position(*table_p_, "position");             
  		
        ScalarColumn<int> numLines(*table_p_, "numLines");             
  		
        ArrayColumn<String> transition(*table_p_, "transition");             
  		
        ArrayColumn<double> restFrequency(*table_p_, "restFrequency");             
  		
        ArrayColumn<double> sysVel(*table_p_, "sysVel");             
  		
        ArrayColumn<double> rangeVel(*table_p_, "rangeVel");             
  		
        ScalarColumn<String> sourceModel(*table_p_, "sourceModel");             
  		
        ScalarColumn<String> frequencyRefCode(*table_p_, "frequencyRefCode");             
  		
        ScalarColumn<int> numFreq(*table_p_, "numFreq");             
  		
        ScalarColumn<int> numStokes(*table_p_, "numStokes");             
  		
        ArrayColumn<double> frequency(*table_p_, "frequency");             
  		
        ArrayColumn<double> frequencyInterval(*table_p_, "frequencyInterval");             
  		
        ArrayColumn<String> stokesParameter(*table_p_, "stokesParameter");             
  		
        ArrayColumn<double> flux(*table_p_, "flux");             
  		
        ArrayColumn<double> fluxErr(*table_p_, "fluxErr");             
  		
        ArrayColumn<double> positionAngle(*table_p_, "positionAngle");             
  		
        ArrayColumn<double> positionAngleErr(*table_p_, "positionAngleErr");             
  		
        ArrayColumn<double> size(*table_p_, "size");             
  		
        ArrayColumn<double> sizeErr(*table_p_, "sizeErr");             
  		
        ScalarColumn<String> velRefCode(*table_p_, "velRefCode");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            sourceId.put(rowIndex, rows.at(i)->getSourceId());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

		
	
            code.put(rowIndex, rows.at(i)->getCode());
	

	
            direction.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getDirection()));
	

	
            properMotion.put(rowIndex, ext2CASA1D<AngularRate,double>(rows.at(i)->getProperMotion()));
	

	
            sourceName.put(rowIndex, rows.at(i)->getSourceName());
	

		
	
            if (rows.at(i)->isDirectionCodeExists())
                directionCode.put(rowIndex, CDirectionReferenceCode::name(rows.at(i)->getDirectionCode()));
	

	
            if (rows.at(i)->isDirectionEquinoxExists())
                directionEquinox.put(rowIndex, rows.at(i)->getDirectionEquinox().get()/(1.0e9));
	

	
            if (rows.at(i)->isCalibrationGroupExists())
                calibrationGroup.put(rowIndex, rows.at(i)->getCalibrationGroup());
	

	
            if (rows.at(i)->isCatalogExists())
                catalog.put(rowIndex, rows.at(i)->getCatalog());
	

	
            if (rows.at(i)->isDeltaVelExists())
                deltaVel.put(rowIndex, rows.at(i)->getDeltaVel().get());
	

	
            if (rows.at(i)->isPositionExists())
                position.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getPosition()));
	

	
            if (rows.at(i)->isNumLinesExists())
                numLines.put(rowIndex, rows.at(i)->getNumLines());
	

	
            if (rows.at(i)->isTransitionExists())
                transition.put(rowIndex, basic2CASA1D<string,String>(rows.at(i)->getTransition()));
	

	
            if (rows.at(i)->isRestFrequencyExists())
                restFrequency.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getRestFrequency()));
	

	
            if (rows.at(i)->isSysVelExists())
                sysVel.put(rowIndex, ext2CASA1D<Speed,double>(rows.at(i)->getSysVel()));
	

	
            if (rows.at(i)->isRangeVelExists())
                rangeVel.put(rowIndex, ext2CASA1D<Speed,double>(rows.at(i)->getRangeVel()));
	

	
            if (rows.at(i)->isSourceModelExists())
                sourceModel.put(rowIndex, CSourceModel::name(rows.at(i)->getSourceModel()));
	

	
            if (rows.at(i)->isFrequencyRefCodeExists())
                frequencyRefCode.put(rowIndex, CFrequencyReferenceCode::name(rows.at(i)->getFrequencyRefCode()));
	

	
            if (rows.at(i)->isNumFreqExists())
                numFreq.put(rowIndex, rows.at(i)->getNumFreq());
	

	
            if (rows.at(i)->isNumStokesExists())
                numStokes.put(rowIndex, rows.at(i)->getNumStokes());
	

	
            if (rows.at(i)->isFrequencyExists())
                frequency.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequency()));
	

	
            if (rows.at(i)->isFrequencyIntervalExists())
                frequencyInterval.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFrequencyInterval()));
	

	
            if (rows.at(i)->isStokesParameterExists())
                stokesParameter.put(rowIndex, enum2CASA1D<StokesParameter,CStokesParameter>(rows.at(i)->getStokesParameter()));
	

	
            if (rows.at(i)->isFluxExists())
                flux.put(rowIndex, ext2CASA2D<Flux,double>(rows.at(i)->getFlux()));
	

	
            if (rows.at(i)->isFluxErrExists())
                fluxErr.put(rowIndex, ext2CASA2D<Flux,double>(rows.at(i)->getFluxErr()));
	

	
            if (rows.at(i)->isPositionAngleExists())
                positionAngle.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getPositionAngle()));
	

	
            if (rows.at(i)->isPositionAngleErrExists())
                positionAngleErr.put(rowIndex, ext2CASA1D<Angle,double>(rows.at(i)->getPositionAngleErr()));
	

	
            if (rows.at(i)->isSizeExists())
                size.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getSize()));
	

	
            if (rows.at(i)->isSizeErrExists())
                sizeErr.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getSizeErr()));
	

	
            if (rows.at(i)->isVelRefCodeExists())
                velRefCode.put(rowIndex, CRadialVelocityReferenceCode::name(rows.at(i)->getVelRefCode()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SPECTRALWINDOW::ASDM_SPECTRALWINDOW() {
        name_ = "ASDM_SPECTRALWINDOW";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset SpectralWindow table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("basebandName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("netSideband", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numChan", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("refFreq", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("sidebandProcessingMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("totBandwidth", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("windowFunction", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("chanFreqStart", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("chanFreqStep", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("chanFreqArray", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("chanWidth", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("chanWidthArray", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("correlationBit", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("effectiveBw", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("effectiveBwArray", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("freqGroup", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("freqGroupName", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<bool>("lineArray", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("measFreqRef", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("name", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("oversampling", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("quantization", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("refChan", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("resolution", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("resolutionArray", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numAssocValues", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("assocNature", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("assocSpectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("imageSpectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("dopplerId", "blabla"));
  		  		
    }

    ASDM_SPECTRALWINDOW::~ASDM_SPECTRALWINDOW() {
    }

    const TableDesc& ASDM_SPECTRALWINDOW::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
			
		
			
		
			
		
			
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
	
    void ASDM_SPECTRALWINDOW::fill(const ASDM& asdm) {
        vector<SpectralWindowRow*> rows = asdm.getSpectralWindow().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
  		
        ScalarColumn<String> basebandName(*table_p_, "basebandName");             
  		
        ScalarColumn<String> netSideband(*table_p_, "netSideband");             
  		
        ScalarColumn<int> numChan(*table_p_, "numChan");             
  		
        ScalarColumn<double> refFreq(*table_p_, "refFreq");             
  		
        ScalarColumn<String> sidebandProcessingMode(*table_p_, "sidebandProcessingMode");             
  		
        ScalarColumn<double> totBandwidth(*table_p_, "totBandwidth");             
  		
        ScalarColumn<String> windowFunction(*table_p_, "windowFunction");             
  		
  		
        ScalarColumn<double> chanFreqStart(*table_p_, "chanFreqStart");             
  		
        ScalarColumn<double> chanFreqStep(*table_p_, "chanFreqStep");             
  		
        ArrayColumn<double> chanFreqArray(*table_p_, "chanFreqArray");             
  		
        ScalarColumn<double> chanWidth(*table_p_, "chanWidth");             
  		
        ArrayColumn<double> chanWidthArray(*table_p_, "chanWidthArray");             
  		
        ScalarColumn<String> correlationBit(*table_p_, "correlationBit");             
  		
        ScalarColumn<double> effectiveBw(*table_p_, "effectiveBw");             
  		
        ArrayColumn<double> effectiveBwArray(*table_p_, "effectiveBwArray");             
  		
        ScalarColumn<int> freqGroup(*table_p_, "freqGroup");             
  		
        ScalarColumn<String> freqGroupName(*table_p_, "freqGroupName");             
  		
        ArrayColumn<bool> lineArray(*table_p_, "lineArray");             
  		
        ScalarColumn<String> measFreqRef(*table_p_, "measFreqRef");             
  		
        ScalarColumn<String> name(*table_p_, "name");             
  		
        ScalarColumn<bool> oversampling(*table_p_, "oversampling");             
  		
        ScalarColumn<bool> quantization(*table_p_, "quantization");             
  		
        ScalarColumn<double> refChan(*table_p_, "refChan");             
  		
        ScalarColumn<double> resolution(*table_p_, "resolution");             
  		
        ArrayColumn<double> resolutionArray(*table_p_, "resolutionArray");             
  		
        ScalarColumn<int> numAssocValues(*table_p_, "numAssocValues");             
  		
        ArrayColumn<String> assocNature(*table_p_, "assocNature");             
  		
        ArrayColumn<String> assocSpectralWindowId(*table_p_, "assocSpectralWindowId");             
  		
        ScalarColumn<String> imageSpectralWindowId(*table_p_, "imageSpectralWindowId");             
  		
        ScalarColumn<int> dopplerId(*table_p_, "dopplerId");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

		
	
            basebandName.put(rowIndex, CBasebandName::name(rows.at(i)->getBasebandName()));
	

	
            netSideband.put(rowIndex, CNetSideband::name(rows.at(i)->getNetSideband()));
	

	
            numChan.put(rowIndex, rows.at(i)->getNumChan());
	

	
            refFreq.put(rowIndex, rows.at(i)->getRefFreq().get());
	

	
            sidebandProcessingMode.put(rowIndex, CSidebandProcessingMode::name(rows.at(i)->getSidebandProcessingMode()));
	

	
            totBandwidth.put(rowIndex, rows.at(i)->getTotBandwidth().get());
	

	
            windowFunction.put(rowIndex, CWindowFunction::name(rows.at(i)->getWindowFunction()));
	

		
	
            if (rows.at(i)->isChanFreqStartExists())
                chanFreqStart.put(rowIndex, rows.at(i)->getChanFreqStart().get());
	

	
            if (rows.at(i)->isChanFreqStepExists())
                chanFreqStep.put(rowIndex, rows.at(i)->getChanFreqStep().get());
	

	
            if (rows.at(i)->isChanFreqArrayExists())
                chanFreqArray.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getChanFreqArray()));
	

	
            if (rows.at(i)->isChanWidthExists())
                chanWidth.put(rowIndex, rows.at(i)->getChanWidth().get());
	

	
            if (rows.at(i)->isChanWidthArrayExists())
                chanWidthArray.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getChanWidthArray()));
	

	
            if (rows.at(i)->isCorrelationBitExists())
                correlationBit.put(rowIndex, CCorrelationBit::name(rows.at(i)->getCorrelationBit()));
	

	
            if (rows.at(i)->isEffectiveBwExists())
                effectiveBw.put(rowIndex, rows.at(i)->getEffectiveBw().get());
	

	
            if (rows.at(i)->isEffectiveBwArrayExists())
                effectiveBwArray.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getEffectiveBwArray()));
	

	
            if (rows.at(i)->isFreqGroupExists())
                freqGroup.put(rowIndex, rows.at(i)->getFreqGroup());
	

	
            if (rows.at(i)->isFreqGroupNameExists())
                freqGroupName.put(rowIndex, rows.at(i)->getFreqGroupName());
	

	
            if (rows.at(i)->isLineArrayExists())
                lineArray.put(rowIndex, basic2CASA1D<bool,bool>(rows.at(i)->getLineArray()));
	

	
            if (rows.at(i)->isMeasFreqRefExists())
                measFreqRef.put(rowIndex, CFrequencyReferenceCode::name(rows.at(i)->getMeasFreqRef()));
	

	
            if (rows.at(i)->isNameExists())
                name.put(rowIndex, rows.at(i)->getName());
	

	
            if (rows.at(i)->isOversamplingExists())
                oversampling.put(rowIndex, rows.at(i)->getOversampling());
	

	
            if (rows.at(i)->isQuantizationExists())
                quantization.put(rowIndex, rows.at(i)->getQuantization());
	

	
            if (rows.at(i)->isRefChanExists())
                refChan.put(rowIndex, rows.at(i)->getRefChan());
	

	
            if (rows.at(i)->isResolutionExists())
                resolution.put(rowIndex, rows.at(i)->getResolution().get());
	

	
            if (rows.at(i)->isResolutionArrayExists())
                resolutionArray.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getResolutionArray()));
	

	
            if (rows.at(i)->isNumAssocValuesExists())
                numAssocValues.put(rowIndex, rows.at(i)->getNumAssocValues());
	

	
            if (rows.at(i)->isAssocNatureExists())
                assocNature.put(rowIndex, enum2CASA1D<SpectralResolutionType,CSpectralResolutionType>(rows.at(i)->getAssocNature()));
	

	
            if (rows.at(i)->isAssocSpectralWindowIdExists())
                assocSpectralWindowId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getAssocSpectralWindowId()));
	

	
            if (rows.at(i)->isImageSpectralWindowIdExists())
                imageSpectralWindowId.put(rowIndex, rows.at(i)->getImageSpectralWindowId().toString());
	

	
            if (rows.at(i)->isDopplerIdExists())
                dopplerId.put(rowIndex, rows.at(i)->getDopplerId());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SQUARELAWDETECTOR::ASDM_SQUARELAWDETECTOR() {
        name_ = "ASDM_SQUARELAWDETECTOR";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset SquareLawDetector table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("squareLawDetectorId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numBand", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("bandType", "blabla"));
  		
  		  		
    }

    ASDM_SQUARELAWDETECTOR::~ASDM_SQUARELAWDETECTOR() {
    }

    const TableDesc& ASDM_SQUARELAWDETECTOR::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
			
		
	
    void ASDM_SQUARELAWDETECTOR::fill(const ASDM& asdm) {
        vector<SquareLawDetectorRow*> rows = asdm.getSquareLawDetector().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> squareLawDetectorId(*table_p_, "squareLawDetectorId");             
  		
  		
        ScalarColumn<int> numBand(*table_p_, "numBand");             
  		
        ScalarColumn<String> bandType(*table_p_, "bandType");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            squareLawDetectorId.put(rowIndex, rows.at(i)->getSquareLawDetectorId().toString());
	

		
	
            numBand.put(rowIndex, rows.at(i)->getNumBand());
	

	
            bandType.put(rowIndex, CDetectorBandType::name(rows.at(i)->getBandType()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_STATE::ASDM_STATE() {
        name_ = "ASDM_STATE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset State table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("stateId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("calDeviceName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("sig", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("ref", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("onSky", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<float>("weight", "blabla"));
  		  		
    }

    ASDM_STATE::~ASDM_STATE() {
    }

    const TableDesc& ASDM_STATE::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_STATE::fill(const ASDM& asdm) {
        vector<StateRow*> rows = asdm.getState().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> stateId(*table_p_, "stateId");             
  		
  		
        ScalarColumn<String> calDeviceName(*table_p_, "calDeviceName");             
  		
        ScalarColumn<bool> sig(*table_p_, "sig");             
  		
        ScalarColumn<bool> ref(*table_p_, "ref");             
  		
        ScalarColumn<bool> onSky(*table_p_, "onSky");             
  		
  		
        ScalarColumn<float> weight(*table_p_, "weight");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            stateId.put(rowIndex, rows.at(i)->getStateId().toString());
	

		
	
            calDeviceName.put(rowIndex, CCalibrationDevice::name(rows.at(i)->getCalDeviceName()));
	

	
            sig.put(rowIndex, rows.at(i)->getSig());
	

	
            ref.put(rowIndex, rows.at(i)->getRef());
	

	
            onSky.put(rowIndex, rows.at(i)->getOnSky());
	

		
	
            if (rows.at(i)->isWeightExists())
                weight.put(rowIndex, rows.at(i)->getWeight());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_STATION::ASDM_STATION() {
        name_ = "ASDM_STATION";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Station table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("stationId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("name", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("position", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("type", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("time", "blabla"));
  		  		
    }

    ASDM_STATION::~ASDM_STATION() {
    }

    const TableDesc& ASDM_STATION::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
			
		
			
		
	
    void ASDM_STATION::fill(const ASDM& asdm) {
        vector<StationRow*> rows = asdm.getStation().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> stationId(*table_p_, "stationId");             
  		
  		
        ScalarColumn<String> name(*table_p_, "name");             
  		
        ArrayColumn<double> position(*table_p_, "position");             
  		
        ScalarColumn<String> type(*table_p_, "type");             
  		
  		
        ScalarColumn<double> time(*table_p_, "time");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            stationId.put(rowIndex, rows.at(i)->getStationId().toString());
	

		
	
            name.put(rowIndex, rows.at(i)->getName());
	

	
            position.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getPosition()));
	

	
            type.put(rowIndex, CStationType::name(rows.at(i)->getType()));
	

		
	
            if (rows.at(i)->isTimeExists())
                time.put(rowIndex, rows.at(i)->getTime().get()/(1.0e9));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SUBSCAN::ASDM_SUBSCAN() {
        name_ = "ASDM_SUBSCAN";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Subscan table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("scanNumber", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("subscanNumber", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("startTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("endTime", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("fieldName", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("subscanIntent", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numIntegration", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("numSubintegration", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("subscanMode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("correlatorCalibration", "blabla"));
  		  		
    }

    ASDM_SUBSCAN::~ASDM_SUBSCAN() {
    }

    const TableDesc& ASDM_SUBSCAN::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
			
		
			
		
			
		
			
			
		
	
    void ASDM_SUBSCAN::fill(const ASDM& asdm) {
        vector<SubscanRow*> rows = asdm.getSubscan().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> execBlockId(*table_p_, "execBlockId");             
  		
        ScalarColumn<int> scanNumber(*table_p_, "scanNumber");             
  		
        ScalarColumn<int> subscanNumber(*table_p_, "subscanNumber");             
  		
  		
        ScalarColumn<double> startTime(*table_p_, "startTime");             
  		
        ScalarColumn<double> endTime(*table_p_, "endTime");             
  		
        ScalarColumn<String> fieldName(*table_p_, "fieldName");             
  		
        ScalarColumn<String> subscanIntent(*table_p_, "subscanIntent");             
  		
        ScalarColumn<int> numIntegration(*table_p_, "numIntegration");             
  		
        ArrayColumn<int> numSubintegration(*table_p_, "numSubintegration");             
  		
  		
        ScalarColumn<String> subscanMode(*table_p_, "subscanMode");             
  		
        ScalarColumn<String> correlatorCalibration(*table_p_, "correlatorCalibration");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            execBlockId.put(rowIndex, rows.at(i)->getExecBlockId().toString());
	

	
            scanNumber.put(rowIndex, rows.at(i)->getScanNumber());
	

	
            subscanNumber.put(rowIndex, rows.at(i)->getSubscanNumber());
	

		
	
            startTime.put(rowIndex, rows.at(i)->getStartTime().get()/(1.0e9));
	

	
            endTime.put(rowIndex, rows.at(i)->getEndTime().get()/(1.0e9));
	

	
            fieldName.put(rowIndex, rows.at(i)->getFieldName());
	

	
            subscanIntent.put(rowIndex, CSubscanIntent::name(rows.at(i)->getSubscanIntent()));
	

	
            numIntegration.put(rowIndex, rows.at(i)->getNumIntegration());
	

	
            numSubintegration.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getNumSubintegration()));
	

		
	
            if (rows.at(i)->isSubscanModeExists())
                subscanMode.put(rowIndex, CSwitchingMode::name(rows.at(i)->getSubscanMode()));
	

	
            if (rows.at(i)->isCorrelatorCalibrationExists())
                correlatorCalibration.put(rowIndex, CCorrelatorCalibration::name(rows.at(i)->getCorrelatorCalibration()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SWITCHCYCLE::ASDM_SWITCHCYCLE() {
        name_ = "ASDM_SWITCHCYCLE";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset SwitchCycle table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("switchCycleId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numStep", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("weightArray", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("dirOffsetArray", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("freqOffsetArray", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("stepDurationArray", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("directionCode", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("directionEquinox", "blabla"));
  		  		
    }

    ASDM_SWITCHCYCLE::~ASDM_SWITCHCYCLE() {
    }

    const TableDesc& ASDM_SWITCHCYCLE::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
			
		
			
		
	
    void ASDM_SWITCHCYCLE::fill(const ASDM& asdm) {
        vector<SwitchCycleRow*> rows = asdm.getSwitchCycle().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> switchCycleId(*table_p_, "switchCycleId");             
  		
  		
        ScalarColumn<int> numStep(*table_p_, "numStep");             
  		
        ArrayColumn<float> weightArray(*table_p_, "weightArray");             
  		
        ArrayColumn<double> dirOffsetArray(*table_p_, "dirOffsetArray");             
  		
        ArrayColumn<double> freqOffsetArray(*table_p_, "freqOffsetArray");             
  		
        ArrayColumn<double> stepDurationArray(*table_p_, "stepDurationArray");             
  		
  		
        ScalarColumn<String> directionCode(*table_p_, "directionCode");             
  		
        ScalarColumn<double> directionEquinox(*table_p_, "directionEquinox");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            switchCycleId.put(rowIndex, rows.at(i)->getSwitchCycleId().toString());
	

		
	
            numStep.put(rowIndex, rows.at(i)->getNumStep());
	

	
            weightArray.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getWeightArray()));
	

	
            dirOffsetArray.put(rowIndex, ext2CASA2D<Angle,double>(rows.at(i)->getDirOffsetArray()));
	

	
            freqOffsetArray.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getFreqOffsetArray()));
	

	
            stepDurationArray.put(rowIndex, interval2CASA1D<Interval,double>(rows.at(i)->getStepDurationArray()));
	

		
	
            if (rows.at(i)->isDirectionCodeExists())
                directionCode.put(rowIndex, CDirectionReferenceCode::name(rows.at(i)->getDirectionCode()));
	

	
            if (rows.at(i)->isDirectionEquinoxExists())
                directionEquinox.put(rowIndex, rows.at(i)->getDirectionEquinox().get()/(1.0e9));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SYSCAL::ASDM_SYSCAL() {
        name_ = "ASDM_SYSCAL";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset SysCal table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("feedId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numChan", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("tcalFlag", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tcalSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("trxFlag", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("trxSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("tskyFlag", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tskySpectrum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("tsysFlag", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("tsysSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("tantFlag", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("tantSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("tantTsysFlag", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("tantTsysSpectrum", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<bool>("phaseDiffFlag", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("phaseDiffSpectrum", "blabla"));
  		  		
    }

    ASDM_SYSCAL::~ASDM_SYSCAL() {
    }

    const TableDesc& ASDM_SYSCAL::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_SYSCAL::fill(const ASDM& asdm) {
        vector<SysCalRow*> rows = asdm.getSysCal().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
        ScalarColumn<int> feedId(*table_p_, "feedId");             
  		
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
        ScalarColumn<int> numChan(*table_p_, "numChan");             
  		
  		
        ScalarColumn<bool> tcalFlag(*table_p_, "tcalFlag");             
  		
        ArrayColumn<double> tcalSpectrum(*table_p_, "tcalSpectrum");             
  		
        ScalarColumn<bool> trxFlag(*table_p_, "trxFlag");             
  		
        ArrayColumn<double> trxSpectrum(*table_p_, "trxSpectrum");             
  		
        ScalarColumn<bool> tskyFlag(*table_p_, "tskyFlag");             
  		
        ArrayColumn<double> tskySpectrum(*table_p_, "tskySpectrum");             
  		
        ScalarColumn<bool> tsysFlag(*table_p_, "tsysFlag");             
  		
        ArrayColumn<double> tsysSpectrum(*table_p_, "tsysSpectrum");             
  		
        ScalarColumn<bool> tantFlag(*table_p_, "tantFlag");             
  		
        ArrayColumn<float> tantSpectrum(*table_p_, "tantSpectrum");             
  		
        ScalarColumn<bool> tantTsysFlag(*table_p_, "tantTsysFlag");             
  		
        ArrayColumn<float> tantTsysSpectrum(*table_p_, "tantTsysSpectrum");             
  		
        ScalarColumn<bool> phaseDiffFlag(*table_p_, "phaseDiffFlag");             
  		
        ArrayColumn<float> phaseDiffSpectrum(*table_p_, "phaseDiffSpectrum");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

	
            feedId.put(rowIndex, rows.at(i)->getFeedId());
	

		
	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

	
            numChan.put(rowIndex, rows.at(i)->getNumChan());
	

		
	
            if (rows.at(i)->isTcalFlagExists())
                tcalFlag.put(rowIndex, rows.at(i)->getTcalFlag());
	

	
            if (rows.at(i)->isTcalSpectrumExists())
                tcalSpectrum.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getTcalSpectrum()));
	

	
            if (rows.at(i)->isTrxFlagExists())
                trxFlag.put(rowIndex, rows.at(i)->getTrxFlag());
	

	
            if (rows.at(i)->isTrxSpectrumExists())
                trxSpectrum.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getTrxSpectrum()));
	

	
            if (rows.at(i)->isTskyFlagExists())
                tskyFlag.put(rowIndex, rows.at(i)->getTskyFlag());
	

	
            if (rows.at(i)->isTskySpectrumExists())
                tskySpectrum.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getTskySpectrum()));
	

	
            if (rows.at(i)->isTsysFlagExists())
                tsysFlag.put(rowIndex, rows.at(i)->getTsysFlag());
	

	
            if (rows.at(i)->isTsysSpectrumExists())
                tsysSpectrum.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getTsysSpectrum()));
	

	
            if (rows.at(i)->isTantFlagExists())
                tantFlag.put(rowIndex, rows.at(i)->getTantFlag());
	

	
            if (rows.at(i)->isTantSpectrumExists())
                tantSpectrum.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getTantSpectrum()));
	

	
            if (rows.at(i)->isTantTsysFlagExists())
                tantTsysFlag.put(rowIndex, rows.at(i)->getTantTsysFlag());
	

	
            if (rows.at(i)->isTantTsysSpectrumExists())
                tantTsysSpectrum.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getTantTsysSpectrum()));
	

	
            if (rows.at(i)->isPhaseDiffFlagExists())
                phaseDiffFlag.put(rowIndex, rows.at(i)->getPhaseDiffFlag());
	

	
            if (rows.at(i)->isPhaseDiffSpectrumExists())
                phaseDiffSpectrum.put(rowIndex, basic2CASA2D<float,float>(rows.at(i)->getPhaseDiffSpectrum()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_SYSPOWER::ASDM_SYSPOWER() {
        name_ = "ASDM_SYSPOWER";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset SysPower table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("feedId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numReceptor", "blabla"));
  		
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("switchedPowerDifference", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("switchedPowerSum", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("requantizerGain", "blabla"));
  		  		
    }

    ASDM_SYSPOWER::~ASDM_SYSPOWER() {
    }

    const TableDesc& ASDM_SYSPOWER::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_SYSPOWER::fill(const ASDM& asdm) {
        vector<SysPowerRow*> rows = asdm.getSysPower().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ScalarColumn<int> feedId(*table_p_, "feedId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
        ScalarColumn<int> numReceptor(*table_p_, "numReceptor");             
  		
  		
        ArrayColumn<float> switchedPowerDifference(*table_p_, "switchedPowerDifference");             
  		
        ArrayColumn<float> switchedPowerSum(*table_p_, "switchedPowerSum");             
  		
        ArrayColumn<float> requantizerGain(*table_p_, "requantizerGain");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            feedId.put(rowIndex, rows.at(i)->getFeedId());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
	
            numReceptor.put(rowIndex, rows.at(i)->getNumReceptor());
	

		
	
            if (rows.at(i)->isSwitchedPowerDifferenceExists())
                switchedPowerDifference.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getSwitchedPowerDifference()));
	

	
            if (rows.at(i)->isSwitchedPowerSumExists())
                switchedPowerSum.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getSwitchedPowerSum()));
	

	
            if (rows.at(i)->isRequantizerGainExists())
                requantizerGain.put(rowIndex, basic2CASA1D<float,float>(rows.at(i)->getRequantizerGain()));
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_TOTALPOWER::ASDM_TOTALPOWER() {
        name_ = "ASDM_TOTALPOWER";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset TotalPower table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("time", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("configDescriptionId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("fieldId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("scanNumber", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("subscanNumber", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("integrationNumber", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("uvw", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("exposure", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeCentroid", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("floatData", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("flagAnt", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<int>("flagPol", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("interval", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("stateId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("execBlockId", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("subintegrationNumber", "blabla"));
  		  		
    }

    ASDM_TOTALPOWER::~ASDM_TOTALPOWER() {
    }

    const TableDesc& ASDM_TOTALPOWER::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_TOTALPOWER::fill(const ASDM& asdm) {
        vector<TotalPowerRow*> rows = asdm.getTotalPower().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<double> time(*table_p_, "time");             
  		
        ScalarColumn<String> configDescriptionId(*table_p_, "configDescriptionId");             
  		
        ScalarColumn<String> fieldId(*table_p_, "fieldId");             
  		
  		
        ScalarColumn<int> scanNumber(*table_p_, "scanNumber");             
  		
        ScalarColumn<int> subscanNumber(*table_p_, "subscanNumber");             
  		
        ScalarColumn<int> integrationNumber(*table_p_, "integrationNumber");             
  		
        ArrayColumn<double> uvw(*table_p_, "uvw");             
  		
        ArrayColumn<double> exposure(*table_p_, "exposure");             
  		
        ArrayColumn<double> timeCentroid(*table_p_, "timeCentroid");             
  		
        ArrayColumn<float> floatData(*table_p_, "floatData");             
  		
        ArrayColumn<int> flagAnt(*table_p_, "flagAnt");             
  		
        ArrayColumn<int> flagPol(*table_p_, "flagPol");             
  		
        ScalarColumn<double> interval(*table_p_, "interval");             
  		
        ArrayColumn<String> stateId(*table_p_, "stateId");             
  		
        ScalarColumn<String> execBlockId(*table_p_, "execBlockId");             
  		
  		
        ScalarColumn<int> subintegrationNumber(*table_p_, "subintegrationNumber");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            time.put(rowIndex, rows.at(i)->getTime().get()/(1.0e9));
	

	
            configDescriptionId.put(rowIndex, rows.at(i)->getConfigDescriptionId().toString());
	

	
            fieldId.put(rowIndex, rows.at(i)->getFieldId().toString());
	

		
	
            scanNumber.put(rowIndex, rows.at(i)->getScanNumber());
	

	
            subscanNumber.put(rowIndex, rows.at(i)->getSubscanNumber());
	

	
            integrationNumber.put(rowIndex, rows.at(i)->getIntegrationNumber());
	

	
            uvw.put(rowIndex, ext2CASA2D<Length,double>(rows.at(i)->getUvw()));
	

	
            exposure.put(rowIndex, interval2CASA2D<Interval,double>(rows.at(i)->getExposure()));
	

	
            timeCentroid.put(rowIndex, at2CASA2D<double>(rows.at(i)->getTimeCentroid()));
	

	
            floatData.put(rowIndex, basic2CASA3D<float,float>(rows.at(i)->getFloatData()));
	

	
            flagAnt.put(rowIndex, basic2CASA1D<int,int>(rows.at(i)->getFlagAnt()));
	

	
            flagPol.put(rowIndex, basic2CASA2D<int,int>(rows.at(i)->getFlagPol()));
	

	
            interval.put(rowIndex, rows.at(i)->getInterval().get()/(1.0e9));
	

	
            stateId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getStateId()));
	

	
            execBlockId.put(rowIndex, rows.at(i)->getExecBlockId().toString());
	

		
	
            if (rows.at(i)->isSubintegrationNumberExists())
                subintegrationNumber.put(rowIndex, rows.at(i)->getSubintegrationNumber());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_WVMCAL::ASDM_WVMCAL() {
        name_ = "ASDM_WVMCAL";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset WVMCal table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("antennaId", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("spectralWindowId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("wvrMethod", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("polyFreqLimits", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numInputAntenna", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numChan", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numPoly", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<float>("pathCoeff", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("refTemp", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<String>("inputAntennaId", "blabla"));
  		
  		  		
    }

    ASDM_WVMCAL::~ASDM_WVMCAL() {
    }

    const TableDesc& ASDM_WVMCAL::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_WVMCAL::fill(const ASDM& asdm) {
        vector<WVMCalRow*> rows = asdm.getWVMCal().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> antennaId(*table_p_, "antennaId");             
  		
        ScalarColumn<String> spectralWindowId(*table_p_, "spectralWindowId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
        ScalarColumn<String> wvrMethod(*table_p_, "wvrMethod");             
  		
        ArrayColumn<double> polyFreqLimits(*table_p_, "polyFreqLimits");             
  		
        ScalarColumn<int> numInputAntenna(*table_p_, "numInputAntenna");             
  		
        ScalarColumn<int> numChan(*table_p_, "numChan");             
  		
        ScalarColumn<int> numPoly(*table_p_, "numPoly");             
  		
        ArrayColumn<float> pathCoeff(*table_p_, "pathCoeff");             
  		
        ArrayColumn<double> refTemp(*table_p_, "refTemp");             
  		
        ArrayColumn<String> inputAntennaId(*table_p_, "inputAntennaId");             
  		
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            antennaId.put(rowIndex, rows.at(i)->getAntennaId().toString());
	

	
            spectralWindowId.put(rowIndex, rows.at(i)->getSpectralWindowId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
	
            wvrMethod.put(rowIndex, CWVRMethod::name(rows.at(i)->getWvrMethod()));
	

	
            polyFreqLimits.put(rowIndex, ext2CASA1D<Frequency,double>(rows.at(i)->getPolyFreqLimits()));
	

	
            numInputAntenna.put(rowIndex, rows.at(i)->getNumInputAntenna());
	

	
            numChan.put(rowIndex, rows.at(i)->getNumChan());
	

	
            numPoly.put(rowIndex, rows.at(i)->getNumPoly());
	

	
            pathCoeff.put(rowIndex, basic2CASA3D<float,float>(rows.at(i)->getPathCoeff()));
	

	
            refTemp.put(rowIndex, ext2CASA2D<Temperature,double>(rows.at(i)->getRefTemp()));
	

	
            inputAntennaId.put(rowIndex, _2CASAString1D<Tag,String>(rows.at(i)->getInputAntennaId()));
	

		
            rowIndex++;		
        }
        table_p_->flush();
    }
	 
    ASDM_WEATHER::ASDM_WEATHER() {
        name_ = "ASDM_WEATHER";
        tableDesc_.comment() = "The verbatim copy of the ASDM's dataset Weather table";
  		
        tableDesc_.addColumn(ScalarColumnDesc<String>("stationId", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("timeInterval", "blabla"));
  		
  		
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("pressure", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("relHumidity", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("temperature", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("windDirection", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("windSpeed", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("windMax", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("dewPoint", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numLayer", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("layerHeight", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("temperatureProfile", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("cloudMonitor", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<int>("numWVR", "blabla"));
  		
        tableDesc_.addColumn(ArrayColumnDesc<double>("wvrTemp", "blabla"));
  		
        tableDesc_.addColumn(ScalarColumnDesc<double>("water", "blabla"));
  		  		
    }

    ASDM_WEATHER::~ASDM_WEATHER() {
    }

    const TableDesc& ASDM_WEATHER::tableDesc() const {
        return tableDesc_;
    }


		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
	
    void ASDM_WEATHER::fill(const ASDM& asdm) {
        vector<WeatherRow*> rows = asdm.getWeather().get();
        unsigned int rowIndex = table_p_->nrow();
        table_p_->addRow(rows.size());
  		
        ScalarColumn<String> stationId(*table_p_, "stationId");             
  		
        ArrayColumn<double> timeInterval(*table_p_, "timeInterval");             
  		
  		
  		
        ScalarColumn<double> pressure(*table_p_, "pressure");             
  		
        ScalarColumn<double> relHumidity(*table_p_, "relHumidity");             
  		
        ScalarColumn<double> temperature(*table_p_, "temperature");             
  		
        ScalarColumn<double> windDirection(*table_p_, "windDirection");             
  		
        ScalarColumn<double> windSpeed(*table_p_, "windSpeed");             
  		
        ScalarColumn<double> windMax(*table_p_, "windMax");             
  		
        ScalarColumn<double> dewPoint(*table_p_, "dewPoint");             
  		
        ScalarColumn<int> numLayer(*table_p_, "numLayer");             
  		
        ArrayColumn<double> layerHeight(*table_p_, "layerHeight");             
  		
        ArrayColumn<double> temperatureProfile(*table_p_, "temperatureProfile");             
  		
        ScalarColumn<double> cloudMonitor(*table_p_, "cloudMonitor");             
  		
        ScalarColumn<int> numWVR(*table_p_, "numWVR");             
  		
        ArrayColumn<double> wvrTemp(*table_p_, "wvrTemp");             
  		
        ScalarColumn<double> water(*table_p_, "water");             
  		  	

        for (unsigned int i = 0; i < rows.size(); i++) {
		
	
            stationId.put(rowIndex, rows.at(i)->getStationId().toString());
	

	
            timeInterval.put(rowIndex, ati2CASA1D<double>(rows.at(i)->getTimeInterval()));
	

		
		
	
            if (rows.at(i)->isPressureExists())
                pressure.put(rowIndex, rows.at(i)->getPressure().get());
	

	
            if (rows.at(i)->isRelHumidityExists())
                relHumidity.put(rowIndex, rows.at(i)->getRelHumidity().get());
	

	
            if (rows.at(i)->isTemperatureExists())
                temperature.put(rowIndex, rows.at(i)->getTemperature().get());
	

	
            if (rows.at(i)->isWindDirectionExists())
                windDirection.put(rowIndex, rows.at(i)->getWindDirection().get());
	

	
            if (rows.at(i)->isWindSpeedExists())
                windSpeed.put(rowIndex, rows.at(i)->getWindSpeed().get());
	

	
            if (rows.at(i)->isWindMaxExists())
                windMax.put(rowIndex, rows.at(i)->getWindMax().get());
	

	
            if (rows.at(i)->isDewPointExists())
                dewPoint.put(rowIndex, rows.at(i)->getDewPoint().get());
	

	
            if (rows.at(i)->isNumLayerExists())
                numLayer.put(rowIndex, rows.at(i)->getNumLayer());
	

	
            if (rows.at(i)->isLayerHeightExists())
                layerHeight.put(rowIndex, ext2CASA1D<Length,double>(rows.at(i)->getLayerHeight()));
	

	
            if (rows.at(i)->isTemperatureProfileExists())
                temperatureProfile.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getTemperatureProfile()));
	

	
            if (rows.at(i)->isCloudMonitorExists())
                cloudMonitor.put(rowIndex, rows.at(i)->getCloudMonitor().get());
	

	
            if (rows.at(i)->isNumWVRExists())
                numWVR.put(rowIndex, rows.at(i)->getNumWVR());
	

	
            if (rows.at(i)->isWvrTempExists())
                wvrTemp.put(rowIndex, ext2CASA1D<Temperature,double>(rows.at(i)->getWvrTemp()));
	

	
            if (rows.at(i)->isWaterExists())
                water.put(rowIndex, rows.at(i)->getWater());
	

            rowIndex++;		
        }
        table_p_->flush();
    }
}

