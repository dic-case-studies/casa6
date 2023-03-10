#
# CASA - Common Astronomy Software Applications
# Copyright (C) 2010 by ESO (in the framework of the ALMA collaboration)
#
# This file is part of CASA.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

casa_add_library( msvis
  MSVis/AsynchronousInterface.cc
  MSVis/AsynchronousTools.cc
  MSVis/AveragingTvi2.cc
  MSVis/AveragingVi2Factory.cc
  MSVis/CalVisBuffer.cc
  MSVis/CalibratingVi2FactoryI.cc
  MSVis/FinalTvi2.cc
  MSVis/GroupProcessor.cc
  MSVis/GroupWorker.cc
  MSVis/IteratingParameters.cc
  MSVis/LayeredVi2Factory.cc
  MSVis/MSCalEnums.cc
  MSVis/MSChecker.cc
  MSVis/MSContinuumSubtractor.cc
  MSVis/MSIter2.cc
  MSVis/MSUVWGenerator.cc
  MSVis/MSUtil.cc
  MSVis/MeasurementSet2.cc
  MSVis/MsAverager.cc
  MSVis/MsRows.cc
  MSVis/Partition.cc
  MSVis/PassThroughTVI.cc
  MSVis/PhaseShiftTvi.cc
  MSVis/PointingDirectionCache.cc
  MSVis/Reweighter.cc
  MSVis/SelectAverageSpw.cc
  MSVis/SimpleSubMS.cc
  MSVis/SimpleSimVi2.cc
  MSVis/SpectralWindow.cc
  MSVis/StatWT.cc
  MSVis/StokesVector.cc
  MSVis/SubMS.cc
  MSVis/TransformingVi2.cc
  MSVis/UtilJ.cc
  MSVis/VBContinuumSubtractor.cc
  MSVis/VBGContinuumSubtractor.cc
  MSVis/VBRemapper.cc
  MSVis/VLAT.cc
  MSVis/VWBT.cc
  MSVis/Vbi2MsRow.cc
  MSVis/ViColumns2.cc
  MSVis/ViFrequencySelection.cc
  MSVis/ViImplementation2.cc
  MSVis/ViiLayerFactory.cc
  MSVis/VisBuffAccumulator.cc
  MSVis/VisBuffGroup.cc
  MSVis/VisBuffGroupAcc.cc
  MSVis/VisBuffer.cc
  MSVis/VisBuffer2.cc
  MSVis/VisBufferAsync.cc
  MSVis/VisBufferAsyncWrapper.cc
  MSVis/VisBufferComponents2.cc
  MSVis/VisBufferImpl2.cc
  MSVis/VisBufferUtil.cc
  MSVis/VisChunkAverager.cc
  MSVis/VisImagingWeight.cc
  MSVis/VisIterator.cc
  MSVis/VisModelDataI.cc
  MSVis/VisSet.cc
  MSVis/VisSetUtil.cc
  MSVis/VisibilityIterator.cc
  MSVis/VisibilityIterator2.cc
  MSVis/VisibilityIteratorImpl.cc
  MSVis/VisibilityIteratorImpl2.cc
  MSVis/VisibilityIteratorImplAsync.cc
  MSVis/VisibilityProcessing.cc 
####   MSVis/AsynchronousInterface2.cc
####  MSVis/VisBufferAsyncWrapper2.cc
####  MSVis/VisBufferImplAsync2.cc
####  MSVis/VisibilityIteratorImplAsync2.cc
  MSVis/statistics/Vi2AntennaDataProvider.cc
  MSVis/statistics/Vi2ArrayIdDataProvider.cc
  MSVis/statistics/Vi2DataDescriptionIdsDataProvider.cc
  MSVis/statistics/Vi2FeedDataProvider.cc
  MSVis/statistics/Vi2FieldIdDataProvider.cc
  MSVis/statistics/Vi2FlagCubeDataProvider.cc
  MSVis/statistics/Vi2FlagRowDataProvider.cc
  MSVis/statistics/Vi2FloatVisDataProvider.cc
  MSVis/statistics/Vi2IntervalDataProvider.cc
  MSVis/statistics/Vi2ScanDataProvider.cc
  MSVis/statistics/Vi2TimeDataProvider.cc
  MSVis/statistics/Vi2UVRangeDataProvider.cc
  MSVis/statistics/Vi2WeightSpectrumDataProvider.cc
  MSVis/statistics/Vi2StatsAmplitudeIterator.cc
  MSVis/statistics/Vi2StatsBoolIterator.cc
  MSVis/statistics/Vi2StatsFlagsIterator.cc
  MSVis/statistics/Vi2StatsImaginaryIterator.cc
  MSVis/statistics/Vi2StatsPhaseIterator.cc
  MSVis/statistics/Vi2StatsRealIterator.cc
  MSVis/statistics/Vi2StatsSigmasIterator.cc
  MSVis/statistics/Vi2StatsUVRangeIterator.cc
  MSVis/statistics/Vi2StatsWeightsIterator.cc)

install (FILES
	MSVis/AsynchronousInterface.h
	MSVis/AsynchronousInterface2.h
	MSVis/AsynchronousTools.h
	MSVis/AveragingTvi2.h
	MSVis/AveragingVi2Factory.h
	MSVis/CalVisBuffer.h
	MSVis/CalibratingVi2FactoryI.h
	MSVis/FinalTvi2.h
	MSVis/GroupProcessor.h
	MSVis/GroupWorker.h
	MSVis/IteratingParameters.h
	MSVis/LayeredVi2Factory.h
	MSVis/MSCalEnums.h
	MSVis/MSChecker.h
	MSVis/MSContinuumSubtractor.h
	MSVis/MSIter2.h
	MSVis/MSUVWGenerator.h
	MSVis/MSUtil.h
	MSVis/MsAverager.h
	MSVis/Partition.h
	MSVis/PassThroughTVI.h
	MSVis/PhaseShiftTvi.h
	MSVis/Reweighter.h
	MSVis/PointingDirectionCache.h
	MSVis/SelectAverageSpw.h
	MSVis/SimpleSimVi2.h
	MSVis/SimpleSubMS.h
	MSVis/StatWT.h
	MSVis/StokesVector.h
	MSVis/SubMS.h
	MSVis/UtilJ.h
	MSVis/VBContinuumSubtractor.h
	MSVis/VBGContinuumSubtractor.h
	MSVis/VBRemapper.h
	MSVis/VLAT.h
	MSVis/VWBT.h
	MSVis/ViiLayerFactory.h
	MSVis/VisBuffAccumulator.h
	MSVis/VisBuffGroup.h
	MSVis/VisBuffGroupAcc.h
	MSVis/VisBuffer.h
	MSVis/VisBuffer.tcc
	MSVis/VisBufferAsync.h
	MSVis/VisBufferAsyncWrapper.h
	MSVis/VisBufferComponents.h
	MSVis/VisBufferUtil.h
	MSVis/VisChunkAverager.h
	MSVis/VisImagingWeight.h
	MSVis/VisIterator.h
	MSVis/VisSet.h
	MSVis/VisSetUtil.h
	MSVis/VisibilityIterator.h
	MSVis/VisibilityIterator.tcc
	MSVis/VisibilityIteratorImpl.h
	MSVis/VisibilityIteratorImplAsync.h
	MSVis/VisibilityProcessing.h
        MSVis/AsynchronousInterface2.h
        MSVis/MeasurementSet2.h
        MSVis/SpectralWindow.h
        MSVis/TransformingVi2.h
        MSVis/ViColumns2.h
        MSVis/ViFrequencySelection.h
        MSVis/ViImplementation2.h
        MSVis/VisBuffer2.h
        MSVis/VisBuffer2Adapter.h
        MSVis/VisBufferAsyncWrapper2.h
        MSVis/VisBufferComponents2.h
        MSVis/VisBufferImpl2.h
        MSVis/VisBufferImplAsync2.h
        MSVis/VisModelDataI.h
        MSVis/VisibilityIterator2.h
        MSVis/VisibilityIteratorImpl2.h
        MSVis/VisibilityIteratorImplAsync2.h
	DESTINATION include/casacode/msvis/MSVis
	)

install (FILES
	MSVis/test/TestUtilsTVI.h 
	DESTINATION include/casacode/msvis/MSVis/test
	)       

install (FILES
	MSVis/statistics/Vi2AntennaDataProvider.h
	MSVis/statistics/Vi2ArrayIdDataProvider.h
	MSVis/statistics/Vi2CorrectedVisDataProvider.h
	MSVis/statistics/Vi2DataDescriptionIdsDataProvider.h
	MSVis/statistics/Vi2DataProvider.h
	MSVis/statistics/Vi2FeedDataProvider.h
	MSVis/statistics/Vi2FieldIdDataProvider.h
	MSVis/statistics/Vi2FlagCubeDataProvider.h
	MSVis/statistics/Vi2FlagRowDataProvider.h
	MSVis/statistics/Vi2FloatVisDataProvider.h
	MSVis/statistics/Vi2IntervalDataProvider.h
	MSVis/statistics/Vi2ModelVisDataProvider.h
	MSVis/statistics/Vi2ObservedVisDataProvider.h
	MSVis/statistics/Vi2ScanDataProvider.h
        MSVis/statistics/Vi2StatisticsIteratee.h
	MSVis/statistics/Vi2TimeDataProvider.h
	MSVis/statistics/Vi2UVRangeDataProvider.h
	MSVis/statistics/Vi2VisAmplitudeProvider.h
	MSVis/statistics/Vi2VisImaginaryProvider.h
	MSVis/statistics/Vi2VisPhaseProvider.h
	MSVis/statistics/Vi2VisRealProvider.h
	MSVis/statistics/Vi2WeightSpectrumDataProvider.h
	MSVis/statistics/Vi2StatsAmplitudeIterator.h
	MSVis/statistics/Vi2StatsBoolIterator.h
	MSVis/statistics/Vi2StatsDataIterator.h
	MSVis/statistics/Vi2StatsFlagsIterator.h
	MSVis/statistics/Vi2StatsImaginaryIterator.h
	MSVis/statistics/Vi2StatsPhaseIterator.h
	MSVis/statistics/Vi2StatsRealIterator.h
	MSVis/statistics/Vi2StatsSigmasIterator.h
	MSVis/statistics/Vi2StatsUVRangeIterator.h
	MSVis/statistics/Vi2StatsWeightsIterator.h
	DESTINATION include/casacode/msvis/MSVis/statistics
	)

install (FILES
	MSVis/util/Future.h
	MSVis/util/Try.h
	DESTINATION include/casacode/msvis/MSVis/util
	)

casa_add_assay( msvis MSVis/test/tVisibilityIterator.cc )
casa_add_assay( msvis MSVis/test/tVisibilityIteratorAsync.cc )
casa_add_unit_test ( MODULES msvis SOURCES MSVis/test/tHanningSmooth.cc )
casa_add_unit_test ( MODULES msvis SOURCES MSVis/test/tMSCalEnums.cc )
casa_add_unit_test ( MODULES msvis SOURCES MSVis/test/tUVSub.cc )
casa_add_unit_test ( MODULES msvis SOURCES MSVis/test/tVisIter.cc )
casa_add_unit_test ( MODULES msvis SOURCES MSVis/test/VisibilityIterator_Test.cc MsFactory.cc )

# Following line commented out until a problem related to the older cmake used on Jenkins (dev machine using 2.8.12.2)
# Uncomment this when that problem is resolved
#

casa_add_google_test (MODULES msvis SOURCES MSVis/test/AveragingTvi2_Test.cc MSVis/test/MsFactory.cc MSVis/test/TestUtilsTVI.cc)

casa_add_google_test (MODULES msvis SOURCES MSVis/test/PointingDirectionCache_Test.cc) 
casa_add_google_test (MODULES msvis SOURCES MSVis/test/VisibilityIterator_Gtest.cc MSVis/test/MsFactory.cc) 
casa_add_google_test (NOT_ON_APPLE MODULES msvis SOURCES MSVis/test/MSChecker_Gtest.cc MSVis/test/MSCheckerTest.cc) 

casa_add_google_test (MODULES msvis SOURCES MSVis/test/tSimpleSimVi2_GT.cc) 
casa_add_google_test (MODULES msvis SOURCES MSVis/test/tViiLayerFactory_GT.cc MSVis/test/MsFactory.cc MSVis/test/TestUtilsTVI.cc) 

casa_add_google_test(MODULES msvis SOURCES MSVis/test/Future_Test.cc)
casa_add_google_test(MODULES msvis SOURCES MSVis/test/Try_Test.cc)

#casa_add_unit_test (msvis MSVis/test/VisibilityIterator_Test.cc MSVis/test/MsFactory.cc) 
