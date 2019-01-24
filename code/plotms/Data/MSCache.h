//# MSCache.h: casacore::MS-specific casacore::Data cache for plotms.
//# Copyright (C) 2009
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: $
#ifndef MSCACHE_H_
#define MSCACHE_H_

#include <plotms/Data/PlotMSCacheBase.h>
#include <plotms/Data/PageHeaderCache.h>


#include <plotms/PlotMS/PlotMSAveraging.h>
#include <plotms/PlotMS/PlotMSConstants.h>
#include <plotms/PlotMS/PlotMSFlagging.h>
//#include <plotms/Threads/PlotMSCacheThread.qo.h>
#include <plotms/Threads/ThreadCommunication.h>
#include <plotms/Data/PlotMSVBAverager.h>
#include <plotms/Data/MSCacheVolMeter.h>

#include <casa/aips.h>
#include <casa/Arrays.h>
#include <casa/Containers/Block.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisBufferUtil.h>
#include <msvis/MSVis/ViFrequencySelection.h>

namespace casa {

//# Forward declarations.
class PlotMSApp;
class PlotMSIndexer;

class MSCache : public PlotMSCacheBase {
    
  // Friend class declarations.
  friend class PlotMSIndexer;

public:    
  
  // Constructor which takes parent PlotMS.
  MSCache(PlotMSApp* parent);
  
  // Destructor
  virtual ~MSCache();

  // Identify myself
  PlotMSCacheBase::Type cacheType() const { return PlotMSCacheBase::MS; };

  // ...not yet casacore::MS-specific... (or ever?)
  // Set up indexing for the plot
  //  void setUpIndexer(PMS::Axis iteraxis=PMS::SCAN,
  //		    casacore::Bool globalXRange=false, casacore::Bool globalYRange=false);

  void setFilename(casacore::String filename) { filename_ = filename; };
  virtual casacore::String polname(casacore::Int ipol);
  const PageHeaderCache& pageHeaderCache() const { return pageHeaderCache_; };

protected:

  // casacore::MS-specific loadIt method
  virtual void loadIt(std::vector<PMS::Axis>& loadAxes,
		      std::vector<PMS::DataColumn>& loadData,
		      /*PlotMSCacheThread**/ThreadCommunication* thread = NULL);

  //Returns whether or not the ephemeris data has been
  //attached to a field - radial velocity and rho.
  virtual bool isEphemeris();
private:
    
  // Forbid copy for now
  MSCache(const MSCache&);

  // Set up:
  // DataColumn
  casacore::String getDataColumn(std::vector<PMS::Axis>& loadAxes, 
                       std::vector<PMS::DataColumn>& loadData);
  PMS::DataColumn checkReqDataColumn(PMS::DataColumn reqDataCol);
  void adjustCurrentAxes(PMS::Axis axis, PMS::DataColumn olddata, 
    PMS::DataColumn newdata);
  casacore::String checkLoadedAxesDatacol();
  casacore::String normalizeColumnName(casacore::String plotmscol);
  // MS String names
  void getNamesFromMS(casacore::MeasurementSet& ms);
  // VisIter
  void setUpVisIter(PlotMSSelection& selection,
		    PlotMSCalibration& calibration,
		    casacore::String dataColumn,
            std::vector<PMS::Axis>& loadAxes,
		    std::vector<PMS::DataColumn>& loadData, 
            casacore::Bool interactive=false,
		    casacore::Bool estimateMemory=false,
            ThreadCommunication* thread=NULL);
  vi::VisibilityIterator2* setUpVisIter(casacore::MeasurementSet& selectedMS,
	casacore::Vector<casacore::Vector<casacore::Slice> > chansel, 
    casacore::Vector<casacore::Vector<casacore::Slice> > corrsel);
  void setUpFrequencySelectionChannels(vi::FrequencySelectionUsingChannels fs,
	casacore::Vector<casacore::Vector<casacore::Slice> > chansel);

  // clean up
  void deleteVi();
  void deleteVm();
  void loadError(casacore::String mesg);

  // Estimate cache size for averaging
  bool countChunks(vi::VisibilityIterator2& vi, 
    casacore::Vector<casacore::Int>& nIterPerAve, 
    std::vector<PMS::Axis>& loadAxes,
	std::vector<PMS::DataColumn>& loadData, 
    /*PlotMSCacheThread**/ThreadCommunication* thread);
  void updateEstimateProgress(ThreadCommunication* thread);

  // Trap attempt to use to much memory (too many points)
  void trapExcessVolume(map<PMS::Axis,casacore::Bool> pendingLoadAxes);
  std::vector<casacore::IPosition> visBufferShapes_;

  // Loop over VisIter, filling the cache
  void loadChunks(vi::VisibilityIterator2& vi,
		  const std::vector<PMS::Axis> loadAxes,
		  const std::vector<PMS::DataColumn> loadData,
		  /*PlotMSCacheThread**/ThreadCommunication* thread);
  void loadChunks(vi::VisibilityIterator2& vi,
		  const PlotMSAveraging& averaging,
		  const casacore::Vector<casacore::Int>& nIterPerAve,
		  const std::vector<PMS::Axis> loadAxes,
		  const std::vector<PMS::DataColumn> loadData,
		  /*PlotMSCacheThread**/ThreadCommunication* thread);
  void updateProgress(ThreadCommunication* thread, casacore::Int chunk);

  // Force read on vb for requested axes 
  //   (so pre-cache averaging treats all data it should)
  void forceVBread(vi::VisBuffer2* vb,
		   std::vector<PMS::Axis> loadAxes,
		   std::vector<PMS::DataColumn> loadData);

  // Tell time averager which data column to read
  void discernData(std::vector<PMS::Axis> loadAxes,
		   std::vector<PMS::DataColumn> loadData,
		   PlotMSVBAverager& vba);

  // Loads the specific axis/metadata into the cache using the given VisBuffer.
  void loadAxis(vi::VisBuffer2* vb, casacore::Int vbnum, PMS::Axis axis,
		PMS::DataColumn data = PMS::DEFAULT_DATACOLUMN);

  // Set flags in the MS
  virtual void flagToDisk(const PlotMSFlagging& flagging,
			  casacore::Vector<casacore::Int>& chunks, 
			  casacore::Vector<casacore::Int>& relids,
			  casacore::Bool setFlag,
			  PlotMSIndexer* indexer, int dataIndex);

  casacore::Vector<casacore::Double> calcVelocity(vi::VisBuffer2* vb);

  // For averaging done by PlotMSVBAverager;
  // Some axes need to come from VB2 attached to VI2
  bool useAveragedVisBuffer(PMS::Axis axis);

  // Page Header Cache loading
  void loadPageHeaderCache(const casacore::MeasurementSet& selectedMS);
  void completeLoadPageHeaderCache();

  // Datacolumn to use (requested or "adjusted")
  casacore::String dataColumn_;

  // Create map of intent names to "intent ids" 
  // since state ids can duplicate intents.
  // Then use map to assign intent ids to replace state ids
  // (stateId -> intent string -> intentId)
  map<casacore::String, casacore::Int> intentIds_; 
  void mapIntentNamesToIds();   // create map
  // Use map to assign intent ids
  casacore::Vector<casacore::Int> assignIntentIds(casacore::Vector<casacore::Int>& stateIds);

  // Provisional flagging helpers
  casacore::Vector<casacore::Int> nVBPerAve_;

  // VisIterator pointer
  vi::VisibilityIterator2* vi_p;

  // Volume meter for volume calculation
  MSCacheVolMeter* vm_;

  map<casacore::Int, casacore::Int> chansPerSpw_; 

  bool ephemerisAvailable;

};
typedef casacore::CountedPtr<MSCache> MSCachePtr;


}

#endif /* MSCACHE_H_ */
