//# TablePlot.h: Implement class for the tableplot DO.
//# Copyright (C) 1994,1995,1996,1997,1998,1999,2000,2001,2002,2003-2008
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
//# You should have receied a copy of the GNU Library General Public License
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
//# $Id$
//#
//# TablePlot Change Log
//# =======================
//# casacore::Data   Name      Description
//# Nov. 2006   Shannon Jaeger  TPPlotter was given a parameter to allow
//#                             flagged/non-flagged data to be displayed at the
//#             same time.
//# Dec 22 2006 Urvashi R.V.    Added/updated code documentation.
//#                             Removed two completely unused function declarations.
//#                            
//# Feb 11 2007 Urvashi R.V.    Added the locateData function. Changed code to
//#                             adapt to refactored TPPlotter.
//#
//# May 14 2007 Urvashi R.V.    Refactored code. Added flag version support,
//#                             additional plot options, and parameter checks.
//# Aug 24 2007 Urvashi R.V.    Refactored code. Made TablePlot a singleton class
//#                             and added call-back hooks for the application layers.
//#
//# Oct 29 2007 Shannon J.      Added new message handling.

#ifndef TABLEPLOT_H
#define TABLEPLOT_H

//# Includes
#include <tools/tables/TablePlot/TPPlotter.h> /* TPPlotter.h includes Python.h which must
                                              be included before any system headers */
#include <casa/aips.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/TableError.h>

#include <tools/tables/TablePlot/BasePlot.h>
#include <tools/tables/TablePlot/CrossPlot.h>
#include <tools/tables/TablePlot/PanelParams.h>

#include <flagging/Flagging/SLog.h>


namespace casa { //# NAMESPACE CASA - BEGIN

// <summary>
// Class that implements X-Y plotting of data from casacore::Table objects.
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="" tests="">
// </reviewed>

// <prerequisite>
//#! Classes or concepts you should understand before using this class.
//   <li> BasePlot
//   <li> TPPlotter
//   <li> PanelParams
//   <li> PlotOptions
//   <li> TPCallBackHooks
// </prerequisite>

// <etymology>
// The TablePlot class allows plotting of data from any kind of
// table.
// </etymology>

// <synopsis>
// Class TablePlot allows the plotting of data from any kind of table 
// using generic TaQL expressions to specify the data to be
// plotted. Access to the tables (reading/writing) is controlled via
// the BasePlot class (one basePlot object per table), and all
// plotting operations (plotting/editing) are controlled via the
// TPPlotter class. The TablePlot class controls the use of the
// BasePlot and TPPlotter classes, and provides an interface that
// application programs can directly use. It manages multiple panels
// and plot overlays.
//
// Relation between an application, TablePlot, BasePlot and TPPlotter : 
//
// The application code constructs and holds the TablePlot object, the TPPlotter object,
// and a list of BasePlot objects.
//
// <li> One TPPlotter object is used for one plot window. 
//   Single and multiple panels are supported.
// 
// <li> One BasePlot object attaches to and operates on one single table/subtable.
//   If the TaQL expressions for the attached table, result in TpDouble Arrays.
//   and correspond to multiple data sets to plot, an overlay plot is created.
//   
// <li> TablePlot holds a list of BasePlot objects (one for each table being simultaneously
//   accessed) and a TPPlotter object.
//   and controls the mapping of the BasePlot objects to the supplied TPPlotter
//   objects. After the data is read from tables into all BasePlots, the list is passed
//   to the TPPlotter class, which reads the data to be plotted from all the BasePlots
//   and (currently) creates an overlay plot of data from all the tables in the list.
//   In the case of multi-panel plots, TablePlot accepts a list of BaseBlot lists,
//   and controls the mapping of each BasePlot list to a separate panel.
// 
//  
// This design allows the following.
//
// <li> Access to casacore::Table data and the Plotting package are independant of each other.
// <li> TaQL expressions can be applied to multiple tables and data from multiple tables
//      are kept independant of each other.
// <li> Editing operations on multiple panels of a plot window are handled independantly.
// <li> Design supports the idea of data parallelization and chunking of data to avoid
//      the simultaneous use of large amounts of memory.
// <li> Since the application level holds the BasePlot object lists and TPPlotter objects,
//      control of mapping BasePlot lists to a TPPlotter (list) is at the application
//      level.
// <li> 
// </synopsis>

// <example>
//
// (See tableplot_cmpt.cc)
//
// <srcblock>
// 
//        TablePlot *TP;
//        TP = casa::TablePlot::TablePlotInstance();
//        TP->setTableT(SelTABS,TabNames,SelStr);
//        PlotOptions pop; casacore::Vector<casacore::String> Errors(2);
//        //... fill in pop.
//        Errors = TP->checkInputs(pop,taqlstr,iteraxes);
//        TP->plotData(pop,taqlstr);   //... or  TP->iterMultiPlotStart(pop,taqlstr,iteraxes)
//                
// </srcblock>
// </example>

// <motivation>
// This class was written to provide a high-level interface to an application
// that needs to plot and interact with data from multiple tables at once. 
// </motivation>

// <templating arg=T>
//    <li>
// </templating>

// <thrown>
//    <li>
//    <li>
// </thrown>


// <todo asof="$DATE:$">
//# A casacore::List of bugs, limitations, extensions or planned refinements.
//   <li> Add functions to read out the data being plotted (into python arrays).
//   <li> Remove "INFO" from the TablePlot::locateData() function. It is not required.
// </todo>

class TablePlot;

extern TablePlot *tableplot;

class TablePlot 
{
   public:
      // Get a reference to an instance of the TablePlot class.
      // TablePlot is a Singleton class, and user applications can 
      // obtain a reference to the single instance via this function.
      // The first call instantiates the class, and subsequent calls
      // return pointers to the TablePlot object.
      static TablePlot* TablePlotInstance() {
         if(!tableplot)
            tableplot = new TablePlot;
         return tableplot;
      }

      // Set a list of casacore::Table objects that are to be used for subsequent 
      //plotting.
      // These Tables are sent into BasePlot objects, which hold onto them. 
      // Tables are released in TablePlot::clearPlot, and a callback
      // (defined in TPCallBackHooks.h) signals this to the user application.
      // The casacore::Table objects passed in can be from Tables on disk, reference
      // Tables obtained via sub-selection, or in-memory Tables.
      // The root casacore::Table names are used for handling Flag Versions. 
      // To disable flag versions, the rootTable names can be empty strings.
      // The selection strings are purely for labelling purposes, and will
      // appear in the Locate and FlagHistory logger output.
      // ArrayColumns in the Tables sent in, need to have the same shape
      // across all rows in the Table. For MSs with varying shape columns,
      // the Tables should be broken up and sent in as a vector<table>.
      casacore::Bool setTableT(casacore::Vector<casacore::Table> &inTabObj, casacore::Vector<casacore::String> &rootTabNames, 
                     casacore::Vector<casacore::String> &selection);   

      // Set a list of Tables by specifying their names.
      // This is similar to TablePlot::setTableT(), but explicitly opens the 
      // Tables with the casacore::Table::Update mode, and no explicit lock acquisition.
      casacore::Bool setTableS(casacore::Vector<casacore::String> &inTabName, casacore::Vector<casacore::String> &rootTabNames, 
                     casacore::Vector<casacore::String> &selection);  

      // Set names of casacore::Table columns to be used for Flags.
      // "dataflags" must be the name of an casacore::ArrayColumn<bool>
      // "rowflags" must be the name of a casacore::ScalarColumn<bool>
      // These are also the columns that will be used for
      // Flag Versions.
      casacore::Bool useFlagColumns(casacore::String dataflags=casacore::String("FLAG"), 
                          casacore::String rowflags=casacore::String("FLAG_ROW"));
      
      // Check if at least one casacore::Table has been set.
      casacore::Bool isTableSet();
      
      // Choose between a GUI-based or completely command-line based interface.
      // This function will reset TablePlot, and instantiate the plotter
      // class in no-gui mode. Default is gui=true.
      // Note : This function has not been tested since the latest refactoring.
      casacore::Bool setGui(casacore::Bool guiVisible=true);

      // Save the plot to an image file, using the pl.savefig command.
      casacore::Bool saveFigure(const casacore::String& filename, const casacore::Int dpi, 
                      const casacore::String& orientation, const casacore::String& papertype, 
                      const casacore::String& facecolor, const casacore::String& edgecolor );
      
      casacore::Bool saveData(const casacore::String& filename);

      // Validation check for input parameters.
      // This function returns errors and warnings, based on the
      // validity of the input parameters. The returned vector<string> is of
      // length 2. The first casacore::String is Errors and the second is Warnings.
      // User applications should not proceed to TablePlot::plotData if 
      // checkInputs fails with errors. TablePlot::plotData re-verifies
      // inputs by calling this function internally as well.
      casacore::Vector<casacore::String> checkInputs(PlotOptions &pop, 
          casacore::Vector<casacore::Vector<casacore::String> > &datastrvector, casacore::Vector<casacore::String> &iteraxes);
      casacore::Vector<casacore::String> checkInputs(PlotOptions &pop, casacore::Vector<casacore::String> &datastr, 
          casacore::Vector<casacore::String> &iteraxes);

      // Check that the "FLAG(dataFlagColName_p)" 
      // casacore::ArrayColumn in input Tables have the same
      // shape across all rows.
      casacore::Vector<casacore::String> checkTables(casacore::Vector<casacore::Table> &tables);
      casacore::Bool checkShapes(casacore::Table &intab);

      // Plot data based on input plot options and [X,Y] TaQL string pairs.
      // Plot options are to be sent in via the PlotOptions class.
      // [X,Y] TaQL string pairs are sent in via datastrvector. 
      // datastrvector.nelements() must be the same as nTabObj.nelements(). 
      // i.e. the same number of Tables currently in use.
      // For each casacore::Table (the inner casacore::Vector), there can be any number of
      // XY TaQL pairs. [X1,Y1,X2,Y2,....] and plots from each pair will
      // appear as overlaid plots on the same panel.
      // each table has one set of taql 
      // each taql set has one or more taql pairs
      casacore::Bool plotData(PlotOptions &pop, casacore::Vector<casacore::Vector<casacore::String> > &datastrvector) {
         return plotData(pop,datastrvector,casacore::String("iter_off"));
      };

      // Plot data based on input plot options and [X,Y] TaQL string pairs.
      // Same as above, but uses only one [X,Y] TaQL pair for all 
      // the input Tables.
      casacore::Bool plotData(PlotOptions &pop, casacore::Vector<casacore::String> &datastr) {
         return plotData(pop,datastr,casacore::String("iter_off"));
      };
      
      // Refresh the plot. 
      // Flags are re-read from all active casacore::Table objects, and the
      // plots are re-made with all the current plot options for
      // each panel and layer.
      casacore::Bool rePlot();

      // Mark regions to flag. This function calls TPPlotter::markRegion().
      // These regions will apply to all layers of the chosen panel.
      casacore::Bool markRegions(casacore::Int nrows, casacore::Int ncols, casacore::Int panel, casacore::Vector<casacore::Double> &region);    

      // Flag all data in regions marked
      // using markFlags(). Direction is 1 for FLAG, and 0 for UNFLAG.
      casacore::Bool flagData(casacore::Int direction);     
      casacore::Bool flagData(casacore::Int direction, casacore::Vector<casacore::Int> regs);     
      
      // Print out information about the data selected using markFlags()
      // It takes in a list of table column names to make lists for.
      // INFO is a structure containing all the locate information.
      // This function automatically prints out INFO, and provides a 
      // callback for custom formatting as well, but this 
      // structure is still accessible here. INFO needs to be empty
      // before being passed in, and needs to be deleted after this call.
      // TODO - remove INFO from here, and make it an internal parameter. 
      // Remember to declare it inside, and clean it up at the end of this
      // function, after calling dumpLocate from inside.
      casacore::Bool locateData(casacore::Vector<casacore::String> columnlist, casacore::PtrBlock<casacore::Record*> &INFO, 
                      casacore::Int doFlag = -1);
      casacore::Bool locateData(casacore::Int doFlag = -1);
                
      // Clear all flags from all the tables attached to the vector of BasePlots
      // Dangerous.  Please use this carefully. It was put in mainly for  
      // debugging purposes.
      casacore::Bool clearAllFlags(casacore::Bool forRootTable);   

      
      // FlagVersion control functions.
      // They apply to all currently set Tables. i.e. all Tables set
      // by the previous TablePlot::setTableT call.
      // When multiple version names are supplied as a vector, the
      // operation is applied to each of them in turn.
      casacore::Bool  saveFlagVersion(casacore::String versionname, casacore::String comment, casacore::String merge);
      casacore::Bool  restoreFlagVersion(casacore::Vector<casacore::String> versionname, casacore::String merge);
      casacore::Bool  deleteFlagVersion(casacore::Vector<casacore::String> versionname);
      casacore::Bool  getFlagVersionList(casacore::Vector<casacore::String> &verlist);

      // Clear on or all plot panels.
      // Specify the co-ordinates of the panels to clear.
      // Setting all parameters to 0, indicates "clear all plots".
      // BasePlot objects are destroyed in this call, and this
      // triggers the "releaseTable" call-back for user-applications.
      casacore::Bool clearPlot(casacore::Int nrows=0, casacore::Int ncols=0, casacore::Int panel=0);   

      // Plotting from a series of subtables created by iterating over 
      // a specified iteration axis.       
      // Iterations are applied to all the latest Tables set via setTableT.
      // This mode of plotting supports multiple panels, with the
      // panel configuration set through plotoptions : nrows,ncols. 
      // Multiple iteration plots can be setup to run as overplots (using
      // the regular "overplot" plot option, and two successive 
      // iterMultiPlotStart calls before the first "next", 
      // as well as multipanels (see the "separateiterplot" plot option).
      casacore::Bool iterMultiPlotStart(PlotOptions &pop, 
         casacore::Vector<casacore::Vector<casacore::String> > &datastrvector,casacore::Vector<casacore::String> &iteraxes);
      casacore::Bool iterMultiPlotStart(PlotOptions &pop, casacore::Vector<casacore::String> &datastr,
         casacore::Vector<casacore::String> &iteraxes);

      // Advances to the next iteration. BasePlot objects persist, but 
      // their attached tables are replaced by the new subtables created 
      //by TableIter.
      casacore::Int iterMultiPlotNext(casacore::Vector<casacore::String> &labelcols, 
                            casacore::Vector<casacore::Vector<casacore::Double> > &labelvals);

      // Terminates the iterations. It is called automatically when the end 
      // of the table is reached. It can also be called before the 
      // iterations end, to cleanly terminate the
      // iterations when desired.
      casacore::Int iterMultiPlotStop( casacore::Bool /*rmplotter*/) {
         iterMultiPlotStop();
         return 0;
      };
      casacore::Bool iterMultiPlotStop();

      // Get the current list of Tables from TablePlot.
      casacore::Bool getTabVec(casacore::Vector<casacore::Table> &tabs);

      // Allow the user app to directly run a python command.
      casacore::Bool runPlotCommand(casacore::String command);
      
      //Change the state of a button on the Gui.
      // button : 'iternext','markregion','flag','unflag','locate'
      // state : 'enabled','disabled'
      // defaults -> iternext=disabled, all others are enabled.
      casacore::Bool changeGuiButtonState( casacore::String button, casacore::String state );

      // Clean up the TablePlot object.
      // Clear all plots and release all tables. This triggers the
      // "reset call-back" for user applications.
      // Note : It is dangerous to call this with closewindow=0.
      //        Extra rows of buttons might appear
      //        TODO - do something about this.
      casacore::Bool resetTP(casacore::Int closewindow=1); 

      // Set the "reset" callback function.
      // This callback is triggered when TablePlot::resetTP() is called.
      casacore::Bool setResetCallBack( casacore::String appname, TPResetCallBack * resetter );
      casacore::Bool clearResetCallBack( casacore::String appname );

   private:
      // Default Constructor
      TablePlot();

      // Copy Constructor (reference semantics)
      TablePlot(TablePlot const &);

      // reference =
      TablePlot& operator=(TablePlot const&);
                
      // Destructor
      ~TablePlot();

      // Instantiate a BasePlot/CrossPlot object for each table.
      // To be used in conjunction with upDateBP(). The second argument
      // should contain 'CROSS' for a CrossPlot object to be instantiated.
      // Otherwise (or left as default), BasePlot objects will be created.
      // One way to use this is to have the incoming TAQL pair have its 'X' TaQL
      // contain 'CROSS', since for CrossPlots, the x-axis corresponds to
      // column indices of the array column selected by the Y-TaQL.
      casacore::Bool createBP(casacore::PtrBlock<BasePlot* > &BPS, casacore::Vector<casacore::Vector<casacore::String> > &taqls);
                
      // Attach each table in the list to a BasePlot object. The first time this
      // function is called, it must be preceded by createBP(). Successive calls
      // to this function only change the tables/subtables being accessed. This
      // feature is used while iterating over an
      // iteration axis - subtables created in each iteration are assigned to a
      // fixed vector of BasePlot objects.
      casacore::Bool upDateBP(casacore::PtrBlock<BasePlot* > &BPS); 

      // Clear panels that lie partially or completely underneath
      // the plot currently being plotted. This is to mimic the
      // behaviour of native matplotlib, and to ensure automatic
      // cleanup of plots and associated data structures, when they
      // are no longer needed.
      casacore::Bool clearOverLaps(casacore::Int nrows, casacore::Int ncols, casacore::Int panel);   

      // Helper function to get the panel index for a particular
      // panel. If this panel does not exist in the list, -1 is returned.
      casacore::Int getPanelIndex(casacore::Int nrows, casacore::Int ncols, casacore::Int panel);

      // TP cleanup + initialization.
      casacore::Bool setupTP();

      // Destroy a BasePlot, and trigger the "releasetable" callback function
      // so that user apps know when to delete/close their tables.
      casacore::Bool deleteBasePlot(casacore::Int panel,casacore::Int layer);
     
      // Same as plotdata described above.
      // The last parameter is for internal use, and signals whether
      // TablePlot::plotData is being called in iteration or non-iteration mode.
      // This information is used to decide cleanup strategy.
      casacore::Bool plotData(PlotOptions &pop, casacore::Vector<casacore::Vector<casacore::String> > &datastrvector,
                    const casacore::String iterstate);
      casacore::Bool plotData(PlotOptions &pop, casacore::Vector<casacore::String> &datastr, 
                    const casacore::String iterstate);
                
      // Read out Flagging information from BasePlots. 
      // This is to be called soon after a TablePlot::MarkRegions and
      // TablePlot::flagData sequence. "flaghist" contains information
      // about the table and plot, the regions marked, and the number of
      // points selected for flagging/unflagging.
      casacore::Bool updateFlagHistory(casacore::PtrBlock<casacore::Record*> &flaghist);     
               
      // Print out (formatted) info about the flag history.
      casacore::Bool dumpFlagHistory(casacore::PtrBlock<casacore::Record*> &flaghist);     

      // Print out (formatted) info about the data points selected
      // during a locate operation.
      // This function triggers the call-back for "printlocater" for
      // user apps to perform custom formatting of "locate" output.
      casacore::Bool dumpLocateInfo(casacore::PtrBlock<casacore::Record*> &INFO, casacore::Int flagdata = -1, 
               casacore::Bool ave = 0);

      // Generate an exception.
      void TablePlotError(casacore::String msg);

      // Private class variables.
      // number of currently set Tables
      casacore::Int nTabs_p;                      
      // vector of currently set Tables
      casacore::Vector<casacore::Table> TABS_p;             
      // root casacore::Table names per Table.
      casacore::Vector<casacore::String> rootTabNames_p;    
      // Strings with selection info oer Table.
      casacore::Vector<casacore::String> tableSelections_p; 
      // "FLAG"
      casacore::String dataFlagColName_p;         
      // "FLAG_ROW"
      casacore::String rowFlagColName_p;          
      
      // Plotter class
      TPPlotter *TPLP;                        
      // Holds all active BasePlots.
      casacore::PtrBlock<casacore::PtrBlock<BasePlot* >*> ATBPS;
      // Holds tableiters for each Table.
      casacore::Vector<casacore::TableIterator> Iters_p;          

      // casacore::List of parameters per panel.
      casacore::PtrBlock<PanelParams* > PAN;      

      // resetTP callback pointer.
      casacore::PtrBlock<TPResetCallBack *> resetters_p; 
      // resetTP callback user-app name
      casacore::Vector<casacore::String> appnames_p;               
                
      // casacore::Input TaqlStrings
      casacore::Vector<casacore::String> DataStr_p;     
      // casacore::Input iteration axes strings.
      casacore::Vector<casacore::String> IterAxes_p;    
      
      // Panel co-ordinates of the latest panel.
      casacore::Int NRows_p,NCols_p;   
      // Temporary instance of Plot options.
      PlotOptions Pop_p;     
      
      // Variables for iteration plots.
      // Overplots and parallel iterations are controlled
      // via these variables.
      casacore::PtrBlock<PlotOptions*> Pops_p;
      casacore::PtrBlock<casacore::Vector<casacore::Vector<casacore::String> >* > TaqlStr_p;
      casacore::PtrBlock<casacore::Vector<casacore::TableIterator>* > OvpIters_p;
      casacore::PtrBlock<casacore::Vector<casacore::String>* > OvpRootTabNames_p;
      casacore::PtrBlock<casacore::Vector<casacore::String>* > OvpTabSelections_p;
      casacore::Vector<casacore::String> TitleStrings_p;
      casacore::String Separate_p;
      
      // Current panel number
      casacore::Int CurrentPanel_p;    
      // Flag to signal that setTableT has been called and that
      // BasePlots needs to be created/updated with the new Tables.
      casacore::Int TableTouch_p;
      // Flags to signal iterplot mode.
      casacore::Bool IterPlotOn_p;
      casacore::Bool IterPlotStarted_p;
      // Flag to check if there is a currently displayed plot or not.
      casacore::Bool isPlotDisplayed_p;
      // Flag to indicate if Plot window is currently displayed.
      casacore::Bool isGuiVisible_p;
      
      // Debug message flags.
      casacore::Timer tmr;
      
      // Logger instance.
      SLog* log;
      static casacore::String clname;
};



// <summary>
// Class for GUI binding - Casa C++ to Python and matplotlib. 
// </summary>

// <use visibility=export>

// <reviewed reviewer="" date="" tests="">
// </reviewed>

// <prerequisite>
//#! Classes or concepts you should understand before using this class.
//   <li> PlotterGlobals.cc
// </prerequisite>

// <etymology>
// Class for GUI binding
// </etymology>

// <synopsis>
// This class controls the actions to be performed by Plotter-specific
// buttons on the plotter matplotlib Gui.
// One global instance of this class is maintained, and is accessed
//   from the Python-C++ binding functions defined in PlotterGlobals.cc 
// </synopsis>

// <example>
//
// (See TablePlot.cc, PlotterGlobals.cc)
//
// <srcblock>
// </srcblock>
// </example>

// <motivation>
// This class was written because the python/C++ binding had to be done as
// global functions, and they need access to the current TablePlot.
// </motivation>

// <templating arg=T>
//    <li>
// </templating>

// <thrown>
//    <li>
//    <li>
// </thrown>


// <todo asof="$DATE:$">
//# A casacore::List of bugs, limitations, extensions or planned refinements.
// </todo>


class TPGuiBinder
{
   public :
      // Constructor
      TPGuiBinder( casa::TablePlot* intp ){
         itsTablePlot = intp;
      }
      // Destructor
      ~TPGuiBinder(){};
      // Gui callback to read in regions.
      casacore::Bool markregion(casacore::Int nrows,casacore::Int ncols, casacore::Int panel, 
                      casacore::Vector<casacore::Double> region) {
         return itsTablePlot->markRegions(nrows, ncols, panel, region); 
      }
      // Gui callback for flagdata
      casacore::Bool flagdata() {
         //casacore::Bool ret = itsTablePlot->flagData(1);
         casacore::Bool ret = itsTablePlot->locateData(1);
         return ret;
      }
      // Gui callback for unflagdata
      casacore::Bool unflagdata() { 
         //casacore::Bool ret = itsTablePlot->flagData(0); 
         casacore::Bool ret = itsTablePlot->locateData(0); 
         return ret;
      }
      // Gui callback for locatedata
      casacore::Bool locatedata() { 
         //casacore::Vector<casacore::String> clist;
         //clist.resize(0);
         //casacore::PtrBlock<casacore::Record*> INFO;
         //casacore::Bool ret = itsTablePlot->locateData(clist, INFO);
         /* Clean up the casacore::PtrBlock */
         //for(casacore::Int i=0;i<(casacore::Int)INFO.nelements();i++)
         //   delete INFO[i];
         //INFO.resize(0);

         casacore::Bool ret = itsTablePlot->locateData();
         return ret;
      }
      // Gui callback for iterplotnext
      casacore::Bool iterplotnext() { 
         casacore::Vector<casacore::String> labcol;
         casacore::Vector<casacore::Vector<casacore::Double> > labval;
         if(itsTablePlot->iterMultiPlotNext(labcol,labval) == -1)
            itsTablePlot->changeGuiButtonState("iternext","disabled");
         return true;
      }
      // Gui callback for iterplotstop.
      // Currently not in use.
      casacore::Bool iterplotstop(){ 
         return itsTablePlot->iterMultiPlotStop(); 
      }
      // Gui callback for clearplot.
      // Currently not in use.
      casacore::Bool clearplot(){ 
         return itsTablePlot->clearPlot(0,0,0);
      }
      // Gui callback for "Quit".
      casacore::Bool quit(casacore::Int closewindow){ 
         return itsTablePlot->resetTP(closewindow); 
      }

   private:
      casa::TablePlot* itsTablePlot;

};


} //# NAMESPACE CASA - END 

//#ifndef AIPS_NO_TEMPLATE_SRC
//#include <tables/TablePlot/TablePlot.cc>
//#endif //# AIPS_NO_TEMPLATE_SRC
#endif

