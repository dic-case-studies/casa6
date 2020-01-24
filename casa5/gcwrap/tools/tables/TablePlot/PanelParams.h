//# PanelParams.h: this defines all of the optional plotting parameters.
//# Copyright (C) 2007-2008
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
//#
//#
//# -------------------------------------------------------------------------
//# Change Log
//# Date   Name       Description
//# 10/29/2007  S.Jaeger    Added new messaging
//# 11/14/2007  S.Jaeger    Added PlotRangesSet option
#ifndef PANELPARAMS_H
#define PANELPARAMS_H

#include <casa/aips.h>
#include <casa/iostream.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Containers/Record.h>
#include <casa/OS/Time.h>
#include <casa/IO/AipsIO.h>
#include <casa/Quanta/MVTime.h>
#include <casa/Quanta/QC.h>

#include <tools/tables/TablePlot/TPCallBackHooks.h>
#include <tables/Tables/Table.h>

namespace casa { //# NAMESPACE CASA - BEGIN
// <summary>
// Classes to hold Plot Options.
// </summary>

// <reviewed reviewer="" date="" tests="">
// </reviewed>

// <prerequisite>
//   <li> 
//   <li> 
// </prerequisite>


// <etymology>
// Plot Options ! 
// </etymology>

// <synopsis>
// User applications need to create instances of PlotOptions, fill them in,
// and send one in per TablePlot::plotData call. These options are then stored
// in the PanelParams class, and TablePlot holds a bunch of PanelParams, 
// corresponding to its BasePlots.
// </synopsis>

// <motivation>
// The need to have all plot options originate from one single place, and to
// have their defaults also set in one single place.
// </motivation>

// <thrown>
//    <li>
//    <li>
// </thrown>


// <todo asof="$DATE:$">
//   <li> 
// </todo>


class PlotOptions
{
   public:
      //Constructor
      PlotOptions();
                
      //Destructor
      ~PlotOptions();

      // Reset to defaults 
      casacore::Bool reset();
       
      // Return a string of current values to print (for debugging) 
      casacore::String print();

      // Validate Current entries 
      // Returns 2 Strings, Errors and Warnings.
      casacore::Vector<casacore::String> validateParams();

      // Fill entries from a record 
      casacore::String fillFromRecord(casacore::Record &rpop);
      
      // operator= : copy semantics
      PlotOptions& operator= (PlotOptions &inplop);
      
      
      // Layer Independant Plot Parameters 

      // default : [1,1,1].
      // [nrows, ncols, panel_index] from matplotlib
      casacore::Vector<casacore::Int> PanelMap;      
                
      // default : []. 
      // [xmin,xmax,ymin,ymax]
      // If min=max, default data range is used.
      casacore::Vector<casacore::Double> PlotRange;
      
      // default : []. 
      // [xminSet,xmaxSet,yminSet,ymaxSet]
      // Bitmap to indicate which range values have been set
      // in the PotRange option.
      casacore::Vector<casacore::Bool> PlotRangesSet;  
                
      // default : 'o' 
      // 'o' : off -> no time formatting 
      // 'x' : time formatting for x axis
      // 'y' : time formatting for y axis
      // 'b' : time formatting for x and y axes
      casacore::String TimePlotChar;       
                
      // default : true : 
      // To use with CrossPlot.
      // true : Columns of the casacore::Array are the x axis.
      // false : Rows of the casacore::Array are the x axis.
      // (example : for a casacore::MeasurementSet and "DATA"
      //            true : channels on x-axis
      //            false : correlations on x-axis )
      casacore::Bool ColumnsXaxis;      
                
      // Default : false
      // true : This plot must sit on top of existing plots
      // false : This plot will be a fresh plot
      //         clearing away all existing layers
      casacore::Bool OverPlot;             
                
      // Default : false
      // true : OverPlot=false will replace only the    
      //        top layer, and will not clear away
      //        all existing layers. 
      //        (example : after making an overplot,
      //         the plotsymbol for the top-most layer
      //         is to be changed  without having to
      //         replot all existing layers)
      // false : OverPlot=false will clear up all 
      //         existing layers, before plotting.
      casacore::Bool ReplaceTopPlot;       
                
      // Default : true
      // true : Make TablePlot mimic the matplotlib
      //        behaviour of automatically clearing panels
      //        when a new panel overlaps on it.
      // false : Do not do the above. The user will 
      //         have to explicitly clear a panel. 
      //         (example : To make a tiny plot
      //          that sits in a corner, inside an
      //          existing larger plot.)
      //
      //
      casacore::Bool RemoveOldPanels;      
                
      // Default : 12 . matplotlib font size
      // This is the title size. The x,ylabels are
      // 80% of this size.
      casacore::Double FontSize; 
                
      // X axis label
      casacore::String XLabel;         
      
      // Y axis label
      casacore::String YLabel;         
      
      // Title
      // For multiline labels, 
      // beware of newline characters. These strings
      // must contain "\\n" in them for newlines.
      casacore::String Title;         
      
      // Default : 8.0 ( in cm ). matplotlib convention
      // Must be same for all panels !
      casacore::Double WindowSize;      
      
      // Default : 1.0 . -> height/width
      // Must be same for all panels !
      casacore::Double AspectRatio;      
      
      // Default : NULL
      // If left as NULL, it gets replaced by a
      // pointer to TPConvertVase. This is the base
      // class that does not modify values read
      // from the table.  Modifications to TaQL
      // expressions, that cannot be written as 
      // TaQL strings, can go
      // in here, and will get applied to every
      // value being read from the TaQL result.
      TPConvertBase *Convert;      
      
      // Default : TPGuiCallBackHooks instance.
      // Applications can supply their own custom
      // function for application specific formatting.
      // (example : chan numbers to frequency vals.)
     TPGuiCallBackHooks *CallBackHooks; 
     
     // Default : false
     // true : If a Scalar or casacore::Vector reduction has been
     //        done using MEAN(), the values are rescaled to account for
     //        the number of flagged values being averaged
     //        by TaQL. Note : Try not to use this.
     // false : Don't do any scaling correction. This is if
     //         averaging is via SUM(...)/SUM()
     casacore::Bool DoScalingCorrection;  
     
     // Default : 'none'
     // 'row' : iteration plots from multiple
     //         TaQLs can be run parallely. The
     //         panels are separated into rows.
     // 'col' : iteration plots from multiple
     //         TaQLs can be run parallely. The
     //         panels are separated into cols.
     casacore::String SeparateIter;       
     
     // Default : false  ( only for CrossPlot )
     // false : Compute the average x axis value
     //        as the middle of the range being avgd
     // true : Compute the average x axis value
     //        accounting for flagged cell rows/cols.
     casacore::Bool HonourXFlags;         
     
     // Default : []
     // This is the list of casacore::Table columns that
     // will get sent into the "locate" function
     // when triggered by the GUI.
     // It will be overridden by anything specified in
     // a DumpLocateInfoBase class.
     casacore::Vector<casacore::String> LocateColumns; 
     
     
     // Layer Dependant Plot Parameters 
     
     // Default : ',' . matplotlib plotsymbols
     casacore::String PlotSymbol;      
     
     // Default : '' . matplotlib colour string
     // Can be a predefined pylab colour 'brown', or
     // an html hex string '#7FFF4e', or '(r,g,b)'.
     // If specified (length>0), this takes 
     // precedence over the colour specified via
     // PlotSymbol. If a non-predefined colour is
     // specified, MultiColour is always false.
     casacore::String ColourString;      
     
     // Default [] : If specified, it will apply to 
     // the first N plotted points.
     casacore::Vector<casacore::String> PointLabels;
     
     // Default : 10.0 . matplotlib markersize
     casacore::Double MarkerSize;      
     
     // Default : 2.0 . matplotlib linewidth
     casacore::Double LineWidth;      
     
     // Default : 'none'
     // 'cellrow' : Cell rows get different colours
     // 'cellcol' : Cell cols get different colours
     // 'both': rows and cols get different colours
     // 'none': rows and cols get the same colour
     casacore::String MultiColour;        
     
     // Default : true
     // true : When multiple Tables are sent in simultaneously
     //        into TablePlot, they go to different Layers and
     //        automatically increment colours.
     // false : All layers begin with the same colour
     casacore::Bool TableMultiColour;     
     
     // Default : false
     // true : Plot flagged points in purple
     // false : Plot unflagged points in specified colour
     // (example : Do an overplot with false and then
     //  true, to see flagged and unflagged points in
     //  different colours.)
     casacore::Bool ShowFlags;          
     
     // Default : "main"
     // Name of the flag version to use while plotting
     // If it doesn't find this name, it uses "main".
     casacore::String FlagVersion;      
     
     // Default : 1
     // Start with the first point, and then 
     // plot only if npoints % SkipNRows == 0
     casacore::Int SkipNRows;           
     
     // Default : 1
     // No averaging. 
     // If >1, averages every N points.
     // Users of plotoptions, need to order/sort/select
     // the input casacore::Table, so that such averaging is
     // accurate.
     casacore::Int AverageNRows;        
     
     casacore::String FlagExt;
     // Default : 'none' : no points are connected by lines.
     // If 'tablerow' : number of plots = nchans x ncorrs
     //                 number of points per plot = nrows
     //             -> points along table rows are connected.
     // If 'cellcol'  : number of plots = nrows x ncorrs
     //                 number of points per plot = nchans
     //             -> points along cell cols are connected
     // If 'cellrow'  : number of plots = nrows x nchans
     //                 number of points per plot = ncorrs
     //             -> points along cell rows are connected.
     // Currently, 'cellcol' and 'cellrow' do the same
     // thing. Nplots = nrow, points per plot = nchanxncorr.
     casacore::String Connect;            
     
     // To allow PanelPArams to get access to some private variables 
     inline casacore::Int getParsedParams(casacore::Int &timeplot, casacore::Int &plotcolour, 
                               casacore::String &pyplotsym) {
         timeplot=TimePlot_p; 
         plotcolour=PlotColour_p;
         pyplotsym=PyPlotSymbol_p;
         return 0;
     };
     
     // Allowed, predefined colours 
     casacore::Int NColours;
     casacore::Bool useLayerColor;

     // A long list of supported colour names.
     casacore::Vector<casacore::String> ColourList;
     casacore::String pylabcolourstring;
     
     //for chan/time averaging
     casacore::Bool doAverage;
     casacore::Matrix<casacore::Int> ChanMap;
     casacore::Matrix<casacore::Int> RowMap;
     casacore::String MSName; 
     casacore::String spwExpr; 

   private:
     
     // 0:1:2:3 -> 'o','x','y','b'
     casacore::Int TimePlot_p;       
     // an casacore::Int for the colour part of plotsymbol
     casacore::Int PlotColour_p;      
     // just the plotsymbol (without colour)
     casacore::String PyPlotSymbol_p;      
     
     // Has the 'Convert' TPConvertXXX been
     // created inside PlotOptions, or by the user ?
     // true : created inside here => needs to be
     //        deleted inside here too, during cleanup. 
     // false : created outside and supplied in
     //         => don't try to delete it !
     casacore::Bool ConvertIsLocal_p;       
                
};

// <summary>
// Classes to hold Plot Options plus other parameters for a panel.
// </summary>

// <reviewed reviewer="" date="" tests="">
// </reviewed>

// <prerequisite>
//   <li> 
//   <li> 
// </prerequisite>


// <etymology>
// Panel Parameters.
// </etymology>

// <synopsis>
// Holds PlotOptions for each BasePlot on one panel. Layers are created via
// overplots, as well as by making one plot using multiple Tables.
// There are layer independant and layer dependant parameters.
// Replots after flag editing, need to remember plot options for all
// layers, so this class does it. TablePlot holds a PanelParams object for each panel.
// </synopsis>

// <motivation>
// The need to remember plot options for multiple plot layers, per panel.
// </motivation>

// <thrown>
//    <li>
//    <li>
// </thrown>


// <todo asof="$DATE:$">
//   <li> 
// </todo>


class PanelParams 
{
   public:
      // Constructor
      PanelParams();

      // Destructor
      ~PanelParams();

      // Change the number of layers
      // Usually called only to increase the layer count, 
      // when overplots are done one-by-one.
      casacore::Int changeNlayers(casacore::Int nlayers);

      // Fill in parameters from the most current PlotOptions.
      casacore::Int updateLayerParams();

      // Reset all parameters.
      casacore::Bool reset();
      
      // A PlotOptions instance to read inputs
      // and to hold layer-independant parameters.
      PlotOptions Plop;
      
      // Layer INDEPENDANT parameters 
      // The current plot-range for this panel.
      casacore::Vector<casacore::Double> PanelZrange;         

      // A list of marked regions for this panel.
      // FlagList(panelId, regionCount); 
      casacore::Vector<casacore::Vector<casacore::Double> > FlagList;

      // The layer number for the topmost layer.
      // NOTE : A layer is a plot from one BasePlot object.
      //        Sometimes, a single TablePlot::plotData call,
      //        generates plots from multiple BasePlots.
      //        All layers from one single TablePlot::plotData
      //        call, get the same "layerNumber".
      casacore::Int MaxLayer;

      // Integer form of the TimePlot parameter. 
      casacore::Int TimePlot;

      // Ingeter form of the plot colour.
      // Needed to allow automatic colour incrementing !
      casacore::Int PlotColour;

      // Python plot symbol.
      casacore::String PyPlotSymbol;

      // LAYER DEPENDANT parameters
      // Number of BasePlot-layers.
      casacore::Int nBP;

      // The "layerNumber" for each BasePlot-Plot
      casacore::Vector<casacore::Int> LayerNumbers;    

      // Other parameters.
      casacore::Vector<casacore::Int> LayerColours;
      casacore::Vector<casacore::String> LayerSymbols;
      casacore::Vector<casacore::Vector<casacore::String> > LayerPointLabels;
      casacore::Vector<casacore::Vector<casacore::String> > LayerLocateColumns;
      casacore::Vector<casacore::Double> LayerMarkerSizes;
      casacore::Vector<casacore::Double> LayerLineWidths;
      casacore::Vector<casacore::Bool> LayerShowFlags;     
      casacore::Vector<casacore::String> LayerMultiColours;
      casacore::Vector<casacore::String> LayerFlagVersions;
      casacore::Vector<casacore::Int> LayerSkipNRows;
      casacore::Vector<casacore::Int> LayerAverageNRows;
      casacore::Vector<casacore::String> LayerFlagExt;
      casacore::Vector<casacore::String> LayerConnects;
      casacore::Vector<casacore::Vector<casacore::String> > LayerXYTaqls;


};


} //# NAMESPACE CASA - END 


#endif

