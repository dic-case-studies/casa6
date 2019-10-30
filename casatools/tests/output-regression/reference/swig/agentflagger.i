/******************** generated by xml-casa (v2) from agentflagger.xml **************
********************* 4b7be1d48b1bb20557804b5d31ee6f39 *****************************/
%module agentflagger
%include <casa_typemaps.i>
%feature("kwargs");
%feature("autodoc", "0");

%feature("docstring", "

Summary:
    Destroy the flag tool

Example:

af.done()

--------------------------------------------------------------------------------
") done;

%feature("docstring", "

Summary:
    Open the MS or a calibration table and attach it to the tool.

Input Parameters:
    msname                    Measurement set or calibration table to be processed. Default:
    ntime                     Time interval. If not given, the default will be used. Default:

Example:

af.open(msname,ntime)

--------------------------------------------------------------------------------
") open;

%feature("docstring", "

Summary:
    Select the data based on the given parameters. For unspecified parameters, the full data range is assumed. All data selection parameters follow the MS Selection syntax. 

Input Parameters:
    config                    The record (dictionary) config may be given or not. If it is not given, n and no specific parameter is given either, the whole MS will be selected. The record n may contain any data selection parameters supported by MS Selection such as:n 
    field                     Field indices or source names : example : '2,3C48'
    spw                       Spectral Window Indices or names : example : '1,2'
    array                     Array Indices or names : example : 'VLAA'
    feed                      Feed index or name : example : '1,2' (not supported yet)
    scan                      Scan number : example : '1,2,3'
    antenna                   Baseline number : example : '2,3,4,5'
    uvrange                   UV-distance range, with a unit : example : '2.0-3000.0 m'
    timerange                 Time range, as MJDs or date strings : example : 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations : example : 'RR,LL,RL,LR,XX,YY,XY,YX,Sol1'
    intent                    Scan intent : example : '*CAL*, *BAND*'
    observation               Observation Id : example : '2~4'

Example:

Select the whole MS
af.selectdata()

Select a portion of the MS
myrecord={}
myrecord['scan']='1~3'
myrecord['spw']='0:1~10'
af.selectdata(myrecord)

Another way to select a portion of the MS
af.selectdata(scan='3~5', spw='0')

--------------------------------------------------------------------------------
") selectdata;

%feature("docstring", "

Summary:
    Parse the parameters for the agent (flagging mode).

Description:

The specific data selection parameters for the agent (flagging mode)
are parsed. These parameters are the data selection and mode-specific parameters.
See the example below:

Input Parameters:
    aparams                   It takes a record (dictionary) with the specific parameters for the flagging mode. n The record may contain any data selection parameters supported by MS Selection, as well as n mode-specific parameters such as:n n (1) array,feed,scan,field,spw,intent,correlation,antenna,uvrange,observation n (2) mode (which can be: manual,clip,quack,shadow,elevation,tfcrop,extendflags,unflag or summary) n     For flagging mode=clip, the parameters are: expression, datacolumn, clipminmax, etc. n     See the documentation of the task flagdata for all the available parameters for each mode.n n (3) apply: default is true (true for flagging and false for unflagging) n Example:nn myrecord={}n myrecord['mode']='clip'n myrecord['scan']='1~3'n myrecord['clipminmax']=[0.02,0.3]n myrecord['apply']=Truen af.parseagentparameters(myrecord)n  

Example:

myrecord={}
myrecord['mode']='clip'
myrecord['scan']='1~3'
myrecord['clipminmax']=[0.02,0.3]
myrecord['apply']=True
af.parseagentparameters(myrecord)

--------------------------------------------------------------------------------
") parseagentparameters;

%feature("docstring", "

Summary:
    Initialize the agents

Description:

This method will initialize the agents and create a list of agents
with their specific parameters. It takes no parameters.

Example:

af.init()

--------------------------------------------------------------------------------
") init;

%feature("docstring", "

Summary:
    Execute a list of flagging agents

Description:


Execute a list of flagging agents and write or not to the MS/cal table. The parameter
writeflags controls whether or not to write to the MS.

Input Parameters:
    writeflags                Write flags to MS
    sequential                Run the agents in the order they are inserted in the list or not. Default is True to run in the original order.

Example:

af.run()

--------------------------------------------------------------------------------
") run;

%feature("docstring", "

Summary:
    Print out a list of saved flag_versions.

Description:


Print out the list of flag versions in the MS, unless the parameter printflags=False. The list of names is returned.

Input Parameters:
    printflags                Print flagversions in logger?

--------------------------------------------------------------------------------
") getflagversionlist;

%feature("docstring", "

Summary:
    Print out a list of current flag selections. 

Description:


Print out a list of current flag selections.

--------------------------------------------------------------------------------
") printflagselection;

%feature("docstring", "

Summary:
    Save current flags with a version name. 

Input Parameters:
    versionname               Version name
    comment                   Comment for this flag table
    merge                     merge type

--------------------------------------------------------------------------------
") saveflagversion;

%feature("docstring", "

Summary:
    Restore flags from a saved flag_version. n versionname : name of flag version to restore to main table n merge : Type of operation to perform during restoration. n        merge = replace  : replaces the main table flags. n        merge = and   : logical AND with main table flags n        merge = or    : logical OR with main table flags n        Default : replace. n 

Input Parameters:
    versionname               Version name
    merge                     merge type

--------------------------------------------------------------------------------
") restoreflagversion;

%feature("docstring", "

Summary:
    Delete a saved flag_version.

Input Parameters:
    versionname               Version name

--------------------------------------------------------------------------------
") deleteflagversion;

%feature("docstring", "

Summary:
    Parse data selection parameters and specific parameters for the manual mode. Data selection follows the MS Selection syntax. 

Input Parameters:
    field                     Field indices or source names. Example: '2,3C48'
    spw                       Spectral Window Indices or names. Example: '1,2'
    array                     Array Indices or names. Example: 'VLAA'
    feed                      Feed index or name. Example: '1,2' (not supported yet)
    scan                      Scan number. Example: '1,2,3'
    antenna                   Baseline number. Example: '2,3,4,5,132'
    uvrange                   UV-distance range, with a unit. Example: '2.0-3000.0 m'
    time                      Time range, as MJDs or date strings. Example: 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations. Example: 'RR,LL,RL,LR,XX,YY,XY,YX'
    intent                    Scan intent. Example: '*CAL*, *BAND*'
    observation               Observation Id. Example: '2~4'
    autocorr                  Parameter to flag only auto-correlations. Default:
    apply                     Parameter to flag or unflag the data. Default:

Example:

af.parsemanualparameters(autocorr=True)

--------------------------------------------------------------------------------
") parsemanualparameters;

%feature("docstring", "

Summary:
    Parse data selection parameters and specific parameters for the clip mode. Data selection follows the MS Selection syntax. 

Input Parameters:
    field                     Field indices or source names. Example: '2,3C48'
    spw                       Spectral Window Indices or names. Example: '1,2'
    array                     Array Indices or names. Example: 'VLAA'
    feed                      Feed index or name. Example: '1,2' (not supported yet)
    scan                      Scan number. Example: '1,2,3'
    antenna                   Baseline number. Example: '2,3,4,5'
    uvrange                   UV-distance range, with a unit. Example: '2.0-3000.0 m'
    time                      Time range, as MJDs or date strings. Example: 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations. Example: 'RR,LL,RL,LR,XX,YY,XY,YX'
    intent                    Scan intent. Example: '*CAL*, *BAND*'
    observation               Observation Id. Example: '2~4'
    datacolumn                Data column to use for clipping. Supported columns for cal tables are FPARAM,CPARAM,SNR. Example: 'DATA'. Default: 
    clipminmax                Range to use for clipping. Example: [100.0,200.0] 
    clipoutside               Clip points outside this range? [True/False]. Default: 
    channelavg                Average data over channels before clipping? [True/False]. Default: 
    chanbin                   Width (bin) of input channels to average to form an output channel.
    timeavg                   Average data over time ranges. [True/False]. Default: 
    timebin                   Bin width for time average. Example: '2s'
    clipzeros                 Clip zero-value data. [True/False]. Default: 
    apply                     Parameter to flag or unflag data. Default:

Example:

The correlation parameter can be used with an operator for the clip mode.
The operator should be written only once and it will affect all the polarizations
given in the string. See the example below:

af.parseclipparameters(clipzeros=True, clipminmax=[0.,4.], correlation='ABS_XX,XY')

or for a calibration table:
af.parseclipparameters(clipzeros=True, clipminmax=[0.,4.], correlation='Sol1')

--------------------------------------------------------------------------------
") parseclipparameters;

%feature("docstring", "

Summary:
    Parse data selection parameters and specific parameters for the quack mode. Data selection follows the MS Selection syntax. 

Input Parameters:
    field                     Field indices or source names. Example: '2,3C48'
    spw                       Spectral Window Indices or names. Example: '1,2'
    array                     Array Indices or names. Example: 'VLAA'
    feed                      Feed index or name. Example: '1,2' (not supported yet)
    scan                      Scan number. Example: '1,2,3'
    antenna                   Baseline number. Example: '2,3,4,5'
    uvrange                   UV-distance range, with a unit. Example: '2.0-3000.0 m'
    time                      Time range, as MJDs or date strings. Example: 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations. Example: 'RR,LL,RL,LR,XX,YY,XY,YX'
    intent                    Scan intent. Example: '*CAL*, *BAND*'
    observation               Observation Id. Example: '2~4'
    quackmode                 Quack mode. Default:
    quackinterval             Quack length in seconds. Default:
    quackincrement            Flag incrementally in time. Default:
    apply                     Parameter to flag or unflag the data. Default:

Example:

af.parsequackparameters(scan='1~3', quackmode='beg', quackinterval=1)

--------------------------------------------------------------------------------
") parsequackparameters;

%feature("docstring", "

Summary:
    Parse data selection parameters and specific parameters for the elevation mode. Data selection follows the MS Selection syntax. 

Input Parameters:
    field                     Field indices or source names. Example: '2,3C48'
    spw                       Spectral Window Indices or names. Example: '1,2'
    array                     Array Indices or names. Example: 'VLAA'
    feed                      Feed index or name. Example: '1,2' (not supported yet)
    scan                      Scan number. Example: '1,2,3'
    antenna                   Baseline number. Example: '2,3,4,5'
    uvrange                   UV-distance range, with a unit. Example: '2.0-3000.0 m'
    time                      Time range, as MJDs or date strings. Example: 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations. Example: 'RR,LL,RL,LR,XX,YY,XY,YX'
    intent                    Scan intent. Example: '*CAL*, *BAND*'
    observation               Observation Id. Example: '2~4'
    lowerlimit                The limiting elevation in degrees. Data obtained at lower antenna elevations will get flagged. Default: 
    upperlimit                The limiting elevation in degrees. Data obtained at higher antenna elevations will get flagged. Default: 
    apply                     Parameter to flag or unflag the data. Default:

Example:

To unflag, set the apply parameter.
af.parseelevationparameters(upperlimit=50.,lowerlimit=10.0, apply=False)

--------------------------------------------------------------------------------
") parseelevationparameters;

%feature("docstring", "

Summary:
    Parse data selection parameters and specific parameters for the time and frequency mode. Data selection follows the MS Selection syntax. 

Input Parameters:
    field                     Field indices or source names. Example: '2,3C48'
    spw                       Spectral Window Indices or names. Example: '1,2'
    array                     Array Indices or names. Example: 'VLAA'
    feed                      Feed index or name. Example: '1,2' (not supported yet)
    scan                      Scan number. Example: '1,2,3'
    antenna                   Baseline number. Example: '2,3,4,5'
    uvrange                   UV-distance range, with a unit. Example: '2.0-3000.0 m'
    time                      Time range, as MJDs or date strings. Example: 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations. Example: 'RR,LL,RL,LR,XX,YY,XY,YX'
    intent                    Scan intent. Example: '*CAL*, *BAND*'
    observation               Observation Id. Example: '2~4'
    ntime                     Time-range to use for each chunk (in seconds or minutes). Default:
    combinescans              Accumulate data across scans depending on the value of ntime. Default:
    datacolumn                Data column to use for clipping. Example: 'DATA'. Default: 
    timecutoff                Flagging thresholds in units of deviation from the fit. Default:
    freqcutoff                Flagging thresholds in units of deviation from the fit. Default:
    timefit                   Fitting function for the time direction (poly/line). Default:
    freqfit                   Fitting function for the frequency direction (poly/line). Default:
    maxnpieces                Number of pieces in the polynomial-fits (for 'freqfit' or 'timefit' = 'poly'). Default:
    flagdimension             Dimensions along which to calculate fits (freq/time/freqtime/timefreq). Default:
    usewindowstats            Calculate additional flags using sliding window statistics (none,sum,std,both). Default:
    halfwin                   Half-width of sliding window to use with 'usewindowstats' (1,2,3). Default:
    extendflags               Extend the flags in time, frequency and correlations. Default:
    apply                     Parameter to flag or unflag the data. Default:
    channelavg                Average data over channels before clipping? [True/False]. Default: 
    chanbin                   Width (bin) of input channels to average to form an output channel.
    timeavg                   Average data over time ranges. [True/False]. Default: 
    timebin                   Bin width for time average. Example: '2s'

Example:

The correlation parameter can be used with an operator for the tfcrop mode.
The operator should be written only once and it will affect all the polarizations
given in the string. Note that if ntime='scan' and combinescans=True, all the scans will
be loaded at once, thus requesting a lot of memory depending on the available spws.
The parameter combinescans should be set to True only when ntime is specified as a
time-interval (not 'scan'). By default, the flags will be extended in time, if
more than 50% of the timeranges are flagged, 80% of the channels are flagged and
it will extend to other polarizations too. This is similar to running the extend
mode after running tfcrop on the MS.

af.parsetfcropparameters(spw='9', ntime=10.0, combinescans=True, correlation='ABS_XX,XY',
extendflags=True)

--------------------------------------------------------------------------------
") parsetfcropparameters;

%feature("docstring", "

Summary:
    Parse data selection parameters and specific parameters for the antint mode. Data selection follows the MS Selection syntax. 

Input Parameters:
    field                     Field indices or source names. Example: '2,3C48'
    spw                       Spectral Window Indices or names. Example: '1,2'
    array                     Array Indices or names. Example: 'VLAA'
    feed                      Feed index or name. Example: '1,2' (not supported yet)
    scan                      Scan number. Example: '1,2,3'
    antenna                   Baseline number. Example: '2,3,4,5'
    uvrange                   UV-distance range, with a unit. Example: '2.0-3000.0 m'
    time                      Time range, as MJDs or date strings. Example: 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations. Example: 'RR,LL,RL,LR,XX,YY,XY,YX'
    intent                    Scan intent. Example: '*CAL*, *BAND*'
    observation               Observation Id. Example: '2~4'
    antint_ref_antenna        Antenna for which the fractions of channels flagged will be checked.
    minchanfrac               Minimum fraction of flagged channels required for a baseline to be deemed as flagged. Default:
    verbose                   Print timestamps of flagged integrations to the log. Default:
    apply                     Parameter to flag or unflag the data. Default:

Example:

af.parseantintparameters(antint_ref_antenna='ea10', minchanfrac=0.45, verbose=True)

--------------------------------------------------------------------------------
") parseantintparameters;

%feature("docstring", "

Summary:
    Parse data selection parameters and specific parameters for the extend mode. Data selection follows the MS Selection syntax. 

Input Parameters:
    field                     Field indices or source names. Example: '2,3C48'
    spw                       Spectral Window Indices or names. Example: '1,2'
    array                     Array Indices or names. Example: 'VLAA'
    feed                      Feed index or name. Example: '1,2' (not supported yet)
    scan                      Scan number. Example: '1,2,3'
    antenna                   Baseline number. Example: '2,3,4,5'
    uvrange                   UV-distance range, with a unit. Example: '2.0-3000.0 m'
    time                      Time range, as MJDs or date strings. Example: 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations. Example: 'RR,LL,RL,LR,XX,YY,XY,YX'
    intent                    Scan intent. Example: '*CAL*, *BAND*'
    observation               Observation Id. Example: '2~4'
    ntime                     Time-range to use for each chunk (in seconds or minutes). Default:
    combinescans              Accumulate data across scans.. Default:
    extendpols                If any correlation is flagged, flag all correlations. Default:
    growtime                  Flag all 'ntime' integrations if more than X% of the timerange is flagged (0-100). Default:
    growfreq                  Flag all selected channels if more than X% of the frequency range is flagged(0-100). Default:
    growaround                Flag data based on surrounding flags. Default:
    flagneartime              Flag one timestep before and after a flagged one. Default:
    flagnearfreq              Flag one channel before and after a flagged one. Default:
    apply                     Parameter to flag or unflag the data. Default:

Example:

af.parseextendparameters(extendpols=True)

--------------------------------------------------------------------------------
") parseextendparameters;

%feature("docstring", "

Summary:
    Parse data selection parameters and specific parameters for the summary mode. Data selection follows the MS Selection syntax. 

Input Parameters:
    field                     Field indices or source names. Example: '2,3C48'
    spw                       Spectral Window Indices or names. Example: '1,2'
    array                     Array Indices or names. Example: 'VLAA'
    feed                      Feed index or name. Example: '1,2' (not supported yet)
    scan                      Scan number. Example: '1,2,3'
    antenna                   Baseline number. Example: '2,3,4,5'
    uvrange                   UV-distance range, with a unit. Example: '2.0-3000.0 m'
    time                      Time range, as MJDs or date strings. Example: 'xx.x.x.x.x~yy.y.y.y.y'
    correlation               Correlations/polarizations. Example: 'RR,LL,RL,LR,XX,YY,XY,YX'
    intent                    Scan intent. Example: '*CAL*, *BAND*'
    observation               Observation Id. Example: '2~4'
    spwchan                   List the number of flags per spw and per channel. Default:
    spwcorr                   List the number of flags per spw and per correlation. Default:
    basecnt                   List the number of flags per baseline. Default:
    fieldcnt                  List the number of flags per field. Default:
    name                      Name of this summary report. Default: summary

Example:

af.parsesummaryparameters(spwchan=True, basecnt=True)

--------------------------------------------------------------------------------
") parsesummaryparameters;

%exception {
   try {
      $action
      } catch (const casacore::AipsError &ae) {
         PyErr_SetString(PyExc_RuntimeError, ae.what());
	 //PyErr_Print();
         return NULL;
      }
}
%include "agentflagger_cmpt.h"

%{
#include <exception>
#include <agentflagger_cmpt.h>
%}
