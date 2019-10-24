// The Private declarations of the PlotMS component.

// Static //

// Constants for how long to wait for the casapy application to launch, in
// microseconds.
// <group>
static const unsigned int LAUNCH_WAIT_INTERVAL_US;
static const unsigned int LAUNCH_TOTAL_WAIT_US;
// </group>

// Implementation of casapy_watcher for use with plotms component.
class plotms_watcher : public CasapyWatcher {
public:
    // Constructor which takes parent.
    plotms_watcher(plotms *pl) : p(pl) { CasapyWatcher::registerWatcher(this); }
    
    // Destructor.
    ~plotms_watcher() { }
    
    // Overrides casapy_watcher::logChanged().
    // <group>
    void logChanged(const casacore::String& sinkLocation) {
        p->setLogFilename(sinkLocation); }
    void logChanged(casacore::LogMessage::Priority filterPriority) {
        p->setLogFilter(casacore::LogMessage::toString(filterPriority).c_str()); }
    // </group>
    
    // Overrides casapy_watcher::casapyClosing().
    void casapyClosing() { p->closeApp(); }
    
private:
    // Plotms parent.
    plotms *p;

};

class plotms_app  {
public:
  plotms_app( ) { }
  const casacore::String &dbusName( ) const { return itsDBusName_; }
  casacore::String &dbusName( ) { return itsDBusName_; }
  // DBus name of the plotms application we're communicating with.
  const QString &getName( ) const { return casa::PlotMSDBusApp::name( ); }
 private:
  casacore::String itsDBusName_;
};

// Non-Static //

// must forward declare & use a pointer because build system
// does not allow extra include in plotms_cmpt.h...
plotms_app app;

// Casapy watcher.
plotms_watcher itsWatcher_;

// PID of the plotms_app
pid_t app_pid;

// Log parameters that are set before the application is launched.
// <group>
casacore::String itsLogFilename_, itsLogFilter_;
// </group>


// Returns true if the DISPLAY environment variable is set, false otherwise.
// If not set, an error message is printed.
bool displaySet();

// Launches the dbus plotms application IF it is not already launched.
// Returns:
//   true = already launched, or launched in launchCasaplotms
//   false = launch failed or plotms (task) has crashed 
bool launchApp();
bool launchCasaplotms();

// Closes the launched dbus plotms application if needed.
void closeApp();

// Helper method for calling an async method.
void callAsync(const casacore::String& methodName);

// Helper method for realizing synchronous behavior
void waitUntilIdle();

// Helper method for setting the casacore::MS selection.
void setPlotMSSelection_(const casa::PlotMSSelection& selection,
        const bool updateImmediately, const int plotIndex);

// Helper method for setting the casacore::MS averaging.
void setPlotMSAveraging_(const casa::PlotMSAveraging& averaging,
        const bool updateImmediately, const int plotIndex);

// Helper method for setting the casacore::MS transformations.
void setPlotMSTransformations_(const casa::PlotMSTransformations& trans,
        const bool updateImmediately, const int plotIndex);

// Helper method for setting the casacore::MS calibration.
void setPlotMSCalibration_(const casa::PlotMSCalibration& calib,
        const bool updateImmediately, const int plotIndex);

// Helper method for setting the iteration parameters
void setPlotMSIterate_(const casa::PlotMSIterParam& iter,
        const bool updateImmediately, const int plotIndex);

// Helper method for setting the flag extension.
void setFlagging_(const casa::PlotMSFlagging& flagging);

bool showGui;
bool asyncCall;
bool isTask;
