
void installLogFilter();
casacore::LogMessage::Priority getLogLevel(const std::string &level);

casacore::LogSinkInterface *thelogsink;
casacore::LogOrigin *itsorigin;
std::string taskname;
std::string processor_name;
std::string logname;
casacore::LogMessage::Priority logLevel;
std::vector<std::string> filterMsgList;
bool globalsink;
