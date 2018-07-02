import os
import fnmatch
import subprocess
import tarfile
from casac import casac
import datetime
import time
import urllib2
import __casac__
import TelemetryLogMonitor

class telemetry:

    def __init__(self, casa):

        self.setCasaVersion()
        self.setHostId()
        self.logdir = casa['dirs']['rc']
        self.logpattern = 'casastats-' + self.casaver + '-' + self.hostid + '*.log'
        self.sendlogpattern = 'casastats-*'+ self.hostid + '*.log'
        self.stampfile = self.logdir + "/telemetry.stamp"
        self.casa = casa

        logfiles = []

        for file in os.listdir(self.logdir):
            if fnmatch.fnmatch(file, self.logpattern):
                 #print "Matched: " + file
                 logfiles.append(file)

        logfiles.sort(reverse=True)
        # Size of the existing (non-active) logfiles
        inactiveTLogSize = 0

        if (logfiles and logfiles[0] != None):
            print "Found an existing telemetry logfile: " + casa['dirs']['rc'] + "/" + logfiles[0]
            casa['files']['telemetry-logfile'] = casa['dirs']['rc'] + "/" + logfiles[0]
            for i in range(1, len(logfiles)):
                inactiveTLogSize = inactiveTLogSize + os.path.getsize(casa['dirs']['rc'] + "/" + logfiles[i])/1024
                #print "Inactive log size: " + str(inactiveTLogSize)
        else :
             print "Creating a new telemetry file"
             self.setNewTelemetryFile()

        # Setup Telemetry log size monitoring
        casa_util = __casac__.utils.utils()

        # Size limit for the telemetry logs
        tLogSizeLimit = 10000
        # File size check interval
        tLogSizeInterval = 60
        try:
            tLogSizeLimit = int(casa_util.getrc("TelemetryLogLimit"))
            tLogSizeInterval = int(casa_util.getrc("TelemetryLogSizeInterval"))
        except:
            pass
        # Subtract the inactive log sizes from the total log file size limit
        tLogSizeLimit = tLogSizeLimit - inactiveTLogSize
        if (tLogSizeLimit <= 0):
            print "Telemetry log size limit exceeded. Disabling telemetry."
            casa['state']['telemetry-enabled'] = False
        else :
            tLogMonitor = TelemetryLogMonitor.TelemetryLogMonitor()
            tLogMonitor.start(casa['files']['telemetry-logfile'],tLogSizeLimit, tLogSizeInterval, casa)
            print "Telemetry initialized."

    def setNewTelemetryFile(self):
        self.casa['files']['telemetry-logfile'] = self.casa['dirs']['rc'] + '/casastats-' + self.casaver +'-'  + self.hostid + "-" + time.strftime("%Y%m%d-%H%M%S", time.gmtime()) + '.log'

    def setCasaVersion(self):
        myUtils = casac.utils()
        ver = myUtils.version()
        self.casaver = str(ver[0])+ str(ver[1]) + str(ver[2])+ "-" + str(ver[3])

    def setHostId(self):
        telemetryhelper = casac.telemetryhelper()
        self.hostid = telemetryhelper.getUniqueId()

    def setCasaLog(self, logger):
        self.logger = logger

    def submitStatistics(self):
        if (self.casa['state']['telemetry-enabled'] == True):
            self.logger.post("Checking telemetry submission interval")
            self.createStampFile()
            if (self.isSubmitInterval()):
                self.send('https://casa.nrao.edu/cgi-bin/crash-report.pl')
                self.refreshStampFile()
                self.setNewTelemetryFile()


    def isSubmitInterval(self):
        currentTime = time.time()
        lastUpdateTime = time.time()
        if (os.path.isfile(self.stampfile)):
            lastUpdateTime = os.path.getmtime(self.stampfile)

        # Check update checkSubmitInterval
        interval = 604800
        utils = casac.utils()
        if (utils.getrc("TelemetrySubmitInterval") != 'Unknown value'):
            interval = float(utils.getrc("TelemetrySubmitInterval"))
        if ((currentTime - lastUpdateTime)> interval):
            self.logger.post("Telemetry submit interval reached, submitting telemetry data.")
            return True
        else:
            self.logger.post("Telemetry submit interval not reached. Not submitting data.")
            #print "lastUpdateTime" +str(lastUpdateTime)
            #print "currentTime" +str(currentTime)
            self.logger.post("Next telemetry data submission in: " + str(datetime.timedelta(  \
                    seconds=(interval-(currentTime-lastUpdateTime)))))
            return False

    def createStampFile(self):
        #print "Checking for stampfile " + self.stampfile
        if not os.path.isfile(self.stampfile):
            self.logger.post("Creating a new telemetry time stamp file." + self.stampfile)
            open(self.stampfile, 'a').close()

    def refreshStampFile(self):
        os.utime(self.stampfile, None)

    def send(self, telemetry_url):

        telemetryhelper = casac.telemetryhelper()
        logfiles = []

        # Test if internet connection is available.
        try:
            urllib2.urlopen('https://casa.nrao.edu/', timeout=2)
        except urllib2.URLError as err:
            return

        # Find logfiles
        for file in os.listdir(self.logdir):
            if fnmatch.fnmatch(file, self.sendlogpattern):
                #print "Matched: " + file
                logfiles.append(file)

        if (len(logfiles) > 0):
            #Tar logfiles
            current_date = datetime.datetime.today().strftime('%Y%m%d%H%M%S')
            tarfileid = self.logdir + "/telemetry-" \
                        + telemetryhelper.getUniqueId() + "-" \
                        + current_date + ".tar.gz"
            tar = tarfile.open(tarfileid, "w:gz")

            for logfile in logfiles:
                tar.add(self.logdir + "/" + logfile,
                        arcname='telemetry/'+logfile)
            tar.close()

            file_param = 'file=@' + tarfileid #+ '\"'
            # Submit tarfile
            #print ['curl', '-F', file_param , telemetry_url]
            proc = subprocess.Popen(['curl', '-F', file_param , telemetry_url],stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            cmd_out, cmd_err = proc.communicate()
            if cmd_out != None:
                self.logger.post(cmd_out, 'DEBUG1')
            if cmd_err != None:
                self.logger.post(cmd_err, 'DEBUG1')

            # Remove files
            for logfile in logfiles:
                os.remove(self.logdir + "/" + logfile)
                #print "Removed " + self.logdir + "/" + logfile
            os.remove(tarfileid)
            self.logger.post("Removed" + tarfileid)
        else:
            self.logger.post("No telemetry files to submit.")
