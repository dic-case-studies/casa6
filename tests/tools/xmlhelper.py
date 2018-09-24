import os
from casatools import quanta

qa_ = quanta( )

def readXML(sdmfile, mytbuff):
    '''
#   readflagxml: reads Antenna.xml and Flag.xml SDM tables and parses
#                into returned dictionary as flag command strings
#      sdmfile (string)  path to SDM containing Antenna.xml and Flag.xml
#      mytbuff (float)   time interval (start and end) padding (seconds)
#
#   Usage: myflags = readflagxml(sdmfile,tbuff)
#
#   Dictionary structure:
#   fid : 'id' (string)
#         'mode' (string)         flag mode ('online')
#         'antenna' (string)
#         'timerange' (string)
#         'reason' (string)
#         'time' (float)          in mjd seconds
#         'interval' (float)      in mjd seconds
#         'command' (string)          string (for COMMAND col in FLAG_CMD)
#         'type' (string)         'FLAG' / 'UNFLAG'
#         'applied' (bool)        set to True here on read-in
#         'level' (int)           set to 0 here on read-in
#         'severity' (int)        set to 0 here on read-in
#
#   Updated STM 2011-11-02 handle new SDM Flag.xml format from ALMA
#   Updated STM 2012-02-14 handle spectral window indices, names, IDs
#   Updated STM 2012-02-21 handle polarization types
#
#   Mode to use for spectral window selection in commands:
#   spwmode =  0 none (flag all spw)
#   spwmode =  1 use name
#   spwmode = -1 use index (counting rows in SpectralWindow.xml)
#
#   Mode to use for polarization selection in commands:
#   polmode =  0 none (flag all pols/corrs)
#   polmode =  1 use polarization type
#
#   CURRENT DEFAULT: Use spw names, flag pols
    '''

    spwmode = 1
    polmode = 1

#
    try:
        from xml.dom import minidom
    except ImportError as e:
        print('failed to load xml.dom.minidom:\n%s' % e)
        exit(1)

    if type(mytbuff) != float:
        print('Found incorrect type for tbuff')
        exit(1)
#        mytbuff = 1.0

    # make sure Flag.xml and Antenna.xml are available (SpectralWindow.xml depends)
    flagexist = os.access(sdmfile + '/Flag.xml', os.F_OK)
    antexist = os.access(sdmfile + '/Antenna.xml', os.F_OK)
    spwexist = os.access(sdmfile + '/SpectralWindow.xml', os.F_OK)
    if not flagexist:
        print('Cannot open ' + sdmfile + '/Flag.xml')
        exit(1)
    if not antexist:
        print('Cannot open ' + sdmfile + '/Antenna.xml')
        exit(1)
    if not spwexist:
        print('Cannot open ' + sdmfile + '/SpectralWindow.xml')

    # construct look-up dictionary of name vs. id from Antenna.xml
    xmlants = minidom.parse(sdmfile + '/Antenna.xml')
    antdict = {}
    rowlist = xmlants.getElementsByTagName('row')
    for rownode in rowlist:
        rowname = rownode.getElementsByTagName('name')
        ant = str(rowname[0].childNodes[0].nodeValue)
        rowid = rownode.getElementsByTagName('antennaId')
        # CAS-4532: remove spaces between content and tags
        antid = str(rowid[0].childNodes[0].nodeValue).strip()
        antdict[antid] = ant
    print('Found ' + str(rowlist.length) + ' antennas in Antenna.xml')

    # construct look-up dictionary of name vs. id from SpectralWindow.xml
    if spwexist:
        xmlspws = minidom.parse(sdmfile + '/SpectralWindow.xml')
        spwdict = {}
        rowlist = xmlspws.getElementsByTagName('row')
        ispw = 0
        for rownode in rowlist:
            rowid = rownode.getElementsByTagName('spectralWindowId')
            # CAS-4532: remove spaces between content and tags
            spwid = str(rowid[0].childNodes[0].nodeValue).strip()
            spwdict[spwid] = {}
            spwdict[spwid]['index'] = ispw
            # SMC: 6/3/2012 ALMA SDM does not have name
            if rownode.getElementsByTagName('name'):
                rowname = rownode.getElementsByTagName('name')
                spw = str(rowname[0].childNodes[0].nodeValue)
                spwdict[spwid]['name'] = spw
            else:
                spwmode = -1
                
#            rowid = rownode.getElementsByTagName('spectralWindowId')
#            spwid = str(rowid[0].childNodes[0].nodeValue)
#            spwdict[spwid] = {}
#            spwdict[spwid]['index'] = ispw
            ispw += 1
        print('Found ' + str(rowlist.length) + ' spw in SpectralWindow.xml')

    # report chosen spw and pol modes
    if spwmode > 0:
        print('Will construct spw flags using names')
    elif spwmode < 0:
        print('Will construct spw flags using table indices')
    else:
        print('')
    #
    if polmode == 0:
        print('Will not set polarization dependent flags (flag all corrs)')
    else:
        print('Will construct polarization flags using polarizationType')

    # now read Flag.xml into dictionary row by row
    xmlflags = minidom.parse(sdmfile + '/Flag.xml')
    flagdict = {}
    rowlist = xmlflags.getElementsByTagName('row')
    nrows = rowlist.length
    newsdm = -1
    newspw = -1
    newpol = -1
    for fid in range(nrows):
        rownode = rowlist[fid]
        rowfid = rownode.getElementsByTagName('flagId')
        fidstr = str(rowfid[0].childNodes[0].nodeValue)
        flagdict[fid] = {}
        flagdict[fid]['id'] = fidstr
        rowid = rownode.getElementsByTagName('antennaId')
        antid = rowid[0].childNodes[0].nodeValue
        # check if there is a numAntenna specified (new format)
        rownant = rownode.getElementsByTagName('numAntenna')
        antname = ''
        if rownant.__len__() > 0:
            xid = antid.split()
            nant = int(rownant[0].childNodes[0].nodeValue)
            if newsdm < 0:
                print('Found numAntenna=' + str(nant) + ' must be a new style SDM')
            newsdm = 1
            # CAS-4698. Flag auto-correlations when there is
            # only one antenna
            if nant == 1:
                aid = xid[2]
                ana = antdict[aid]
                if antname == '':
                    antname = ana+'&&*'
                else:
                    antname += ',' + ana
                    
            elif nant > 1:
                for ia in range(nant):
                    aid = xid[2 + ia]
                    ana = antdict[aid]
                    if antname == '':
                        antname = ana
                    else:
                        antname += ',' + ana
            else:
            # numAntenna = 0 means flag all antennas
                antname = ''
        else:
            if newsdm < 0:
                print('No numAntenna entry found, must be an old style SDM')
            newsdm = 0
            nant = 1
            aid = antid
            ana = antdict[aid]
            antname = ana
        # start and end times in mjd ns
        rowstart = rownode.getElementsByTagName('startTime')
        start = int(rowstart[0].childNodes[0].nodeValue)
        startmjds = float(start) * 1.0E-9 - mytbuff
        t = qa_.quantity(startmjds, 's')
        starttime = qa_.time(t, form='ymd', prec=9)[0]
        rowend = rownode.getElementsByTagName('endTime')
        end = int(rowend[0].childNodes[0].nodeValue)
        endmjds = float(end) * 1.0E-9 + mytbuff
        t = qa_.quantity(endmjds, 's')
        endtime = qa_.time(t, form='ymd', prec=9)[0]
    # time and interval for FLAG_CMD use
        times = 0.5 * (startmjds + endmjds)
        intervs = endmjds - startmjds
        flagdict[fid]['time'] = times
        flagdict[fid]['interval'] = intervs
        # reasons
        rowreason = rownode.getElementsByTagName('reason')
        reas = str(rowreason[0].childNodes[0].nodeValue)
        # Replace any white space with underscores
        reason = reas.replace(' ','_')
    # NEW SDM ADDITIONS 2011-11-01
        rownspw = rownode.getElementsByTagName('numSpectralWindow')
        spwstring = ''
        if spwmode != 0 and rownspw.__len__() > 0:
            nspw = int(rownspw[0].childNodes[0].nodeValue)
        # has a new-style spw specification
            if newspw < 0:
                if not spwexist:
                    print('Cannot open ' + sdmfile + '/SpectralWindow.xml')
                    exit(1)
                print('Found SpectralWindow=' + str(nspw) + ' must be a new style SDM')
            newspw = 1
            if nspw > 0:
                rowspwid = \
                    rownode.getElementsByTagName('spectralWindowId')
                spwids = rowspwid[0].childNodes[0].nodeValue
                xspw = spwids.split()
                for isp in range(nspw):
                    spid = str(xspw[2 + isp])
                    if spwmode > 0:
                        spstr = spwdict[spid]['name']
                    else:
                        spstr = str(spwdict[spid]['index'])
                    if spwstring == '':
                        spwstring = spstr
                    else:
                        spwstring += ',' + spstr
        polstring = ''
        rownpol = rownode.getElementsByTagName('numPolarizationType')
        if polmode != 0 and rownpol.__len__() > 0:
            npol = int(rownpol[0].childNodes[0].nodeValue)
        # has a new-style pol specification
            if newpol < 0:
                print('Found numPolarizationType=' + str(npol) + ' must be a new style SDM')
            newpol = 1
            if npol > 0:
                rowpolid = \
                    rownode.getElementsByTagName('polarizationType')
                polids = rowpolid[0].childNodes[0].nodeValue
                xpol = polids.split()
                for ipol in range(npol):
                    polid = str(xpol[2 + ipol])
                    if polstring == '':
                        polstring = polid
                    else:
                        polstring += ',' + polid
    #
        # Construct antenna name and timerange and reason strings
        flagdict[fid]['antenna'] = antname
        timestr = starttime + '~' + endtime
        flagdict[fid]['timerange'] = timestr
        flagdict[fid]['reason'] = reason
        # Construct command strings (per input flag)
        cmd = "antenna='" + antname + "' timerange='" + timestr + "'"
        if spwstring != '':
            cmd += " spw='" + spwstring + "'"
            flagdict[fid]['spw'] = spwstring
#        if polstring != '':
#            cmd += " poln='" + polstring + "'"
#            flagdict[fid]['poln'] = polstring
        if polstring != '':
            # Write the poln translation in correlation
            if polstring.count('R')>0:
                if polstring.count('L')>0:
                    corr = 'RR,RL,LR,LL'
                else:
                    corr = 'RR,RL,LR'
            elif polstring.count('L')>0:
                corr = 'LL,LR,RL'
            elif polstring.count('X')>0:
                if polstring.count('Y')>0:
                    corr = 'XX,XY,YX,YY'
                else:
                    corr = 'XX,XY,YX'
            elif polstring.count('Y')>0:
                corr = 'YY,YX,XY'

            cmd += " correlation='" + corr + "'"
#            flagdict[fid]['poln'] = polstring
        flagdict[fid]['command'] = cmd
    #
        flagdict[fid]['type'] = 'FLAG'
        flagdict[fid]['applied'] = False
        flagdict[fid]['level'] = 0
        flagdict[fid]['severity'] = 0
        flagdict[fid]['mode'] = 'online'

    flags = {}
    if rowlist.length > 0:
        flags = flagdict
        print('Found ' + str(rowlist.length) + ' flags in Flag.xml')
    else:
        print('No valid flags found in Flag.xml')

    # return the dictionary for later use
    return flags


