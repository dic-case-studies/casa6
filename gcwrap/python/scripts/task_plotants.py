import os
import numpy as np
import pylab as pl
from taskinit import gentools, qatool

def plotants(vis=None, antindex=None, logpos=None, figfile=None):
    """Plot the antenna distribution in the local reference frame:

    The location of the antennas in the MS will be plotted with
    X-toward local east; Y-toward local north.  The name of each 
    antenna is shown next to its respective location.

    Keyword arguments:
    vis -- Name of input visibility file.
            default: none. example: vis='ngc5921.ms'

    antindex -- Label antennas with name and antenna ID
            default: False. example: antindex=True

    logpos -- Produce a logarithmic position plot.
            default: False. example: logpos=True

    figfile -- Save the plotted figure in this file.
            default: ''. example: figfile='myFigure.png'

    You can zoom in by pressing the magnifier button (bottom,
    third from right) and making a rectangular region with
    the mouse.  Press the home button (left most button) to
    remove zoom.

    A hard-copy of this plot can be obtained by pressing the
    button on the right at the bottom of the display. A file 
    dialog will allow you to choose the directory, filename,
    and format of the export.
    """

    try:
        telescope, antNames, antXs, antYs = getAntennaInfo(vis)
        pl.clf()
        if logpos:
            plotAntennasLog(antXs, antYs, antNames, antindex, telescope)
        else:
            plotAntennas(antXs, antYs, antNames, antindex, telescope)
        pl.title("Antenna positions for " + os.path.basename(vis))
        if figfile:
            pl.savefig(figfile)
    except Exception, instance:
          print '*** Error ***',instance

def getAntennaInfo(msname):

    msmd, me, tb = gentools(['msmd', 'me', 'tb'])
    qa = qatool()

    msmd.open(msname)
    telescope = msmd.observatorynames()[0]
    arrayPos = msmd.observatoryposition()
    """
    allAntennas = np.array([], dtype=int)
    for scan in msmd.scannumbers():
        allAntennas = np.append(allAntennas, msmd.antennasforscan(scan))
    antIdsUsed = set(allAntennas)
    """
    msmd.close()

    arrayWgs84 = me.measure(arrayPos, 'WGS84')
    arrayLon, arrayLat, arrayAlt = [arrayWgs84[i]['value'] 
        for i in ['m0','m1','m2']]

    # Open the ANTENNA subtable to get the names of the antennas in this MS and
    # their positions.  Note that the entries in the ANTENNA subtable are pretty
    # much in random order, so antNames translates between their index and name
    # (e.g., index 11 = STD155).  We'll need these indices for later, since the
    # main data table refers to the antennas by their indices, not names.
    tb.open(msname)
    ants1 = tb.getcol('ANTENNA1')
    ants2 = tb.getcol('ANTENNA2')
    tb.close()
    allantennas = np.append(ants1, ants2)
    antIdsUsed = set(allantennas)
    print "Number of points being plotted:", len(antIdsUsed)

    anttabname = msname + '/ANTENNA'
    tb.open(anttabname)
    if telescope == 'VLBA':
        antNames = np.array(tb.getcol('STATION'))
    else:
        antNames = np.array(tb.getcol('NAME'))
    antNames = [antNames[i] for i in antIdsUsed]

    antPositions = np.array([me.position('ITRF', qa.quantity(x, 'm'),
        qa.quantity(y, 'm'), qa.quantity(z, 'm'))
        for (x, y, z) in tb.getcol('POSITION').transpose()])
    tb.close()
    antPositions = [antPositions[i] for i in antIdsUsed]
    
    # Get the names, indices, and lat/lon/alt coords of our "good" antennas.
    antWgs84s = np.array([me.measure(pos, 'WGS84') for pos in antPositions])
    antLons, antLats, antAlts = [np.array( [pos[i]['value'] 
        for pos in antWgs84s]) for i in ['m0','m1','m2']]
    
    # Convert from lat, lon, alt to X, Y, Z (unless VLBA)
    # where X is east, Y is north, Z is up,
    # and 0, 0, 0 is the center of the LWA1.  
    # Note: this conversion is NOT exact, since it doesn't take into account
    # Earth's ellipticity!  But it's close enough.
    radE = 6370000.
    antXs = (antLons - arrayLon) * radE * np.cos(arrayLat)
    antYs = (antLats - arrayLat) * radE
    antZs = antAlts - arrayAlt
    
    return telescope, antNames, antXs, antYs

def plotAntennasLog(antXs, antYs, antNames, antIndex, telescope):
    fig = pl.figure(1)
    # Set up subplot.
    ax = fig.add_subplot(1, 1, 1, polar=True, projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    # Do not show azimuth labels.
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # Do not show grid.
    ax.grid(False)

    # code from pipeline summary.py 
    # PlotAntsChart draw_polarlog_ant_map_in_subplot
    if telescope in ('VLA', 'EVLA'):
        # For (E)VLA, set a fixed local center position that has been
        # tuned to work well for its array configurations (CAS-7479).
        xcenter, ycenter = -32, 0
        rmin_min, rmin_max = 12.5, 350
    else:
        # For non-(E)VLA, take the median of antenna offsets as the
        # center for the plot.
        xcenter = np.median(antXs)
        ycenter = np.median(antYs)
        rmin_min, rmin_max = 3, 350

    # Derive radial offset w.r.t. center position.
    r = ((antXs-xcenter)**2 + (antYs-ycenter)**2)**0.5
    # Set rmin, clamp between a min and max value, ignore station
    # at r=0 if one is there.
    rmin = min(rmin_max, max(rmin_min, 0.8*np.min(r[r > 0])))
    # Update r to move any points below rmin to r=rmin.
    r[r <= rmin] = rmin
    rmin = np.log(rmin)
    # Set rmax.
    rmax = np.log(1.5*np.max(r))
    # Derive angle of offset w.r.t. center position.
    theta = np.arctan2(antXs-xcenter, antYs-ycenter)

    # Draw circles at specific distances from the center.
    angles = np.arange(0, 2.01*np.pi, 0.01*np.pi)
    show_circle = True
    circles = [30, 100, 300, 1000, 3000, 10000]
    if telescope == "VLBA":
        circles = [1e5, 3e5, 1e6, 3e6, 1e7]
    for cr in circles:

        # Only draw circles outside rmin.
        if cr > np.min(r) and show_circle:

            # Draw the circle.
            radius = np.ones(len(angles))*np.log(cr)
            ax.plot(angles, radius, 'k:')

            # Draw tick marks on the circle at 1 km intervals.
            inc = 0.1*10000/cr
            if telescope == "VLBA":
                inc = 0.1*1e8/cr
            if cr > 100:
                for angle in np.arange(inc/2., 2*np.pi+0.05, inc):
                    ax.plot([angle, angle],
                               [np.log(0.95*cr), np.log(1.05*cr)], 'k-')

            # Add text label to circle to denote distance from center.
            va = 'top'
            circle_label_angle = -20.0 * np.pi / 180.
            if cr >= 1000:
                if np.log(cr) < rmax:
                    ax.text(circle_label_angle, np.log(cr),
                               '%d km' % (cr/1000), size=8, va=va)
                    ax.text(circle_label_angle + np.pi, np.log(cr),
                               '%d km' % (cr / 1000), size=8, va=va)
            else:
                ax.text(circle_label_angle, np.log(cr), '%dm' % (cr), 
                    size=8, va=va)
                ax.text(circle_label_angle + np.pi, np.log(cr), 
                        '%dm' % (cr), size=8, va=va)

        # Find out if most recently drawn circle was outside all antennas, 
        # if so, no more circles will be drawn.
        if np.log(cr) > rmax:
            show_circle = False

    # plot points and antenna names/ids
    for i, antenna in enumerate(antNames):
        ax.plot(theta[i], np.log(r[i]), 'ko', ms=5, mfc='r')
        if antIndex:
            antenna += ' (' + str(i) + ')'
        ax.text(theta[i], np.log(r[i]), '   '+antenna, size=8, va='top')
    # Set minimum and maximum radius.
    ax.set_rmax(rmax)
    ax.set_rmin(rmin)

def plotAntennas(antXs, antYs, antNames, antIndex, telescope):
    fig = pl.figure(1)
    ax = fig.add_subplot(111)

    # use m or km units
    units = ' (m)'
    if np.median(antXs) > 1e6 or np.median(antYs) > 1e6:
        antXs /= 1e3
        antYs /= 1e3
        units = ' (km)'

    # plot points and antenna names/ids
    for i, (x, y, name) in enumerate(zip(antXs, antYs, antNames)):
        ax.plot(x, y, 'ro')
        if antIndex:
            name += ' (' + str(i) + ')'
        ax.text(x, y+2, '  ' + name, size=8, va='top')
        fig.show()

    pl.xlabel('X' + units)
    pl.ylabel('Y' + units)
    pl.margins(0.1, 0.1)
