from __future__ import absolute_import
from __future__ import print_function

import datetime
import re

from casatasks.private.casa_transition import is_CASA6
if is_CASA6:
    from urllib.request import urlopen
    from urllib.error import URLError
    from casatools import table, quanta
    from casatasks import casalog
    from casatools.platform import bytes2str

    _tb = table( )
    _qa = quanta( )
else:
    from urllib2 import urlopen
    from urllib2 import URLError

    from taskinit import *

    (_tb,)=gentools(['tb'])

    # to make the following code the same as the CASA6 version
    _qa = qa

######################################################################
def correct_ant_posns_evla (vis_name, print_offsets=False):
    '''
    Given an input visibility MS name (vis_name), find the antenna
    position offsets that should be applied.  This application should
    be via the gencal task, using caltype='antpos'.

    If the print_offsets parameter is True, will print out each of
    the found offsets (or indicate that none were found), otherwise
    runs silently.

    A list is returned where the first element is the returned error
    code, the second element is a string of the antennas, and the 
    third element is a list of antenna Bx,By,Bz offsets.  An example 
    return list might look like:
    [ 0, 'ea01,ea19', [0.0184, -0.0065, 0.005, 0.0365, -0.0435, 0.0543] ]

    Usage examples:

       CASA <1>: antenna_offsets = correct_ant_posns('test.ms')
       CASA <2>: if (antenna_offsets[0] == 0):
       CASA <3>:     gencal(vis='test.ms', caltable='cal.G', \
                     caltype='antpos', antenna=antenna_offsets[1], \
                     parameter=antenna_offsets[2])

    This function does NOT work for VLA datasets, only EVLA.  If an
    attempt is made to use the function for VLA data (prior to 2010),
    an error code of 1 is returned.

    The offsets are retrieved over the internet.  A description and the
    ability to manually examine and retrieve offsets is at:
    http://www.vla.nrao.edu/astro/archive/baselines/
    If the attempt to establish the internet connection fails, an error
    code of 2 is returned.

    Uses the same algorithm that the AIPS task VLANT does.


    BJB
    NRAO
    Spring 2020 (fixed version)
    '''

    MONTHS = [ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' ]
    URL_BASE = 'http://www.vla.nrao.edu/cgi-bin/evlais_blines.cgi?Year='

#
# get start date+time of observation
#
    observation = _tb.open(vis_name+'/OBSERVATION')
    time_range = _tb.getcol('TIME_RANGE')
    _tb.close()
    MJD_start_time = time_range[0][0] / 86400
    q1 = _qa.quantity(time_range[0][0],'s')
    date_time = _qa.time(q1,form='ymd')
# date_time looks like: '2011/08/10/06:56:49'
    [obs_year,obs_month,obs_day,obs_time_string] = date_time[0].split('/')
    if int(obs_year) < 2010:
        if print_offsets:
            print('Does not work for VLA observations')
        return [1, '', []]
    [obs_hour,obs_minute,obs_second] = obs_time_string.split(':')
    obs_time = 10000*int(obs_year) + 100*int(obs_month) + int(obs_day) + \
               int(obs_hour)/24.0 + int(obs_minute)/1440.0 + \
               int(obs_second)/86400.0

#
# get antenna to station mappings
#
    observation = _tb.open(vis_name+'/ANTENNA')
    ant_names = _tb.getcol('NAME')
    ant_stations = _tb.getcol('STATION')
    _tb.close()
    ant_num_stas = []
    for ii in range(len(ant_names)):
        ant_num_stas.append([int(ant_names[ii][2:]), ant_names[ii], \
                            ant_stations[ii], 0.0, 0.0, 0.0, False])

    correction_lines = []
    current_year = datetime.datetime.now().year
# first, see if the internet connection is possible
    try:
        response = urlopen(URL_BASE + '2010')
    except URLError as err:
        if print_offsets:
            print('No internet connection to antenna position correction URL ', err.reason)
        return [2, '', []]
    response.close()
    for year in range(2010,current_year+1):
        response = urlopen(URL_BASE + str(year))
        if is_CASA6:
            html = bytes2str(response.read())
        else:
            html = response.read()
        response.close()
        html_lines = html.split('\n')

        for correction_line in html_lines:
            if len(correction_line) and correction_line[0] != '<' and correction_line[0] != ';':
                for month in MONTHS:
                    if month in correction_line:
                        correction_lines.append(str(year)+' '+correction_line)
                        break

    corrections_list = []
    for correction_line in correction_lines:
        correction_line_fields = correction_line.split()
        if len(correction_line_fields) > 9:
            [c_year, moved_date, obs_date, put_date, put_time_str, ant, pad, Bx, By, Bz] = correction_line_fields
            s_moved = moved_date[:3]
            i_month = 1
            for month in MONTHS:
                if moved_date.find(month) >= 0:
                    break
                i_month = i_month + 1
            moved_time = 10000 * int(c_year) + 100 * i_month + \
                         int(moved_date[3:])
        else:
            [c_year, obs_date, put_date, put_time_str, ant, pad, Bx, By, Bz] = correction_line_fields
            moved_date = '     '
            moved_time = 0
        s_obs = obs_date[:3]
        i_month = 1
        for month in MONTHS:
            if s_obs.find(month) >= 0:
                break
            i_month = i_month + 1
        obs_time_2 = 10000 * int(c_year) + 100 * i_month + int(obs_date[3:])
        s_put = put_date[:3]
        i_month = 1
        for month in MONTHS:
            if s_put.find(month) >= 0:
                break
            i_month = i_month + 1
        put_time = 10000 * int(c_year) + 100 * i_month + int(put_date[3:])
        [put_hr, put_min] = put_time_str.split(':')
        put_time += (int(put_hr)/24.0 + int(put_min)/1440.0)
        corrections_list.append([c_year, moved_date, moved_time, obs_date, obs_time_2, put_date, put_time, int(ant), pad, float(Bx), float(By), float(Bz)])

    for correction_list in corrections_list:
        [c_year, moved_date, moved_time, obs_date, obs_time_2, put_date, put_time, ant, pad, Bx, By, Bz] = correction_list
        ant_ind = -1
        for ii in range(len(ant_num_stas)):
            if ant_num_stas[ii][0] == ant:
                ant_ind = ii
                break
# make sure the antenna in this correction is in the observation,
# and is not done
        if ant_ind != -1 and not ant_num_stas[ant_ind][6]:
            ant_num_sta = ant_num_stas[ant_ind]
            if moved_time:
# the antenna moved
                if moved_time > obs_time:
# we are done considering this antenna
                    ant_num_stas[ant_ind][6] = True
                else:
# otherwise, it moved, so the offsets should be reset
                    ant_num_stas[ant_ind][3] = 0.0
                    ant_num_stas[ant_ind][4] = 0.0
                    ant_num_stas[ant_ind][5] = 0.0
            if put_time > obs_time and not ant_num_stas[ant_ind][6] and pad == ant_num_stas[ant_ind][2]:
# it's the right antenna/pad; add the offsets to those already accumulated
                ant_num_stas[ant_ind][3] += Bx
                ant_num_stas[ant_ind][4] += By
                ant_num_stas[ant_ind][5] += Bz

    ants = []
    parms = []
    for ant_num_sta in ant_num_stas:
        if ant_num_sta[3] != 0.0 or ant_num_sta[4] != 0.0 or ant_num_sta[3] != 0.0:
            if print_offsets:
                print("Offsets for antenna %4s on pad %3s: %8.5f  %8.5f  %8.5f" % \
                      (ant_num_sta[1], ant_num_sta[2], ant_num_sta[3], ant_num_sta[4], ant_num_sta[5]))
            else:
                casalog.post("offsets for antenna %4s : %8.5f  %8.5f  %8.5f" % \
                      (ant_num_sta[1], ant_num_sta[3], ant_num_sta[4], ant_num_sta[5]))
            ants.append(ant_num_sta[1])
            parms.append(ant_num_sta[3])
            parms.append(ant_num_sta[4])
            parms.append(ant_num_sta[5])
    if len(parms) == 0 and print_offsets:
        print("No offsets found for this MS")
    ant_string = ','.join(["%s" % ii for ii in ants])
    return [ 0, ant_string, parms ]

######################################################################
