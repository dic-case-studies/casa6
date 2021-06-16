import certifi
import collections
import csv
import datetime
import json
import os
import re
import ssl
import string
import urllib

import numpy as np

from casatools import ms as mstool
from casatools import msmetadata 
from casatools import quanta 
from casatools import table 

# get is_CASA6 and is_python3
if is_CASA6:
    from casatasks import casalog
else:
    from taskinit import casalog


### web api part

QueryStruct = collections.namedtuple('QueryStruct', ['param', 'subparam'])
ResponseStruct = collections.namedtuple('ResponseStruct', ['response', 'subparam'])


class ALMAJyPerKDatabaseAccessBase(object):
    BASE_URL = 'https://asa.alma.cl/science/jy-kelvins'
    ENDPOINT_TYPE = None

    @property
    def url(self):
        assert self.ENDPOINT_TYPE is not None, \
            '{} cannot be instantiated. Please use subclasses.'.format(self.__class__.__name__)

        s = '/'.join([self.BASE_URL, self.ENDPOINT_TYPE])
        if not s.endswith('/'):
            s += '/'
        return s

    def __init__(self, context=None):
        """
        ALMAJyPerKDatabaseAccessBase is a base class for accessing Jy/K
        DB to retrieve conversion factor for ALMA TP data.
        ALMAJyPerKDatabaseAccessBase is kind of a template class that
        only provides a standard workflow to get a list of conversion
        factors. Each access class must inherit this class and
        implement/override some methods according to the target API.
        Subclasses must implement properties and methods listed below:

            ENDPOINT_TYPE (property): Must be a string representing the API
            access (method): Receive a list of queries as a generator,
                             access the DB through the generator, and
                             return the formatted response. Return value
                             should be a dictionary with 'query', 'data',
                             and 'total' fields. The 'query' field holds
                             the query whil the 'data' field stores the
                             response. The 'total' fields is the number of
                             response. Each item of the 'data' field should
                             consist of single conversion factor ('Factor')
                             with the meta-data, 'MS', 'Antenna', 'Spwid',
                             'Polarization'.
            get_params (method): Receive a name of the MS and generate
                                 a dictionary containing a list of query
                                 parameters. Required parameters depend on
                                 the API.

        Keyword Arguments:
            context {Context} -- Pipeline Context object (default: {None})
        """
        self.context = context

    def _get_observing_band(self, ms):
        if self.context is None:
            return 'Unknown'

        spws = ms.get_spectral_windows(science_windows_only=True)
        bands = [spw.band for spw in spws]
        return np.unique(bands)

    def _generate_query(self, url, params):
        try:
            for p in params:
                # encode params
                encoded = urllib.parse.urlencode(p.param)

                # try opening url
                query = '?'.join([url, encoded])
                casalog.post('Accessing Jy/K DB: query is "{}"'.format(query))
                ssl_context = ssl.create_default_context(cafile=certifi.where())
                # set timeout to 3min (=180sec)
                response = urllib.request.urlopen(query, context=ssl_context, timeout=180)
                retval = json.load(response)
                if not retval['success']:
                    msg = 'Failed to get a Jy/K factor from DB: {}'.format(retval['error'])
                    casalog.post(msg)
                    raise RuntimeError(msg)
                yield ResponseStruct(response=retval, subparam=p.subparam)
        except urllib.error.HTTPError as e:
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: HTTPError(code={0}, Reason="{1}")\n'.format(e.code, e.reason)
            casalog.post(msg)
            raise e
        except urllib.error.URLError as e:
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: URLError(Reason="{0}")\n'.format(e.reason)
            casalog.post(msg)
            raise e

    def validate_is_alma_data(self, vis):
        msmd = msmetadata()
        msmd.open(vis)
        try:
            array_name = msmd.observatorynames()[0]
        except KeyError:
            LOG.error('{} is not registered to observatorynames'.format(vis))
        msmd.close()
        if array_name != 'ALMA':
            raise RuntimeError('{} is not ALMA data'.format(basename))

    def getJyPerK(self, vis):
        """
        getJyPerK returns list of Jy/K conversion factors with their
        meta data (MS name, antenna name, spwid, and pol string).

        Arguments:
            vis {str} -- Name of MS

        Returns:
            [list] -- List of Jy/K conversion factors with meta data
        """
        # sanity check
        self.validate_is_alma_data(vis)

        # get Jy/K value from DB
        jyperk = self.get(vis)

        # convert to pipeline-friendly format
        formatted = self.format_jyperk(vis, jyperk)
        #casalog.post('formatted = {}'.format(formatted))
        filtered = self.filter_jyperk(vis, formatted)
        #casalog.post('filtered = {}'.format(filtered))

        return filtered

    def get_params(self, vis):
        raise NotImplementedError

    def access(self, queries):
        raise NotImplementedError

    def get(self, vis):
        """
        Access Jy/K DB and return its response.

        Arguments:
            vis {str} -- Name of MS

        Raises:
            urllib2.HTTPError
            urllib2.URLError

        Returns:
            [dict] -- Response from the DB as a dictionary. It should contain
                      the following keys:
                          'query' -- query data
                          'total' -- number of data
                          'data'  -- data
        """
        # set URL
        url = self.url

        params = self.get_params(vis)

        queries = self._generate_query(url, params)

        retval = self.access(queries)
        # retval should be a dict that consists of
        # 'query': query data
        # 'total': number of data
        # 'data': response data
        return retval

    def format_jyperk(self, vis, jyperk):
        """
        Format given dictionary to the formatted list as below.

            [['MS_name', 'antenna_name', 'spwid', 'pol string', 'factor'],
             ['MS_name', 'antenna_name', 'spwid', 'pol string', 'factor'],
             ...
             ['MS_name', 'antenna_name', 'spwid', 'pol string', 'factor']]

        Arguments:
            vis {str} -- Name of MS
            jyperk {dict} -- Dictionary containing Jy/K factors with meta data

        Returns:
            [list] -- Formatted list of Jy/K factors
        """
        template = string.Template('$vis $Antenna $Spwid I $factor')
        data = jyperk['data']
        basename = os.path.basename(vis.rstrip('/'))
        factors = [list(map(str, template.safe_substitute(vis=basename, **d).split())) for d in data]
        return factors

    def filter_jyperk(self, vis, factors):
        ms = self.context.observing_run.get_ms(vis)
        science_windows = [x.id for x in ms.get_spectral_windows(science_windows_only=True)]
        filtered = [i for i in factors if (len(i) == 5) and (i[0] == ms.basename) and (int(i[2]) in science_windows)]
        return filtered


class JyPerKAbstractEndPoint(ALMAJyPerKDatabaseAccessBase):
    def get_params(self, vis):
        ms = self.context.observing_run.get_ms(vis)

        # parameter dictionary
        params = {}

        # date
        params['date'] = mjd_to_datestring(ms.start_time)

        # temperature
        params['temperature'] = get_mean_temperature(vis)

        # other
        params.update(self._aux_params())

        # loop over antennas and spws
        for ant in ms.antennas:
            # antenna name
            params['antenna'] = ant.name

            # elevation
            params['elevation'] = get_mean_elevation(self.context, vis, ant.id)

            for spw in ms.get_spectral_windows(science_windows_only=True):
                # observing band is taken from the string spw.band
                # whose format should be "ALMA Band X"
                params['band'] = int(spw.band.split()[-1])

                # baseband
                params['baseband'] = int(spw.baseband)

                # mean frequency
                params['frequency'] = get_mean_frequency(spw)

                # subparam is dictionary holding vis and spw id
                subparam = {'vis': vis, 'spwid': spw.id}
                yield QueryStruct(param=params, subparam=subparam)

    def access(self, queries):
        data = []
        for result in queries:
            # response from DB
            response = result.response

            # subparam is dictionary holding vis and spw id
            subparam = result.subparam
            assert isinstance(subparam, dict)
            assert ('vis' in subparam) and ('spwid' in subparam)
            spwid = subparam['spwid']
            assert isinstance(spwid, int)
            vis = subparam['vis']
            assert isinstance(vis, str)
            basename = os.path.basename(vis.rstrip('/'))

            factor = self._extract_factor(response)
            polarization = 'I'
            antenna = response['query']['antenna']
            data.append({'MS': basename, 'Antenna': antenna, 'Spwid': spwid,
                         'Polarization': polarization, 'factor': factor})

        return {'query': '', 'data': data, 'total': len(data)}

    def _aux_params(self):
        return {}

    def _extract_factor(self, response):
        raise NotImplementedError


class JyPerKAsdmEndPoint(ALMAJyPerKDatabaseAccessBase):
    ENDPOINT_TYPE = 'asdm'

    def get_params(self, vis):
        # subparam is vis
        yield QueryStruct(param={'uid': vis_to_uid(vis)}, subparam=vis)

    def access(self, queries):
        responses = list(queries)

        # there should be only one query
        assert len(responses) == 1

        response = responses[0].response
        response['total'] = response['data']['length']
        response['data'] = response['data']['factors']
        return response


class JyPerKModelFitEndPoint(JyPerKAbstractEndPoint):
    ENDPOINT_TYPE = 'model-fit'

    def _extract_factor(self, response):
        return response['data']['factor']


class JyPerKInterpolationEndPoint(JyPerKAbstractEndPoint):
    ENDPOINT_TYPE = 'interpolation'

    def _aux_params(self):
        return {'delta_days': 1000}

    def _extract_factor(self, response):
        return response['data']['factor']['mean']


def vis_to_uid(vis):
    """
    Convert MS name like uid___A002_Xabcd_X012 into uid://A002/Xabcd/X012

    Arguments:
        vis {str} -- Name of MS

    Raises:
        RuntimeError:

    Returns:
        str -- Corresponding ASDM uid
    """
    basename = os.path.basename(vis.rstrip('/'))
    pattern = '^uid___A[0-9][0-9][0-9]_X[0-9a-f]+_X[0-9a-f]+\.ms$'
    if re.match(pattern, basename):
        return basename.replace('___', '://').replace('_', '/').replace('.ms', '')
    else:
        raise RuntimeError('MS name is not appropriate for DB query: {}'.format(basename))


def mjd_to_datestring(epoch):
    # casa_tools
    me = casa_tools.measures
    qa = casa_tools.quanta

    if epoch['refer'] != 'UTC':
        try:
            epoch = me.measure(epoch, 'UTC')
        finally:
            me.done()

    t = qa.splitdate(epoch['m0'])
    dd = datetime.datetime(t['year'], t['month'], t['monthday'], t['hour'], t['min'], t['sec'], t['usec'])
    #datestring = dd.strftime('%Y-%m-%dT%H:%M:%S.%f')
    datestring = dd.strftime('%Y-%m-%dT%H:%M:%S')
    return datestring


def get_mean_frequency(spw):
    return float(spw.mean_frequency.convert_to(measures.FrequencyUnits.HERTZ).value)


def get_mean_temperature(vis):
    with casa_tools.TableReader(os.path.join(vis, 'WEATHER')) as tb:
        valid_temperatures = np.ma.masked_array(
            tb.getcol('TEMPERATURE'),
            tb.getcol('TEMPERATURE_FLAG')
        )

    return valid_temperatures.mean()


def get_mean_elevation(vis, antenna_id):
    ms = mstool()
    ms.open(vis)
    ms.msselect({'spw': spw, 'scanintent': 'OBSERVE_TARGET#ON_SOURCE'})
    selected = ms.msselectedindices()
    ms.close()
    stateid = selected['stateid']
    ddid = selected['spwdd'][0]

    tb = table()
    tb.open(vis)
    query = f'ANTENNA1=={antenna_id}&&ANTENNA2=={antenna_id}&&DATA_DESC_ID=={ddid}&&STATE_ID IN {list(stateid)}' 
    tsel = tb.query(query)
    rows = tsel.rownumbers()
    tsel.close()
    tb.close()

    qa = quanta() 
    msmd = msmetadata()
    msmd.open(vis) 
    elevations = [] 
    for row in rows: 
        p = msmd.pointingdirection(row, initialrow=row) 
        assert p['antenna1']['pointingdirection']['refer'].startswith('AZEL') 
        el_deg = qa.convert(p['antenna1']['pointingdirection']['m1'], 'deg') 
        elevations.append(el_deg['value']) 
    msmd.close()
    elevations = np.asarray(elevations)


def translate_spw(data, ms):
    vis = ms.name
    science_windows = np.asarray(ms.get_spectral_windows(science_windows_only=True))
    with casa_tools.TableReader(os.path.join(vis, 'ASDM_SPECTRALWINDOW')) as tb:
        idcol = tb.getcol('spectralWindowId')
        namecol = tb.getcol('name')

    translated = []
    science_window_names = np.asarray([x.name for x in science_windows])
    casalog.post('Translate ASDM Spws to MS Spws:')
    for d in data:
        asdm_spw_id = d['Spwid']
        asdm_spw_names = namecol[np.where(idcol == 'SpectralWindow_{}'.format(asdm_spw_id))]
        assert len(asdm_spw_names) == 1
        asdm_spw_name = asdm_spw_names[0]
        if asdm_spw_name.endswith('CH_AVG'):
            chan_avg_name = asdm_spw_name
            full_res_name = asdm_spw_name.replace('CH_AVG', 'FULL_RES')
        elif asdm_spw_name.endswith('FULL_RES'):
            chan_avg_name = asdm_spw_name.replace('FULL_RES', 'CH_AVG')
            full_res_name = asdm_spw_name
        else:
            chan_avg_name = asdm_spw_name
            full_res_name = asdm_spw_name
        i = np.where(science_window_names == full_res_name)
        if len(i[0]) == 0:
            i = np.where(science_window_names == chan_avg_name)
        if len(i[0]) > 0:
            spws = science_windows[i]
            assert len(spws) == 1
            spw = spws[0]
            t = d.copy()
            t['Spwid'] = '{}'.format(spw.id)
            translated.append(t)
            casalog.post('   * ASDM Spw {} (name {})'.format(asdm_spw_id, asdm_spw_name))
            casalog.post('    -> MS Spw {} (name {})'.format(spw.id, spw.name))
    return translated


# file part
class JyPerKReader4File():
    def __init__(self, filename):
        self.filename = filename
        
    def get(self):
        """
        Reads jyperk factors from a file and returns a string list
        of [['MS','ant','spwid','polid','factor'], ...]
        """
        with open(self.filename, 'r') as f:
            return list(self._extract_jyperk_from_csv(f))

    def _extract_jyperk_from_csv(self, stream):
        reader = csv.reader(stream)
        # Check if first line is header or not
        line = next(reader)
        if len(line) == 0 or line[0].strip().upper() == 'MS' or line[0].strip()[0] == '#':
            # must be a header, commented line, or empty line
            pass
        elif len(line) == 5:
            # may be a data record
            yield line
        else:
            casalog.post('Jy/K factor file is invalid format')
        for line in reader:
            if len(line) == 0 or len(line[0]) == 0 or line[0][0] == '#':
                continue
            elif len(line) == 5:
                yield line
            else:
                casalog.post('Jy/K factor file is invalid format')
