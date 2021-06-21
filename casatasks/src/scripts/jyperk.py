import collections
import csv
import datetime
import json
import os
import re
import ssl
import string
from urllib.parse import urlencode
from urllib.request import urlopen
from urllib.error import HTTPError, URLError

import numpy as np

from casatasks import casalog
from casatools import ms as mstool
from casatools import msmetadata 
from casatools import measures
from casatools import quanta 
from casatools import table 


### web api part

QueryStruct = collections.namedtuple('QueryStruct', ['param', 'subparam'])
ResponseStruct = collections.namedtuple('ResponseStruct', ['response', 'subparam'])


class ASDMParamsGenerator():
    """
    Generate required parameters for Jy/K Web API.
    This class has a public class method to generate the parameters.

    Usage:
        vis = "./uid___A002_X85c183_X36f.ms"
        params = ASDMParamsGenerator.get_params(vis)
    """

    @classmethod
    def get_params(cls, vis):
    """
    Generate required parameters for Jy/K Web API.

    Arguments:
        vis {str} -- File path of MS

    Returns:
        Generator Object -- yield QueryStruct() object. A sample like below. 
            QueryStruct(
                param={'uid': 'uid://A002/X85c183/X36f'}, 
                subparam='./uid___A002_X85c183_X36f.ms'
            )
    """  
        yield QueryStruct(param={'uid': cls._vis_to_uid(vis)}, subparam=vis)

    @staticmethod
    def _vis_to_uid(vis):
        """
        Convert MS name like uid___A002_Xabcd_X012 into uid://A002/Xabcd/X012

        Arguments:
            vis {str} -- File path of MS

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

class InterpolationParamsGenerator():
    """
    Usage:
        vis = './uid___A002_X85c183_X36f.ms'
        params = InterpolationParamsGenerator.get_params(vis)
    """
    @classmethod
    def get_params(cls, vis, spw=''):
        params = {}
        
        if spw == '':
            spw = '*'

        ms = mstool()
        selected = ms.msseltoindex(vis, spw=spw)
        science_windows = selected['spw']

        msmd = msmetadata()
        msmd.open(vis) 
        timerange = msmd.timerangeforobs(0)
        antenna_names = msmd.antennanames()
        basebands = dict((i, msmd.baseband(i)) for i in science_windows)
        mean_freqs = dict((i, msmd.meanfreq(i)) for i in science_windows)
        spwnames = msmd.namesforspws(science_windows)
        msmd.close()
        
        bands = {}
        for i, n in zip(science_windows, spwnames):
            if '#' in n and '_' in n:
                bands[i] = int(n.split('#')[0].split('_')[-1])
                
        # bands = dict((i, int(n.split('#')[0].split('_')[-1])) for i, n in zip(science_windows, spwnames))
        params['date'] = cls._mjd_to_datestring(timerange['begin'])

        tb = table()
        tb.open(os.path.join(vis, 'SPECTRAL_WINDOW'))
        spw_names = [tb.getcell('NAME', i) for i in science_windows]

        mean_freqs = [tb.getcell('CHAN_FREQ', i).mean() for i in science_windows]

        for antenna_id, antenna_name in enumerate(antenna_names):
            params['antenna'] = antenna_name
            params['elevation'] = cls._get_mean_elevation(vis, antenna_id, spw)

            for spw in science_windows:
                params['band'] = bands[spw]
                params['baseband'] = basebands[spw]
                params['frequency'] = mean_freqs[spw]
                subparam = {'vis': vis, 'spwid': spw}
                yield QueryStruct(param=params, subparam=subparam)

    @staticmethod
    def _mjd_to_datestring(epoch):
        me = measures()
        qa = quanta()

        if epoch['refer'] != 'UTC':
            try:
                epoch = me.measure(epoch, 'UTC')
            finally:
                me.done()

        datestring = qa.time(epoch['m0'], form='fits')
        return datestring

    @staticmethod
    def _get_mean_elevation(vis, antenna_id, spw):
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

        return elevations.mean()



class ModelFitParamsGenerator(InterpolationParamsGenerator):
    pass


class JyPerKDatabaseClient():
    BASE_URL = 'https://asa.alma.cl/science/jy-kelvins'

    def __init__(self, endpoint_type, id=0):
        assert endpoint_type in ['asdm', 'model-fit', 'interpolation'], \
            'Please set endpoint_type: asdm, model-fit, interpolation'
        self.web_api_url = self._generate_web_api_url(endpoint_type)
        self.id = 0

    def get(self, param):
        request_url = self._generate_query(param)
        body = self._retrieve(request_url)
        retval = self._convert_to_json(body)
        self._check_retval(retval)
        return retval

    def _generate_web_api_url(self, endpoint_type):
        web_api_url = '/'.join([self.BASE_URL, endpoint_type])
        if not web_api_url.endswith('/'):
            web_api_url += '/'
        return web_api_url

    def _generate_query(self, param):
        # encode params
        encoded = urlencode(param)
        query = '?'.join([self.web_api_url, encoded])
        casalog.post('Accessing Jy/K DB: query is "{}"'.format(query))
        return query

    def _retrieve(self, url, timeout=180):
        try:
            with urlopen(url) as resp:
                body = resp.read()
                return body.decode('utf-8')
        except HTTPError as e:
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: HTTPError(code={0}, Reason="{1}")\n'.format(e.code, e.reason)
            casalog.post(msg)
            return {'success': False}
        except URLError as e:
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: URLError(Reason="{0}")\n'.format(e.reason)
            casalog.post(msg)
            return {'success': False}

    def _convert_to_json(self, response):
        try:
            return json.loads(response)
        except ValueError as e:
            msg = 'Failed to get a Jy/K factor from DB: The response is not JSON format'
            casalog.post(msg)
            return {'success': False}

    def _check_retval(self, retval):
        if not retval['success']:
            msg = 'Failed to get a Jy/K factor from DB: {}'.format(retval['error'])
            casalog.post(msg)
            return {'success': False}


class ALMAJyPerKDatabaseAccessBase(object):
    ENDPOINT_TYPE = None

    @property
    def url(self):
        assert self.ENDPOINT_TYPE is not None, \
            '{} cannot be instantiated. Please use subclasses.'.format(self.__class__.__name__)

    def __init__(self):
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
        """

    def _generate_query(self, params):
        for p in params:
            client = JyPerKDatabaseClient(self.ENDPOINT_TYPE)
            retval = client.get(p)
            yield ResponseStruct(response=retval, subparam=p.subparam)

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

        Returns:
            [dict] -- Response from the DB as a dictionary. It should contain
                      the following keys:
                          'query' -- query data
                          'total' -- number of data
                          'data'  -- data
        """
        params = self.get_params(vis)

        queries = self._generate_query(params)

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
        ms = mstool()
        selected = ms.msseltoindex(vis=vis, spw=spw)
        science_windows = selected['spw']
        filtered = [
            i for i in factors if (len(i) == 5) 
            and (
                i[0] == os.path.basename(vis.rstrip('/')) and (int(i[2]) in science_windows)
            )
        ]
        return filtered


class JyPerKAbstractEndPoint(ALMAJyPerKDatabaseAccessBase):


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
