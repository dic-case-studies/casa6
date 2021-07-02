import certifi
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
from casatasks.private.sdutil import table_selector, tbmanager, toolmanager


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
    def get_params(cls, vis, spw='*'):
        if spw == '':
            spw='*'

        if spw == '*':
            spw = cls._get_available_spw(vis, spw)
        
        params = {}

        science_windows = cls._get_science_windows(vis, spw)
        timerange, antenna_names, basebands, mean_freqs, spwnames = cls._extract_msmetadata(science_windows, vis)

        mean_freqs = cls._get_mean_freqs(vis, science_windows)
        bands = Bands.get(science_windows, spwnames, mean_freqs, vis)

        params['date'] = cls._mjd_to_datestring(timerange['begin'])
        params['temperature'] = cls._get_mean_temperature(vis)
        params.update(cls._get_aux_params())

        for antenna_id, antenna_name in enumerate(antenna_names):
            params['antenna'] = antenna_name
            params['elevation'] = MeanElevation.get(vis, antenna_id)
 
            for sw_id in science_windows:
                params['band'] = bands[sw_id]
                params['baseband'] = basebands[sw_id]
                params['frequency'] = mean_freqs[sw_id]
                subparam = {'vis': vis, 'spwid': sw_id}
                yield QueryStruct(param=params, subparam=subparam)

    @staticmethod
    def _get_science_windows(vis, spw):
        ms = mstool()
        selected = ms.msseltoindex(vis, spw=spw)
        science_windows = selected['spw']
        return science_windows

    @staticmethod
    def _extract_msmetadata(science_windows, vis):
        with toolmanager(vis, msmetadata) as msmd:  
            timerange = msmd.timerangeforobs(0)
            antenna_names = msmd.antennanames()
            basebands = dict((i, msmd.baseband(i)) for i in science_windows)
            mean_freqs = dict((i, msmd.meanfreq(i)) for i in science_windows)
            spwnames = msmd.namesforspws(science_windows)

        return timerange, antenna_names, basebands, mean_freqs, spwnames
   
    @staticmethod
    def _get_available_spw(vis, spw):
        science_windows = InterpolationParamsGenerator._get_science_windows(vis, spw=spw)
        with toolmanager(vis, msmetadata) as msmd:  
            spwnames = msmd.namesforspws(science_windows)

        spw = ','.join(map(str, [i for i, name in enumerate(spwnames) if not name.startswith('WVR')]))
        return spw

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
        return datestring[0]

    @staticmethod
    def _get_mean_temperature(vis):
        with tbmanager(os.path.join(vis, 'WEATHER')) as tb:
            valid_temperatures = np.ma.masked_array(
                tb.getcol("TEMPERATURE"),
                tb.getcol("TEMPERATURE_FLAG")
            )       
        return valid_temperatures.mean()
    
    @staticmethod
    def _get_mean_freqs(vis, science_windows):
        with tbmanager(os.path.join(vis, 'SPECTRAL_WINDOW')) as tb:
            mean_freqs = dict((i, tb.getcell('CHAN_FREQ', i).mean()) for i in science_windows)
        return mean_freqs

    @staticmethod
    def _get_aux_params(self):
        return {'delta_days': 1000}


class ModelFitParamsGenerator(InterpolationParamsGenerator):
    @staticmethod
    def _get_aux_params(self):
        return {}


class Bands():
    @classmethod
    def get(cls, science_windows, spwnames, mean_freqs, vis):
        """ Return all bands corresponding to the 'science_window' given in the input.
        First the method scan 'spwnames', if the band can be detect, the method will
        adopt this value. In other case, the method compare the freq with the 'mean_freqs'
        at which the band was detect, the method detect the band from the frequencies 
        that are closest to the result.
        """
        bands = cls._extract_bands_from_spwnames(science_windows, spwnames)
        mean_freqs_with_undetected_band = cls._filter_mean_freqs_with_undetected_band(
                                            science_windows, spwnames, mean_freqs)
        if len(mean_freqs_with_undetected_band) > 0:
            bands.update(
                cls._detect_bands_from_mean_freqs(mean_freqs_with_undetected_band, vis)
            )
        return bands

    @staticmethod
    def _extract_bands_from_spwnames(science_windows, spwnames):
        """ Extract bands that contain band information in the spwname.
        The spwnames is like 'X835577456#ALMA_RB_06#BB_2#SW-01#CH_AVG'.
        """
        bands = {}
        for sw, spwname in zip(science_windows, spwnames):
            if 'ALMA_RB_' in spwname:
                bands[sw] = int(re.findall(r'^.*?ALMA_RB_(\d+)#.*', spwname)[0])
        return bands

    @staticmethod
    def _filter_mean_freqs_with_undetected_band(science_windows, spwnames, mean_freqs):
        """ Filter mean freqs without 'ALMA_RB_'.
        """
        filtered_mean_freqs = {}
        for sw, spwname in zip(science_windows, spwnames):
            if not 'ALMA_RB_' in spwname:
                filtered_mean_freqs[sw] = mean_freqs[sw]
        return filtered_mean_freqs

    @staticmethod
    def _detect_bands_from_mean_freqs(target_mean_freqs, vis):
        """ Extract bands using the mean freqs.
        
        Params:
            target_mean_freqs {dict}: The mean freqs which does not been detected the bands.
            vis {char}: The file path of vis.
        """
        known_bands = Bands._get_known_bands(vis)
        science_windows = list(known_bands.keys())
        mean_freqs = Bands._extract_mean_freqs(science_windows, vis)

        extracted_bands = {}
        for spw, target_mean_freq in target_mean_freqs.items():
            nearest_spw = Bands._calc_nearest_spw(mean_freqs, target_mean_freq)
            extracted_bands[spw] = known_bands[nearest_spw]
        return extracted_bands

    @staticmethod
    def _calc_nearest_spw(available_mean_freqs, mean_freq):
        available_mean_freqs_list = list(available_mean_freqs.values())
        available_spw = list(available_mean_freqs.keys())
        nearest_i = np.argmin(np.abs(np.array(available_mean_freqs_list) - mean_freq))
        return available_spw[nearest_i]


    @staticmethod
    def _get_spwnames(vis, science_windows):
        with toolmanager(vis, msmetadata) as msmd:
            spwnames = msmd.namesforspws(science_windows)
        return spwnames

    @staticmethod
    def _get_known_bands(vis):
        science_windows = Bands._get_all_science_windows(vis)
        with toolmanager(vis, msmetadata) as msmd:
            spwnames = msmd.namesforspws(science_windows)    
        bands = Bands._extract_bands_from_spwnames(science_windows, spwnames)
        return bands
    
    @staticmethod
    def _get_all_science_windows(vis):
        ms = mstool()
        selected = ms.msseltoindex(vis, spw='*')
        science_windows = selected['spw']
        return science_windows

    @staticmethod
    def _extract_mean_freqs(science_windows, vis):
        with toolmanager(vis, msmetadata) as msmd:  
            mean_freqs = dict((i, msmd.meanfreq(i)) for i in science_windows)
        return mean_freqs


class MeanElevation():
    @classmethod
    def get(cls, vis, antenna_id):
        stateid = cls._get_stateid(vis)
        science_dd = cls._get_science_dd(vis)
        rows = cls._query_rows(vis, science_dd, stateid, antenna_id)

        return cls._calc_elevation_mean(rows)

    @staticmethod
    def _get_stateid(vis): ###
        with toolmanager(vis, mstool) as ms:
            ms.msselect({'scanintent': 'OBSERVE_TARGET#ON_SOURCE'})
            selected = ms.msselectedindices()

        stateid = selected['stateid']
        return stateid

    @staticmethod
    def _get_science_dd(vis):
        with toolmanager(vis, msmetadata) as msmd:  
            science_spw = list(np.intersect1d(
                msmd.almaspws(tdm=True, fdm=True),
                msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
            ))
            science_dd = [msmd.datadescids(spw=i)[0] for i in science_spw]

        return science_dd
        
    @staticmethod
    def _query_rows(vis, science_dd, stateid, antenna_id):
        query = f'ANTENNA1=={antenna_id}&&ANTENNA2=={antenna_id}&&DATA_DESC_ID=={science_dd[0]}&&STATE_ID IN {list(stateid)}'
        with table_selector(vis, query) as tb:
            rows = tb.rownumbers()

        return rows

    @staticmethod
    def _calc_elevation_mean(rows):
        elevations = []       
        qa = quanta()

        with toolmanager(vis, msmetadata) as msmd:  
            for row in rows: 
                p = msmd.pointingdirection(row, initialrow=row) 
                assert p['antenna1']['pointingdirection']['refer'].startswith('AZEL') 
                el_deg = qa.convert(p['antenna1']['pointingdirection']['m1'], 'deg') 
                elevations.append(el_deg['value']) 

        elevations = np.asarray(elevations)
        return elevations.mean()


class JyPerKDatabaseClient():
    """
    Get values from Jy/K Web API (https://asa.alma.cl/science/jy-kelvins).
    
    Arguments:
        endpoint_type {str} -- Endpoint of Jy/K Web API.
            The value to be entered must be one of asdm, model-fit or interpolation.
        timeout {int} --- Maximum waiting time when accessing the web API.
        retry {int} -- Number of times to retry when the web API access fails.
    """
    BASE_URL = 'https://asa.alma.cl/science/jy-kelvins'

    def __init__(self, endpoint_type, timeout=180, retry=3):
        assert endpoint_type in ['asdm', 'model-fit', 'interpolation'], \
            'Please set endpoint_type: asdm, model-fit, interpolation'
        self.web_api_url = self._generate_web_api_url(endpoint_type)
        self.timeout = timeout
        self.retry = retry

    def get(self, param):
        request_url = self._generate_query(param)
        body = self._try_to_get_resonse(request_url)
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

    def _retrieve(self, url):
        """
        Arguments:
            url {str} -- url to retrieve in the Jy/K Web API.

        Returns:
            dict -- If the request is successfull, dictionary contain below information.
                success {bool} -- Boolean stating whether the request succeeded or not.
                timestamp {str} -- Timestamp of when the request was received by the Jy/K Service.
                elapsed {bool} -- Boolean stating whether the request succeeded or not.
                error {dict} -- Error field which has a message describing the problem.
                query {dict} -- Dictionary with all the parameters passed to the request.
                data {dict} -- 	Data field which varies for each endpoint.

                Please check the following source for details.
                https://confluence.alma.cl/pages/viewpage.action?pageId=35258466
        """
        try:
            ssl_context = ssl.create_default_context(cafile=certifi.where())
            with urlopen(url, context=ssl_context, timeout=self.timeout) as resp:
                body = resp.read()
                return {'status': 'Success', 'code': 200, 'body': body}
        except HTTPError as e: # 4xx, 5xx
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: HTTPError(code={0}, Reason="{1}")\n'.format(e.code, e.reason)
            casalog.post(msg)
            return {'HTTPError': False, 'code': e.code, 'error_reason': e.reason}
        except URLError as e: # not connect
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: URLError(Reason="{0}")\n'.format(e.reason)
            casalog.post(msg)
            return {'URLError': False, 'code': None, 'error_reason': e.reason}

    def _try_to_get_resonse(self, url):
        for i in range(self.retry):
            response_with_tag = self._retrieve(url)
            if response_with_tag['status'] == 'Success':
                return response_with_tag['body']

        if response_with_tag['status'] == HTTPError:
            self._raise_http_error(url, response_with_tag)
            
    def _raise_http_error(self, url, response_with_tag):
        msg = 'Failed to load URL: {0}\n'.format(url) \
            + 'Error Message: HTTPError(code={0}, Reason="{1}")\n'.format(e.code, e.reason)
        casalog.post(msg)
        raise RuntimeError(msg)

    def _raise_url_error(self, response_with_tag):
        msg = 'Failed to load URL: {0}\n'.format(url) \
            + 'Error Message: URLError(Reason="{0}")\n'.format(e.reason)
        casalog.post(msg)
        raise RuntimeError(msg)

    def _convert_to_json(self, response):
        try:
            return json.loads(response)

        except json.JSONDecodeError as e:
            msg = 'Failed to get a Jy/K factor from DB: JSON Syntax error. {}'.format(e)
            casalog.post(msg)
            raise RuntimeError(msg)

    def _check_retval(self, retval):
        """
        This method only checks if the api was able to complete the process successfully or not.
        """
        if not retval['success']:
            msg = 'Failed to get a Jy/K factor from DB: {}'.format(retval['error'])
            casalog.post(msg)
            raise RuntimeError(msg)


class ASDMRspTranslator():
    def convert(self, response):
        """
        Arguments:
        Returns:
            [list] -- List of Jy/K conversion factors with meta data
        """
        # convert to pipeline-friendly format
        formatted = self.format_jyperk(vis, jyperk)
        #casalog.post('formatted = {}'.format(formatted))
        filtered = self.filter_jyperk(vis, formatted)
        #casalog.post('filtered = {}'.format(filtered))

        return filtered

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

    def format(self, queries):
        responses = list(queries)

        # there should be only one query
        assert len(responses) == 1

        response = responses[0].response
        response['total'] = response['data']['length']
        response['data'] = response['data']['factors']
        return response


class InterpolationRspTranslator(ASDMRspTranslator):
    def format(self, queries):
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

    def _extract_factor(self, response):
        return response['data']['factor']['mean']


class ModelFitRspTranslator(InterpolationRspTranslator):
    def _extract_factor(self, response):
        return response['data']['factor']


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
