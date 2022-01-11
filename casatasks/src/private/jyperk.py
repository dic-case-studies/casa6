import collections
import csv
import json
import os
import re
from socket import timeout as socket_timeout
import ssl
import string
from time import sleep
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import certifi
import numpy as np

from casatasks import casalog
from casatasks.private.sdutil import table_selector, table_manager, tool_manager
from casatools import measures
from casatools import ms as mstool
from casatools import msmetadata, quanta


# Jy/K DB part
def gen_factor_via_web_api(vis, spw='*',
                           endpoint='asdm',
                           timeout=180, retry=3, retry_wait_time=5):
    """Generate factors via Jy/K Web API.

    This function will be used task_gencal.

    Arguments:
        vis {str} -- The file path of the visibility data.
        spw {str} -- Spectral windows.
        endpoint {str} -- The endpoint of Jy/K DB Web API to access. Options are
            'asdm' (default), 'model-fit', 'interpolation'.
        timeout {int} -- Maximum waiting time [sec] for the Web API access, defaults
            to 180 sec.
        retry {int} -- Number of retry when the Web API access fails, defaults to 3
            times.
        retry_wait_time {int} -- Waiting time [sec] until next query when the Web API
            access fails, defaults to 5 sec.
    """
    if spw == '':
        spw = '*'

    assert endpoint in ['asdm', 'model-fit', 'interpolation'], \
        'The JyPerKDatabaseClient class requires one of endpoint: asdm, model-fit or interpolation'

    return __factor_creator_via_jy_per_k_db(endpoint=endpoint, vis=vis, spw=spw,
                                            factory=__jyperk_factory[endpoint],
                                            timeout=timeout, retry=retry,
                                            retry_wait_time=retry_wait_time)


def __factor_creator_via_jy_per_k_db(endpoint='', vis=None, spw='*',
                                     factory=None,
                                     timeout=180, retry=3, retry_wait_time=5):
    params_generator = factory[0]
    response_translator = factory[1]

    params = params_generator.get_params(vis, spw=spw)
    client = JyPerKDatabaseClient(endpoint,
                                  timeout=timeout, retry=retry,
                                  retry_wait_time=retry_wait_time)
    manager = RequestsManager(client)
    resps = manager.get(params)
    return response_translator.convert(resps, vis, spw=spw)


QueryStruct = collections.namedtuple('QueryStruct', ['param', 'subparam'])
ResponseStruct = collections.namedtuple('ResponseStruct', ['response', 'subparam'])


class ASDMParamsGenerator():
    """A class to generate required parameters for Jy/K Web API with asdm.

    Usage:
        vis = "./uid___A002_X85c183_X36f.ms"
        params = ASDMParamsGenerator.get_params(vis)
    """

    @classmethod
    def get_params(cls, vis, spw=None):
        """Generate required parameters for Jy/K Web API.

        Arguments:
            vis {str} -- The file path of the visibility data.
            spw {str} -- This parameter is not used. It is provided to align the
                calling method with other classes (InterpolationParamsGenerator and
                ModelFitParamsGenerator).

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
        """Convert MS name like 'uid___A002_Xabcd_X012 into uid://A002/Xabcd/X012'.

        Arguments:
            vis {str} -- The file path of the visibility data.

        Returns:
            str -- Corresponding ASDM uid.
        """
        basename = os.path.basename(os.path.abspath(vis))
        pattern = '^uid___A[0-9][0-9][0-9]_X[0-9a-f]+_X[0-9a-f]+\.ms$'
        if re.match(pattern, basename):
            return basename.replace('___', '://').replace('_', '/').replace('.ms', '')
        else:
            raise RuntimeError('MS name is not appropriate for DB query: {}'.format(basename))


class InterpolationParamsGenerator():
    """A class to generate required parameters for Jy/K Web API with interpolation.

    Usage:
        vis = './uid___A002_X85c183_X36f.ms'
        params = InterpolationParamsGenerator.get_params(vis)
    """

    @classmethod
    def get_params(cls, vis, spw='*'):
        if spw == '':
            spw = '*'

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
                subparam = {'vis': vis, 'spwid': int(sw_id)}
                yield QueryStruct(param=params, subparam=subparam)

    @staticmethod
    def _get_science_windows(vis, spw):
        ms = mstool()
        selected = ms.msseltoindex(vis, spw=spw)
        science_windows = selected['spw']
        return science_windows

    @staticmethod
    def _extract_msmetadata(science_windows, vis):
        with tool_manager(vis, msmetadata) as msmd:
            timerange = msmd.timerangeforobs(0)
            antenna_names = msmd.antennanames()
            basebands = dict((i, msmd.baseband(i)) for i in science_windows)
            mean_freqs = dict((i, msmd.meanfreq(i)) for i in science_windows)
            spwnames = msmd.namesforspws(science_windows)

        return timerange, antenna_names, basebands, mean_freqs, spwnames

    @staticmethod
    def _get_available_spw(vis, spw):
        science_windows = InterpolationParamsGenerator._get_science_windows(vis, spw=spw)
        with tool_manager(vis, msmetadata) as msmd:
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
        with table_manager(os.path.join(vis, 'WEATHER')) as tb:
            valid_temperatures = np.ma.masked_array(
                tb.getcol("TEMPERATURE"),
                tb.getcol("TEMPERATURE_FLAG")
            )
        return valid_temperatures.mean()

    @staticmethod
    def _get_mean_freqs(vis, science_windows):
        with table_manager(os.path.join(vis, 'SPECTRAL_WINDOW')) as tb:
            mean_freqs = dict((i, tb.getcell('CHAN_FREQ', i).mean()) for i in science_windows)
        return mean_freqs

    @staticmethod
    def _get_aux_params():
        return {'delta_days': 1000}


class ModelFitParamsGenerator(InterpolationParamsGenerator):
    """A class to generate required parameters for Jy/K Web API with model-fit."""

    @staticmethod
    def _get_aux_params():
        return {}


class Bands():
    """A class to extract all bands corresponding from VIS file.

    Usage:
        bands = Bands.get(science_windows, spwnames, mean_freqs, vis)
    """

    @classmethod
    def get(cls, science_windows, spwnames, mean_freqs, vis):
        """Return all bands corresponding to the 'science_window' given in the input.

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
        """Extract bands that contain band information in the spwname.

        The spwnames is like 'X835577456#ALMA_RB_06#BB_2#SW-01#CH_AVG'.
        """
        bands = {}
        for sw, spwname in zip(science_windows, spwnames):
            if 'ALMA_RB_' in spwname:
                bands[sw] = int(re.findall(r'^.*?ALMA_RB_(\d+)#.*', spwname)[0])
        return bands

    @staticmethod
    def _filter_mean_freqs_with_undetected_band(science_windows, spwnames, mean_freqs):
        """Filter mean freqs without 'ALMA_RB_'."""
        filtered_mean_freqs = {}
        for sw, spwname in zip(science_windows, spwnames):
            if 'ALMA_RB_' not in spwname:
                filtered_mean_freqs[sw] = mean_freqs[sw]
        return filtered_mean_freqs

    @staticmethod
    def _detect_bands_from_mean_freqs(target_mean_freqs, vis):
        """Extract bands using the mean freqs.

        Params:
            target_mean_freqs {dict} -- The mean freqs which does not been detected the bands.
            vis {str} -- The file path of the visibility data.
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
        with tool_manager(vis, msmetadata) as msmd:
            spwnames = msmd.namesforspws(science_windows)
        return spwnames

    @staticmethod
    def _get_known_bands(vis):
        science_windows = Bands._get_all_science_windows(vis)
        with tool_manager(vis, msmetadata) as msmd:
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
        with tool_manager(vis, msmetadata) as msmd:
            mean_freqs = dict((i, msmd.meanfreq(i)) for i in science_windows)
        return mean_freqs


class MeanElevation():
    """A class to extract elevations from the VIS file and calcurate elevations average.

    Usage:
        mean_elevation = MeanElevation.get(vis, antenna_id)
    """

    @classmethod
    def get(cls, vis, antenna_id):
        """Get elevations aveage."""
        stateid = cls._get_stateid(vis)
        science_dd = cls._get_science_dd(vis)
        rows = cls._query_rows(vis, science_dd, stateid, antenna_id)

        return cls._calc_elevation_mean(rows, vis)

    @staticmethod
    def _get_stateid(vis):
        with tool_manager(vis, mstool) as ms:
            ms.msselect({'scanintent': 'OBSERVE_TARGET#ON_SOURCE'})
            selected = ms.msselectedindices()

        stateid = selected['stateid']
        return stateid

    @staticmethod
    def _get_science_dd(vis):
        with tool_manager(vis, msmetadata) as msmd:
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
    def _calc_elevation_mean(rows, vis):
        elevations = []
        qa = quanta()

        with tool_manager(vis, msmetadata) as msmd:
            for row in rows:
                dirs = msmd.pointingdirection(row, initialrow=row)
                assert dirs['antenna1']['pointingdirection']['refer'].startswith('AZEL')
                el_deg = qa.convert(dirs['antenna1']['pointingdirection']['m1'], 'deg')
                elevations.append(el_deg['value'])

        elevations = np.asarray(elevations)

        if len(elevations) == 0:
            elevation_mean = np.nan
        else:
            elevation_mean = elevations.mean()

        return elevation_mean


class RequestsManager():
    """A class to manage the Jy/K Database access by the param.

    Usage:
        vis = "./uid___A002_Xb32033_X9067.ms"
        client = JyPerKDatabaseClient('asdm')
        params = ASDMParamsGenerator.get_params(vis)
        manager = RequestsManager(client)
        manager.get(params)
    """

    def __init__(self, client):
        """Set client."""
        self.client = client

    def get(self, params):
        """Get the responses of the Jy/K DB."""
        dataset = [{'response': self.client.get(param.param), 'aux': param.subparam} for param in params]
        return self._filter_success_is_true(dataset)

    def _filter_success_is_true(self, dataset):
        return [data for data in dataset if data['response']['success']]


class JyPerKDatabaseClient():
    """A class to get values from Jy/K Web API.

    The Jy/K Web API address is 'https://asa.alma.cl/science/jy-kelvins'. The address
    can be changed with the environment variable 'JYPERKDB_URL'.
    """

    BASE_URL = os.getenv('JYPERKDB_URL', 'https://asa.alma.cl/science/jy-kelvins')

    def __init__(self, endpoint, timeout=180, retry=3, retry_wait_time=5):
        """Set the parameters to be used when accessing the Web API.

        Arguments:
            endpoint {str} -- The endpoint of Jy/K DB Web API to access. Options are
                'asdm' (default), 'model-fit', 'interpolation'.
            timeout {int} --- Maximum waiting time [sec] for the Web API access, defaults
                to 180 sec.
            retry {int} -- Number of retry when the Web API access fails, defaults to 3
                times.
            retry_wait_time {int} -- Waiting time [sec] until next query when the Web API
                access fails, defaults to 5 sec.
        """
        assert endpoint in ['asdm', 'model-fit', 'interpolation'], \
            'The JyPerKDatabaseClient class requires one of endpoint: asdm, model-fit or interpolation'
        self.web_api_url = self._generate_web_api_url(endpoint)
        self.timeout = timeout
        self.retry = retry
        self.retry_wait_time = retry_wait_time

    def get(self, param):
        """Get the Web API response.

        Arguments:
            param {dict} -- The parameters used in the Web API.
        """
        request_url = self._generate_url(param)
        body = self._try_to_get_response(request_url)
        retval = self._convert_to_json(body)
        self._check_retval(retval)
        return retval

    def _generate_web_api_url(self, endpoint_type):
        web_api_url = '/'.join([self.BASE_URL, endpoint_type])
        if not web_api_url.endswith('/'):
            web_api_url += '/'
        return web_api_url

    def _generate_url(self, param):
        # encode params
        encoded = urlencode(param)
        query = '?'.join([self.web_api_url, encoded])
        return query

    def _retrieve(self, url):
        """Access to Jy/K DB and return response.

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
                return {'status': 'Success', 'err_msg': None, 'body': body}
        except HTTPError as e:  # 4xx, 5xx
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: HTTPError(code={0}, Reason="{1}")\n'.format(e.code, e.reason)
            casalog.post(msg)
            return {'status': 'HTTPError', 'err_msg': msg}
        except URLError as e:  # not connect
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: URLError(Reason="{0}")\n'.format(e.reason)
            casalog.post(msg)
            return {'status': 'URLError', 'err_msg': msg}
        except socket_timeout as e:  # not connect
            msg = 'Failed to load URL: {0}\n'.format(url) \
                + 'Error Message: URLError(Reason="{0}")\n'.format(e)
            casalog.post(msg)
            return {'status': 'URLError', 'err_msg': msg}


    def _try_to_get_response(self, url):
        casalog.post(f'Accessing Jy/K DB: request URL is "{url}"')
        for i in range(self.retry):
            response_with_tag = self._retrieve(url)
            if response_with_tag['status'] == 'Success':
                casalog.post(f'Got a response successfully')
                return response_with_tag['body']
            
            if i < self.retry - 1:
                casalog.post(response_with_tag['err_msg'])
                casalog.post(f'Sleeping for {str(self.retry_wait_time)} seconds because the request failed')
                sleep(self.retry_wait_time)
                casalog.post(f'Retry to access Jy/K DB ({str(i + 2)}/{str(self.retry)})')

        if response_with_tag['status'] != 'Success':
            raise RuntimeError(response_with_tag['err_msg'])

    def _convert_to_json(self, response):
        try:
            return json.loads(response)

        except json.JSONDecodeError as e:
            msg = 'Failed to get a Jy/K factor from DB: JSON Syntax error. {}'.format(e)
            casalog.post(msg)
            raise RuntimeError(msg)

    def _check_retval(self, retval):
        """ Check if 'success' of retval dict is True.

        This method only checks if the api was able to complete the process successfully or not.
        It is expected that 'success' will be False as a response, so the mothod does not raise
        RuntimeError. If the 'success' is False, the *Transelator classes will reject the factor
        value.
        """
        if not retval['success']:
            msg = 'Failed to get a Jy/K factor from DB: {}'.format(retval['error'])
            casalog.post(msg)


class Translator():
    """A class containing the methods required to convert Jy/K DB responses into factors."""

    @staticmethod
    def format_cal_table_format(factors):  # _format_jyperk
        """Format given dictionary to the formatted list.

        Sample formated list:
            [['MS_name', 'antenna_name', 'spwid', 'pol string', 'factor'],
             ['MS_name', 'antenna_name', 'spwid', 'pol string', 'factor'],
             ...
             ['MS_name', 'antenna_name', 'spwid', 'pol string', 'factor']]

        Arguments:
            factors {dict} -- Dictionary containing Jy/K factors with meta data.

        Returns:
            list -- Formatted list of Jy/K factors.
        """
        template = string.Template('$MS $Antenna $Spwid I $factor')
        factors = [list(map(str, template.safe_substitute(**factor).split())) for factor in factors]
        return factors

    @staticmethod
    def filter_jyperk(vis, factors, spw):
        """Filter factors using spw."""
        ms = mstool()
        selected = ms.msseltoindex(vis=vis, spw=spw)
        science_windows = selected['spw']
        filtered = [
            i for i in factors if (len(i) == 5)
            and (
                i[0] == os.path.basename(os.path.abspath(vis)) and (int(i[2]) in science_windows)
            )
        ]
        return filtered


class ASDMRspTranslator():
    """A class to convert the response for asdm from the Jy/K DB to factors."""

    @classmethod
    def convert(cls, data, vis, spw='*'):
        """Convert from the response to list with factor.

        Arguments:
            spw {str} -- This parameter is not used. It is provided to align the
                    calling method with other classes (InterpolationRspTranslator and
                    ModelFitRspTranslator).

        Returns:
            list -- List of Jy/K conversion factors with meta data.
        """
        assert len(data) == 1
        data = data[0]

        factors = cls._extract_factors(data)
        formatted = Translator.format_cal_table_format(factors)

        return Translator.filter_jyperk(vis, formatted, spw)

    @staticmethod
    def _extract_factors(data):
        return data['response']['data']['factors']


class InterpolationRspTranslator():
    """A class to convert the responses for interpolation from the Jy/K DB to factors."""

    @classmethod
    def convert(cls, data_set, vis, spw='*'):
        """Convert from the response to list with factor.

        Arguments:
            data_set {dict} -- The result of the Web API.
            vis {str} -- The file path of the visibility data.
            spw {str} -- Spectral windows.

        Returns:
            list -- List of Jy/K conversion factors with meta data.
        """
        cal_dict = cls._dataset_to_cal_dict(data_set, cls._extract_factor)
        formatted = Translator.format_cal_table_format(cal_dict)

        return Translator.filter_jyperk(vis, formatted, spw)

    @staticmethod
    def _extract_factor(data):
        if data['response']['query']['elevation'] == 'nan':
            casalog.post('The elevation mean is np.nan, so the factor is set 1.',
                         priority='WARN',
                         origin='jyperk.MeanElevation._calc_elevation_mean')
            return 1

        return data['response']['data']['factor']['mean']

    @staticmethod
    def _dataset_to_cal_dict(dataset, _extract_factor):
        return_data = []

        for data in dataset:
            # aux is dictionary holding vis and spw id
            aux = data['aux']
            if not isinstance(aux, dict):
                raise TypeError('The response.aux in the JSON obtained from Jy/K db must be dict.')

            if 'vis' not in aux:
                raise KeyError('The response.aux in the JSON obtained from Jy/K db must contain vis.')

            if 'spwid' not in aux:
                raise KeyError('The response.aux in the JSON obtained from Jy/K db must contain spwid.')

            spwid = aux['spwid']
            if not isinstance(spwid, int):
                raise TypeError('The response.aux.spwid in the JSON obtained from Jy/K db must be int.')

            vis = aux['vis']
            if not isinstance(vis, str):
                raise TypeError('The response.aux.vis in the JSON obtained from Jy/K db must be str.')

            basename = os.path.basename(os.path.abspath(vis))

            factor = _extract_factor(data)
            polarization = 'I'
            antenna = data['response']['query']['antenna']

            return_data.append({'MS': basename, 'Antenna': antenna, 'Spwid': spwid,
                                'Polarization': polarization, 'factor': factor})
        return return_data


class ModelFitRspTranslator(InterpolationRspTranslator):
    """A class to convert the responses for model-fit from the Jy/K DB to factors."""

    @staticmethod
    def _extract_factor(data):
        return data['response']['data']['factor']


__jyperk_factory = {
    'asdm': (ASDMParamsGenerator, ASDMRspTranslator),
    'interpolation': (InterpolationParamsGenerator, InterpolationRspTranslator),
    'model-fit': (ModelFitParamsGenerator, ModelFitRspTranslator),
}


# file part
class JyPerKReader4File():
    """A class to read from CSV file and format factors.

    This function will be used task_gencal.
    """

    def __init__(self, infile):
        """Set a parameter.

        Arguments:
            infile {str} -- The file path of CSV which the factors are stored.
        """
        self.infile = infile

    def get(self):
        """Read jyperk factors from a file and returns a string list.

        Returns:
            list -- [['MS','ant','spwid','polid','factor'], ...]
        """
        if not os.path.isfile(self.infile):
            raise OSError(f'There is no Jy/K db-derived factor file: {self.infile}')

        with open(self.infile, 'r') as f:
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
