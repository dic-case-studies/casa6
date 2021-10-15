from copy import deepcopy
import json
import math
import os
import shutil
import unittest
from unittest.mock import patch, MagicMock
from urllib.error import HTTPError, URLError

from casatasks.private import jyperk
from casatasks.private.casa_transition import is_CASA6

if is_CASA6:
    from casatools import ctsys


class TestASDMParamsGenerator(unittest.TestCase):
    """test ASDMParamsGenerator class.
    """
    vis = "./uid___A002_X85c183_X36f.ms"
        
    def test_get_params(self):
        params = jyperk.ASDMParamsGenerator.get_params(self.vis)
        for param in params:
            self.assertEqual(param.param, {'uid': 'uid://A002/X85c183/X36f'})
            self.assertEqual(param.subparam, './uid___A002_X85c183_X36f.ms')


class JyPerKWithVisTestCase(unittest.TestCase):
    """This is a test case for Jy/K with VIS data.
    """
    working_directory = 'working_directory_for_jyperk'
    data_path = 'measurementset/almasd'
    vis = 'uid___A002_X85c183_X36f.ms'
    original = f'{vis}.sel'

    @classmethod
    def setUpClass(cls):
        cls.casa_cwd_path = os.getcwd()

        if os.path.isdir(cls.working_directory):
            shutil.rmtree(cls.working_directory)

        os.mkdir(cls.working_directory)
        os.chdir(cls.working_directory)

        ms_datapath = ctsys.resolve(cls.data_path)
        original_vis = os.path.join(ms_datapath, cls.original)
        shutil.copytree(original_vis, cls.vis, symlinks=False)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.casa_cwd_path)
        shutil.rmtree(cls.working_directory)


class TestCollection4InterpolationParamsGenerator(JyPerKWithVisTestCase):
    """test collection for interpolation logic classes.

    This test target classes are:
    * InterpolationParamsGenerator
    * Bands
    * MeanElevation
    """

    def test_get_params_in_InterpolationParamsGenerator(self):
        params = jyperk.InterpolationParamsGenerator.get_params(self.vis, spw='1')
        
        param = params.__next__()
        self.assertEqual(param.param, {'date': '2014-07-01T21:49:32', 'temperature': 266.50347483801465, 
                                       'delta_days': 1000, 'antenna': 'DA61', 'elevation': 51.11212932686397, 
                                       'band': 3, 'baseband': 1, 'frequency': 90994575000.0})
        self.assertEqual(param.subparam, {'vis': 'uid___A002_X85c183_X36f.ms', 'spwid': 1})
      
    def _generate_spwnames(self):
        spwnames = []
        for i in range(3):
            for bb in range(1, 5):
                for resolution in ['FULL_RES', 'CH_AVG']:
                    spwnames.append(f'ALMA_RB_03#BB_{str(bb)}#SW-01#{resolution}')
        return spwnames

    def _generate_ref_bands(self):
        ref_bands = {}
        for i in range(1, 25):
            ref_bands[i] = 3
        return ref_bands

    def test_get_params_in_Bands(self):
        science_windows = np.array(range(1, 25))
        spwnames = self._generate_spwnames()
        mean_freqs = {1: 90994575000.0, 2: 90978950000.0, 3: 92932075000.0, 
                      4: 92924262500.0, 5: 102994575000.0, 6: 102986762500.0, 
                      7: 104994575000.0, 8: 104986762500.0, 9: 100949999999.89998, 
                      10: 100926562499.9, 11: 102765150000.0, 12: 102741712500.0, 
                      13: 112807150000.0, 14: 112783712500.0, 15: 114682150000.0, 
                      16: 114658712500.0, 17: 100949999999.89996, 18: 100949755859.275, 
                      19: 102765150000.0, 20: 102764905859.375, 21: 112807150000.0, 
                      22: 112806905859.375, 23: 114682150000.0, 24: 114681905859.375}
        bands = jyperk.Bands.get(science_windows, spwnames, mean_freqs, self.vis)
        ref_bands = self._generate_ref_bands()
        
        self.assertEqual(bands, ref_bands)

    def test_get_params_in_MeanElevation(self):
        mean_elevation = jyperk.MeanElevation.get(self.vis, 0)
        assert math.isclose(mean_elevation, 51.11212932686397, rel_tol=1e-8)


class TestInterpolationRspTranslator(JyPerKWithVisTestCase):
    response = ({'response': {'success': True, 'timestamp': '2021-10-05T05:51:20', 'elapsed': '0:00:00.461', 
               'error': None, 'query': {'elevation': '51.11212932686397', 'temperature': '266.50347483801465',
               'antenna': 'DA61', 'band': '3', 'frequency': '90994575000.0', 'date': '2014-07-01T21:49:32',
               'baseband': '1', 'delta_days': '1000'}, 'data': {'factor': {'std': 2.004, 'n': 18, 'mean':
               44.672}}}, 'aux': {'vis': './uid___A002_X85c183_X36f.ms', 'spwid': 1}}, {'response': {'success':
               True, 'timestamp': '2021-10-05T05:51:21', 'elapsed': '0:00:00.442', 'error': None, 'query':
               {'elevation': 'nan', 'temperature': '266.50347483801465', 'antenna': 'PM03', 'band': '3',
               'frequency': '90994575000.0', 'date': '2014-07-01T21:49:32', 'baseband': '1', 'delta_days':
               '1000'}, 'data': {'factor': {'std': 2.004, 'n': 18, 'mean': 44.672}}}, 'aux': {'vis':
               './uid___A002_X85c183_X36f.ms', 'spwid': 1}})
    response2 = [deepcopy(response1[0])]
    response2[0]['aux'] = 'some string'
    response3 = [deepcopy(response1[0])]
    response3[0]['aux'].pop('vis')
    response4 = [deepcopy(response1[0])]
    response4[0]['aux'].pop('spwid')
    response5 = [deepcopy(response1[0])]
    response5[0]['aux']['spwid'] = 'na'
    response6 = [deepcopy(response1[0])]
    response6[0]['aux']['vis'] = None
    factors = (['uid___A002_X85c183_X36f.ms', 'DA61', '1', 'I', '44.672'],
               ['uid___A002_X85c183_X36f.ms', 'PM03', '1', 'I', '44.672'])

    def test_convert(self):
        converted = jyperk.InterpolationRspTranslator.convert(self.response, self.vis)
        self.assertEqual(converted, self.factors)

    def test_convert_as_aux1(self):
        """Check the response.aux in the JSON obtained from Jy/K is dict."""
        with self.assertRaises(TypeError) as cm:
            converted = jyperk.InterpolationRspTranslator.convert(self.response2, self.vis)

        self.assertEqual(cm.exception.args[0],
                         'The response.aux in the JSON obtained from Jy/K db must be dict.')

    def test_convert_as_aux2(self):
        """Check the response.aux in the JSON obtained from Jy/K contain vis."""
        with self.assertRaises(KeyError) as cm:
            converted = jyperk.InterpolationRspTranslator.convert(self.response3, self.vis)

        self.assertEqual(cm.exception.args[0],
                         'The response.aux in the JSON obtained from Jy/K db must contain vis.')

    def test_convert_as_aux3(self):
        """Check the response.aux in the JSON obtained from Jy/K contain spwid."""
        with self.assertRaises(KeyError) as cm:
            converted = jyperk.InterpolationRspTranslator.convert(self.response4, self.vis)

        self.assertEqual(cm.exception.args[0],
                         'The response.aux in the JSON obtained from Jy/K db must contain spwid.')

    def test_convert_as_spwid(self):
        """Check the response.aux.spwid in the JSON obtained from Jy/K is int."""
        with self.assertRaises(TypeError) as cm:
            converted = jyperk.InterpolationRspTranslator.convert(self.response5, self.vis)

        self.assertEqual(cm.exception.args[0],
                         'The response.aux.spwid in the JSON obtained from Jy/K db must be int.')

    def test_convert_as_vis(self):
        """Check the response.aux.vis in the JSON obtained from Jy/K is str."""
        with self.assertRaises(TypeError) as cm:
            converted = jyperk.InterpolationRspTranslator.convert(self.response6, self.vis)

        self.assertEqual(cm.exception.args[0],
                         'The response.aux.vis in the JSON obtained from Jy/K db must be str.')


class TestModelFitRspTranslator(JyPerKWithVisTestCase):
    response = [{'response': {'success': True, 'timestamp': '2021-10-05T06:53:58', 'elapsed': '0:00:00.015',
               'error': None, 'query': {'elevation': '51.11212932686397', 'temperature': '266.50347483801465',
               'antenna': 'DA61', 'band': '3', 'frequency': '90994575000.0', 'date': '2014-07-01T21:49:32',
               'baseband': '1'}, 'data': {'factor': 43.719}}, 'aux': {'vis': './uid___A002_X85c183_X36f.ms',
               'spwid': 1}}, {'response': {'success': True, 'timestamp': '2021-10-05T06:53:59', 'elapsed':
               '0:00:00.004', 'error': None, 'query': {'elevation': 'nan', 'temperature': '266.50347483801465',
               'antenna': 'PM03', 'band': '3', 'frequency': '90994575000.0', 'date': '2014-07-01T21:49:32',
               'baseband': '1'}, 'data': {'factor': 43.719}}, 'aux': {'vis': './uid___A002_X85c183_X36f.ms',
               'spwid': 1}}]
    factors = [['uid___A002_X85c183_X36f.ms', 'DA61', '1', 'I', '43.719'],
               ['uid___A002_X85c183_X36f.ms', 'PM03', '1', 'I', '43.719']]

    def test_convert(self):
        converted = jyperk.ModelFitRspTranslator.convert(self.response, self.vis)
        self.assertEqual(converted, self.factors)


class TestModelFitParamsGenerator(JyPerKWithVisTestCase):
    """test ModelFitParamsGenerator class.
    """

    def test_get_params(self):
        params = jyperk.ModelFitParamsGenerator.get_params(self.vis, spw='1')
        
        param = params.__next__()
        self.assertEqual(param.param, {'date': '2014-07-01T21:49:32', 'temperature': 266.50347483801465,
                                       'antenna': 'DA61', 'elevation': 51.11212932686397, 'band': 3, 
                                       'baseband': 1, 'frequency': 90994575000.0})
        self.assertEqual(param.subparam, {'vis': 'uid___A002_X85c183_X36f.ms', 'spwid': 1})


class TestRequestsManager(unittest.TestCase):
    """test RequestsManager class.
    """

    vis = 'uid___A002_Xb32033_X9067.ms'
    
    params = MagicMock()
    params.param.return_value = {'uid': 'uid://A002/Xb32033/X9067'}

    @patch('casatasks.private.jyperk.urlopen')
    def test_get_with_success(self, urlopen_patch):
        content_body = '''
        {"success": true, "data": "dummy"}
        '''

        mock = MagicMock()
        mock.__enter__.return_value.read.return_value = content_body
        urlopen_patch.return_value = mock
        
        client = jyperk.JyPerKDatabaseClient('asdm')
        manager = jyperk.RequestsManager(client)
        
        params = jyperk.ASDMParamsGenerator.get_params(self.vis)
        result = manager.get(params)
        
        reference = [{'response': {'success': True, 'data': 'dummy'}, 'aux': 'uid___A002_Xb32033_X9067.ms'}]
        self.assertEqual(result, reference)
        self.assertTrue(urlopen_patch.called)
        self.assertTrue(urlopen_patch.call_count == 1)

    @patch('casatasks.private.jyperk.urlopen')
    def test_get_without_success(self, urlopen_patch):
        content_body = '''
        {"success": false, "data": "dummy", "error": "dummy"}
        '''

        mock = MagicMock()
        mock.__enter__.return_value.read.return_value = content_body
        urlopen_patch.return_value = mock
        
        client = jyperk.JyPerKDatabaseClient('asdm')
        manager = jyperk.RequestsManager(client)
        
        params = jyperk.ASDMParamsGenerator.get_params(self.vis)
        result = manager.get(params)
        
        reference = []
        self.assertEqual(result, reference)
        self.assertTrue(urlopen_patch.called)
        self.assertTrue(urlopen_patch.call_count == 1)


class TestJyPerKDatabaseClient(unittest.TestCase):
    """test JyPerKDatabaseClient class.
    """
    content_body = '''{"success": true}'''
    param = {'uid': 'uid://A002/X85c183/X36f'}
    
    @patch('casatasks.private.jyperk.urlopen')
    def test_get_as_success(self, urlopen_patch):
        mock = MagicMock()
        mock.__enter__.return_value.read.return_value = self.content_body
        urlopen_patch.return_value = mock
        
        client = jyperk.JyPerKDatabaseClient('asdm')
        json_obj = client.get(self.param)
              
        self.assertEqual(json_obj, json.loads(self.content_body))
        self.assertTrue(urlopen_patch.called)
        self.assertTrue(urlopen_patch.call_count == 1)

    @patch('casatasks.private.jyperk.urlopen')
    def test_get_as_httperror(self, urlopen_patch):
        urlopen_patch.side_effect = HTTPError('', 500, '', {}, None)
        
        client = jyperk.JyPerKDatabaseClient('asdm', retry=1, retry_wait_time=0.1)
        
        with self.assertRaises(RuntimeError) as cm:
            client.get(self.param)
       
        msg = 'Failed to load URL: https://asa.alma.cl/science/jy-kelvins/asdm/?uid=uid%3A%2F%2FA002%2FX85c183%2FX36f'
        self.assertEqual(cm.exception.args[0].split('\n')[0], msg)
        self.assertTrue(urlopen_patch.called)
        self.assertTrue(urlopen_patch.call_count == 1)

    @patch("casatasks.private.jyperk.urlopen")
    def test_get_as_500_500_200(self, urlopen_patch):
        error = HTTPError('', 500, '', {}, None)
        
        mock = MagicMock()
        mock.__enter__.return_value.read.side_effect = [error, error, self.content_body]
        urlopen_patch.return_value = mock
        
        client = jyperk.JyPerKDatabaseClient('asdm', retry=3, retry_wait_time=0.1)
        json_obj = client.get(self.param)
              
        self.assertEqual(json_obj, json.loads(self.content_body))
        self.assertTrue(urlopen_patch.called)
        self.assertTrue(urlopen_patch.call_count == 3)

    @patch('casatasks.private.jyperk.urlopen')
    def test_get_as_urlerror(self, urlopen_patch):
        urlopen_patch.side_effect = URLError('')
        
        client = jyperk.JyPerKDatabaseClient('asdm', retry=1, retry_wait_time=0.1)
        
        with self.assertRaises(RuntimeError) as cm:
            client.get(self.param)
       
        msg = 'Failed to load URL: https://asa.alma.cl/science/jy-kelvins/asdm/?uid=uid%3A%2F%2FA002%2FX85c183%2FX36f'
        self.assertEqual(cm.exception.args[0].split('\n')[0], msg)
        self.assertTrue(urlopen_patch.called)
        self.assertTrue(urlopen_patch.call_count == 1)


class TestTranslator(JyPerKWithVisTestCase):
    responsed_factors = [{'Antenna': 'DA61', 'Spwid': 17, 'origSpwid': 20, 'Polarization': 
                'Polarization_0', 'MS': 'uid___A002_X85c183_X36f.ms', 'factor': 43.768}, 
                {'Antenna': 'DA61', 'Spwid': 19, 'origSpwid': 22, 'Polarization': 
                'Polarization_0', 'MS': 'uid___A002_X85c183_X36fa.ms', 'factor': 43.776},
                {'Antenna': 'DA61', 'Spwid': 21, 'origSpwid': 24, 'Polarization':
                'Polarization_0', 'MS': 'uid___A002_X85c183_X36f.ms', 'factor': 43.824}]

    factors = [['uid___A002_X85c183_X36f.ms', 'DA61', '17', 'I', '43.768'],
              ['uid___A002_X85c183_X36fa.ms', 'DA61', '19', 'I', '43.776'],
              ['uid___A002_X85c183_X36f.ms', 'DA61', '21', 'I', '43.824']]

    def test_format_cal_table_format(self):
        factors = jyperk.Translator.format_cal_table_format(self.responsed_factors)
        self.assertEqual(factors, self.factors)
        
    def test_filter_jyperk_by_vis(self):
        selected = jyperk.Translator.filter_jyperk(self.vis, self.factors, '*')
        self.assertEqual(selected, [self.factors[0], self.factors[2]])

    def test_filter_jyperk_by_vis_and_spw(self):
        selected = jyperk.Translator.filter_jyperk(self.vis, self.factors, '17')
        self.assertEqual(selected, [self.factors[0]])


class TestASDMRspTranslator(JyPerKWithVisTestCase):
    response = [{'response': {'success': True, 'timestamp': '2021-10-05T03:46:16',
                'elapsed': '0:00:09.795','error': None, 'query': {'uid': 'uid://A002/X85c183/X36f'},
                'data': {'length': 12, 'factors': [{'Antenna': 'DA61', 'Spwid': 17, 'origSpwid': 20,
                'Polarization': 'Polarization_0', 'MS': 'uid___A002_X85c183_X36f.ms','factor': 43.768},
                {'Antenna': 'DA61','Spwid': 19, 'origSpwid': 22, 'Polarization': 'Polarization_0',
                'MS': 'uid___A002_X85c183_X36f.ms', 'factor': 43.776}]}}}]
    factors = [['uid___A002_X85c183_X36f.ms', 'DA61', '17', 'I', '43.768'],
               ['uid___A002_X85c183_X36f.ms', 'DA61', '19', 'I', '43.776']]

    def test_convert(self):
        converted = jyperk.ASDMRspTranslator.convert(self.response, self.vis)
        self.assertEqual(converted, self.factors)


class TestJyPerKReader4File(unittest.TestCase):
    """test TestJyPerKReader class.
    """
    working_directory = 'working_directory_for_jyperk'
    jyperk_factor_path = 'jyperk_factor.csv'

    @classmethod
    def setUpClass(cls):
        cls.casa_cwd_path = os.getcwd()

        if os.path.isdir(cls.working_directory):
            shutil.rmtree(cls.working_directory)

        os.mkdir(cls.working_directory)
        os.chdir(cls.working_directory)

        cls._generate_jyperk_factor_csv()

        ms_datapath = ctsys.resolve('measurementset/almasd')

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.casa_cwd_path)
        shutil.rmtree(cls.working_directory)

    @classmethod
    def _generate_jyperk_factor_csv(cls):
        jyperk_factor_content = """MS,Antenna,Spwid,Polarization,Factor
uid___A002_X85c183_X36f.ms,DA61,17,I,51.890035198
uid___A002_X85c183_X36f.ms,PM03,17,I,51.890035198
uid___A002_X85c183_X36f.ms,PM04,17,I,51.890035198
uid___A002_X85c183_X60b.ms,DA61,17,I,51.890035198
uid___A002_X85c183_X60b.ms,PM03,17,I,51.890035198
uid___A002_X85c183_X60b.ms,PM04,17,I,51.890035198
uid___A002_X8602fa_X2ab.ms,PM02,17,I,51.4634859397
uid___A002_X8602fa_X2ab.ms,PM03,17,I,51.4634859397
uid___A002_X8602fa_X2ab.ms,PM04,17,I,51.4634859397"""

        with open(cls.jyperk_factor_path, 'w') as fp:
            fp.write(jyperk_factor_content)
        
    def test_get(self):
        jyperk_factor_list = [['uid___A002_X85c183_X36f.ms', 'DA61', '17', 'I', '51.890035198'],
                             ['uid___A002_X85c183_X36f.ms', 'PM03', '17', 'I', '51.890035198'],
                             ['uid___A002_X85c183_X36f.ms', 'PM04', '17', 'I', '51.890035198'],
                             ['uid___A002_X85c183_X60b.ms', 'DA61', '17', 'I', '51.890035198'],
                             ['uid___A002_X85c183_X60b.ms', 'PM03', '17', 'I', '51.890035198'],
                             ['uid___A002_X85c183_X60b.ms', 'PM04', '17', 'I', '51.890035198'],
                             ['uid___A002_X8602fa_X2ab.ms', 'PM02', '17', 'I', '51.4634859397'],
                             ['uid___A002_X8602fa_X2ab.ms', 'PM03', '17', 'I', '51.4634859397'],
                             ['uid___A002_X8602fa_X2ab.ms', 'PM04', '17', 'I', '51.4634859397']]
        f = jyperk.JyPerKReader4File(self.jyperk_factor_path)

        self.assertEqual(f.get(), jyperk_factor_list)


if __name__ == '__main__':
    unittest.main()