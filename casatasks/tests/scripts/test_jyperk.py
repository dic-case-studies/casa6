import os
import shutil
import unittest

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
    vis = 'uid___A002_X85c183_X36f.ms'

    @classmethod
    def setUpClass(cls):
        cls.casa_cwd_path = os.getcwd()

        if os.path.exists(cls.working_directory):
            shutil.rmtree(cls.working_directory)

        os.mkdir(cls.working_directory)
        os.chdir(cls.working_directory)

        ms_datapath = ctsys.resolve('measurementset/almasd')
        original_vis = os.path.join(ms_datapath, f'{cls.vis}.sel')
        shutil.copytree(original_vis, cls.vis, symlinks=False)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.casa_cwd_path)
        shutil.rmtree(cls.working_directory)


class TestCollection4InterpolationParamsGenerator(JyPerKWithVisTestCase):
    """test collection for interpolation logic classes.
    """        
    def test_get_params_in_InterpolationParamsGenerator(self):
        params = jyperk.InterpolationParamsGenerator.get_params(self.vis, spw='1')
        
        param = params.__next__()
        self.assertEqual(param.param, {'date': '2014-07-01T21:49:32', 'temperature': 266.50347483801465, 
                                       'delta_days': 1000, 'antenna': 'DA61', 'elevation': 51.11212932686397, 
                                       'band': 3, 'baseband': 1, 'frequency': 90994575000.0})
        self.assertEqual(param.subparam, {'vis': 'uid___A002_X85c183_X36f.ms', 'spwid': 1})
      
    def test_get_params_in_Bands(self):
        science_windows = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                    21, 22, 23, 24])
        spwnames = ['ALMA_RB_03#BB_1#SW-01#FULL_RES', 'ALMA_RB_03#BB_1#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_2#SW-01#FULL_RES', 'ALMA_RB_03#BB_2#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_3#SW-01#FULL_RES', 'ALMA_RB_03#BB_3#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_4#SW-01#FULL_RES', 'ALMA_RB_03#BB_4#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_1#SW-01#FULL_RES', 'ALMA_RB_03#BB_1#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_2#SW-01#FULL_RES', 'ALMA_RB_03#BB_2#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_3#SW-01#FULL_RES', 'ALMA_RB_03#BB_3#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_4#SW-01#FULL_RES', 'ALMA_RB_03#BB_4#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_1#SW-01#FULL_RES', 'ALMA_RB_03#BB_1#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_2#SW-01#FULL_RES', 'ALMA_RB_03#BB_2#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_3#SW-01#FULL_RES', 'ALMA_RB_03#BB_3#SW-01#CH_AVG',
                    'ALMA_RB_03#BB_4#SW-01#FULL_RES', 'ALMA_RB_03#BB_4#SW-01#CH_AVG']
        mean_freqs = {1: 90994575000.0, 2: 90978950000.0, 3: 92932075000.0, 
                      4: 92924262500.0, 5: 102994575000.0, 6: 102986762500.0, 
                      7: 104994575000.0, 8: 104986762500.0, 9: 100949999999.89998, 
                      10: 100926562499.9, 11: 102765150000.0, 12: 102741712500.0, 
                      13: 112807150000.0, 14: 112783712500.0, 15: 114682150000.0, 
                      16: 114658712500.0, 17: 100949999999.89996, 18: 100949755859.275, 
                      19: 102765150000.0, 20: 102764905859.375, 21: 112807150000.0, 
                      22: 112806905859.375, 23: 114682150000.0, 24: 114681905859.375}

        bands = jyperk.Bands.get(science_windows, spwnames, mean_freqs, self.vis)
        
        ref_bands = {1: 3, 2: 3, 3: 3, 4: 3, 5: 3, 6: 3, 7: 3, 8: 3, 9: 3, 10: 
                 3, 11: 3, 12: 3, 13: 3, 14: 3, 15: 3, 16: 3, 17: 3, 18: 3,
                 19: 3, 20: 3, 21: 3, 22: 3, 23: 3, 24: 3}
        self.assertEqual(bands, ref_bands)

    def test_get_params_in_MeanElevation(self):
        mean_elevation = jyperk.MeanElevation.get(self.vis, 0)
        self.assertEqual(bands, ref_bands)


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


class TestJyPerKReader4File(unittest.TestCase):
    """test TestJyPerKReader class.
    """
    working_directory = 'working_directory_for_jyperk'
    jyperk_factor_path = 'jyperk_factor.csv'

    @classmethod
    def setUpClass(cls):
        cls.casa_cwd_path = os.getcwd()

        if os.path.exists(cls.working_directory):
            shutil.rmtree(cls.working_directory)

        os.mkdir(cls.working_directory)
        os.chdir(cls.working_directory)

        cls._generate_jyperk_factor_csv()

        ms_datapath = ctsys.resolve('measurementset/almasd')

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.casa_cwd_path)
        shutil.rmtree(cls.working_directory)

    def _delete_dir(self, path):
        if os.path.exists(path):
            shutil.rmtree(path)

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