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


class TestInterpolationParamsGenerator(JyPerKWithVisTestCase):
    """test ASDMParamsGenerator class.
    """        
    def test_get_params(self):
        params = jyperk.InterpolationParamsGenerator.get_params(self.vis, spw='1')
        
        param = params.__next__()
        self.assertEqual(param.param, {'date': '2014-07-01T21:49:32', 'temperature': 266.50347483801465, 
                                       'delta_days': 1000, 'antenna': 'DA61', 'elevation': 51.11212932686397, 
                                       'band': 3, 'baseband': 1, 'frequency': 90994575000.0})
        self.assertEqual(param.subparam, {'vis': 'uid___A002_X85c183_X36f.ms', 'spwid': 1})


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