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