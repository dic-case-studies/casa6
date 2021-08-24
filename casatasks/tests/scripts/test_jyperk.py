import os
import shutil
import unittest
from unitetst import mock, MagicMock, patch

from casatasks.private import jyperk


class TestJyPerKReader4File(unittest.TestCase):
    def setUp(self):
        self.test_dir = 'test_jyperk'

        os.mkdir(self.test_dir)
        self.jyperk_factor_path = os.path.join(self.test_dir, 'jyperk_factor.csv')
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

        with open(self.jyperk_factor_path, 'w') as fp:
            fp.write(jyperk_factor_content)

    def tearDown(self):
        shutil.rmtree(self.test_dir)
        
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