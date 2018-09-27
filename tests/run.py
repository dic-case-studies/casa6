###
### use like:
###
###      PYTHONPATH=build/lib.macosx-10.12-x86_64-3.6 python tests/run.py tests/tools/image/test_ia_imageconcat.py tests/tools/image/test_ia_makecomplex.py
###
import os
import sys
import unittest

###
### collect paths and test module names
###
test_paths = set( )
test_modules = [ ]
for i in sys.argv[1:]:
    ## note directory
    suite_path = os.path.dirname(i)
    test_paths.add(suite_path)
    ## discover module name
    suite_module, xxx = os.path.splitext(os.path.basename(i))
    test_modules.append(suite_module)

sys.path = list(test_paths) + sys.path

suite = unittest.defaultTestLoader.loadTestsFromNames(test_modules)
unittest.TextTestRunner().run(suite)
