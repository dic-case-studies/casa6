import os
import shutil
from taskinit import gentools

# A function object that can be passed to ignore parameter
# of shutil.copytree. It will ignore subversion directory
# when data are copied to working directory.
ignore_subversion = shutil.ignore_patterns('.svn')

def copytree_ignore_subversion(datadir, name, outname=None):
    if outname is None:
        outname = name
    if not os.path.exists(name):
        shutil.copytree(os.path.join(datadir, name), outname,
                        ignore=ignore_subversion)


def get_table_cache():
    (mytb,) = gentools(['tb'])
    cache = mytb.showcache()
    #print 'cache = {}'.format(cache)
    return cache


class TableCacheValidator(object):
    def __init__(self):
        self.original_cache = get_table_cache()
        
    def validate(self):
        cache = get_table_cache()
        #print 'original {} current {}'.format(self.original_cache, cache)
        return len(cache) == 0 or cache == self.original_cache
