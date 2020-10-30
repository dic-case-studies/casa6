import os
import shutil
from taskinit import gentools


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
