from __future__ import absolute_import

def as_list( v ):
    if type(v) is list:
        return v
    else:
        try:
            iter(v)
            return [v] if type(v) is str else list(v)
        except:
            return [v]
