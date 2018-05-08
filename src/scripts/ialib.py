from CASAtools import image
from CASAtools import ctsys
from .. import casalog
from . import cvt

def write_image_history(myia, tname, param_names, param_vals, myclog=None):
    """
    Update image attached to image tool with the parameters that task tname was called with.

    myia - attached image tool or image name (string)
    tname - name of the calling task.
    param_names - list of parameter names.
    param_vals - list of parameter values (in the same order as param_names).
    myclog - a casalog instance (optional)
    """
    param_names = cvt.as_list(param_names)
    param_vals = cvt.as_list(param_vals)
    myia_is_string = type(myia) == str
    if myia_is_string:
        if not myia:
            # empty string
            return
        _ia = image()
        _ia.open(myia)
    elif not hasattr(myia, 'sethistory'):
        return False
    else:
        _ia = myia
    try:
        if not myclog and hasattr(casalog, 'post'):
            myclog = casalog
    except Exception as instance:
        # There's no logger to complain to, and I don't want to exit
        # just because of that.
        pass
    try:
        vestr = 'version: ' + ctsys.version_info( )
        _ia.sethistory(tname, vestr)
        # Write the arguments.
        s = tname + "("
        n = len(param_names)
        for argnum in range(n):
            s += str(param_names[argnum]) + "="
            val = str(param_vals[argnum])
            if type(val) == str:
                s += '"'
            s += str(val)
            if type(val) == str:
                s += '"'
            if argnum < n-1:
                s += ", "
        s += ")" 
        _ia.sethistory(tname, s)
    except Exception as instance:
        if hasattr(myclog, 'post'):
            myclog.post("*** Error \"%s\" updating HISTORY of " % (instance),
                        'SEVERE')
        return False
    finally:
        if myia_is_string:
            _ia.done()
    return True

def get_created_images(outfile, target_time):
    dirpath = os.path.dirname(outfile)
    if not dirpath:
        dirpath = "."
    base = os.path.basename(outfile)
    # get all entries in the directory w/ stats
    entries = []
    for fn in os.listdir(dirpath):
        if os.path.basename(fn).startswith(base):
            entries.append((os.stat(fn), fn))
    # leave only directories, insert creation date
    entries = ((stat.st_mtime, path)
        for stat, path in entries if S_ISDIR(stat[ST_MODE]))
    # reverse sort by time
    zz = sorted(entries)
    zz.reverse()
    created_images = []
    for mdate, path in zz:
        # kludge because all of a sudden, some mdates are less than time.time() value
        # that was gotten before these files were created on OSX. Weird.
        if mdate < target_time - 1:
            break
        if os.path.basename(path).startswith(base):
            created_images.append(path)
    return created_images
