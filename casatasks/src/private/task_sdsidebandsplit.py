from casatasks.private.casa_transition import is_CASA6

if is_CASA6:
    from casatools import quanta, sidebandseparator

    from . import sdutil
else:
    import sdutil
    from casac import casac
    from taskinit import qatool as quanta

    def sidebandseparator():
        return casac.sidebandseparator()


@sdutil.sdtask_decorator
def sdsidebandsplit(imagename, outfile, overwrite, signalshift, imageshift,
                    getbothside, refchan, refval, useother, threshold):

    separator = sidebandseparator()
    try:
        separator.open(imagename)
        separator.setshift(signalshift, True)
        if len(imageshift) > 0:
            separator.setshift(imageshift, False)
        separator.setlimit(threshold)
        separator.setboth(getbothside)
        if getbothside:
            if refval == '':
                qrefval = -1.0
            else:
                myqa = quanta()
                qrefval = myqa.quantity(refval)
            separator.set_imageband_frequency(refchan, qrefval)
        separator.setsolveother(useother)
        separator.separate(outfile, overwrite)
    finally:
        separator.close()
