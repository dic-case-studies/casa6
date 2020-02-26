import collections
import glob
import itertools
import numpy
import os

try:
    from casatools import msmetadata
    from casatools import table
    from casatools import measures
    from casatools import quanta
    from casatools import ms as mstool
    from casatools import image
except Exception:
    from __casac__.msmetadata import msmetadata
    from __casac__.table import table
    from __casac__.measures import measures
    from __casac__.quanta import quanta
    from __casac__.ms import ms as mstool
    from __casac__.image import image

MetaDataSet = collections.namedtuple(
    'MetaDataSet',
    ['msmeta', 'ephemmeta']
)

MsMeta = collections.namedtuple(
    'MsMeta',
    ['positions', 'times', 'freqmin', 'freqmax']
)

EphemMeta = collections.namedtuple(
    'EphemMeta',
    ['table', 'times', 'unit_time', 'velocities', 'unit_vel', 'frame_vel']
)

FrequencyRange = collections.namedtuple(
    'FrequencyRange',
    ['min', 'max', 'ref']
)

COMMON_VELOCITY_UNIT = 'km/s'

DEBUG = False


def debug_print(msg):
    if DEBUG:
        for m in msg.split('\n'):
            print('DEBUG: {}'.format(m))


def inspect_ms(vis, fieldid, spwid, chanstart=0, nchan=-1):
    """Inspect MS

    Arguments:
        vis {str} -- name of MS
        fieldid {int} -- FIELD_ID for target
        spwid {int} -- SPW_ID for target spw
        chanstart {int} -- start channel (default: 0)
        nchan {int} -- number of channels
                       (default: -1 => from chanstart to end channel of spw)

    Returns:
        MsMeta -- namedtuple containing metadata
                      positions: list of antenna positions
                      times: list of timestamps for target field
                      freqmin: minimum frequency of given spw as measure
                      freqmax: maximum frequency of given spw as measure
    """
    msmd = msmetadata()
    msmd.open(vis)
    try:
        positions = list(map(msmd.antennaposition, msmd.antennaids()))
        chanfreqs = msmd.chanfreqs(spwid)
        chanwidths = msmd.chanwidths(spwid)
    finally:
        msmd.close()

    cw = chanwidths.mean()
    start = chanstart
    end = min(chanstart + nchan, len(chanfreqs)) if nchan >= 0 else len(chanfreqs)
    debug_print('start {}, end {}'.format(start, end))
    assert 0 <= start
    assert 0 < end
    assert start < end
    freqmin = chanfreqs[start:end].min() - cw / 2
    freqmax = chanfreqs[start:end].max() + cw / 2
    debug_print('freqmin {}, freqmax {}'.format(freqmin, freqmax))

    ms = mstool()
    ms.open(vis)
    try:
        ms.msselect({'spw': str(spwid), 'field': str(fieldid), 'scanintent': 'OBSERVE_TARGET#ON_SOURCE'})
        data = ms.getdata(['time'])
        times = data['time']
    finally:
        ms.close()

    tb = table()
    try:
        tb.open(os.path.join(vis, 'SPECTRAL_WINDOW'))
        freq_ref_id = tb.getcell('MEAS_FREQ_REF', spwid)
    finally:
        tb.close()

    me = measures()
    codes = me.listcodes(me.frequency())
    refmap = codes['normal']
    assert 0 <= freq_ref_id and freq_ref_id < len(refmap)
    freq_ref_str = refmap[freq_ref_id]

    qa = quanta()

    metadata = MsMeta(
        positions,
        times,
        me.frequency(rf=freq_ref_str, v0=qa.quantity(freqmin, 'Hz')),
        me.frequency(rf=freq_ref_str, v0=qa.quantity(freqmax, 'Hz'))
    )

    me.done()

    return metadata


def get_ephem_table(vis, fieldid):
    """Return a name of the Ephemeris table corresponding to given FIELD_ID

    Arguments:
        vis {str} -- name of the MS
        fieldid {int} -- FIELD_ID

    Returns:
        str -- name of the Ephemeris table
    """
    tb = table()
    field_table = os.path.join(vis, 'FIELD')
    tb.open(field_table)
    try:
        ephem_id = tb.getcell('EPHEMERIS_ID', fieldid)
    finally:
        tb.close()

    pattern = os.path.join(field_table, 'EPHEM{}*'.format(ephem_id))
    candidates = glob.glob(pattern)
    assert len(candidates) == 1
    ephem_table = candidates[0]

    return ephem_table


def inspect_ephem(name):
    """inspect Ephemeris table

    Arguments:
        name {str} -- name of Ephemeris table

    Returns:
        EphemMeta -- data of Ephemeris table
                         time: time list
                         unit_time: unit of the time
                         velocities: velocity list
                         unit_vel: unit of the velocity
                         frame_vel: reference frame of the velocity
    """
    tb = table()
    tb.open(name)
    try:
        eph_time = tb.getcol('MJD')
        eph_time_unit = tb.getcolkeyword('MJD', 'UNIT')
        eph_radvel = tb.getcol('RadVel')
        eph_radvel_unit = tb.getcolkeyword('RadVel', 'UNIT')
        eph_geo_dist = tb.getkeyword('GeoDist')
    finally:
        tb.close()

    # Logic borrowed from FTMachine::initSourceFreqConv
    # The eph_geo_dist is a distrance from the GEOCENTER in km.
    # If eph_geo_dist > 1e-3 km (=1m), velocity reference
    # frame is regarded as TOPO (TOPOCENTRIC).
    # Otherwise, the frame should be GEO (GEOCENTRIC).
    qa = quanta()
    eph_geo_dist = qa.quantity(eph_geo_dist, 'km')
    eph_geo_threshold = qa.quantity(1.0e-3, 'km')
    if qa.gt(eph_geo_dist, eph_geo_threshold):
        eph_vel_frame = 'TOPO'
    else:
        eph_vel_frame = 'GEO'

    data = EphemMeta(
        table=name,
        times=eph_time,
        unit_time=eph_time_unit,
        velocities=eph_radvel,
        unit_vel=eph_radvel_unit,
        frame_vel=eph_vel_frame
    )

    return data


def update_measure(measures_instance, position=None, epoch=None, direction=None, comet_table=None):
    if position is not None:
        measures_instance.doframe(position)

    if epoch is not None:
        measures_instance.doframe(epoch)

    if direction is not None:
        measures_instance.doframe(direction)

    if comet_table is not None:
        assert isinstance(comet_table, str)
        assert os.path.exists(comet_table)
        measures_instance.framecomet(comet_table)
        measures_instance.doframe(measures_instance.direction('COMET'))

    return measures_instance


def get_doppler(measure_instance, radial_velocity, velocity_unit, velocity_frame):
    qa = quanta()
    vel = qa.convert(qa.quantity(radial_velocity, velocity_unit), COMMON_VELOCITY_UNIT)
    if velocity_frame == 'GEO':
        # relative velocity between GEO and TOPO must be subtracted
        radvel_zero = measure_instance.measure(v=measure_instance.radialvelocity(rf='TOPO', v0=qa.quantity(0, COMMON_VELOCITY_UNIT)), rf='GEO')
        qzero = qa.convert(radvel_zero['m0'], COMMON_VELOCITY_UNIT)
        debug_print('velocity in GEO frame. Require conversion to TOPO.')
        debug_print('Original Velocity: {value} {unit}'.format(**vel))
        debug_print('Delta Velocity: {value} {unit}'.format(**qzero))
        vel = qa.sub(vel, radvel_zero['m0'])
    debug_print('TOPO velocity: {value} {unit}'.format(**vel))
    doppler = measure_instance.doppler(rf='RELATIVISTIC', v0=vel)

    return doppler


def ms_freq_range(metadataset):
    msmeta = metadataset.msmeta
    ephem_data = metadataset.ephemmeta
    ephem_table = ephem_data.table

    qa = quanta()
    min_frequency = None
    max_frequency = None
    eph_time = qa.convert(qa.quantity(ephem_data.times, ephem_data.unit_time), 's')['value']
    eph_vel = ephem_data.velocities
    interpolated_velocities = numpy.interp(msmeta.times, eph_time, eph_vel)
    for item in itertools.product(msmeta.positions, zip(msmeta.times, interpolated_velocities)):
        position = item[0]
        timestamp = item[1][0]
        velocity = item[1][1]
        me = measures()
        me.done()
        epoch = me.epoch('UTC', qa.quantity(timestamp, 's'))
        me = update_measure(me, epoch=epoch, position=position, comet_table=ephem_table)
        doppler = get_doppler(me, velocity, ephem_data.unit_vel, ephem_data.frame_vel)
        fmin, fmax = map(lambda x: me.torestfrequency(x, doppler), [msmeta.freqmin, msmeta.freqmax])
        assert fmin['refer'] == 'REST'
        assert fmax['refer'] == 'REST'
        if min_frequency is None or qa.lt(fmin['m0'], min_frequency) is True:
            min_frequency = fmin['m0']
        if max_frequency is None or qa.gt(fmax['m0'], max_frequency) is True:
            max_frequency = fmax['m0']
        debug_print('min freq: {}'.format(qa.tos(min_frequency)))
        debug_print('max freq: {}'.format(qa.tos(max_frequency)))

    return FrequencyRange(min_frequency, max_frequency, 'REST')


def image_freq_range(imagename):
    ia = image()
    ia.open(imagename)
    csys = ia.coordsys()
    try:
        imshape = ia.shape()
        chmin = -0.5
        chmax = imshape[3] - 1 + 0.5
        refpix = [0, 0, 0, chmin]
        wmin = csys.toworld(refpix, format='m')
        refpix[3] = chmax
        wmax = csys.toworld(refpix, format='m')
    finally:
        csys.done()
        ia.close()

    refmin = wmin['measure']['spectral']['frequency']['refer']
    refmax = wmax['measure']['spectral']['frequency']['refer']
    assert refmin == refmax
    fmin = wmin['measure']['spectral']['frequency']['m0']
    fmax = wmax['measure']['spectral']['frequency']['m0']

    return FrequencyRange(fmin, fmax, refmax)


def get_lorentz_factor(metadataset):
    msmeta = metadataset.msmeta
    ephem_data = metadataset.ephemmeta
    ephem_table = ephem_data.table

    me = measures()
    me.done()
    qa = quanta()
    epoch = me.epoch('UTC', qa.quantity(msmeta.times[0], 's'))
    position = me.observatory('ALMA')
    me = update_measure(me, epoch=epoch, position=position, comet_table=ephem_table)
    velocity = ephem_data.velocities[0]
    unit = ephem_data.unit_vel
    frame = ephem_data.frame_vel
    doppler = get_doppler(me, velocity, unit, frame)
    refvel = qa.convert(doppler['m0'], COMMON_VELOCITY_UNIT)
    speed_of_light = qa.convert(qa.constants('c'), COMMON_VELOCITY_UNIT)
    return qa.div(refvel, speed_of_light)['value']


def frequency_value(freq):
    if isinstance(freq, (int, float)):
        val = freq
    elif isinstance(freq, dict) and 'value' in freq:
        val = freq['value']
    else:
        val = None

    return val


def is_frequency_close(freq1, freq2, lorentz_factor, rtol=1e-1):
    val1 = frequency_value(freq1)
    val2 = frequency_value(freq2)
    tolerance = abs(lorentz_factor * rtol)
    reldiff = abs((val2 - val1) / val1)
    debug_print('values: {} {}'.format(val1, val2))
    debug_print('tolerance: {} (factor {})'.format(tolerance, lorentz_factor))
    debug_print('relative diff: {}'.format(reldiff))
    return reldiff <= tolerance


def get_metadataset(vis, fieldid, spwid, chanstart=0, nchan=-1):
    msmeta = inspect_ms(vis, fieldid, spwid)

    ephem_table = get_ephem_table(vis, fieldid)

    ephem_data = inspect_ephem(ephem_table)

    return MetaDataSet(msmeta, ephem_data)
