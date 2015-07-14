#!/usr/bin/python
"""Calculate phase centre directions.

Usage:
    casapy --nogui --nologger --log2term -c calculate_pointing_directions.py
"""

import numpy as np
import os
import numpy.random as rand


if __name__ == '__main__':
    # -------------------------------------------------------------------------
    # Telescope position
    # VLA
    lon0 = -107.6184
    lat0 = 34.0790

    alt = 0.0
    # Initial start datetime (note: hour and minute is randomly generated)
    t0_year = 2015
    t0_month = 7
    t0_day = 10
    t0_sec = 0
    # Minumum elevation at the higest point.
    min_el = 75
    # Number of pointings to generate.
    num_pointings = 10
    outfile = os.path.join('pointings.txt')
    seed = 1
    # -------------------------------------------------------------------------

    fh = file(outfile, 'w')
    fh.write('# %s\n' % ('-' * 70))
    fh.write('# + %i pointings for telescope at:\n' % num_pointings)
    fh.write('#   - lon0 = %.7f deg.\n' % lon0)
    fh.write('#   - lat0 = %.7f deg.\n' % lat0)
    fh.write('# + Minimum elevation = %.2f deg.\n' % min_el)
    fh.write('# + Seed = %i\n' % seed)
    fh.write('# + Note: Times, az, and el are for the mid point of the '
             'observation.\n')
    fh.write('# %s\n' % ('-' * 70))
    fh.write('# RA, Dec, CASA date-time, MJD, az, el\n')

    rand.seed(seed)
    telescope_position = me.position('WGS84',
                                     qa.quantity(lon0, 'deg'),
                                     qa.quantity(lat0, 'deg'),
                                     qa.quantity(alt, 'm'))

    for p in range(0, num_pointings):

        # Random time for the centre of the observation.
        t0_hour = rand.randint(0, 24)
        t0_min = rand.randint(0, 60)
        t0_utc = '%04i/%02i/%02i/%02i:%02i:%05.2f' % \
                 (t0_year, t0_month, t0_day, t0_hour, t0_min, t0_sec)
        epoch0 = me.epoch('UTC', t0_utc)
        t0_mjd_utc = epoch0['m0']['value']
        time0 = qa.quantity(t0_mjd_utc, unitname='d')

        # Az,El of the phase centre at highest point.
        az0 = rand.randint(0, 2) * 180.0
        el0 = rand.randint(min_el, 90)
        direction0 = me.direction('AZEL', qa.quantity(az0, 'deg'),
                                  qa.quantity(el0, 'deg'))

        # Attach telescope position and time to measure
        me.doframe(epoch0)
        me.doframe(telescope_position)

        # Convert position to celestial coordinates.
        direction_radec = me.measure(direction0, 'J2000')
        ra = direction_radec['m0']['value'] * (180. / np.pi)
        dec = direction_radec['m1']['value'] * (180. / np.pi)

        fh.write('% 15.10f % 15.10f %25s %.11f % 9.2f % 9.2f\n' %
                 (ra, dec, t0_utc, t0_mjd_utc, az0, el0))
    fh.close()
