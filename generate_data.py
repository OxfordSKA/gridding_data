#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import math
import os
import sys
import time

import numpy

import ephem
import oskar


def obs_params(longitude_deg, latitude_deg, mid_azimuth_deg, mid_elevation_deg,
               mid_date_string, obs_length_sec):
    """Evaluate the observation parameters (RA, Dec, MJD start) for the
    requested azimuth, elevation and time.

    Args:
        longitude_deg (float):     Telescope longitude, in degrees.
        latitude_deg (float):      Telescope latitude, in degrees.
        mid_azimuth_deg (float):   Azimuth at mid-point, in degrees.
        mid_elevation_deg (float): Elevation at mid-point, in degrees.
        mid_date_string (str):     Date and time of mid-point, as a string.
        obs_length_sec (float):    Target observation length, in seconds.

    Returns:
        Tuple containing RA, Dec and MJD start for the given parameters.
    """
    obs = ephem.Observer()
    obs.lon, obs.lat, obs.elevation = \
        math.radians(longitude_deg), math.radians(latitude_deg), 0.0
    obs.date = mid_date_string
    ra, dec = obs.radec_of(math.radians(mid_azimuth_deg),
                           math.radians(mid_elevation_deg))
    ra, dec = math.degrees(ra), math.degrees(dec)
    mjd_mid = ephem.julian_date(obs.date) - 2400000.5
    mjd_start = mjd_mid - obs_length_sec / (2 * 86400.0)
    return ra, dec, mjd_start


def main():
    # Global options.
    precision = b'double'
    telescope_model = os.path.join('models', 'SKA1_LOW_v6.tm')
    global_freq_start_hz = 42.5e6
    bandwidth_hz = 5e3
    dump_time_sec = 0.9
    obs_length_sec = 45.0
    num_times = 50

    # Get channel start ID and number of channels from command line arguments.
    # These are set in the job submission script, which runs this repeatedly.
    if len(sys.argv) != 3:
        raise RuntimeError('Usage: ./generate_data.py '
                           '<start channel> <num channels>')
    start_channel = int(sys.argv[1])
    num_channels = int(sys.argv[2])
    freq_start_hz = bandwidth_hz * start_channel + global_freq_start_hz
    output_root = b'channel_data_%04d-%04d' % (start_channel,
                                               start_channel + num_channels - 1)

    # Calculate RA, Dec and start MJD for zenith pointing.
    position = numpy.loadtxt(
        os.path.join(telescope_model, 'position.txt'), delimiter=b', ')
    ra, dec, mjd_start = obs_params(
        position[0], position[1], 0.0, 90.0, b'2016/10/01 00:00',
        obs_length_sec)

    # Set up sky model.
    sky = oskar.Sky.generate_grid(ra, dec, side_length=80, fov_deg=40.0,
                                  precision=precision)

    # Set up telescope model.
    tel = oskar.Telescope(precision)
    tel.set_channel_bandwidth(bandwidth_hz)
    tel.set_time_average(dump_time_sec)
    tel.set_pol_mode(b'Scalar')
    tel.load(telescope_model)
    tel.set_phase_centre(ra, dec)

    # Set up imagers.
    # imagers = []
    # for i in range(1):
    #     imagers.append(oskar.Imager(precision))
    #     imagers[i].set(fov_deg=8.0, image_size=8192, algorithm='W-projection')
    #     imagers[i].set(channel_snapshots=False, time_snapshots=False)
    #     imagers[i].set(output_root=output_root)

    # Set up the simulator.
    # simulator = oskar.ImagingSimulator(imagers, precision)
    simulator = oskar.Simulator(precision)
    simulator.set_sky_model(sky)
    simulator.set_telescope_model(tel)
    simulator.set_observation_frequency(freq_start_hz, bandwidth_hz,
                                        num_channels)
    simulator.set_observation_time(mjd_start, obs_length_sec, num_times)
    simulator.set_output_vis_file(output_root + '.vis')
    simulator.set_max_times_per_block(2)
    simulator.set_horizon_clip(False)

    # Run.
    print('Running simulation...')
    start = time.time()
    simulator.run()
    print('Completed after %.3f seconds' % (time.time() - start))

if __name__ == '__main__':
    main()
