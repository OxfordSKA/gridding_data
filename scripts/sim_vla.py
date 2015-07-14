# -*- coding: utf-8 -*-
"""Simulation of gridding test data using the VLA-A layout."""

import collections
import os
import numpy as np
import subprocess


def create_sky_model(file_name, ra, dec, stokes_i):
    """."""
    fh = open(file_name, 'w')
    for ra_, dec_, flux in zip(ra, dec, stokes_i):
        fh.write('%.14f, %.14f, %.3f\n' % (ra_, dec_, flux))
    fh.close()


def dict_to_ini(settings_dict, ini):
    """Convert a dictionary of settings to and OSKAR settings ini file."""
    ini_dir = os.path.dirname(ini)
    if not ini_dir == "" and not os.path.isdir(ini_dir):
        os.makedirs(ini_dir)
    for group in sorted(settings_dict):
        for key in sorted(settings_dict[group]):
            key_ = group + key
            value_ = settings_dict[group][key]
            subprocess.call(["oskar_settings_set", "-q", ini,
                            key_, str(value_)])


def create_settings(sky_file, freq_hz, start_time_mjd, t_acc, num_times, ra0,
                    dec0, vis_name):
    """."""
    s = collections.OrderedDict()
    s['simulator/'] = {
        'max_sources_per_chunk': 20,
        'double_precision': 'true',
        'keep_log_file': 'false'
    }
    s['sky/'] = {
        'oskar_sky_model/file': sky_file
    }
    s['observation/'] = {
        'start_frequency_hz': freq_hz,
        'num_channels': 1,
        'start_time_utc': start_time_mjd,
        'length': num_times * t_acc,
        'num_time_steps': num_times,
        'phase_centre_ra_deg': ra0,
        'phase_centre_dec_deg': dec0
    }
    s['telescope/'] = {
        'longitude_deg': -107.6184,
        'latitude_deg': 34.0790,
        'input_directory': os.path.join('models', 'VLA_A.tm'),
        'pol_mode': 'Scalar',
        'station_type': 'Isotropic beam'
    }
    s['interferometer/'] = {
        'time_average_sec': 0.0,
        'channel_bandwidth_hz': 0.0,
        'ms_filename': vis_name + '.ms',
        'oskar_vis_filename': vis_name + '.vis'
    }
    return s


def main():
    """."""
    dtype = [('ra', 'f8'), ('dec', 'f8'), ('date-time', 'a25'),
             ('mjd', 'f8'), ('az', 'f8'), ('el', 'f8')]
    pointing = np.loadtxt(os.path.join('models', 'vla_pointings.txt'),
                          dtype=dtype)

    # ----------------------------------------
    pointing_idx = 0
    ra0 = pointing['ra'][pointing_idx]
    dec0 = pointing['dec'][pointing_idx]
    sky_file = 'test.osm'
    freq_hz = 500.e6
    start_time_mjd = pointing['mjd'][pointing_idx]
    # t_acc = 5.0  # seconds
    num_times = 43
    t_acc = (6.0 * 3600.0) / float(num_times)  # seconds
    vis_name = 'test_vla'
    ini = 'test.ini'
    # ----------------------------------------

    ra = [ra0, ra0 + 0.05, ra0 + 0.18, ra0 + 0.3]
    dec = [dec0, dec0 + 0.05, dec0, dec0 - 0.1]
    flux = np.ones((len(ra),), dtype='f8')
    create_sky_model(sky_file, ra, dec, flux)
    s = create_settings(sky_file, freq_hz, start_time_mjd, t_acc, num_times,
                        ra0, dec0, vis_name)
    dict_to_ini(s, ini)

    subprocess.call(["oskar_sim_interferometer", ini])

if __name__ == '__main__':
    main()
