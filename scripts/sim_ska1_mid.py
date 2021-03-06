# -*- coding: utf-8 -*-
"""Simulation of gridding test data using the VLA-A layout."""

import collections
import os
import numpy as np
import subprocess
import sys


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
        'max_sources_per_chunk': 1000,
        'double_precision': 'true',
        'keep_log_file': 'false'
    }
    s['sky/'] = {
        #'oskar_sky_model/file': sky_file
        'generator/grid/side_length': 11,
        'generator/grid/fov_deg': 3,
        'generator/grid/mean_flux_jy': 1
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
        'input_directory': os.path.join('models',
                                        'ska1_meerkat_mid_combined_july_2015.tm'),
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


def main(output_path):
    """."""
    dtype = [('ra', 'f8'), ('dec', 'f8'), ('date-time', 'a25'),
             ('mjd', 'f8'), ('az', 'f8'), ('el', 'f8')]
    pointing = np.loadtxt(os.path.join('models', 'vla_pointings.txt'),
                          dtype=dtype)

    # ----------------------------------------
    pointing_idx = 0
    ra0 = pointing['ra'][pointing_idx]
    dec0 = pointing['dec'][pointing_idx]
    sky_file = os.path.join(output_path, 'test.osm')
    freq_hz = 500.e6
    start_time_mjd = pointing['mjd'][pointing_idx]
    num_times = 10
    t_acc = (6. * 3600.) / float(num_times)  # seconds
    vis_name = os.path.join(output_path, 'test_ska')
    ini = os.path.join(output_path, 'test.ini')
    # ----------------------------------------

    s = create_settings(sky_file, freq_hz, start_time_mjd, t_acc, num_times,
                        ra0, dec0, vis_name)
    dict_to_ini(s, ini)

    subprocess.call(["oskar_sim_interferometer", ini])

if __name__ == '__main__':

    output_path = os.path.abspath(sys.argv[1])
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    main(output_path)
