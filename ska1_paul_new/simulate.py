#!../venv/bin/python
import os
from os.path import join
import collections
import numpy
import subprocess


class OskarSettings(object):
    """Interface between OSKAR settings ini files with Python dictionaries"""

    def __init__(self):
        self.ini_separator = '/'
        self.ini = ''

    def dict_to_ini(self, settings, ini_file, overwrite=False):
        """Convert dictionary of settings to ini file."""
        self.ini = os.path.abspath(ini_file)
        if overwrite and os.path.exists(self.ini):
            print "Removing ini file", self.ini
            os.remove(self.ini)
        ini_dir = os.path.dirname(self.ini)
        if not ini_dir == "" and not os.path.isdir(ini_dir):
            os.makedirs(ini_dir)
        for key in sorted(settings):
            if isinstance(settings[key], dict):
                self.__iterate_settings(settings[key],
                                        key + self.ini_separator)

    def __iterate_settings(self, node, group):
        for key in sorted(node):
            if isinstance(node[key], dict):
                self.__iterate_settings(node[key],
                                        group + key + self.ini_separator)
            else:
                self.__set_setting(group + key, node[key])

    def __set_setting(self, key, value):
        subprocess.call(["oskar_settings_set", "-q",
                         self.ini, key, str(value)])


def create_settings(pointing_index, ini_file, snapshot_id):
    """Create settings dictionary."""
    dtype = [('ra', 'f8'), ('dec', 'f8'),
             ('date-time', 'a25'), ('mjd', 'f8'),
             ('az', 'f8'), ('el', 'f8')]
    pointing = numpy.loadtxt(join('data', 'pointings_ska1_mid.txt'),
                             dtype=dtype)
    ra = pointing['ra'][pointing_index]
    dec = pointing['dec'][pointing_index]
    start_time_mjd = pointing['mjd'][pointing_index]
    start_offset = (-3. + (20. / 60.) * snapshot_id) / 24.
    start_time_mjd += start_offset

    num_times = 200
    dump_time_s = 0.1
    obs_length = num_times * dump_time_s
    vis_path_root = 'p%02i_s%02i' % (pointing_index, snapshot_id)

    s = collections.OrderedDict()
    s['simulator'] = {
        'double_precision': 'true',
        'keep_log_file': 'false'
    }
    s['sky'] = {
        # 'oskar_sky_model/file': 'test.osm',
        'generator': {
            'random_power_law': {
                'num_sources': int(25.0e6),
                'flux_min': 1.0e-4,  # Jy
                'flux_max': 10.0,  # Jy
                'power': -2.0,
                'filter': {
                    'radius_outer_deg': 1.5
                },
            },
        },
        'output_text_file': 'sky_rpl.osm'
    }
    s['observation'] = {
        'start_frequency_hz': 700.0e6,
        'num_channels': 1,
        'start_time_utc': start_time_mjd,
        'length': obs_length,
        'num_time_steps': num_times,
        'phase_centre_ra_deg': '%.10f' % ra,
        'phase_centre_dec_deg': '%.10f' % dec
    }
    s['telescope'] = {
        'longitude_deg': 21.4429090,
        'latitude_deg': -30.7394750,
        'input_directory': join('data',
                                'ska1_meerkat_mid_combined_july_2015.tm'),
        'pol_mode': 'Scalar',
        # 'station_type': 'Isotropic beam'
        'station_type': 'Gaussian beam',
        'gaussian_beam': {
            'fwhm_deg': 1.64,
            'ref_freq_hz': 700.0e6
        }
    }
    s['interferometer'] = {
        'time_average_sec': dump_time_s,
        'channel_bandwidth_hz': 700.0e6 / (2**16),  # ~65k channels, ~10kHz
        # 'time_average_sec': 0.0,
        # 'channel_bandwidth_hz': 0.0,
        'ms_filename': vis_path_root + '.ms',
        'oskar_vis_filename': vis_path_root + '.vis'
    }
    settings = OskarSettings()
    settings.dict_to_ini(s, os.path.abspath(ini_file), overwrite=True)

if __name__ == '__main__':

    for i in range(0, 9):
        ini = 'test_%02i.ini' % i
        create_settings(pointing_index=0, ini_file=ini, snapshot_id=i)
        subprocess.call(["oskar_sim_interferometer", ini])
