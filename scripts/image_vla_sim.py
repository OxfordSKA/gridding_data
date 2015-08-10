# -*- coding: utf-8 -*-
"""Imaging script for the VLA simulated data."""

import argparse
import numpy as np
import time
import os
import cProfile
import pstats
import sys
from parula import parula_map
import matplotlib.pyplot as plt
from read_oskar_vis import OskarVis
from make_image_standard import *
import pyfits


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_usage()
        sys.exit(2)


def grid_cell_size(cell_size_lm_arcsec, im_size):
    return (180. * 3600.) / (im_size * cell_size_lm_arcsec * np.pi)


def load_vis(vis_file):
    oskar_vis = OskarVis(vis_file)
    uu, vv, ww = oskar_vis.uvw(flatten=True)
    freq_hz = oskar_vis.frequency()
    wave_length_m = 299792458. / freq_hz
    uvw = np.transpose(np.array([uu, vv, ww]))
    uvw /= wave_length_m
    vis = oskar_vis.amplitudes(flatten=True)
    return uvw, vis


def main(vis_file, support, over_sample, fov, im_size):
    cell_size_lm_arcsec = fov_to_cell_size(fov, im_size)
    print 'cell size lm (arcsec):', cell_size_lm_arcsec

    uvw, vis = load_vis(vis_file)
    # num_vis = 1
    # uvw = np.empty((num_vis, 3), dtype='f8')
    # vis = np.ones((num_vis,), dtype='c16')
    # # uvw[0, :] = [500, 500, 0]
    # uvw[0, :] = [0, 0, 0]

    t0 = time.time()
    gcf = GriddingConvolutionFunction(support, over_sample)
    # gcf.pillbox_2d(1.0)
    gcf.exponential_2d(width=1.0)
    # gcf.exponential_sinc_2d()
    # gcf.sinc_2d()
    print 'Generate GCF = %.3f s' % (time.time() - t0)

    t0 = time.time()
    cell_size_uv_m = grid_cell_size(cell_size_lm_arcsec, im_size)
    grid = Grid(im_size, cell_size_uv_m)
    grid_data(uvw, vis, gcf, grid)
    print 'Gridding took = %.3f s' % (time.time() - t0)

    t0 = time.time()
    image = Image(grid)
    print 'FFT took = %.3f s' % (time.time() - t0)

    t0 = time.time()
    # image.grid_correct(gcf)
    print 'Grid correct took = %.3f s' % (time.time() - t0)

    return grid, image, gcf


def run_profile(vis_file, support, over_sample, fov, im_size):
    prof = cProfile.Profile()
    prof.enable()
    # Code to profile goes here
    # ------------------------------------------------
    grid, image, gcf = main(vis_file, support, over_sample, fov, im_size)
    # ------------------------------------------------
    prof.disable()
    prof_stats = pstats.Stats(prof, stream=sys.stdout)
    prof_stats.print_stats()
    return grid, image, gcf

if __name__ == '__main__':
    parser = MyParser(description='Test imaging script.')
    parser.format_help()
    plot_modes = ['none', 'grid', 'image', 'gcf']
    parser.add_argument('-p', metavar='MODE',
                        help='Space separated list of plot mode(s). Allowed '\
                             'values are: ' + ', '.join(plot_modes) + '.',
                        choices=plot_modes,
                        default=plot_modes[2],
                        nargs='+',
                        type=str)
    parser.add_argument('--prof', help='Run profile.',
                        action='store_true')
    parser.add_argument('--fits', help='Write FITS image.',
                        action='store_true')

    args = parser.parse_args()
    plot_modes = args.p if isinstance(args.p, list) else [args.p]

    # ------------------------------------------------------------------------
    vis_file = os.path.abspath(os.path.join('ska1_6h_scalar', 'test_ska.vis'))
    support = 5
    over_sample = 63
    fov = 3.5  # deg
    im_size = 4096
    # ------------------------------------------------------------------------

    if args.prof:
        grid, image, gcf = run_profile(vis_file, support, over_sample, fov,
                                       im_size)
    else:
        t0 = time.time()
        grid, image, gcf = main(vis_file, support, over_sample, fov, im_size)
        print 'Time taken: %.3f seconds' % (time.time() - t0)

    if args.fits:
        hdr = pyfits.header.Header()
        fits_image = os.path.join(os.path.dirname(vis_file),
                                  'test_vla_image.fits')
        if os.path.exists(fits_image):
            os.remove(fits_image)
        print 'Writing FITS image: %s' % fits_image
        pyfits.writeto(fits_image, np.real(image.data), hdr, clobber=True)

    for plot_mode in plot_modes:
        print 'Plotting:', plot_mode
        if plot_mode == 'grid':
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, aspect='equal')
            im = np.real(grid.data)
            im_h = ax.imshow(im, interpolation='nearest', cmap=parula_map)
            plt.colorbar(im_h)
            plt.show()

        if plot_mode == 'image':
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, aspect='equal')
            im = np.real(image.data)
            im = np.flipud(im)
            im = np.fliplr(im)
            im_h = ax.imshow(im, interpolation='nearest', cmap=parula_map)
            plt.colorbar(im_h)
            plt.show()

        if plot_mode == 'gcf':
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, aspect='equal')
            im_h = ax.imshow(gcf.data, interpolation='nearest',
                             cmap=parula_map)
            plt.colorbar(im_h)
            plt.show()
