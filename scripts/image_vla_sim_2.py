# -*- coding: utf-8 -*-
"""Imaging script for the VLA simulated data."""

import argparse
import numpy as np
import math
import time
import os
import cProfile
import pstats
import sys
from parula import parula_map
import matplotlib.pyplot as plt
from read_oskar_vis import OskarVis
from make_image_standard import *


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_usage()
        sys.exit(2)


def grid_cell_size(cell_size_lm_arcsec, im_size):
    return 0.5 / ((im_size * cell_size_lm_arcsec * np.pi) / (3600. * 180.))


def load_vis(vis_file):
    oskar_vis = OskarVis(vis_file)
    uu, vv, ww = oskar_vis.uvw(flatten=True)
    uvw = np.array([uu, vv, ww])
    vis = oskar_vis.amplitudes(flatten=True)
    return uvw, vis


def w_gcf(im_size, over_sample, cell_size_arcsec, w, taper_width=0.15):
    """Generate w-projection GCF.

    1. Generate in the image plane.
        a. Need to make a symmetric with the zero in the correct place for
           the FFT.
        b. The function will be eventually padded to generate the required
           over-sample.
    2. FFT the padded (by over_sample) image plane function to generate the uv
       plane GCF.

    Open questions.
        - is there a better way to generate the GCF than starting with a
          image plane kernel the size of the FoV?
    """
    # FIXME(BM) must be a better way to get the gcf cell size here!
    t0 = time.time()
    fov = cell_size_to_fov(cell_size_arcsec, im_size)
    cell_size = fov_to_cell_size(fov, im_size * over_sample)
    cell_size *= np.pi / (3600. * 180.)
    inc = cell_size
    gcf_x = np.arange(-im_size / 2, im_size / 2, dtype='f8')
    x, y = np.meshgrid(gcf_x, -gcf_x)
    l = x * inc
    m = y * inc
    r2 = l**2 + m**2
    # Image plane w-projection kernel.
    gcf_lm = np.exp(-1.j * 2. * w * np.pi * (np.sqrt(1. - r2) - 1.))
    # Add a taper.
    sigma = taper_width * im_size
    gcf_lm *= np.exp(-(x**2 + y**2) / (2.0 * sigma**2))

    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(111, aspect='equal')
    # plt.imshow(np.real(gcf_lm), interpolation='nearest', cmap=parula_map)
    # plt.title('GCF lm, w=%.3f' % w)
    # plt.colorbar()
    # plt.show()

    # Pad and FFT.
    pad_size = im_size
    # FIXME pad_width adds gcf_size pixels to each side of the image.
    # so pad size == im_size over sample of 3?
    gcf_uv = np.pad(gcf_lm, pad_width=(pad_size, ), mode='constant',
                    constant_values=(0.0,))
    gcf_uv = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(gcf_uv)))
    gcf_uv /= np.max(np.abs(gcf_uv))
    print 'over sample = %f' % (gcf_uv.shape[0] / float(im_size))

    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(111, aspect='equal')
    # plt.imshow(np.real(gcf_uv), interpolation='nearest', cmap=parula_map)
    # plt.title('GCF uv, w=%.3f' % w)
    # plt.colorbar()
    # plt.show()

    size = gcf_uv.shape[0]
    centre = size/2
    print size, centre

    plot_gcf_uv = False
    if plot_gcf_uv:
        a = 50
        fig = plt.figure(figsize=(15, 5))
        fig.add_subplot(121, aspect='equal')
        im = np.real(gcf_uv)
        im = im[centre-a+1:centre+a, centre-a+1:centre+a]
        plt.imshow(im, interpolation='nearest', cmap=parula_map)
        plt.title('Real. GCF uv, w=%.3f' % w)
        plt.colorbar()
        fig.add_subplot(122, aspect='equal')
        im = np.imag(gcf_uv)
        im = im[centre-a+1:centre+a, centre-a+1:centre+a]
        plt.imshow(im, interpolation='nearest', cmap=parula_map)
        plt.title('Imag. GCF uv, w=%.3f' % w)
        plt.colorbar()
        plt.show()

    # Work out GCF clipping level.
    cut_level = 1.e-1
    print gcf_uv.shape
    print gcf_uv[centre, :].shape
    # TODO(BM) round offset from centre to whole number of support cells.
    cut_idx0 = np.argmax(np.abs(gcf_uv[centre + 1, :]) > cut_level)
    cut_idx1 = centre + (centre - cut_idx0)

    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(111)
    # ax.semilogy(range(0, size), np.abs(gcf_uv[centre + 1, :]), '+--')
    # ax.semilogy([cut_idx0, cut_idx0], [1.e-5, 1], 'r-', lw=1)
    # ax.semilogy([cut_idx1, cut_idx1], [1.e-5, 1], 'r-', lw=1)
    # ax.semilogy([centre, centre], [1.e-5, 1], 'g--', lw=0.5)
    # ax.semilogy([0, size], [cut_level, cut_level], 'r--', lw=0.5)
    # plt.show()

    # Cut and recentre.
    gcf_uv = gcf_uv[cut_idx0+1:cut_idx1, cut_idx0+1:cut_idx1]

    centre = int(math.floor(gcf_uv.shape[0] / 2.))
    width = centre
    fig = plt.figure(figsize=(15, 5))
    fig.add_subplot(121, aspect='equal')
    im = np.real(gcf_uv)
    im = im[centre-width:centre+width+1, centre-width:centre+width+1]
    plt.imshow(im, interpolation='nearest', cmap=parula_map)
    plt.title('Real. GCF uv, w=%.3f' % w)
    plt.colorbar()
    fig.add_subplot(122, aspect='equal')
    im = np.imag(gcf_uv)
    im = im[centre-width:centre+width+1, centre-width:centre+width+1]
    plt.imshow(im, interpolation='nearest', cmap=parula_map)
    plt.title('Imag. GCF uv, w=%.3f' % w)
    plt.colorbar()
    plt.show()

    print 'Time taken to generate w-projection GCF = %3f s' % (time.time() - t0)


def bin_vis(uvw, vis, w_planes = 32):
    """."""
    # Bin uvw data data into w_plane bins in sqrt(abs(w)) space.
    sqrt_ww = np.sqrt(np.abs(uvw[2]))
    # ww bin edges in sqrt(abs(ww)) space.
    ww_bin_limits_sqrt = np.linspace(min(sqrt_ww), max(sqrt_ww), w_planes + 1)
    # ww bin edges in linear abs(ww) space.
    ww_bin_limits = ww_bin_limits_sqrt**2

    # Sort coordinates into abs(ww) order
    sort_indices = np.abs(uvw[2, :]).argsort()
    uvw = uvw[:, sort_indices]
    vis = vis[sort_indices]

    # Get indices and counts of visibilities in each w bin.
    bin_indices = [0]
    for i in range(1, w_planes):
        w_upper = ww_bin_limits[i]
        idx = np.argmin(uvw[2] <= w_upper)
        bin_indices.append(idx)
    bin_indices.append(len(uvw[2]))

    visibility_bin = {}
    total_vis = 0
    for i in range(0, w_planes):
        i0 = bin_indices[i]
        i1 = bin_indices[i + 1] - 1
        uvw_bin = uvw[:, i0:i1 + 1]
        vis_bin = vis[i0:i1 + 1]
        total_vis += uvw_bin.shape[1]
        mean_w = np.mean(np.abs(uvw_bin[2]))
        min_w = np.min(np.abs(uvw_bin[2]))
        max_w = np.max(np.abs(uvw_bin[2]))
        centre_w = ((max_w - min_w) / 2.) + min_w
        visibility_bin[i] = {'i0': i0,
                             'i1': i1,
                             'min_w': min_w,
                             'max_w': max_w,
                             'mean_w': mean_w,
                             'centre_w': centre_w,
                             'uvw': uvw_bin,
                             'amp': vis_bin}
        assert uvw_bin.shape[1] == vis_bin.shape[0]
    assert total_vis == len(uvw[2])
    return visibility_bin


def main():
    # As the taper applied to the w-kernel has an analytical image plane
    # function use this for exact grid correction.

    vis_file = 'test_data_3/test_vla.vis'
    w_planes = 32
    uvw, vis = load_vis(vis_file)
    t0 = time.time()
    w_bin = bin_vis(uvw, vis, w_planes=w_planes)
    print 'Time taken to pre-process %i w planes = %.3fs' % (w_planes,
                                                             time.time() - t0)
    im_size = 512
    fov = 5.0  # deg
    over_sample = 1
    cell_size_arcsec = fov_to_cell_size(fov, im_size)

    # w_gcf(im_size, over_sample, cell_size_arcsec, w=w_bin[w_planes-1]['mean_w'])
    w_gcf(im_size, over_sample, cell_size_arcsec, w=w_bin[w_planes/2]['mean_w'])

    # for b in w_bin:
    #     # print b, w_bin[b]['mean_w']
    #     # Generate w-projection GCF.
    #     w_gcf(im_size, over_sample, gcf_cell_size, w_bin[b]['mean_w'])
    #     # Grid


if __name__ == '__main__':
    main()
