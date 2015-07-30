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


def grid_cell_size(cell_size_lm_arcsec, im_size):
    """Get grid cell size from image cell size and number of image pixels."""
    return (180. * 3600.) / (im_size * cell_size_lm_arcsec * np.pi)


def load_vis(vis_file):
    """Load the oskar visibility file."""
    oskar_vis = OskarVis(vis_file)
    uu, vv, ww = oskar_vis.uvw(flatten=True)
    # uvw = np.transpose(np.array([uu, vv, ww]))
    uvw = np.array([uu, vv, ww])
    freq_hz = oskar_vis.frequency()
    wave_length_m = 299792458. / freq_hz
    uvw /= wave_length_m  # Scale into wavelengths
    vis = oskar_vis.amplitudes(flatten=True)
    return uvw, vis


def w_gcf(w, size, fov, over_sample, taper_width=0.15, cut_level=0.1):
    """Generate w-projection GCF.

    Args:
        w (float) : W-coordinate value, in wavelengths.
        size (int) : Initial kernel size (pixels).
        fov (float) : Image field of view, in degrees.
        over_sample (int) : Number of GCF pixels per grid cell to generate.
        taper_width (float) : Width of aliasing function in the image plane.
        cut_level (float) : Amplitude cutoff when generating final GCF.

    Computation steps:
        1. Generate image plane w-function to the size of the FOV.
        2. Pad to build in the over sample.
        3. FT to the uv-plane.
        4. Cut to given level to generate a smaller kernel.
    """
    plot_w_lm = False
    plot_w_uv_raw = False
    plot_cut_debug = False
    plot_w_uv_final = True

    t0 = time.time()

    # Image cell size for the GCF.
    cell_size = fov_to_cell_size(fov, size)
    x = np.arange(-size / 2, size / 2, dtype='f8')
    x, y = np.meshgrid(x, -x)
    l = x * math.radians(cell_size / 3600.)
    m = y * math.radians(cell_size / 3600.)
    r2 = l**2 + m**2

    # Image plane w-projection kernel. TODO(BM) use better function for this.
    w_lm = np.exp(-1.j * 2. * w * np.pi * (np.sqrt(1. - r2) - 1.))

    # Add a taper. (TODO replace with better anti-aliasing function)
    sigma = taper_width * size
    w_lm *= np.exp(-(x**2 + y**2) / (2.0 * sigma**2))

    # w_lm = np.zeros_like(w_lm)
    # w_lm[size/2, size/2] = 1.0 + 0.0j

    # Plot the image plane kernel.
    if plot_w_lm:
        fig = plt.figure(figsize=(20, 8))
        fig.add_subplot(121, aspect='equal')
        plt.imshow(np.real(w_lm), interpolation='nearest', cmap=parula_map)
        plt.colorbar()
        plt.title('Real(Image plane w-kernel), w=%.3f' % w, fontsize='x-small')
        fig.add_subplot(122, aspect='equal')
        plt.imshow(np.imag(w_lm), interpolation='nearest', cmap=parula_map)
        plt.title('Imag(Image plane w-kernel), w=%.3f' % w, fontsize='x-small')
        plt.colorbar()
        plt.show()

    # Pad and FFT to obtain the desired over sample.
    t_ft = time.time()
    new_size = size * over_sample
    new_centre = new_size / 2
    pad_size = (new_size - size) / 2
    gcf_uv = np.pad(w_lm, pad_width=(pad_size, ), mode='constant',
                    constant_values=(0.0,))

    # # Plot the image plane kernel.
    # fig = plt.figure(figsize=(20, 8))
    # fig.add_subplot(121, aspect='equal')
    # plt.imshow(np.real(w_lm[size/2-10:size/2+10, size/2-10:size/2+10]), interpolation='nearest', cmap=parula_map)
    # plt.colorbar()
    # plt.title('old', fontsize='x-small')
    # fig.add_subplot(122, aspect='equal')
    # plt.imshow(np.real(gcf_uv[new_centre-10:new_centre+10, new_centre-10:new_centre+10]), interpolation='nearest', cmap=parula_map)
    # plt.title('new', fontsize='x-small')
    # plt.colorbar()
    # plt.show()

    # gcf_uv = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(gcf_uv)))
    gcf_uv[::2, :] *= -1
    gcf_uv[:, ::2] *= -1
    gcf_uv = np.fft.fft2(gcf_uv)
    gcf_uv[::2, :] *= -1
    gcf_uv[:, ::2] *= -1
    gcf_uv /= np.max(np.abs(gcf_uv))
    print 'T FT = %.3f s' % (time.time() - t_ft)

    # Plot raw uv plane kernel (full size)
    if plot_w_uv_raw:
        fig = plt.figure(figsize=(10, 10))
        fig.add_subplot(111, aspect='equal')
        plt.imshow(np.real(gcf_uv), interpolation='nearest', cmap=parula_map)
        plt.title('Raw (uncut) GCF uv, w=%.3f' % w, fontsize='x-small')
        plt.colorbar()
        plt.show()

    # fig = plt.figure(figsize=(16, 7))
    # fig.add_subplot(121, aspect='equal')
    # plt.imshow(np.real(gcf_uv[new_centre-5:new_centre+6, new_centre-5:new_centre+6]), interpolation='nearest', cmap=parula_map)
    # plt.colorbar()
    # fig.add_subplot(122, aspect='equal')
    # plt.imshow(np.imag(gcf_uv[new_centre-5:new_centre+6, new_centre-5:new_centre+6]), interpolation='nearest', cmap=parula_map)
    # plt.colorbar()
    # plt.show()

    # Work out GCF clipping level.
    size = gcf_uv.shape[0]
    centre = size / 2
    cut_idx0 = np.argmax(np.abs(gcf_uv[centre, :]) > cut_level)
    cut_idx1 = centre + (centre - cut_idx0)
    w_uv_width = cut_idx1 - cut_idx0 - 1
    w_uv_width = math.ceil(float(w_uv_width) / over_sample)
    w_uv_support = (w_uv_width - 1) / 2.0
    c0 = centre - w_uv_support * over_sample - over_sample / 2
    c1 = centre + w_uv_support * over_sample + over_sample / 2

    # print 'cut_idx0      : %i' % cut_idx0
    # print 'cut_idx1      : %i' % cut_idx1
    # print 'rounded width : %i' % w_uv_width
    # print 'new_support   : %i' % w_uv_support
    # print 'over sample   : %i' % over_sample
    # print 'GCF size      : %i' % (over_sample * (w_uv_support * 2 + 1))
    # print 'centre        : %i' % centre
    # print 'c0            : %i' % c0
    # print 'c1            : %i' % c1

    if plot_cut_debug:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        ax.semilogy(range(0, size), np.abs(gcf_uv[centre + 1, :]), 'o--')
        ax.semilogy([c0, c0], [1.e-5, 1], 'c-', lw=1)
        ax.semilogy([c1, c1], [1.e-5, 1], 'c-', lw=1)
        ax.semilogy([cut_idx0, cut_idx0], [1.e-5, 1], 'r--', lw=1)
        ax.semilogy([cut_idx1, cut_idx1], [1.e-5, 1], 'r--', lw=1)
        ax.semilogy([centre, centre], [1.e-5, 1], 'g--', lw=0.5)
        ax.semilogy([0, size], [cut_level, cut_level], 'r--', lw=0.5)
        plt.show()

    # Cut and re-centre.
    gcf_uv = gcf_uv[c0:c1+1, c0:c1+1]

    if plot_w_uv_final:
        centre = int(math.floor(gcf_uv.shape[0] / 2.))
        width = centre
        # plt.close()
        fig = plt.figure(figsize=(15, 5))
        ax =fig.add_subplot(121, aspect='equal')
        im = np.real(gcf_uv)
        im = im[centre-width:centre+width+1, centre-width:centre+width+1]
        w_uv_width = w_uv_support * 2 + 1
        extent = [-w_uv_width / 2, w_uv_width / 2,
                  -w_uv_width / 2, w_uv_width / 2]
        plt.imshow(im, interpolation='nearest', cmap=parula_map, extent=extent)
        # ticks = np.linspace(-w_uv_width / 2, w_uv_width / 2,
        #                     w_uv_support * 2 + 2)
        # ax.xaxis.set_ticks(ticks)
        # ax.yaxis.set_ticks(ticks)
        plt.grid()
        plt.title('Real. GCF uv, w=%.3f, s%i, o%i' % (w, w_uv_support,
                                                      over_sample))
        plt.colorbar()
        fig.add_subplot(122, aspect='equal')
        im = np.imag(gcf_uv)
        im = im[centre-width:centre+width+1, centre-width:centre+width+1]
        plt.imshow(im, interpolation='nearest', cmap=parula_map, extent=extent)
        plt.title('Imag. GCF uv, w=%.3f, s%i, o%i' % (w, w_uv_support,
                                                      over_sample))
        plt.grid()
        plt.colorbar()
        plt.show(block=True)

    print ('  w-projection GCF = %3f s [support = %i, over-sample = %i]' %
           (time.time() - t0, w_uv_support, over_sample))
    return gcf_uv, w_uv_support


def bin_vis(uvw, vis, w_planes=32):
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


def grid_data(uvw, vis, gcf, support, over_sample, cell_size_uv_m,
              grid):
    assert uvw.shape[0] == vis.shape[0]
    assert uvw.shape[1] == 3
    assert isinstance(support, int)
    assert isinstance(over_sample, int)
    print '  Gridding %i data points' % uvw.shape[0]
    grid_size = grid.shape[0]
    g_centre = grid_size / 2
    k_centre = gcf.shape[0] / 2
    g_sum = 0.0
    kxy = np.arange(-support, support + 1) * over_sample
    kxy_x, kxy_y = np.meshgrid(kxy, kxy)

    # TODO(BM) memory reordering of kernel.

    for i, (coord, amp) in enumerate(zip(uvw, vis)):
        # uv coord scaled to grid space.
        x = -coord[0] / cell_size_uv_m
        y = coord[1] / cell_size_uv_m

        # Closest grid cell.
        xg = int(round(x))
        yg = int(round(y))

        # Grid matrix index.
        ix = xg + g_centre
        iy = yg + g_centre

        # Check point is inside the grid.
        if ix + support >= grid_size or ix - support < 0:
            continue
        if iy + support >= grid_size or iy - support < 0:
            continue

        # Scaled distance from the nearest grid point.
        x_offset = xg - x
        y_offset = yg - y

        # gcf offset
        x_delta = x_offset * over_sample
        y_delta = y_offset * over_sample
        ikx = int(round(x_delta)) + k_centre
        iky = int(round(y_delta)) + k_centre
        # print y_offset, over_sample, y_offset * over_sample
        # print iky
        # print kxy_y + iky
        # print gcf.shape

        # Convolve onto the grid.
        k = gcf[kxy_y + iky, kxy_x + ikx]
        grid[iy-support:iy+support+1, ix-support:ix+support+1] += k * amp
        g_sum += np.sum(k)

    return grid, g_sum


def main():
    # Inputs ------------------------------------------------------------------
    vis_file = 'test_data_3/test_vla.vis'
    w_planes = 32
    im_size = 4096
    fov = 2.0  # deg
    gcf_size = 256
    over_sample = 3
    # -------------------------------------------------------------------------

    print '-' * 80
    print 'w planes    : %i' % w_planes
    print 'im size     : %i' % im_size
    print 'fov         : %.2f deg' % fov
    print 'over sample : %i' % over_sample
    print '-' * 80

    # Load visibility data
    uvw, vis = load_vis(vis_file)
    t0 = time.time()
    w_bin = bin_vis(uvw, vis, w_planes=w_planes)
    print 'Time taken to pre-process %i w planes = %.3fs' % (w_planes,
                                                             time.time() - t0)
    print ''

    grid = np.zeros((im_size, im_size), dtype='c16')
    cell_size_lm_arcsec = fov_to_cell_size(fov, im_size)
    cell_size_uv_m = grid_cell_size(cell_size_lm_arcsec, im_size)
    print 'cell size uv : %f [m]' % cell_size_uv_m
    #
    # b = w_planes / 2
    # b_gcf, support = w_gcf(w_bin[b]['mean_w'], gcf_size, fov, over_sample)

    # Grid one w-plane at a time.
    print 'WARNING:'
    print 'WARNING: There is a bug BUG due to symmetry of the w-kernel?'
    print 'WARNING'
    # TODO(BM) check for symmetry by putting single pixel at fft origin and
    # checking for DC result.
    g_sum = 0.0
    for b in w_bin:
        print '(%2i) mean bin abs(w) : %.3f' % (b, w_bin[b]['mean_w'])
        # Generate w-projection GCF.
        b_gcf, support = w_gcf(w_bin[b]['mean_w'], gcf_size, fov, over_sample)
        tg = time.time()
        grid, b_g_sum = grid_data(np.transpose(w_bin[b]['uvw']), w_bin[b]['amp'],
                                  b_gcf, int(support), over_sample,
                                  cell_size_uv_m, grid)
        print '  Gridding took %.3f s' % (time.time() - tg)
        g_sum += b_g_sum
    # plt.imshow(np.real(grid), interpolation='nearest', cmap=parula_map)
    # plt.colorbar()
    # plt.show()
    grid /= g_sum

    t_im = time.time()
    # image = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(grid)))
    image = np.copy(grid)
    image[::2, :] *= -1
    image[:, ::2] *= -1
    image = np.fft.fft2(image)
    image[::2, :] *= -1
    image[:, ::2] *= -1
    print 'T IM = %.3f s' % (time.time() - t_im)

    # plt.imshow(np.real(image), interpolation='nearest', cmap=parula_map)
    # plt.colorbar()
    # plt.show()

    # Compute grid correction.
    sigma = 0.2 * im_size
    x = np.arange(-im_size / 2, im_size / 2, dtype='f8')
    x, y = np.meshgrid(x, -x)
    correction = np.exp(-(x**2 + y**2) / (2.0 * sigma**2))

    image /= correction

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, aspect='equal')
    im = np.real(image)
    im = np.flipud(im)
    im = np.fliplr(im)
    extent = [-im.shape[0]/2, im.shape[0]/2,
              -im.shape[0]/2, im.shape[0]/2]
    im_h = ax.imshow(im, interpolation='nearest', cmap=parula_map, extent=extent)
    plt.colorbar(im_h)
    plt.show()


if __name__ == '__main__':
    main()
