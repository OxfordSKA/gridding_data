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
    return (180. * 3600.) / (im_size * cell_size_lm_arcsec * np.pi)


def sphfn(eta, iflag, ialf):
    """Port of AIPS 31DEC15/APL/SUB/NOTST/SPHFN.FOR.

    eta goes from 0 at the centre to 1 at the edge

    NOTE: current limited to support 5
    """
    p5 = [3.722238e-3, -4.991683e-2,  1.658905e-1, -2.387240e-1,
          1.877469e-1, -8.159855e-2,  3.051959e-2,  8.182649e-3,
          -7.325459e-2, 1.945697e-1, -2.396387e-1,  1.667832e-1,
          -6.620786e-2, 2.224041e-2,  1.466325e-2, -9.858686e-2,
          2.180684e-1, -2.347118e-1,  1.464354e-1, -5.350728e-2,
          1.624782e-2,  2.314317e-2, -1.246383e-1,  2.362036e-1,
          -2.257366e-1, 1.275895e-1, -4.317874e-2,  1.193168e-2,
          3.346886e-2, -1.503778e-1,  2.492826e-1, -2.142055e-1,
          1.106482e-1, -3.486024e-2,  8.821107e-3]
    p5 = np.array(p5, dtype='f8', order='F')
    p5 = p5.reshape((7, 5), order='F')
    q5 = np.array([2.418820e-1, 2.291233e-1, 2.177793e-1, 2.075784e-1,
                   1.983358e-1], dtype='f8')
    alpha = np.array([0.0, 0.5, 1.0, 1.5, 2.0], dtype='f8')

    eta2 = eta**2
    x = eta2 - 1.0
    j = ialf  # weighting exponent (1, 2, 3, 4, 5)
    psi = (p5[0, j] + x * (p5[1, j] + x * (p5[2, j] + x * (p5[3, j] +\
        x * (p5[4, j] + x * (p5[5, j] + x * p5[6, j]))))))
    psi /= (1.0 + x * q5[j])

    if iflag > 0 or ialf == 1 or eta == 0.0:
        return psi

    if math.fabs(eta) == 1.0:
        return 0.0

    print 'here', alpha[ialf], (1.0 - eta2)**alpha[ialf]
    return (1.0 - eta2)**alpha[ialf] * psi


def load_vis(vis_file):
    oskar_vis = OskarVis(vis_file)
    uu, vv, ww = oskar_vis.uvw(flatten=True)
    uvw = np.array([uu, vv, ww])
    freq_hz = oskar_vis.frequency()
    wave_length_m = 299792458. / freq_hz
    uvw /= wave_length_m  # Scale into wavelengths
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
    t0 = time.time()

    # FIXME(BM) must be a better way to get the gcf cell size here!
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

    # Pad and FFT.4,
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
    uvw, vis = load_vis(vis_file) # FIXME(BM) scale into wavelengths!!!
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
    # main()

    do_sphfn = False
    if do_sphfn:
        # This computes half of a 1d spheroidal function.
        # support = 5 only?!
        num_points = 20
        gcf_1 = np.empty((num_points, ), dtype='f8')
        gcf_2 = np.empty((num_points, ), dtype='f8')
        x = np.linspace(0., 1., num_points)

        for i, eta in enumerate(x):
            gcf_1[i] = sphfn(eta, iflag=1, ialf=3)
            gcf_2[i] = sphfn(eta, iflag=-1, ialf=3)

        plt.plot(x, gcf_1, 'rx-')
        plt.plot(x, gcf_2, 'b+--')
        plt.show()

    # CONVFN.FOR
    param = [3]
    nrow = max(param[0], 1)
    nrow = nrow * 2 + 1
    linc = 100.0
    linc = (linc / 2) * 2
    lim = nrow * linc
    bias = (linc / 2.0) * nrow
    umax = param[0]
    xinc = 1.0 / linc

    print 'param[0]', param[0]
    print 'nrow    ', nrow
    print 'linc    ', linc
    print 'lim     ', lim
    print 'bias    ', bias
    print 'umax    ', umax
    print 'xinc    ', xinc

    p1 = 1.0 / 1.0
    p1 = 2.0
    p2 = 3.0
    buffer_size = (2 * param[0] + 1) * int(linc)
    buffer = np.zeros((buffer_size,), dtype='f8')
    buffer_u = np.zeros((buffer_size,), dtype='f8')
    umax = param[0]
    for i in range(0, int(lim)):
        u = (i - bias) * xinc
        buffer_u[i] = u
        absu = math.fabs(u)
        if absu <= umax:
            buffer[i] = math.exp(-((p1*absu)**p2))
        # print i, buffer_u[i], buffer[i]
    # plt.plot(buffer_u, buffer)
    # plt.show()
    buffer /= np.sum(buffer)  # Normalise by area.

    # Combine to make 2d image of the kernel.
    im_buffer = np.empty((buffer_size, buffer_size), dtype='f8')
    for i in range(0, buffer_size):
        im_buffer[:, i] = buffer
    for j in range(0, buffer_size):
        im_buffer[j, :] *= buffer
    # plt.imshow(im_buffer, interpolation='nearest')
    # plt.show()

    # Compute grid correction.
    im_size = 64
    c_buffer = np.empty(shape=(im_size), dtype='f8')
    image_l = np.arange(-im_size/2, im_size/2, dtype='f8')
    fov_deg = 10.0
    l_max = math.sin(math.radians(fov_deg) / 2.)
    l_inc = l_max / (0.5 * im_size)
    image_l *= l_inc

    # DFT of kernel modulated by sinc of over-sample top hats
    print buffer.shape
    print im_size
    for i, l in enumerate(image_l):
        if l == 0.0:
            c_buffer[i] = 1.0
        else:
            c_buffer[i] = 0.0
            for j, k in enumerate(buffer):
                c_buffer[i] += k * math.cos(2.0 * math.pi * l * buffer_u[j])
            arg = math.pi * xinc * l
            c_buffer[i] *= math.sin(arg) / arg

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(211)
    ax.plot(buffer_u, buffer)
    ax.set_title('GCF')
    ax.grid()
    ax = fig.add_subplot(212)
    ax.plot(image_l, c_buffer)
    ax.set_title('correction')
    ax.grid()
    plt.ylim([0.9, 1.0])
    plt.show()
