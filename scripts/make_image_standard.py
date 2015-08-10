# -*- coding: utf-8 -*-
"""Demonstration of gridding & imaging."""


import numpy as np
import math
import sphfn


def fov_to_cell_size(fov, im_size):
    """Evaluate image pixel size for given FoV and number of pixels.

    Args:
        fov (float): Field of view, in degrees
        im_size (int): Image size (1D) in pixels.

    Returns:
        double: pixel size, in arcsec.
    """
    r_max = math.sin(math.radians(fov) / 2.)
    inc = r_max / (0.5 * im_size)
    return math.degrees(math.asin(inc)) * 3600.0


def cell_size_to_fov(cell_size, im_size):
    """Obtain image fov from cell size and image size."""
    inc = math.sin(math.radians(cell_size / 3600.0))
    r_max = inc * (0.5 * im_size)
    return math.degrees(2.0 * math.asin(r_max))


class GriddingConvolutionFunction:

    """Class holding a gridding convolution function.

    Aliasing & choice of GCF: p142 Synthesis imaging.
    """

    def __init__(self, support, over_sample):
        self.support = support
        self.over_sample = over_sample
        self.num_cells = support * 2 + 1
        self.num_pixels = self.num_cells * over_sample
        self.centre = self.num_pixels / 2
        self.inc = 1.0 / over_sample
        self.data = np.ones((self.num_pixels, self.num_pixels), dtype='f8')

    def exponential_2d(self, width=1.0):
        """Generate an exponential (Gaussian) GCF."""
        p1 = 1.0 / width
        p2 = 2.0
        x, y = self.coord_grid()
        arg_x = (np.abs(x) * p1)**p2
        arg_y = (np.abs(y) * p1)**p2
        self.data *= np.exp(-(arg_x + arg_y))

    def pillbox_2d(self, width=1.0):
        """Generate a pillbox GCF."""
        x, y = self.coord_grid()
        self.data[(np.abs(x) >= width) | (np.abs(y) >= width)] = 0.0

    def sinc_2d(self):
        """Generate a sin(x)/x GCF."""
        p1 = 1.0
        x, y = self.coord_grid()
        ax = np.pi * x / p1
        ay = np.pi * y / p1
        k = np.ones(self.data.shape, dtype='f8')
        with np.errstate(divide='ignore', invalid='ignore'):
            kx = np.sin(ax) / ax
            kx[ax == 0.0] = 1.0
            ky = np.sin(ay) / ay
            ky[ay == 0.0] = 1.0
            k = kx * ky
        self.data *= k

    def exponential_sinc_2d(self):
        """."""
        self.sinc_2d()
        self.exponential_2d(2.0)  # FIXME(BM) exp_2d cpar1, cpar2 etc...

    def sphfn(self, exponent_id):
        """."""
        support = self.support
        over_sample = self.over_sample


    def coord_grid(self):
        """Evaluate the GCF coordinate grid in pixel space."""
        n = self.num_pixels
        # Note: this linspace can be used as the kernel is odd sized.
        # FIXME(BM) double check this!
        c = np.linspace(-(n/2), n/2, n) * self.inc
        x, y = np.meshgrid(c, c)
        return x, y


class Grid:

    """Class holding a UV grid."""

    def __init__(self, size, cell_size_m):
        self.size = size
        self.cell_size_m = cell_size_m  # m
        self.centre = size / 2  # zero based grid centre index (even grid only!)
        self.data = np.zeros((size, size), dtype='c16')

    def is_hermitian(self):
        """Returns True is the grid is Hermitian."""
        return np.max(np.abs(np.conj(self.data.T) - self.data)) == 0.0


class Image:

    """."""

    def __init__(self, grid, pad=1.2):
        self.size = grid.size
        self.pad = pad
        self.cell_size_rad = 1.0 / (grid.cell_size_m * grid.size)
        self.correction = np.empty(grid.data.shape, dtype='f8')
        image = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(grid.data)))
        self.data = image

    def cell_size(self):
        return self.cell_size_rad * (180./np.pi) * 3600.0

    def fov(self):
        return cell_size_to_fov(self.cell_size(), self.size)

    def lm(self):
        inc = np.sin(self.cell_size_rad)
        centre = self.size / 2
        x0 = -centre - 0.5
        x1 = centre - 0.5
        lm = np.linspace(x0 * inc, x1 * inc, self.size)
        return lm

    def grid_correct(self, gcf):
        """Correct for the image taper introduced by the GCF.

        p138 Synthesis imaging in radio astronomy. 'Grid correction is not
        an exact correction, except in the limit of large numbers of well
        distributed visibility measurements.'
        """
        # TODO(BM) use padding rather than sampling of the GCF onto the grid.??
        # TODO(BM) use analytical FT pairs!
        pad_width = self.data.shape[0] - gcf.num_cells
        print 'xxx', pad_width, pad_width/2
        print np.arange(0, gcf.num_cells)
        print 63/2
        k_idx = np.arange(0, gcf.num_cells) * gcf.over_sample + gcf.over_sample / 2
        k = gcf.data[k_idx, k_idx]
        print '--', k.shape
        print gcf.num_cells, gcf.over_sample, gcf.data.shape
        # k = np.pad(gcf.data, pad_width=(pad_size, ), mode='constant',
        #             constant_values=(0.0,))

        k = np.zeros(self.data.shape, dtype='f8')
        for j in range(0, gcf.num_cells):
            ky = j * gcf.over_sample + gcf.over_sample / 2
            jk = (self.data.shape[0] / 2) + j - (gcf.num_cells / 2)
            for i in range(0, gcf.num_cells):
                kx = i * gcf.over_sample + gcf.over_sample / 2
                ik = self.data.shape[0] / 2 + i - gcf.num_cells / 2
                k[jk, ik] = gcf.data[kx, ky]
        correction = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(k)))
        correction /= np.max(correction)
        correction = 1.0 / correction
        self.correction = np.real(correction)
        self.data *= self.correction


def grid_data(uvw, vis, gcf, grid, normalise=True, do_hermitian_copy=False):
    g_centre = grid.size / 2
    k_centre = gcf.num_pixels / 2
    g_sum = 0.0
    support = gcf.support
    cell_size_m = grid.cell_size_m
    over_sample = gcf.over_sample
    kxy = np.arange(-support, support + 1) * gcf.over_sample
    kxy_x, kxy_y = np.meshgrid(kxy, kxy)
    assert uvw.shape[0] == vis.shape[0]
    assert uvw.shape[1] == 3

    for i, (coord, amp) in enumerate(zip(uvw, vis)):
        # uv coord scaled to grid space.
        x = -coord[0] / cell_size_m
        y = coord[1] / cell_size_m

        # Closest grid cell.
        xg = int(round(x))
        yg = int(round(y))

        # Grid matrix index.
        ix = xg + g_centre
        iy = yg + g_centre

        # Check point is inside the grid.
        if ix + support >= grid.size or ix - support < 0:
            continue
        if iy + support >= grid.size or iy - support < 0:
            continue

        # Scaled distance from the nearest grid point.
        x_offset = xg - x
        y_offset = yg - y

        # gcf offset
        x_delta = x_offset * over_sample
        y_delta = y_offset * over_sample
        ikx = int(round(x_delta)) + k_centre
        iky = int(round(y_delta)) + k_centre

        # Convolve onto the grid.
        k = gcf.data[kxy_y + iky, kxy_x + ikx]
        grid.data[iy-support:iy+support+1, ix-support:ix+support+1] += k * amp
        g_sum += np.sum(k)

        # Convolve hermitian copy onto the grid, not needed if taking the real
        # part of the image.
        if do_hermitian_copy:
            ix = g_centre - xg
            iy = g_centre - yg
            ikx = k_centre - int(round(x_delta))
            iky = k_centre - int(round(y_delta))
            k = gcf.data[kxy_y + iky, kxy_x + ikx]
            grid.data[iy-support:iy+support+1, ix-support:ix+support+1] += \
                k * np.conj(amp)
            g_sum += np.sum(k)

    if normalise:
        grid.data /= g_sum
