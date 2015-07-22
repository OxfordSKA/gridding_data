# -*- coding: utf-8 -*-
"""Demonstration of gridding & imaging."""


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker


def fov_to_cell_size(fov, im_size):
    """Obtain cell_size from a fov and image size.

    Args:
        fov (float, list): Field of view, in degrees
        im_size (int, list): Image size (1D) in pixels.

    Returns:
        double, list: Cell size(s) in arcsec.
    """
    r_max = np.sin(np.array(fov, 'f8') / 2. * (np.pi / 180.))
    inc = r_max / (0.5 * np.array(im_size, 'i8'))
    cell = np.arcsin(inc) * ((180. * 3600.) / np.pi)
    return cell.tolist()


def cell_size_to_fov(cell_size, im_size):
    """Obtain image fov from cell size and image size."""
    cell = (cell_size / 3600.0) * (np.pi / 180.0)
    inc = np.sin(cell)
    r_max = inc * (0.5 * im_size)
    fov = 2 * np.arcsin(r_max) * 180.0 / np.pi
    return fov


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
        arg_x = np.power(np.abs(x) * p1, p2)
        arg_y = np.power(np.abs(y) * p1, p2)
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

    def kaiser(self, dk, W, alpha):
        beta = np.pi * np.sqrt((W**2 / alpha / alpha) * (alpha - 0.5)**2 - 0.8)
        # temp3 =

    def coord_grid(self):
        """Evaluate the GCF coordinate grid in pixel space."""
        n = self.num_pixels
        c = np.linspace(-(n/2), n/2, n) * self.inc
        x, y = np.meshgrid(c, c)
        return x, y

    def plot(self):
        centre = self.num_cells / 2.0
        extent_ = [-centre, +centre, -centre, +centre]
        plt.imshow(self.data, interpolation='nearest', cmap='jet',
                   extent=extent_)
        plt.title('support %i, over-sample %i, cells %i, '
                  'shape (%i, %i), centre %.2f' %
                  (self.support, self.over_sample, self.num_cells,
                   self.data.shape[0],
                   self.data.shape[1], self.centre), fontsize='x-small')
        plt.colorbar()
        plt.grid(which='both', color='w', linestyle='--')
        plt.show()


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

    def plot(self, grid_coords=None):
        centre = self.size / 2
        extent = [-centre - 0.5, centre - 0.5,
                  -centre - 0.5, centre - 0.5]
        extent = np.array(extent, dtype='f8') * self.cell_size_m

        data = np.flipud(self.data)
        if grid_coords:
            data = data[grid_coords[0]:grid_coords[1],
                        grid_coords[2]:grid_coords[3]]

        fig = plt.figure(figsize=(12, 6))

        ax = fig.add_subplot(121)
        plt.imshow(np.real(data), interpolation='nearest', cmap='jet',
                   extent=extent)
        plt.title('real(grid)')
        plt.grid(which='both', color='w', linestyle='--')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        plt.colorbar(cax=cax)

        ax = fig.add_subplot(122)
        plt.imshow(np.imag(data), interpolation='nearest', cmap='jet',
                   extent=extent)
        plt.title('imag(grid)')
        plt.grid(which='both', color='w', linestyle='--')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        plt.colorbar(cax=cax)
        plt.show()


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

    def plot(self):
        size = self.size
        centre = size / 2
        extent = [-centre - 0.5, centre - 0.5,
                  -centre - 0.5, centre - 0.5]
        extent = np.array(extent, dtype='f8')
        extent *= np.sin(self.cell_size_rad)

        data = np.flipud(self.data)

        plt.figure(figsize=(12, 6))

        # http://matplotlib.org/users/gridspec.html
        gs1 = gridspec.GridSpec(4, 4)
        gs1.update(left=0.05, right=0.45, wspace=0.3, hspace=0.3)

        # Image - Real part
        plt.subplot(gs1[0:3, 0:3])
        im = plt.imshow(np.real(data), extent=extent, interpolation='nearest')
        cbar = plt.colorbar(im, shrink=0.85, pad=0.01)
        cbar.ax.tick_params(labelsize=6)
        plt.title('real(image)')
        plt.grid(which='both', color='w', linestyle='--')
        plt.tick_params(axis='both', which='major', labelsize=6)

        # Y-axis slice
        ax2 = plt.subplot(gs1[0:3, 3])  # gs1[y, x]
        x2 = np.real(self.data[:, centre])
        y2 = np.linspace(extent[2], extent[3], self.size)
        ax2.plot(x2, y2)
        start, end = ax2.get_xlim()
        ax2.set_ylim(y2[0], y2[-1])
        ax2.xaxis.set_ticks(np.linspace(start, end, 4))
        ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax2.tick_params(axis='both', which='both', labelsize=6)
        ax2.grid()

        # X-axis slice
        ax3 = plt.subplot(gs1[3, 0:3])
        x3 = np.linspace(extent[0], extent[1], self.size)
        y3 = np.real(data[centre, :])
        ax3.plot(x3, y3)
        start, end = ax3.get_ylim()
        ax3.set_xlim(x3[0], x3[-1])
        ax3.yaxis.set_ticks(np.linspace(start, end, 4))
        ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax3.grid()
        plt.tick_params(axis='both', which='major', labelsize=6)

        # Image - Imaginary part
        gs2 = gridspec.GridSpec(4, 4)
        gs2.update(left=0.55, right=0.95, wspace=0.3, hspace=0.3)

        # Image plot
        plt.subplot(gs2[0:3, 0:3])
        im = plt.imshow(np.imag(data), extent=extent, interpolation='nearest')
        cbar = plt.colorbar(im, shrink=0.85, pad=0.01)
        cbar.ax.tick_params(labelsize=6)
        plt.title('imag(image)')
        plt.grid(which='both', color='w', linestyle='--')
        plt.tick_params(axis='both', which='major', labelsize=6)

        # Y-axis slice
        ax2 = plt.subplot(gs2[0:3, 3])  # gs1[y, x]
        x2 = np.imag(self.data[:, centre])
        y2 = np.linspace(extent[2], extent[3], self.size)
        ax2.plot(x2, y2)
        start, end = ax2.get_xlim()
        ax2.set_ylim(y2[0], y2[-1])
        ax2.xaxis.set_ticks(np.linspace(start, end, 4))
        ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
        ax2.tick_params(axis='both', which='both', labelsize=6)
        ax2.grid()

        # X-axis slice
        ax3 = plt.subplot(gs2[3, 0:3])
        x3 = np.linspace(extent[0], extent[1], self.size)
        y3 = np.imag(data[centre, :])
        ax3.plot(x3, y3)
        start, end = ax3.get_ylim()
        ax3.set_xlim(x3[0], x3[-1])
        ax3.yaxis.set_ticks(np.linspace(start, end, 4))
        ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
        ax3.grid()
        plt.tick_params(axis='both', which='major', labelsize=6)

        plt.show()

    def plot_correction(self):
        plt.imshow(self.correction, interpolation='nearest')
        plt.colorbar()
        plt.title('correction')
        plt.show()


def add_to_grid(amp, grid_i0, gcf_i0, gcf, grid):
    # TODO(BM) better memory order of gcf
    # TODO(BM) padding of grid to avoid edge issues?
    # TODO(BM) need catch for invalid cell size
    size = 2 * gcf.support + 1
    gx0 = grid_i0[0] - gcf.support
    gy0 = grid_i0[1] - gcf.support
    kxy = np.arange(-gcf.support, gcf.support + 1) * gcf.over_sample
    kxy_x, kxy_y = np.meshgrid(kxy + gcf_i0[0], kxy + gcf_i0[1])
    k = gcf.data[kxy_x, kxy_y]
    grid.data[gx0:gx0+size, gy0:gy0+size] += k * amp
    # k = gcf.data[kxy_y, kxy_x]
    # grid.data[gy0:gy0+size, gx0:gx0+size] += k * amp
    return np.sum(k)



def grid_data(uvw, vis, gcf, grid, normalise=True):
    g_centre = grid.size / 2
    k_centre = gcf.num_pixels / 2
    g_sum = 0.0
    assert uvw.shape[0] == vis.shape[0]
    assert uvw.shape[1] == 3

    for coord, amp in zip(uvw, vis):
        # print coord, amp

        # TODO(BM) function to get grid coordinate and gcf offset.
        # uv coord scaled to grid space.
        x = coord[0] / grid.cell_size_m
        y = coord[1] / grid.cell_size_m

        # Closest grid cell.
        xg = int(round(x))
        yg = int(round(y))

        # Grid matrix index.
        ix = xg + g_centre
        iy = yg + g_centre

        # Scaled distance from the nearest grid point.
        x_offset = xg - x
        y_offset = yg - y

        # gcf offset
        x_delta = x_offset * gcf.over_sample
        y_delta = y_offset * gcf.over_sample
        ikx = int(round(x_delta)) + k_centre
        iky = int(round(y_delta)) + k_centre

        # Convolve onto the grid.
        g_sum += add_to_grid(amp, [ix, iy], [ikx, iky], gcf, grid)

        # # NOTE: Don't need to do this... just take the real part of the FFT.
        # ix = g_centre - xg
        # iy = g_centre - yg
        # ikx = k_centre - int(round(x_delta))
        # iky = k_centre - int(round(y_delta))
        #
        # # NOTE this function only ever uses values in the gcf with a
        # # certain stride (offset + over-sample * n)... Generate GCF
        # # with stride packing and the adds to grid is then a
        # # element multiply of two matrices.
        # g_sum += add_to_grid(amp, [ix, iy], [ikx, iky], gcf, grid)

    if normalise:
        grid.data /= g_sum
