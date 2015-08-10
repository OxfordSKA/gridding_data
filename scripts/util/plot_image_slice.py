# -*- coding: utf-8 -*-

import pyfits
import numpy
import matplotlib.pyplot as plt
import math
import os
from matplotlib.patches import Rectangle

# Load the fits image.
# image_file = os.path.join('ska1_test_grid_6h', 'test_ska_dirty_04096_w256.fits')
# image_file = os.path.join('ska1_test_grid_6h', 'test_ska_dirty_08192_w256.fits')
# image_file = os.path.join('ska1_test_grid_60s', 'test_ska_dirty_08192_w256.fits')
image_file = os.path.join('ska1_test_grid_60s', 'test_ska_dirty_08192_w000.fits')

image_file = os.path.abspath(image_file)
image_file_root = os.path.splitext(image_file)[0]
image = numpy.squeeze(pyfits.getdata(image_file))

y_step = int(math.ceil(float(image.shape[0]) / 11))
# fig = plt.figure(figsize=(10, 6))
# ax = fig.add_subplot(111, aspect='equal')
# plt.imshow(image, interpolation='nearest')
# plt.colorbar()

for iy, y in enumerate(range(0, image.shape[0], y_step)):
    image_strip = image[y:y+y_step, :]
    strip_max_y = numpy.max(image_strip, axis=0)

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(211)
    plt.imshow(image_strip,
               interpolation='nearest')
    plt.title('y slice %i' % iy)
    plt.colorbar()
    ax = fig.add_subplot(212)
    plt.plot(strip_max_y, '-')
    plt.ylim(-0.1, 1.1)
    plt.grid()
    plt.savefig('%s_ystrip_%02i.png' % (image_file_root, iy))
    # ax.add_patch(Rectangle(
    #     (0, y), image.shape[0], y_step,
    #     edgecolor='w', fill=False, linewidth=1))

# plt.savefig('%s_strips.png' % (image_file_root))
