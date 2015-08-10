# -*- coding: utf-8 -*-

import pyfits
import numpy
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# Load the fits image.
# image_file = 'ska1_test_3/test_ska_dirty_w256.fits'
image_file = 'ska1_test_3/test_ska_dirty_16384_w256.fits'
image = numpy.squeeze(pyfits.getdata(image_file))


# Find list of pixels > level
level = 0.8
pixels = numpy.where(image > level)
pixels_x = numpy.argsort(pixels[0])
pixels = (pixels[0][pixels_x], pixels[1][pixels_x])

fig = plt.figure(figsize=(16, 6))
plt.show(block=False)

width = 40
for i, loc in enumerate(zip(pixels[0], pixels[1])):

    plt.clf()

    ax = fig.add_subplot(121, aspect='equal')
    plt.imshow(image, interpolation='nearest')
    plt.colorbar()


    ax = fig.add_subplot(122, aspect='equal')
    im_block = image[loc[0]-width:loc[0]+width, loc[1]-width:loc[1]+width]
    plt.imshow(im_block, interpolation='nearest')
    plt.title('%i, %i: max = %.3f' % (loc[0], loc[1], numpy.max(im_block)))
    plt.colorbar()

    print i, loc, numpy.max(im_block)

    plt.draw()

