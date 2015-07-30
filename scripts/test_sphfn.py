# -*- coding: utf-8 -*-
"""Test of the sphfn.py module."""

import sphfn
import numpy
import matplotlib.pyplot as plt


def test1():
    """."""
    support = 5
    over_sample = 100
    exponent_id = 4

    width = support * 2 + 1
    size = width * over_sample
    centre = (over_sample / 2) * width
    print 'Support     : %i' % support
    print 'Over sample : %i' % over_sample
    print 'Width       : %i' % width
    print 'Size        : %i' % size
    print 'Centre      : %i' % centre

    # Declare memory for the GCF and its FT
    gcf = numpy.zeros((size, ), dtype='f8')
    correction = numpy.zeros((size, ), dtype='f8')

    # eta values.
    eta = numpy.linspace(0.0, 1.0, size / 2)
    for i, eta_ in enumerate(eta):
        idx = centre + i
        gcf[idx] = sphfn.sphfn(eta_, exponent_id, support, gridding=True)
        correction[idx] = sphfn.sphfn(eta_, exponent_id, support,
                                      gridding=False)

    gcf[0:centre] = gcf[:centre-1:-1]
    correction[0:centre] = correction[:centre-1:-1]

    x = numpy.linspace(-1.0, 1.0, size)

    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(211)
    ax.plot(x, gcf, 'r-')
    ax.grid()
    ax.set_title('GCF')

    ax = fig.add_subplot(212)
    ax.plot(x, correction, 'b-')
    ax.grid()
    ax.set_title('correction')

    plt.show()

if __name__ == '__main__':
     test1()
