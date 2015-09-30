# -*- coding: utf-8 -*-
"""Simple CASA imaging script.

Run with:
    casapy --nogui --nologger --log2term --nologfile -c casa_image.py <ms>

See:
    http://casa.nrao.edu/docs/CasaRef/CasaRef.html
for documentation of methods on CASA objects 'im' and 'ia'
"""

import os
import time
import sys
import math
import glob


def fov_to_cell_size(fov_deg, im_size):
    """Convert field of view in degrees to cell-size in arcsec"""
    inc = math.sin(math.radians(fov_deg) / 2.) / (0.5 * im_size)
    return math.degrees(math.asin(inc)) * 3600.


if __name__ == '__main__':

    ms_path = os.path.abspath(sys.argv[-1])
    if not os.path.isdir(ms_path):
        raise RuntimeError('Specified MS ({}) not found!'.format(ms_path))

    # ------------------------------------------------------------------------
    size = 4096
    fov = 4.5  # deg
    # ra = -23.464210816
    # dec = -153.231195727
    ra = None
    dec = None
    w_planes = 256
    weighting = 'natural'
    niter = 40000
    gain = 0.1
    # ------------------------------------------------------------------------
    root_name = '%s_%04i_%3.1f_w%03i_%s' % (os.path.splitext(ms_path)[0],
                                         size, fov, w_planes, weighting)
    images = {
        'dirty': '%s.dirty.img' % root_name,
        'restored': '%s.restored.img' % root_name,
        'residual': '%s.residual.img' % root_name,
        'psf': '%s.psf.img' % root_name,
        'model': '%s.cc.img' % root_name
    }
    # ------------------------------------------------------------------------

    im_size = [size, size]
    cell = fov_to_cell_size(fov, size)  # arcsec
    cell_size = ['%.10farcsec' % cell, '%.10farcsec' % cell]

    im.open(ms_path, usescratch=False, compress=False)
    if ra is not None and dec is not None:
        im.defineimage(nx=im_size[0], ny=im_size[1],
                       cellx=cell_size[0], celly=cell_size[1],
                       stokes='I', mode='mfs', step=1, spw=[-1],
                       outframe='', veltype='radio',
                       phasecenter=me.direction('J2000',
                                                '%.14fdeg' % ra,
                                                '%.14fdeg' % dec))
    else:
        im.defineimage(nx=im_size[0], ny=im_size[1],
                       cellx=cell_size[0], celly=cell_size[1],
                       stokes='I', mode='mfs', step=1, spw=[-1], outframe='',
                       veltype='radio')
    im.weight(type=weighting)
    if w_planes > 0:
        im.setoptions(ftmachine='wproject', wprojplanes=w_planes,
                      gridfunction='SF', padding=1.2,
                      dopbgriddingcorrections=True,
                      applypointingoffsets=False)
    else:
        im.setoptions(ftmachine='ft', gridfunction='SF', padding=1.2,
                      dopbgriddingcorrections=True, applypointingoffsets=False)

    t0 = time.time()
    print '*' * 80
    print '* Starting imaging...'
    im.makeimage(image=images['dirty'], type='observed', verbose=True)
    if niter > 0:
        im.clean(niter=niter, gain=gain, displayprogress=True,
                 model=images['model'],
                 residual=images['residual'],
                 image=images['restored'],
                 psfimage=images['psf'])
    im.close()
    for key in images.keys():
        if os.path.isdir(images[key]):
            ia.open(images[key])
            ia.tofits('%s.fits' % images[key], overwrite=True)
            ia.close()
    print '* Time taken to make images = %.3f s' % (time.time() - t0)
    print '*' * 80

    # Remove ipython log files.
    for f in glob.glob("ipython*.log"):
        os.remove(f)
