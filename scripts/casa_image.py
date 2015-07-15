# -*- coding: utf-8 -*-
"""Simple CASA imaging script.

Run with:
    casapy --nogui --nologger --log2term -c casa_image.py <ms>

See:
    http://casa.nrao.edu/docs/CasaRef/CasaRef.html
for documentation of methods on CASA objects 'im' and 'ia'
"""

import os
import time
import sys

# -------------------------------------
ms_path = os.path.abspath(sys.argv[-1])
ms_name = os.path.basename(ms_path)
image_root_name = os.path.splitext(ms_name)[0]
# im_size = [8192, 32768]
im_size = [1024, 1024]
cell_size = ['0.75arcsec', '0.75arcsec']
w_planes = 256
make_psf = False
# -------------------------------------


print '=' * 80
print 'ms = %s' % ms_path
print '=' * 80

im.open(ms_path, usescratch=False, compress=False)
im.defineimage(nx=im_size[0], ny=im_size[1],
               cellx=cell_size[0], celly=cell_size[1],
               stokes='I', mode='mfs', step=1, spw=[-1], outframe='',
               veltype='radio')
im.weight(type='natural')
im.setoptions(ftmachine='wproject', wprojplanes=w_planes, gridfunction='SF',
              padding=1.2, dopbgriddingcorrections=True,
              applypointingoffsets=False)
dirty_image = image_root_name + '_dirty'
t0 = time.time()
print '*' * 80
print '* Starting imaging...'
im.makeimage(image=dirty_image + '.img', type='observed', verbose=True)
print '* Time taken to make dirty image = %.3f s' % (time.time() - t0)
print '*' * 80
if make_psf:
    psf_image = image_root_name + '_psf'
    im.makeimage(image=psf_image + '.img', type='psf', verbose=True)
im.close()
ia.open(dirty_image + '.img')
ia.tofits(dirty_image + '.fits', overwrite=True)
ia.close()
if make_psf:
    ia.open(psf_image + '.img')
    ia.tofits(psf_image + '.fits', overwrite=True)
    ia.close()
