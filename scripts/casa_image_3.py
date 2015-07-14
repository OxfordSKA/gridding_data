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
import shutil


# -------------------------------------
ms = os.path.basename(os.path.abspath(sys.argv[-1]))
image_root_name = os.path.splitext(ms)[0]
im_size = [512, 512]
cell_size = ['2.0arcsec', '2.0arcsec']
w_planes = 0
make_psf = False
grid_function = 'SF'  # SF | BOX
# -------------------------------------

im.open(ms, usescratch=False, compress=False)
im.defineimage(nx=im_size[0], ny=im_size[1],
               cellx=cell_size[0], celly=cell_size[1],
               stokes='I', mode='mfs', step=1, spw=[-1], outframe='',
               veltype='radio')
im.weight(type='natural')
if w_planes > 0:
    im.setoptions(ftmachine='wproject', wprojplanes=w_planes,
                  gridfunction=grid_function,
                  padding=1.2, dopbgriddingcorrections=True,
                  applypointingoffsets=False)
else:
    im.setoptions(ftmachine='ft', gridfunction=grid_function, padding=1.2,
                  dopbgriddingcorrections=True, applypointingoffsets=False)
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
shutil.rmtree(dirty_image + '.img')

if make_psf:
    ia.open(psf_image + '.img')
    ia.tofits(psf_image + '.fits', overwrite=True)
    ia.close()
    shutil.rmtree(psf_image + '.img')
