#!/usr/bin/env python
#Why not /usr/bin/python ... see: http://goo.gl/KFjTo

#
# casapy --nogui --nologger --log2term -c casa_make_image.py
#

default(clean)

name='test_1'
vis=name + '.ms'
imagename=name
niter = 0
imsize=4096
#! cell=['9.0arcsec']
cell=['12.0arcsec']
#! Enable w-projection by having the following 3 lines.
gridmode='widefield'
wprojplanes=32
imagename+='_wproj_%03i_%03i' % (wprojplanes, imsize)

t0 = time.time()
go()
print '='*60
print "Time taken to image '%s' = %.3f s" % (vis, time.time()-t0)
print '='*60

exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits')
exportfits(imagename=imagename+'.psf',fitsimage=imagename+'.psf.fits')

import shutil
import os
for ext in ('.image','.psf','.flux','.model','.residual'):
    if os.path.isdir(imagename+ext):
        shutil.rmtree(imagename+ext)
