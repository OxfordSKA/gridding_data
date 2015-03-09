#!/usr/bin/env python

#
# casapy --nogui --nologger --log2term -c casa_make_image.py <vis>
#

default(clean)

vis=sys.argv[-1]
imagename=name
niter = 0
imsize=2048
#! cell=['9.0arcsec']
cell=['12.0arcsec']
#! Enable w-projection by including the following 3 lines.
#gridmode='widefield'
#wprojplanes=128
#imagename+='_wproj_%03i_%03i' % (wprojplanes, imsize)

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
