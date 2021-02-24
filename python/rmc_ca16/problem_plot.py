#
# First do, for example:
#   radmc3d mctherm
#   radmc3d image lambda 0.4e3 theta 45 phi 45 nostar
#
# or viewed at an angle:
#
#   radmc3d image lambda 0.4e3 theta 45 phi 45 nostar
#
# Then:
#
#   ipython --matplotlib
#   %run problem_plot
#
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import math
import numpy as np
from radmc3dPy.image import *
from plotpoldir import *
import glob
import astropy.io.fits as pyfits
au  = 1.49598e13     # Astronomical Unit       [cm]

def write_fits(fname,array):
    hdu = pyfits.PrimaryHDU(array.swapaxes(0,1))
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(fname,overwrite=True)

for fname in glob.glob('*image*.out'):
    print("Reading",fname)

    a = readImage(fname)
    theta = int(fname.split(".")[0].split("_")[1])
    phi   = int(fname.split(".")[0].split("_")[2])
    print(theta,phi)
    for stokes in ['Q','U','I']:
        index = 'IQU'.index(stokes)
        plotImage(a,cmap=cm.hot,au=True,bunit='inu', stokes=stokes)
        plotpoldir(a.x/au,a.y/au,a.image[:,:,0,0],a.image[:,:,1,0],a.image[:,:,2,0])
        oname = fname.split('.')[0]+"_%s.png"%(stokes)
        oname =  "ca16_cyl_DD%04d_th%04d_ph%04d_%s.png"% (1,theta,phi,stokes)
        onamef = "ca16_cyl_DD%04d_th%04d_ph%04d_%s.fits"%(1,theta,phi,stokes)
        plt.savefig(oname)


        write_fits(onamef, a.image[:,:,index,0])








