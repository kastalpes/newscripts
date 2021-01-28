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
    hdu = pyfits.PrimaryHDU(array)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(fname,overwrite=True)

for fname in glob.glob('image*.out'):
    print("Reading",fname)

    a = readImage(fname)
    angle = int(fname.split("_")[1].split(".")[0])
    print(angle)
    for stokes in ['Q','U','I']:
        index = 'IQU'.index(stokes)
        plotImage(a,cmap=cm.hot,au=True,bunit='inu', stokes=stokes)
        plotpoldir(a.x/au,a.y/au,a.image[:,:,0,0],a.image[:,:,1,0],a.image[:,:,2,0])
        oname = fname.split('.')[0]+"_%s.png"%(stokes)
        oname = "ca02_units_DD%04d_%s_th%04d.png"%(490,stokes,angle)
        plt.savefig(oname)


        write_fits("ca02_units_DD%04d_%s_th%04d.fits"%(490,stokes,angle), a.image[:,:,index,0])








