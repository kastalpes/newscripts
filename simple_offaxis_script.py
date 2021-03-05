"""
Make simple off-axis plot using data
      
"""
from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import yt
reload(queb3)

simdir="half_half"
frame=30
plotdir="/Users/Kye/512reruns/512rerunplots/offaxisplots"

ds = yt.load("/Users/Kye/512reruns/frbs/half_half/DD%04d/data%04d"%(frame,frame))
L = [1,1,0] # vector normal to cutting plane
north_vector = [-1,1,0]
W = [0.02, 0.02, 0.02]
c = [0.5, 0.5, 0.5]
N = 512
image = yt.off_axis_projection(ds, c, L, W, N, "density")
yt.write_image(np.log10(image), "%s/%s/%04d_offaxis_density_projection.png" %(plotdir,simdir,frame))
