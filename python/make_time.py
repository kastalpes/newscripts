"""
Make E plot and slope using 2 methods.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
import queb3
reload(queb3)
frames=list(range(0,500,10)) + list(range(0,10))
frames=[499]
sim_dir = "/scratch1/dcollins/Paper49_EBQU/ca02_turb"
plot_dir =  "/home/dcollins4096/PigPen"
prefix='slope_snag'

pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix, BoxSize=2)
#produce all QUEB products.
#pack.EBall()
for frame in frames:
    proj_dx1 =pack.read_queb(frame=frame,ax='x',bin_style='dx1')
    proj_5deg=pack.read_queb(frame=frame,ax='x',bin_style='5deg')

    proj_dx1.compute_harmonic_products()
    proj_5deg.compute_harmonic_products()
    fitrange_dx1=proj_dx1.determine_fit_range()  #something better 
    fitrange_5deg=proj_5deg.determine_fit_range()  #something better 
    #do the fit.  Not much fancy here
    slopes_dx1=proj_dx1.fit_eb_slopes(fitrange=fitrange_dx1)
    slopes_5deg=proj_5deg.fit_eb_slopes_2(fitrange=fitrange_5deg)
    #this plotting package can be generalized.
    pack.plot_2eb(proj_dx1,proj_5deg, 
                  fname='%s/%s_reb_x_n%04d.png'%(plot_dir,prefix,frame),
                  slopes0=slopes_dx1, slopes1=slopes_5deg)
