"""
Make many queb data products.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
from davetools import *
import queb3
reload(queb3)
frames=[490]
sim_dir = "/scratch1/dcollins/Paper49d/radmc3d_pol_3/"
plot_dir =  "./"

#a thing that describes the simulation
pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix="ca02_radmc", 
                                dataset_name='ca02_units_DD',frbname="./",plotdir=plot_dir)
#produce all QUEB products.
#pack.EBall()
angles=rainbow_map(360)
for frame in frames:
    #read Q&U.  Compute E&B.
    fig,ax=plt.subplots(1,1)
    for theta in range(0,370,10):
        this_proj=pack.read_queb(frame=frame,theta=theta)
        this_proj.compute_harmonic_products()
        pack.image_fields(frame,axis="0",ts=this_proj,field_list=['Q','U','E','B','T'],theta=theta)
        ax.hist(this_proj['E'].flatten(),histtype='step',color=angles(theta),bins=100)
    fig.savefig("%s/angle_hist.pdf"%plot_dir)

    continue

    #here we will make improvements.
    #Also here we can vary the fit range.
    #Further, we can average this_proj.ClEE over several frames to 
    #get a more meaningful result.
    fitrange=this_proj.determine_fit_range()  #something better 
    #do the fit.  Not much fancy here
    slopes=this_proj.fit_eb_slopes(fitrange=fitrange)
    #this plotting package can be generalized.
    pack.plot_eb(this_proj, fname='%s/ca02_reb_x_n%04d.png'%(plot_dir,frame),slopes=slopes)
