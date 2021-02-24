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
frames=[1]
sim_dir = "rmc_ca16"
plot_dir =  "./"

#a thing that describes the simulation
#ca14_cyl_DD0001_th0000_ph0000_I.fits
pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix="ca16_cyl_rmc_pipe",
                                dataset_name='ca16_cyl_DD',frbname="./",plotdir=plot_dir)
#produce all QUEB products.
#pack.EBall()
angles=rainbow_map(3)
for frame in frames:
    #read Q&U.  Compute E&B.
    fig,ax=plt.subplots(1,1)
    theta_phi =  [[0,0,'zish'], [90,-90,'xish'], [90,180,'yish']]
    theta_phi =  [[90,-90,'xish']]
    #theta_phi = [90,180]
    for theta,phi,lab in theta_phi:
        this_proj=pack.read_queb(frame=frame,theta=theta,phi=phi)
        this_proj.compute_harmonic_products()
        pack.image_fields(frame,axis="0",ts=this_proj,field_list=['Q','U','E','B','T'],theta=theta,phi=phi)

        ax.hist(this_proj['E'].flatten(),histtype='step',bins=100,label=lab)
    ax.legend(loc=0)
    fig.savefig("%s/%s_E_hist.pdf"%(plot_dir,pack.prefix))

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
