"""
Make many queb data products.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
import queb3
reload(queb3)
frames=list(range(0,500,10)) + list(range(0,10))
frames=[1]
sim_dir = "/scratch1/dcollins/Paper49_EBQU/ca16_cyl2/"
plot_dir =  "/home/dcollins4096/PigPen"

#a thing that describes the simulation
pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix="ca16_dcpipe_")
#produce all QUEB products.
pack.EBall()
for frame in frames:
    #makes 3d spectra
    #pack.make_spectra(frame=frame)
    #read Q&U.  Compute E&B.
    for axis in 'xyz':
        this_proj=pack.read_queb(frame=frame,ax=axis)
        this_proj.compute_harmonic_products()
        pack.image_fields(frame,ts=this_proj, axis=axis)
    #pack.image_fields(frame,axis="0",ts=this_proj,field_list=['Q','U','E','B','T'],theta=theta)
    #read 3d spectra    
    #this_proj.read_spectra(frame)
    #here we will make improvements.
    #Also here we can vary the fit range.
    #Further, we can average this_proj.ClEE over several frames to 
    #get a more meaningful result.
#    fitrange=this_proj.determine_fit_range()  #something better 
#    #do the fit.  Not much fancy here
#    slopes=this_proj.fit_eb_slopes(fitrange=fitrange)
#    #this plotting package can be generalized.
#    pack.plot_eb(this_proj, fname='%s/ca02_reb_x_n%04d.png'%(plot_dir,frame),slopes=slopes)
