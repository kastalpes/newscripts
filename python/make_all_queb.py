from GL import *
import queb3 as fr
reload(fr)
frames=list(range(0,500,10)) + list(range(0,10))
frames=[499]
sim_dir = "/scratch1/dcollins/Paper49_EBQU/ca02_turb"
plot_dir =  /home/dcollins4096/PigPen
pack = fr.simulation_package( directory=sim_dir,frames=frames,prefix="ca02", BoxSize=2)
#pack.EBall()
for frame in frames:
    this_proj=pack.read_queb(frame=frame,ax='x')
    this_proj.read_spectra(frame)
    this_proj.compute_harmonic_products()
    fitrange=this_proj.determine_fit_range()  #something better 
    slopes=this_proj.fit_eb_slopes(fitrange=fitrange)
    pack.plot_eb(this_proj, fname='%s/ca02_reb_x_n%04d.png'%(plot_dir,frame),slopes=slopes)
