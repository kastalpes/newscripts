"""
Test binning for spectra fitting.

Edit queb3.compute_bins_horse_around
to change the binning.
      
"""
from GL import *
import queb3
reload(queb3)
frames=list(range(0,500,10)) + list(range(0,10))
frames=[200]
sim_dir  = "/scratch1/dcollins/Paper49d/za01_M10_MA1_64_quan"
plot_dir =  "%s/PlotsToMove"%os.environ['HOME']
plot_dir =  "%s/PigPen"%os.environ['HOME']
prefix='za01'
#a thing that describes the simulation
pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix)
#produce all QUEB products.
#pack.EBall()
ax='x'
for frame in frames:
    #read Q&U.
    #dx1 is probably the right thing, horse is for testing
    proj_dx1 =pack.read_queb(frame=frame,ax=ax,bin_style='dx1')
    proj_test=pack.read_queb(frame=frame,ax=ax,bin_style='horse')

    proj_dx1.compute_harmonic_products()
    proj_test.compute_harmonic_products()
    fitrange_dx1=proj_dx1.determine_fit_range()  
    fitrange_test=proj_test.determine_fit_range()

    slopes_dx1=proj_dx1.fit_eb_slopes_2(fitrange=fitrange_dx1)
    slopes_test=proj_test.fit_eb_slopes_2(fitrange=fitrange_test)

    proj_dx1.scatter_clee(
                  fname='%s/%s_dx1_scatter_n%04d_%s.png'%(plot_dir,prefix,frame,ax))
    proj_test.scatter_clee(
                  fname='%s/%s_bin_test_scatter_n%04d_%s.png'%(plot_dir,prefix,frame,ax))

