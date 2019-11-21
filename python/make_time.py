"""
Make E plot and slope using 2 methods.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
import queb3
import davetools as dt
reload(queb3)
frames=range(400,500,20)
sim_dir = "/scratch1/dcollins/Paper49_EBQU/ca02_turb"
plot_dir =  "/home/dcollins4096/PigPen"
prefix='slope_snag'


pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix)

#produce all QUEB products.
#pack.EBall()

fig,ax=plt.subplots(1,1)
projections=[]
avg_clee = 0
for frame in frames:
    proj=pack.read_queb(frame=frame,ax='x',bin_style='dx1')
    projections.append(proj)
    projections[-1].compute_harmonic_products()
    ax.plot(proj.lcent,proj['ClEE'])
    avg_clee += proj['ClEE']
avg_clee /= len(frames)

ax.plot(proj.lcent,avg_clee,c='k')
dt.axbonk(ax,xlabel='k',ylabel='<ClEE>',xscale='log',yscale='log')
fig.savefig(os.environ['HOME']+"/PigPen/clee_average.png")
