"""
Make h power spectra time averaged plots.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
reload(queb3)

plot_dir = "/Users/Kye/512reruns/512rerunplots/general"
prefix="power_h_time"
simlist=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
framelist=[range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(65,85),range(11,31),range(11,31),range(72,91),range(56,75),range(20,40)]
avg_h_list = [[],[],[],[],[],[],[],[],[],[],[],[]]
avg_k_list = [[],[],[],[],[],[],[],[],[],[],[],[]]
avg_slopes_list = [0,0,0,0,0,0,0,0,0,0,0,0]
pos_k=slice(None)
for i in range(12):
  frames=framelist[i]
  simdes=simlist[i]
  sim_dir = "/Users/Kye/512reruns/frbs/%s"%simdes
  h=0
  kval=0

  pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix)

  for frame in frames:
    print(frame,simdes)
    ds=pack.read_queb(frame=frame,ax='x',bin_style='dx1')
    ds.read_spectra(frame)
    data=ds['hspec'][1][pos_k]
    k = ds['hspec'][0]
    this_k=(k*2*np.pi)#[1:]
    this_k = 0.5*(this_k[1:]+this_k[:-1])
    kval += this_k
    h += np.abs(data[1:])
    fitrange=ds.determine_fit_range()
  kval/=len(frames)
  h/=len(frames)
  avg_h_list[i]=h
  avg_k_list[i]=kval

  slopefind=list(plfit(ds.lcent,h,fitrange))  #returns [slope,amp,res]
  avg_slopes_list[i]=slopefind[0]                           
  #ellmask = (fitrange[0] < ds.lcent)*(ds.lcent < fitrange[1])
  #xslope=ds.lcent[ellmask]                             
  #yslope=avg_slopes[1]*xslope**(avg_slopes[0])
  #ax1.plot(xslope,yslope,label=r'$\alpha^{EE}$ %0.2f'%avg_slo
fig,ax=plt.subplots(1,1)
ax.set_prop_cycle(cycler(color=['r','g','b','r','g','b','r','g','b','r','g','b'])+cycler(marker=['.','.','.','x','x','x','*','*','*','1','1','1']))
for i in range(12):
    xs=avg_k_list[i]
    ys=avg_h_list[i]
    ls=simlist[i]
    ax.plot(xs, ys,ms=5,label='{simulation}={slopes}'.format(simulation=ls,slopes=round(avg_slopes_list[i],3)),linestyle='--') 

#axis.set_xlim([1,4.5])
#axis.set_ylim(1e-7,50)
#axis.set_title('Parameter Space of Averages')
    ax.legend(loc=1)
    ax.set_xlabel(r'k/k1')
    ax.set_ylabel(r'Average Power of $h$')
    ax.set_yscale('log')
    ax.set_xscale('log')
    fig.savefig("%s/power_h_avg.png"%(plot_dir))
    plt.clf()
    plt.close()
