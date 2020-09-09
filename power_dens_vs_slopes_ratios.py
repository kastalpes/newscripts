"""
Make dens vs slopes and ratios.

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

for axes in ['x','y','z']:
    plot_dir = "/Users/Kye/512reruns/512rerunplots/general"
    prefix="power_dens_time"
    simlist=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
    framelist=[range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(65,85),range(11,31),range(11,31),range(72,91),range(56,75),range(20,40)]
    avg_dens_list = [[],[],[],[],[],[],[],[],[],[],[],[]]
    avg_k_list = [[],[],[],[],[],[],[],[],[],[],[],[]]
    avg_slopes_list = [0,0,0,0,0,0,0,0,0,0,0,0]
    avg_ee_slopes = [0,0,0,0,0,0,0,0,0,0,0,0]
    avg_bb_slopes = [0,0,0,0,0,0,0,0,0,0,0,0]
    avg_ratio_list = [0,0,0,0,0,0,0,0,0,0,0,0]
    projections=[]
    pos_k=slice(None)

    for i in range(12):
      frames=framelist[i]
      simdes=simlist[i]
      sim_dir = "/Users/Kye/512reruns/frbs/%s"%simdes
      dens=0
      kval=0
      avg_clee=0
      avg_clbb=0

      pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix)

      for frame in frames:
        print(frame,simdes)
        ds=pack.read_queb(frame=frame,ax=axes,bin_style='dx1')
        projections.append(ds)
        projections[-1].compute_harmonic_products()
        ds.read_spectra(frame)
        data=ds['dspec'][1][pos_k]
        k = ds['dspec'][0]
        this_k=(k*2*np.pi)#[1:]
        this_k = 0.5*(this_k[1:]+this_k[:-1])
        kval += this_k
        dens += np.abs(data[1:])
        avg_clee += ds['ClEE']
        avg_clbb += ds['ClBB']
        fitrange=ds.determine_fit_range() 
        slopes=ds.fit_eb_slopes(fitrange=fitrange)
      avg_clee /= len(frames)
      avg_clbb /= len(frames)
      kval/=len(frames)
      dens/=len(frames)
      avg_dens_list[i]=dens
      avg_k_list[i]=kval
      slopefind=list(plfit(ds.lcent,dens,fitrange))  #returns [slope,amp,res]
      avg_slopes_list[i]=slopefind[0]
      eeslopes=list(plfit(ds.lcent,avg_clee,fitrange))  #returns [slope,amp,res]
      bbslopes=list(plfit(ds.lcent,avg_clbb,fitrange))  #returns [slope,amp,res]
      avg_ee_slopes[i]=eeslopes[0]
      avg_bb_slopes[i]=bbslopes[0]
      avg_ratio=avg_clee/avg_clbb
      meanquan=np.mean(avg_ratio[2:30])
      avg_ratio_list[i]=meanquan
    xs=avg_slopes_list
    ys1=avg_ee_slopes
    ys2=avg_bb_slopes
    ys3=avg_ratio_list
    m1,b1=np.polyfit(xs,ys1,1)
    m2,b2=np.polyfit(xs,ys2,1)
    m3,b3=np.polyfit(xs,ys3,1)
    fig1,ax1=plt.subplots(1,1)
    fig2,ax2=plt.subplots(1,1)
    fig3,ax3=plt.subplots(1,1)
    ax1.plot(xs, ys1,'.') 
    y1=[i*m1 for i in xs]
    yy1=[i+b1 for i in y1]
    ax1.plot(xs,yy1,'-',label=np.round(m1,decimals=3))
    ax1.set_xlabel(r'Density Slope')
    ax1.set_ylabel(r'ClEE Slope')
    ax1.set_title(r'Density EE Parameter Space')
    ax1.legend(loc=0)
    fig1.savefig("%s/ee_density_slopes_parspace_%s.png"%(plot_dir,axes))
    ax2.plot(xs, ys2,'.')
    y2=[i*m2 for i in xs]
    yy2=[i+b2 for i in y2]
    ax2.plot(xs,yy2,'-',label=np.round(m2,decimals=3))
    ax2.set_xlabel(r'Density Slope')
    ax2.set_ylabel(r'ClBB Slope')
    ax2.set_title(r'Density BB Parameter Space')
    ax2.legend(loc=0)
    fig2.savefig("%s/bb_density_slopes_parspace_%s.png"%(plot_dir,axes))
    ax3.plot(xs, ys3,'.')
    y3=[i*m3 for i in xs]
    yy3=[i+b3 for i in y3]
    ax3.plot(xs,yy3,'-',label=np.round(m3,decimals=3))
    ax3.set_xlabel(r'Density Slope')
    ax3.set_ylabel(r'EE/BB Ratio')
    ax3.set_title(r'Density E/B Ratio Parameter Space')
    ax3.legend(loc=0)
    fig3.savefig("%s/density_ratio_parspace_%s.png"%(plot_dir,axes))
