"""
Make E plot and slope using 2 methods.

*  queb3.simulation_package points to a simulation.  
   BoxSize is in units of 64 zones.
*  
      
"""
from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
reload(queb3)

for axes in ['x','y','z']:
    shortprefix="time"
    simlist=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
    framelist=[range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(65,85),range(11,31),range(11,31),range(72,91),range(56,75),range(20,40)]
    ee_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    bb_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    eb_ratios=[0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(12):
      frames=framelist[i]
      simdes=simlist[i]
      sim_dir = "/Users/Kye/512reruns/frbs/%s"%simdes
      plot_dir = "/Users/Kye/512reruns/512rerunplots/%s"%simdes
      avg_clee=0
      avg_clbb=0
      fig1,ax1=plt.subplots(1,1)
      fig2,ax2=plt.subplots(1,1)
      projections=[]
      longprefix='%s_%s_%s'%(simdes,shortprefix,axes)
      pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=longprefix)

      for frame in frames:
        print(frame,simdes)
        proj=pack.read_queb(frame=frame,ax=axes,bin_style='dx1')
        projections.append(proj)
        projections[-1].compute_harmonic_products()
#    ax.plot(proj.lcent,proj['ClEE'])
        avg_clee += proj['ClEE']
        avg_clbb += proj['ClBB']
        fitrange=proj.determine_fit_range()  #something better
#        pack.make_spectra(frame=frame)
        proj.read_spectra(frame)
        pack.plot_many_spectra(proj,simdes,fname='%s/%s_spectra_%04d_%s.png'%(plot_dir,simdes,frame,axes))
    #do the fit.  Not much fancy here
        slopes=proj.fit_eb_slopes(fitrange=fitrange)
      avg_clee /= len(frames)
      avg_clbb /= len(frames)
      ##------slope plot
      avg_slopes_ee=list(plfit(proj.lcent,avg_clee,fitrange))  #returns [slope,amp,res]
      ellmask = (fitrange[0] < proj.lcent)*(proj.lcent < fitrange[1])
      xslope=proj.lcent[ellmask]
      yslope=avg_slopes_ee[1]*xslope**(avg_slopes_ee[0])
      ax1.plot(xslope,yslope,label=r'$\alpha^{EE}$ %0.2f'%avg_slopes_ee[0],c='r')
      avg_slopes_bb=list(plfit(proj.lcent,avg_clbb,fitrange))  #returns [slope,amp,res]
      yslope=avg_slopes_bb[1]*xslope**(avg_slopes_bb[0])
      ax1.plot(xslope,yslope,label=r'$\alpha^{BB}$ %0.2f'%avg_slopes_bb[0],c='k')
      ax1.legend(loc=0)
      ##------
      ##------ratio plot
      avg_ratio = avg_clee/avg_clbb
      mean_quantity = np.mean(avg_ratio[2:30]) #was 4:10
      ymean = (np.zeros(np.size(xslope))+mean_quantity)
      ax2.plot(proj.lcent, avg_ratio,label='E/B ratio',c='k')
      ax2.plot(xslope, ymean,label='average ratio %0.2f'%mean_quantity,c='r')
      ax2.set_title('ratio %0.2f'%mean_quantity)
      ax2.legend(loc=0)
      ##------
      ee_slopes[i]=avg_slopes_ee[0]
      bb_slopes[i]=avg_slopes_bb[0]
      eb_ratios[i]=mean_quantity
      ee_slopes=np.array(ee_slopes)
      bb_slopes=np.array(bb_slopes)
      eb_ratios=np.array(eb_ratios)
      print('ee_slopes =')
      print(ee_slopes)
      ##------
      ax1.plot(proj.lcent,avg_clee,label='ClEE',c='g')
      ax1.plot(proj.lcent,avg_clbb,label='ClBB',c='b')
      dt.axbonk(ax1,xlabel='k',ylabel='ClEE and ClBB',xscale='log',yscale='log')
      #ax1.set_ylim((10**(-6)),300)
      fig1.savefig("%s/%s_cl_average_new.png"%(plot_dir,longprefix))
#      dt.axbonk(ax2,xlabel='k',ylabel='ratio',xscale='log',yscale='log')
      dt.axbonk(ax2,xlabel='k',ylabel='ratio',xscale='log')
      ax2.set_ylim(0,10)
      fig2.savefig("%s/%s_ratio_average.png"%(plot_dir,longprefix))
      fig1.clf()
      fig2.clf()
    Fptr = h5py.File("slopes_ratios_dict_%s.h5"%axes,'w')
    Fptr.create_dataset("ee_slopes", data=ee_slopes)
    Fptr.create_dataset("bb_slopes", data=bb_slopes)
    Fptr.create_dataset("eb_ratios", data=eb_ratios)
    Fptr.close()
