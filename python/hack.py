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
from matplotlib.pyplot import cm
import sim_colors
reload(queb3)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)
for axes in ['x','y','z']:
    plt.close('all')
    shortprefix="time"
    simlist=sim_colors.simlist
    framelist=[range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(65,85),range(11,31),range(11,31),range(72,91),range(56,75),range(20,40)]
    ee_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    bb_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    te_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    tb_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    eb_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    eb_ratios=[0,0,0,0,0,0,0,0,0,0,0,0]
    te_ratios=[0,0,0,0,0,0,0,0,0,0,0,0]
    tb_ratios=[0,0,0,0,0,0,0,0,0,0,0,0]
    v_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    h_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    d_slopes=[0,0,0,0,0,0,0,0,0,0,0,0]
    reb=[0,0,0,0,0,0,0,0,0,0,0,0]
    figratio,axratio=plt.subplots(1,1)
    fige,axe=plt.subplots(1,1)
    figb,axb=plt.subplots(1,1)
    figt,axt=plt.subplots(1,1)
    figv,axv=plt.subplots(1,1)
    figh,axh=plt.subplots(1,1)
    figd,axd=plt.subplots(1,1)
    figte,axte=plt.subplots(1,1)
    figtb,axtb=plt.subplots(1,1)
    figeb,axeb=plt.subplots(1,1)
    figreb,axreb=plt.subplots(1,1)

    axratio.set_prop_cycle(cycler(color  =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axe.set_prop_cycle(cycler(color  =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axb.set_prop_cycle(cycler(color  =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axt.set_prop_cycle(cycler(color  =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axv.set_prop_cycle(cycler(color  =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axh.set_prop_cycle(cycler(color  =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axd.set_prop_cycle(cycler(color  =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axte.set_prop_cycle(cycler(color =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axtb.set_prop_cycle(cycler(color =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axeb.set_prop_cycle(cycler(color =sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))
    axreb.set_prop_cycle(cycler(color=sim_colors.color_list)+  cycler(linestyle=sim_colors.line_list))

    for i in range(12):
        frames=framelist[i]
        simdes=simlist[i]
        spectra_fname = "avg_spectra_%s_%s.h5"%(simdes,axes)
        sim_dir = "/archive2/kas14d/512reruns/frbs/%s"%simdes
        plot_dir = "./plots"
        gen_dir = "./plots"
        avg_clee=0
        avg_clbb=0
        avg_cleb=0
        avg_clte=0
        avg_cltb=0
        avg_cltt=0
        avg_v=0
        avg_h=0
        avg_d=0
        avg_reb=0
        fig1,ax1=plt.subplots(1,1)
        fig2,ax2=plt.subplots(1,1)
        fig3,ax3=plt.subplots(1,1)
        projections=[]
        longprefix='%s_%s_%s'%(simdes,shortprefix,axes)
        pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=longprefix)
        nplots=0
        proj=pack.read_queb(frame=frames[0],ax=axes,bin_style='dx1')
        fitrange=proj.determine_fit_range()  #something better
        if not os.path.exists(spectra_fname):

            for frame in frames: 
                
            #print(frame,simdes)
                proj=pack.read_queb(frame=frame,ax=axes,bin_style='dx1')
                if proj is None:
                    continue
                nplots+=1

                proj.read_spectra(frame)
                projections.append(proj)
                projections[-1].compute_harmonic_products()
                #    ax.plot(proj.lcent,proj['ClEE'])
                if verbose:
                    print('proj[ClEE]')
                    print(proj['ClEE'])
                avg_clee += proj['ClEE']
                avg_clbb += proj['ClBB']
                avg_cleb += proj['ClEB']
                avg_clte += proj['ClTE']
                avg_cltb += proj['ClTB']
                avg_cltt += proj['ClTT']
                avg_v1 = proj['vspec'][1]
                avg_v += np.abs(avg_v1[1:])
                avg_h1 = proj['hspec'][1]
                avg_h += np.abs(avg_h1[1:])
                avg_d1 = proj['dspec'][1]
                avg_d += np.abs(avg_d1[1:])
                avg_reb += proj['ClEB']/((proj['ClEE']*proj['ClBB'])**0.5)
                fitrange=proj.determine_fit_range()  #something better
                #        pack.make_spectra(frame=frame)
                #       proj.read_spectra(frame)
                #        pack.plot_many_spectra(proj,simdes,fname='%s/%s_spectra_%04d_%s.png'%(plot_dir,simdes,frame,axes))
                #do the fit.  Not much fancy here
            slopes=proj.fit_eb_slopes(fitrange=fitrange)
            avg_clee /= nplots
            avg_clbb /= nplots
            avg_cleb /= nplots
            avg_clte /= nplots
            avg_cltb /= nplots
            avg_cltt /= nplots
            avg_v    /= nplots
            avg_h    /= nplots
            avg_d    /= nplots
            avg_reb  /= nplots
            fptr = h5py.File(spectra_fname,"w")
            try:
                fptr.create_dataset("avg_clee",data = avg_clee)
                fptr.create_dataset("avg_clbb",data = avg_clbb)
                fptr.create_dataset("avg_cleb",data = avg_cleb)
                fptr.create_dataset("avg_clte",data = avg_clte)
                fptr.create_dataset("avg_cltb",data = avg_cltb)
                fptr.create_dataset("avg_cltt",data = avg_cltt)
                fptr.create_dataset("avg_v"   ,data = avg_v   )
                fptr.create_dataset("avg_h"   ,data = avg_h   )
                fptr.create_dataset("avg_d"   ,data = avg_d   )
                fptr.create_dataset("avg_reb" ,data = avg_reb )
            except:
                raise
            finally:
                fptr.close()
        else:
            fptr=h5py.File(spectra_fname,"r")
            try:
                avg_clee = fptr["avg_clee"][()]
                avg_clbb = fptr["avg_clbb"][()]
                avg_cleb = fptr["avg_cleb"][()]
                avg_clte = fptr["avg_clte"][()]
                avg_cltb = fptr["avg_cltb"][()]
                avg_cltt = fptr["avg_cltt"][()]
                avg_v    = fptr["avg_v"   ][()]
                avg_h    = fptr["avg_h"   ][()]
                avg_d    = fptr["avg_d"   ][()]
                avg_reb  = fptr["avg_reb" ][()]

            except:
                raise
            finally:
                fptr.close()
        avg_slopes_ee=list(plfit(proj.lcent,avg_clee,fitrange))  #returns [slope,amp,res]
        avg_slopes_bb=list(plfit(proj.lcent,avg_clbb,fitrange))  #returns [slope,amp,res]
        avg_slopes_tt=list(plfit(proj.lcent,avg_cltt,fitrange))  #returns [slope,amp,res]
        all_slopes['ee'].append(avg_slopes_ee[0])
        all_slopes['bb'].append(avg_slopes_bb[0])
        all_slopes['tt'].append(avg_slopes_tt[0])
        if 0:
            ##------ee bb tt slope plot
            ellmask = (fitrange[0] < proj.lcent)*(proj.lcent < fitrange[1])
            xslope=proj.lcent[ellmask]
            yslope=avg_slopes_ee[1]*xslope**(avg_slopes_ee[0])
            ax1.plot(xslope,yslope,label=r'$\alpha^{EE}$ %0.2f'%avg_slopes_ee[0],c='m')
            yslope=avg_slopes_bb[1]*xslope**(avg_slopes_bb[0])
            ax1.plot(xslope,yslope,label=r'$\alpha^{BB}$ %0.2f'%avg_slopes_bb[0],c='c')
            yslope=avg_slopes_tt[1]*xslope**(avg_slopes_tt[0])
            ax1.plot(xslope,yslope,label=r'$\alpha^{TT}$ %0.2f'%avg_slopes_tt[0],c='k')
            ax1.plot(proj.lcent,avg_clee,label='ClEE',c='r')
            ax1.plot(proj.lcent,avg_clbb,label='ClBB',c='b')
            ax1.plot(proj.lcent,avg_cltt,label='ClTT',c='y')
            #      ax1.plot(proj.lcent,avg_cltb,label='ClTB',c='c')
            dt.axbonk(ax1,xlabel='k',ylabel='ClXX Power',xscale='log',yscale='log')
            #ax1.set_ylim((10**(-6)),300)
            ax1.legend(loc=0)
            fig1.savefig("%s/%s_Cl_XX_average.pdf"%(plot_dir,longprefix))
            fig1.clf()
        ##-------te tb symlog plot
        avg_slopes_te=list(plfit(proj.lcent,avg_clte,fitrange))  #returns [slope,amp,res]
        avg_slopes_tb=list(plfit(proj.lcent,avg_cltb,fitrange))  #returns [slope,amp,res]
        avg_slopes_eb=list(plfit(proj.lcent,avg_cleb,fitrange))  #returns [slope,amp,res]
        if 0:
            yslope=avg_slopes_te[1]*xslope**(avg_slopes_te[0])
            ax3.plot(xslope,yslope,label=r'$\alpha^{TE}$ %0.2f'%avg_slopes_te[0],c='m')
            ax3.plot(xslope,yslope,label=r'$\alpha^{TB}$ %0.2f'%avg_slopes_tb[0],c='c')
            yslope=avg_slopes_eb[1]*xslope**(avg_slopes_eb[0])
            ax3.plot(xslope,yslope,label=r'$\alpha^{EB}$ %0.2f'%avg_slopes_eb[0],c='k')
            ax3.plot(proj.lcent,avg_clte,label='ClTE',c='r')
            ax3.plot(proj.lcent,avg_cltb,label='ClTB',c='b')
            ax3.plot(proj.lcent,avg_cleb,label='ClEB',c='y')
            dt.axbonk(ax3,xlabel='k',ylabel='ClXY Power',xscale='log',yscale='symlog')
            ax3.set_yscale('symlog',linthreshy=1e-6)
            ax3.set_ylim(-1e4,1e4)
            #ax3.set_ylim((10**(-6)),300)
            ax3.legend(loc=0)
            fig3.savefig("%s/%s_Cl_XY_average.pdf"%(plot_dir,longprefix))
      ##------
      ##------ratio plo
        if 0:
            avg_ratio = avg_clee/avg_clbb
            mean_quantity = np.mean(avg_ratio[2:30]) #was 4:10
            ymean = (np.zeros(np.size(xslope))+mean_quantity)
            ax2.plot(proj.lcent, avg_ratio,label='EE/BB ratio',c='r')
            ax2.plot(xslope, ymean,label='average EB ratio %0.2f'%mean_quantity,c='m')
            avg_ratio2 = avg_clte/avg_clee
            mean_quantity2 = np.mean(avg_ratio2[2:30]) #was 4:10
            ymean2 = (np.zeros(np.size(xslope))+mean_quantity2)
            ax2.plot(proj.lcent, avg_ratio2,label='TE/EE ratio',c='y')
            ax2.plot(xslope, ymean2,label='average TE ratio %0.2f'%mean_quantity2,c='k')
            avg_ratio3=avg_cltb/avg_clbb
            mean_quantity3 = np.mean(avg_ratio3[2:30]) #was 4:10
            ymean3 = (np.zeros(np.size(xslope))+mean_quantity3)
            ax2.plot(proj.lcent, avg_ratio3,label='T/B ratio',c='b')
            ax2.plot(xslope, ymean3,label='average TB ratio %0.2f'%mean_quantity3,c='c')
            ax2.set_title('Cl Power Ratios')
            ax2.legend(loc=0)
            dt.axbonk(ax2,xlabel='k',ylabel='Cl Power Ratios',xscale='log')
            ax2.set_ylim(0,10)
            fig2.savefig("%s/%s_Cl_XY_ratio_average.pdf"%(plot_dir,longprefix))
            fig2.clf()
        ##------12 panel plots
        if 1:
            avg_slopes_v=list(plfit(proj.lcent,avg_v,fitrange))  #returns [slope,amp,res]
            avg_slopes_h=list(plfit(proj.lcent,avg_h,fitrange))  #returns [slope,amp,res]
            avg_slopes_d=list(plfit(proj.lcent,avg_d,fitrange))  #returns [slope,amp,res]
            all_slopes['v'].append(avg_slopes_v[0])
            all_slopes['h'].append(avg_slopes_h[0])
            all_slopes['d'].append(avg_slopes_d[0])
            axratio.plot(proj.lcent,avg_clee/avg_clbb,label='%s,%0.2f'%(simdes,avg_slopes_ee[0]))
            axe.plot(proj.lcent,avg_clee,label='%s,%0.2f'%(simdes,avg_slopes_ee[0]))
            axb.plot(proj.lcent,avg_clbb,label='%s,%0.2f'%(simdes,avg_slopes_bb[0]))
            axt.plot(proj.lcent,avg_cltt,label='%s,%0.2f'%(simdes,avg_slopes_tt[0]))
            axv.plot(proj.lcent,avg_v,label='%s,%0.2f'%(simdes,avg_slopes_v[0]))
            axh.plot(proj.lcent,avg_h,label='%s,%0.2f'%(simdes,avg_slopes_h[0]))
            axd.plot(proj.lcent,avg_d,label='%s,%0.2f'%(simdes,avg_slopes_d[0]))
            axte.plot(proj.lcent,avg_clte,label='%s,%0.2f'%(simdes,avg_slopes_te[0]))
            axte.set_yscale('symlog',linthreshy=1e-5)
            axte.set_xscale('log')
            axte.set_ylim(-1e4,1e4)
            axtb.plot(proj.lcent,avg_cltb,label='%s,%0.2f'%(simdes,avg_slopes_tb[0]))
            axtb.set_yscale('symlog',linthreshy=1e-5)
            axtb.set_xscale('log')
            axtb.set_ylim(-1e4,1e4)
            axeb.plot(proj.lcent,avg_cleb,label='%s,%0.2f'%(simdes,avg_slopes_eb[0]))
            axeb.set_yscale('symlog',linthreshy=1e-5)
            axeb.set_xscale('log')
            axeb.set_ylim(-1e4,1e4)
            axreb.plot(proj.lcent,avg_reb,label='%s,%0.2f'%(simdes,avg_reb[15]))
            axreb.set_yscale('symlog',linthreshy=1e-5)
            axreb.set_xscale('log')
            axreb.set_ylim(-1e4,1e4)
        if 0:
        ##------
            reb[i]=avg_reb[15] #choose whatever makes sense for reb
            ee_slopes[i]=avg_slopes_ee[0]
            bb_slopes[i]=avg_slopes_bb[0]
            eb_slopes[i]=avg_slopes_eb[0]
            te_slopes[i]=avg_slopes_te[0]
            tb_slopes[i]=avg_slopes_tb[0]
            v_slopes[i]=avg_slopes_v[0]
            h_slopes[i]=avg_slopes_h[0]
            d_slopes[i]=avg_slopes_d[0]
            eb_ratios[i]=mean_quantity
            te_ratios[i]=mean_quantity2
            tb_ratios[i]=mean_quantity3
            ee_slopes=np.array(ee_slopes)
            bb_slopes=np.array(bb_slopes)
            eb_slopes=np.array(eb_slopes)
            te_slopes=np.array(te_slopes)
            tb_slopes=np.array(tb_slopes)
            v_slopes=np.array(v_slopes)
            h_slopes=np.array(h_slopes)
            d_slopes=np.array(d_slopes)
            eb_ratios=np.array(eb_ratios)
            te_ratios=np.array(te_ratios)
            tb_ratios=np.array(tb_ratios)
            if verbose:
                print('ee_slopes =')
                print(ee_slopes)
        ##------
#      ax1.plot(proj.lcent,avg_clee,label='ClEE',c='r')
#      ax1.plot(proj.lcent,avg_clbb,label='ClBB',c='b')
#      ax1.plot(proj.lcent,avg_clte,label='ClTT',c='y')
#      ax1.plot(proj.lcent,avg_cltb,label='ClTB',c='c')
#      dt.axbonk(ax1,xlabel='k',ylabel='Cl Power Spectra',xscale='log',yscale='log')
      #ax1.set_ylim((10**(-6)),300)
#      fig1.savefig("%s/%s_cl_average_all.pdf"%(plot_dir,longprefix))
#      dt.axbonk(ax2,xlabel='k',ylabel='ratio',xscale='log',yscale='log')
#      dt.axbonk(ax2,xlabel='k',ylabel='Cl Power Ratios',xscale='log')
#      ax2.set_ylim(0,10)
#      fig2.savefig("%s/%s_ratio_average_all.pdf"%(plot_dir,longprefix))
#      fig2.clf()
    ##------12 panel savefig
    axratio.legend(loc=0)
    axratio.set_title('ClEE Power Spectra')
    dt.axbonk(axratio,xlabel='k',ylabel='ClEE/ClBB',xscale='log',yscale='log')
    figratio.savefig("%s/ClEEClBB_%s_12panel_average.pdf"%(gen_dir,axes))
#
    axe.legend(loc=0)
    axe.set_title('ClEE Power Spectra')
    dt.axbonk(axe,xlabel='k',ylabel='ClEE Power',xscale='log',yscale='log')
    fige.savefig("%s/ClEE_%s_12panel_average.pdf"%(gen_dir,axes))
    axb.legend(loc=0)
    axb.set_title('ClBB Power Spectra')
    dt.axbonk(axb,xlabel='k',ylabel='ClBB Power',xscale='log',yscale='log')
    figb.savefig("%s/ClBB_%s_12panel_average.pdf"%(gen_dir,axes))
    axt.legend(loc=0)
    axt.set_title('ClTT Power Spectra')
    dt.axbonk(axt,xlabel='k',ylabel='ClTT Power',xscale='log',yscale='log')
    figt.savefig("%s/ClTT_%s_12panel_average.pdf"%(gen_dir,axes))
    axv.legend(loc=0)
    axv.set_title('Velocity Power Spectra')
    dt.axbonk(axv,xlabel='k',ylabel='Velocity Power',xscale='log',yscale='log')
    figv.savefig("%s/V_%s_12panel_average.pdf"%(gen_dir,axes))
    axh.legend(loc=0)
    axh.set_title('H-Field Power Spectra')
    dt.axbonk(axh,xlabel='k',ylabel='H-field Power',xscale='log',yscale='log')
    figh.savefig("%s/H_%s_12panel_average.pdf"%(gen_dir,axes))
    axd.legend(loc=0)
    axd.set_title('Density Power Spectra')
    dt.axbonk(axd,xlabel='k',ylabel='Density Power',xscale='log',yscale='log')
    figd.savefig("%s/D_%s_12panel_average.pdf"%(gen_dir,axes))
    axte.set_title('ClTE Power Spectra')
    figte.savefig("%s/ClTE_%s_12panel_average.pdf"%(gen_dir,axes))
    axtb.set_title('ClTB Power Spectra')
    figtb.savefig("%s/ClTB_%s_12panel_average.pdf"%(gen_dir,axes))
    axeb.set_title('ClEB Power Spectra')
    figeb.savefig("%s/ClEB_%s_12panel_average.pdf"%(gen_dir,axes))
    axreb.legend(loc=0)
    axreb.set_title('rEB Spectra')
    dt.axbonk(axreb,xlabel='k',ylabel='Correlation Ratio rEB')
    figreb.savefig("%s/rEB_%s_12panel_average.pdf"%(gen_dir,axes))
    ##------
    Fptr = h5py.File("slopes_ratios_dict_%s.h5"%axes,'w')
    Fptr.create_dataset("ee_slopes", data=ee_slopes)
    Fptr.create_dataset("bb_slopes", data=bb_slopes)
    Fptr.create_dataset("eb_slopes", data=eb_slopes)
    Fptr.create_dataset("te_slopes", data=te_slopes)
    Fptr.create_dataset("tb_slopes", data=tb_slopes)
    Fptr.create_dataset("eb_ratios", data=eb_ratios)
    Fptr.create_dataset("te_ratios", data=te_ratios)
    Fptr.create_dataset("tb_ratios", data=tb_ratios)
    Fptr.create_dataset("v_slopes", data=v_slopes)
    Fptr.create_dataset("h_slopes", data=h_slopes)
    Fptr.create_dataset("d_slopes", data=d_slopes)
    Fptr.close()


