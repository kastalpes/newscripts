from GL import *
import astropy.io.fits as pyfits
import matplotlib.colors as colors
import sim_colors
import davetools as dt
reload(sim_colors)
plt.close('all')
#ms_ma; ms = line, ma = color
#simlist=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
#color=['r','g','b','r','g','b','r','g','b','r','g','b']
#linestyle=['-','-','-','-.','-.','-.','--','--','--',':',':',':']
#simlist=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
#color=['r','g','b','r','g','b','r','g','b','r','g','b']
#linestyle=['-','-','-','-.','-.','-.','--','--','--',':',':',':']
basedir="/archive2/kas14d/512reruns/"
outdir="plots"
hist_density=False
for axis in 'xy': #axis='y'
    density_array={}
    b_array={}
    e_array={}

    for sim in sim_colors.plot_order:
        directory = "%s/frbs/%s"%(basedir,sim)
        #pick the last density frb
        density_frb_list = glob.glob("%s/*density_%s.fits"%(directory,axis))
        last_density = sorted(density_frb_list)[-1]
        density_array[sim]= pyfits.open(last_density)[0].data 

        b_frb_list = glob.glob("%s/*B%s.fits"%(directory,axis))
        last_b = sorted(b_frb_list)[-1]
        b_array[sim] =  pyfits.open(last_b)[0].data 

        e_frb_list = glob.glob("%s/*E%s.fits"%(directory,axis))
        last_e = sorted(e_frb_list)[-1]
        e_array[sim] =  pyfits.open(last_e)[0].data 



    fig, axes = plt.subplots(3,4)
    ax_list=axes.flatten()

    den_max = 1.1# max([b.max() for b in density_array])
    den_min = 0.9 # min([b.min() for b in density_array])
    norm = colors.Normalize(vmin=den_min,vmax=den_max)
    norm = colors.SymLogNorm(linthresh = 0.1, vmin=-1.1,vmax=1.1)
    fig3,ax3 = plt.subplots(1,1)
    for i,sim in enumerate(sim_colors.plot_order):

        this_array= density_array[sim] - 1
        plot=ax_list[i].imshow( this_array,origin='lower',interpolation='nearest',norm=norm)
        ax3.hist( density_array[sim].flatten(),histtype='step',label='%s'%sim,
                 color=sim_colors.color[sim],
                 linestyle=sim_colors.linestyle[sim],density=hist_density)
        ax_list[i].set_xticks([])
        ax_list[i].set_yticks([])
        ax_list[i].text(10,10,sim,c=[1.]*3)
    ax3.legend(loc=0)
    dt.axbonk(ax3,xlabel=r'$T-mode$',ylabel=r'$N$')
    fig3.savefig("%s/Thist_%s.pdf"%(outdir,axis))
    plt.close(fig3)
    cb=fig.colorbar(plot, ax=axes[:,3])
    fig.subplots_adjust(left=0,bottom=0.05,top=0.95,wspace=0, hspace=0)
    fig.savefig("%s/multiplot_density_%s.pdf"%(outdir,axis))
    plt.close(fig)

    fig, axes = plt.subplots(3,4)
    ax_list=axes.flatten()

    fig3,ax3 = plt.subplots(1,1)

    norm = colors.SymLogNorm(linthresh = 0.01, vmin=-1,vmax=1)
    for i,sim in enumerate(sim_colors.plot_order):
        this_array =  b_array[sim]
        plot=ax_list[i].imshow( this_array,origin='lower',interpolation='nearest',norm=norm)
        ax3.hist( b_array[sim].flatten(),histtype='step',label='%s'%sim,
                 color=sim_colors.color[sim],
                 linestyle=sim_colors.linestyle[sim],density=hist_density)
        ax_list[i].set_xticks([])
        ax_list[i].set_yticks([])
        ax_list[i].text(10,10,sim,c=[1.]*3)

    ax3.legend(loc=0)
    dt.axbonk(ax3,ylabel=r'$N$',xlabel=r'$B-mode$')
    #fig3.subplots_adjust(left=0,bottom=0.05,top=0.95,wspace=0, hspace=0)
    fig3.savefig("%s/Bhist_%s.pdf"%(outdir,axis))
    plt.close(fig3)
    cb=fig.colorbar(plot, ax=axes[:,3])
    #fig.subplots_adjust(wspace=0, hspace=0)
    fig.subplots_adjust(left=0,bottom=0.05,top=0.95,wspace=0, hspace=0)
    outname="%s/multiplot_B_%s.pdf"%(outdir,axis)
    fig.savefig(outname)
    print(outname)
    plt.close(fig)
    fig, axes = plt.subplots(3,4)
    ax_list=axes.flatten()

    b_max =  0.2 #max([b.max() for b in b_array])*2
    b_min = -0.2 #min([b.min() for b in b_array])*2
#norm = colors.SymLogNorm(linthresh=.1,vmin=b_min,vmax=b_max)
#norm = colors.Normalize(vmin=b_min,vmax=b_max)
    minmax=2.5
    norm = colors.SymLogNorm(linthresh = 0.02, vmin=-minmax,vmax=minmax)
    rm=dt.rainbow_map(12)
    fig3,ax3 = plt.subplots(1,1)

    for i,sim in enumerate(sim_colors.plot_order):
        ax_list[i].clear()
        this_array = e_array[sim] - e_array[sim].mean()
        ax3.hist( e_array[sim].flatten(),histtype='step',label='%s'%sim,
                 color=sim_colors.color[sim],
                 linestyle=sim_colors.linestyle[sim],density=hist_density)
        ax_list[i].imshow( this_array,origin='lower',interpolation='nearest')#,norm=norm)
        ax_list[i].set_xticks([])
        ax_list[i].set_yticks([])
        ax_list[i].text(10,10,sim,c=[1.]*3)

    ax3.legend(loc=0)
    dt.axbonk(ax3,ylabel=r'$N$',xlabel=r'$E-mode$')
    fig3.savefig("%s/Ehist_%s.pdf"%(outdir,axis))
    plt.close(fig3)
    cb=fig.colorbar(plot, use_gridspec=True, ax=axes[:,3])
    #fig.subplots_adjust(wspace=0, hspace=0)
    fig.subplots_adjust(left=0,bottom=0.05,top=0.95,wspace=0, hspace=0)
#fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    fig.savefig("%s/multiplot_E_%s.pdf"%(outdir,axis))
    plt.close(fig)

