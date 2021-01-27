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
frames=range(10,15)
plot_dir =  "/home/dcollins4096/PigPen"
prefix='p49d_multi'
do_ratio = False

if 1:
    simcolor={}
    simcolor['1_1']='r'
    simcolor['1_2']='r--'
    simcolor['1_half']='r:'
    simcolor['2_1']='g'
    simcolor['2_2']='g--'
    simcolor['2_half']='g:'
    simcolor['3_1']='b'
    simcolor['3_2']='b--'
    simcolor['3_half']="b:"
    simcolor['half_1']="c"
    simcolor['half_2']="c--"
    simcolor['half_half']="c:"
    simlist = ['1_1','1_2','1_half','2_1','2_2','2_half','3_1','3_2','3_half','half_1','half_2','half_half']
    simlist = ['1_1'] #kludge
    simulation_base = "/scratch1/dcollins/Paper49d/512reruns/"

    simlist = ['B02/512']
    simulation_base = "/scratch1/dcollins/Paper08/"

if 1:
    plt.close('all')
    fig_cl,axes_cl=plt.subplots(2, 3,figsize=(12,12))
    fig_po,axes_po=plt.subplots(2, 3,figsize=(12,12))

    ax_ee= axes_cl[0][0]
    ax_bb= axes_cl[0][1]
    ax_tt= axes_cl[0][2]
    ax_te= axes_cl[1][0]
    ax_tb= axes_cl[1][1]
    ax_eb= axes_cl[1][2]

    ax_v = axes_po[0][0]
    ax_d = axes_po[0][1]
    ax_h = axes_po[0][2]
    ax_r = axes_po[1][0]
    ax_a = axes_po[1][1]
    ax_hh= axes_po[1][2]

    positive_definite = ['aspec','vspec','dspec','hspec', 'ClEE','ClBB', 'ClHH']
    fields = ['aspec','vspec','dspec','hspec', 'ClEE','ClBB','ClTE','ClTT', 'ClHH', 'ClTB','ClEB']
    fields = ['vspec'] #kludge
    axes =  [ax_a,   ax_v,    ax_d,   ax_h,   ax_ee, ax_bb, ax_te, ax_tt, ax_hh, ax_tb, ax_eb]
    axdict = dict( zip( fields, axes))
                   
                     


if 1:
    frame_ranges={}
    frame_ranges['2_half'] = range(65,85)
    frame_ranges['3_half'] = range(72,91)
    frame_ranges['3_1'] = range(56,74)
    frame_ranges['3_2'] = range(30,40)
    frame_ranges['B02/512'] = [10]
    projections={}
    spectra_monster = {}
    slope_monster = {}
    spectra_3d = ['aspec','vspec','dspec','hspec']
    spectra_list = spectra_3d + ['ClEE','ClBB','ClTE','ClTT', 'ClHH', 'ClTB','ClEB']
    spectra_list = fields #kludge
    ratio_avg = {}
    avg_ratio = {}
    k_avg = {}
    l_avg = {}
    for  field in spectra_list:
        spectra_monster[field]={}

    for sim in simlist:
        slope_monster[sim]=queb3.slope_package()

    for sim in simlist:
        sim_dir = simulation_base + sim
        product_directory = simulation_base + sim
        prefix = "P49d_"+sim

        pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=prefix, product_directory=product_directory)

        #setup
        projections[sim]=[]
        for field in spectra_monster:
            spectra_monster[field][sim]  = 0
        ratio_avg[sim]=0
        k_avg[sim]  = 0
        l_avg[sim]  = 0
        avg_ratio[sim]=0

        frames=frame_ranges.get(sim,range(11,30))
        for frame in frames:
            proj=pack.read_queb(frame=frame,ax='y',bin_style='dx1')
            proj.read_spectra(frame=frame)
            projections[sim].append(proj)
            projections[sim][-1].compute_harmonic_products()
            for field in spectra_list:
                if field in spectra_3d:
                    spectra_monster[field][sim] += proj[field][1].real
                    k_avg[sim] += proj[field][0]
                else:
                    spectra_monster[field][sim] += proj[field]
                    l_avg[sim] += proj.lcent
            print("WWW sim %s frame %d ksize %d"%(sim,frame, k_avg[sim].size))
            if do_ratio:
                ratio_avg[sim]  += proj['ClEE']/proj['ClBB']


        if 1:
            for field in spectra_monster:
                spectra_monster[field][sim] /= len(frames)
            if do_ratio:
                avg_ratio[sim] = spectra_monster['ClEE'][sim]/spectra_monster['ClBB'][sim]
                ratio_avg[sim]/=len(frames)
            k_avg[sim] /= len(frames)
            for field in spectra_list:
                if field not in positive_definite:
                    continue
                fit_range=projections[sim][0].determine_fit_range()
                fit_range = [2e-2,7e-2]
                if field in spectra_3d:
                    the_x = k_avg[sim]/k_avg[sim].max()
                else:
                    the_x = l_avg[sim]/l_avg[sim].max()
                slope_monster[sim].ingest( field, *queb3.powerlaw_fit( the_x, spectra_monster[field][sim], fit_range))


if 0: #kludge
    signed_extents = dt.extents()
    for field in ['ClTE','ClEB','ClTB']:
        for sim in simlist:
            signed_extents( spectra_monster[field][sim])

if 1:
    def compensate(k, p):
        comp = p*k**(2.5)
        comp *= 1/comp[1]
        return k, comp
    for field in spectra_monster:
        for sim in simlist:
            if field in spectra_3d:
                the_y = spectra_monster[field][sim]
                the_x = k_avg[sim]/k_avg[sim].max()
            else:
                the_y = spectra_monster[field][sim]
                the_x = l_avg[sim]/l_avg[sim].max()
            axdict[field].plot( the_x, the_y,simcolor[sim],label=sim)#,c='g')
            if field in positive_definite:
                slope_monster[sim].plot(axdict[field],field)
            axdict[field].set_title("%s"%(field))
            ax_r.plot(l_avg[sim], ratio_avg[sim],simcolor[sim],label=sim)

        #ax_bb.plot( *compensate(proj.lcent,avg_clbb[sim]),simcolor[sim],label=sim)#,c='b')
        #ax_r.plot(proj.lcent, avg_ratio[sim],simcolor[sim],label=sim)#,c='r',label=sim)
        #ax_ratio.plot(proj.lcent, ratio_avg[sim],c='g')
if 1:
    #ax_ee.legend(loc=0)
    #ax_bb.legend(loc=0)
    dt.axbonk(ax_ee,xlabel='k',ylabel='<ClEE>',xscale='log',yscale='log')
    dt.axbonk(ax_bb,xlabel='k',ylabel='<ClBB>',xscale='log',yscale='log')
    dt.axbonk(ax_tt,xlabel='k',ylabel='<ClTT>',xscale='log',yscale='log')

    dt.axbonk(ax_tb,xlabel='k',ylabel='<ClTB>',xscale='log',yscale='log')
    dt.axbonk(ax_te,xlabel='k',ylabel='<ClTE>',xscale='log',yscale='log')
    dt.axbonk(ax_eb,xlabel='k',ylabel='<ClEB>',xscale='log',yscale='log')

    dt.axbonk(ax_v,xlabel='k',ylabel=r'$P_v$',xscale='log',yscale='log')
    dt.axbonk(ax_d,xlabel='k',ylabel=r'$P_\rho$',xscale='log',yscale='log')
    dt.axbonk(ax_h,xlabel='k',ylabel=r'$P_H$',xscale='log',yscale='log')

    dt.axbonk(ax_a,xlabel='k',ylabel=r'$P_A$',xscale='log',yscale='log')
    dt.axbonk(ax_h,xlabel='k',ylabel=r'$P_H$',xscale='log',yscale='log')
    dt.axbonk(ax_hh,xlabel='k',ylabel=r'$C_\ell^{HH}$',xscale='log',yscale='log')

    dt.axbonk(ax_r,xlabel='k',ylabel=r'$C_\ell^{EE}/C_\ell^{BB}$',xscale='log',yscale='log')

    if 0: #kludge   
        for this_ax in [ax_tb,ax_te,ax_eb]:
            this_ax.set_yscale('symlog',linthreshy=1e-3)
            this_ax.set_ylim( signed_extents.minmax )



    fig_cl.savefig(os.environ['HOME']+"/PigPen/CLxx_average.png")
    fig_po.savefig(os.environ['HOME']+"/PigPen/SPECxx_average.png")
