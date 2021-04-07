


from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(queb3)
reload(sim_colors)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)


shortprefix="time"
def make_spectra_files():
    for axes in ['x','y','z']:
        simlist=sim_colors.simlist
        for i,sim in enumerate(simlist):
            frames=sim_colors.framelist[i]
            simdes=sim
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
            avg_rte=0
            avg_rtb=0
            avg_reb=0
            projections=[]
            longprefix='%s_%s_%s'%(simdes,shortprefix,axes)
            pack = queb3.simulation_package( directory=sim_dir,frames=frames,prefix=longprefix)
            nplots=0
            proj=pack.read_queb(frame=frames[0],ax=axes,bin_style='dx1')
            fitrange=proj.determine_fit_range()  #something better
            if not os.path.exists(spectra_fname) or True:

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
                    avg_rte += proj['ClTE']/((proj['ClTT']*proj['ClEE'])**0.5)
                    avg_rtb += proj['ClTB']/((proj['ClTT']*proj['ClBB'])**0.5)
                    avg_reb += proj['ClEB']/((proj['ClEE']*proj['ClBB'])**0.5)
                    if np.isnan(avg_rtb).any():
                        pdb.set_trace()
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
                avg_rte  /= nplots
                avg_rtb  /= nplots
                avg_reb  /= nplots
                print(avg_rtb)
                if np.isnan(avg_rtb).any():
                    pdb.set_trace()
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
                    fptr.create_dataset("avg_rte" ,data = avg_rte )
                    fptr.create_dataset("avg_rtb" ,data = avg_rtb )
                    fptr.create_dataset("avg_reb" ,data = avg_reb )

                finally:
                    fptr.close()

make_spectra_files()
