
from GL import *
from cycler import cycler
import queb3
from queb3 import powerlaw_fit as plfit
import davetools as dt
import h5py
from matplotlib.pyplot import cm
import sim_colors
reload(sim_colors)
reload(queb3)
verbose=False
from collections import defaultdict
all_slopes=defaultdict(list)

class bin_tool():
    def __init__(self,Nx,Ny):
        self.N = np.array([Nx,Ny])

    def compute_bins_dx1(self):
        #We arrange things so that Delta=Delta x = 1.
        #Thus Delta L = 2pi/N
        #This is the best thing to use.
        self.xsize = self.N[0]

        self.size2d = np.array([self.xsize]*2)
        self.Delta = self.size2d/self.N
        self.Deltal = 2*np.pi/(self.N*self.Delta) #cmbtools.Delta2l(self.Delta,self.N) #2 pi/(N Delta)

        self.lmax = self.Deltal[0]*self.N[0]/2
        self.lbins = np.arange(0,self.N[0]//2) *self.Deltal[0]
        self.lcent = self.lbins[:-1] + np.diff(self.lbins)/2.
        return self.lcent

a=bin_tool(512,512)
lcent=a.compute_bins_dx1()



class slopes():
    def __init__(self,name):
        self.name=name
        self.spectra={}
        self.slopes={}
        self.amps={}
        self.res={}
        self.pack = None
    def fit(self,lcent,fitrange):
        for field in self.spectra:
            slope,amp,res=plfit(lcent,self.spectra[field],fitrange)
            self.amps[field]=amp
            self.slopes[field]=slope
            self.res[field]=res

if 'spectra_dict' not in dir():
    spectra_dict={}
    shortprefix='time'
    for axes in ['x','y','z']:
        spectra_dict[axes]={}
        for i, sim in enumerate(sim_colors.simlist):
            spectra_dict[axes][sim]=slopes(sim) 

            sim_dir = "/archive2/kas14d/512reruns/frbs/%s"%sim
            longprefix='%s_%s_%s'%(sim,shortprefix,axes)
            spectra_fname = "avg_spectra_%s_%s.h5"%(sim,axes)
            pack = queb3.simulation_package( directory=sim_dir,frames=sim_colors.frames[sim],prefix=longprefix)
            proj=pack.read_queb(frame=sim_colors.frames[sim][0],ax=axes,bin_style='dx1')
            fitrange=proj.determine_fit_range()  #something better
            if os.path.exists(spectra_fname):
                try:
                    fptr=h5py.File(spectra_fname,'r')
                    for field in fptr:
                        spectra_dict[axes][sim].spectra[field]=fptr[field][()]
                except:
                    raise
                finally:
                    fptr.close()
            spectra_dict[axes][sim].fit(proj.lcent,fitrange)

#
# Spectra plots
#

if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb']):
        for sim in sim_colors.simlist:
            axlist[nf].plot(proj.lcent, spectra_dict['x'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb']):
        for sim in sim_colors.simlist:
            axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            axlist[nf+3].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
    for a in axlist[3:]:
        a.set_xlabel(r'$k/k_max$')
    for a in [axlist[0],axlist[3]]:
        a.set_ylabel(r'$P(k)$')

    plotdir='/home/dcollins4096/PigPen'
    fig.savefig('%s/TT_EE_BB.pdf'%plotdir)

if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb']):
        for sim in sim_colors.simlist:
            qu = field[-2:].upper()
            label = r'$%s, \parallel$'%qu
            axlist[nf  ].text(2.5,-5,label)
            axlist[nf  ].scatter(sim_colors.Ms[sim], spectra_dict['x'][sim].slopes[field], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            label = r'$%s, \perp$'%qu
            axlist[nf+3].text(2.5,-5,label)
            axlist[nf+3].scatter(sim_colors.Ms[sim], spectra_dict['y'][sim].slopes[field], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            axlist[nf+3].scatter(sim_colors.Ms[sim], spectra_dict['z'][sim].slopes[field], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
    for a in axlist:
        dt.axbonk(a,xscale='linear',yscale='linear',xlabel=None,ylabel=None)
    for a in axlist[3:]:
        a.set_xlabel(r'$M_s$')
    for a in [axlist[0],axlist[3]]:
        a.set_ylabel(r'$\alpha_{XX}$')

    plotdir='/home/dcollins4096/PigPen'
    fig.savefig('%s/Slopes_Ms.pdf'%plotdir)

if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_cltt','avg_clee','avg_clbb']):
        for sim in sim_colors.simlist:
            qu = field[-2:].upper()
            label = r'$%s, \parallel$'%qu
            axlist[nf  ].text(1.5,-5,label)
            axlist[nf  ].scatter(sim_colors.Ma[sim], spectra_dict['x'][sim].slopes[field], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            label = r'$%s, \perp$'%qu
            axlist[nf+3].text(1.5,-5,label)
            axlist[nf+3].scatter(sim_colors.Ma[sim], spectra_dict['y'][sim].slopes[field], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
            axlist[nf+3].scatter(sim_colors.Ma[sim], spectra_dict['z'][sim].slopes[field], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
    for a in axlist:
        dt.axbonk(a,xscale='linear',yscale='linear',xlabel=None,ylabel=None)
    for a in axlist[3:]:
        a.set_xlabel(r'$M_a$')
    for a in [axlist[0],axlist[3]]:
        a.set_ylabel(r'$\alpha_{XX}$')

    plotdir='/home/dcollins4096/PigPen'
    fig.savefig('%s/Slopes_Ma.pdf'%plotdir)

if 0:
    plt.close('all')
    fig,ax = plt.subplots(2,3, sharex=True,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    axlist=ax.flatten()

    for nf,field in enumerate(['avg_reb', 'avg_rtb', 'avg_rte']):
        for sim in sim_colors.simlist:
            qu = field[-2:].upper()
            label = r'$r_{%s} \parallel$'%qu
            axlist[nf  ].text(0.1, -0.8, label)
            label = r'$r_{%s} \perp$'%qu
            axlist[nf+3].text(0.1, -0.8, label)
            axlist[nf  ].plot(proj.lcent, spectra_dict['x'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            axlist[nf+3].plot(proj.lcent, spectra_dict['y'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
            #axlist[nf+6].plot(proj.lcent, spectra_dict['z'][sim].spectra[field], c=sim_colors.color[sim], linestyle=sim_colors.linestyle[sim])
    for a in axlist:
        dt.axbonk(a,xscale='log',yscale='log',xlabel=None,ylabel=None)
        a.set_yscale('symlog',linthresh=0.9)
        #a.set_yscale('linear')
        a.set_ylim([-1,1])
    for a in axlist[3:]:
        a.set_xlabel(r'$k/k_max$')
    for a in [axlist[0],axlist[3]]:
        a.set_ylabel(r'$r_{XY}$')
        a.set_yticks([-1,-0.1,0.1,1])

    plotdir='/home/dcollins4096/PigPen'
    fig.savefig('%s/r_XY.pdf'%plotdir)

if 1:
    plt.close('all')
    fig,ax = plt.subplots(3,2)#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for sim in sim_colors.simlist:
        ax[0][0].scatter(sim_colors.Ma[sim], spectra_dict['y'][sim].slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][0].scatter(sim_colors.Ma[sim], spectra_dict['y'][sim].slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][0].scatter(sim_colors.Ma[sim], spectra_dict['y'][sim].slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        print(sim_colors.Ma[sim], spectra_dict['y'][sim].slopes['avg_d'])
        
        ax[0][0].set_ylabel(r'$\alpha_v$')
        ax[1][0].set_ylabel(r'$\alpha_\rho$')
        ax[2][0].set_ylabel(r'$\alpha_H$')
        ax[2][0].set_xlabel(r'$M_{\rm{A}}$')
        ax[2][1].set_xlabel(r'$M_{\rm{s}}$')
        ax[0][1].scatter(sim_colors.Ms[sim], spectra_dict['y'][sim].slopes['avg_v'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[1][1].scatter(sim_colors.Ms[sim], spectra_dict['y'][sim].slopes['avg_d'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        ax[2][1].scatter(sim_colors.Ms[sim], spectra_dict['y'][sim].slopes['avg_h'], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        for a in [ax[0][1],ax[1][1],ax[2][1]]:
            a.yaxis.tick_right()
    fig.savefig('%s/prim_spectra.pdf'%plotdir)

if 1:
    plt.close('all')
    fig,ax = plt.subplots(1,1)#,sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    for sim in sim_colors.simlist:
        ax.scatter(sim_colors.Ma[sim], sim_colors.Ms[sim], c=sim_colors.color[sim], marker=sim_colors.marker[sim])
        dt.axbonk(ax,xlabel=r'$M_{\rm{A}}$',ylabel=r'$M_{\rm{S}}$',xlim=[0.45,2.1],ylim=[0.45,3.2])

    fig.savefig('%s/point_legend.pdf'%plotdir)
            

            

#dict_keys(['avg_clbb', 'avg_cleb', 'avg_clee', 'avg_cltb', 'avg_clte', 'avg_cltt', 'avg_d', 'avg_h', 'avg_reb', 'avg_rtb', 'avg_rte', 'avg_v'])





