"""
Improved tools for generating QUEB maps, fitting spectra, and plotting stuff.
Two primary  objects: simulation_package and queb_snapshot.
simulation_package holds the location of the simulation and meta data.
queb_snapshot is a container for a single projection (frame and axis).  
Simulation_package has several useful methods:
    EBall computes FRBs of Q and U for all directions and frames
    make_frbs makes frbs.  Can be expensive
    make_spectra makes 3d spectra of primitive quantities. Can be expensive.
    read_queb reads QUEBTH fields from disk
    image_fields makes images of the FRBs
    plot_eb plots EE and BB
    plot_many_spectra plots all the spectra.


queb_snapshot has several useful methods:
    compute_harmonic_products takes Q&U and produces E&B, transforms, and spectra
    determine_fit_range is a stub that will get more physics soon.
    fit_eb_slopes fits EE and BB spectra to powerlaw
Also in this update is "spectra_tools" which creates 3d spectra.

Usage can be found in make_all_queb
"""
from GL import *
import yt
import p49_QU2EB
import p49_fields
import cmbtools
from scipy.optimize import leastsq
import davetools as dt
reload(dt)
import spectra_tools as st
reload(st)
frbname="frbs"
reload(p49_QU2EB)

def powerlaw_fit(x,y,fitrange):  
    mb=nar([0,0])
    ellmask = (fitrange[0] < x)*(x < fitrange[1])
    ellmask = np.logical_and( ellmask, y>0 )
    logx = np.log10(x[ellmask])
    logy = np.log10(y[ellmask])
    fig,ax=plt.subplots(1,1)
    ax.plot(logx,logy)
    res = np.polyfit(logx,logy,1)
    slope = res[0]
    amp = pow(10,res[1])
    ax.plot(logx, slope*logx+res[1],c='k')
    return slope,amp, res

def read_fits(fitname):
    """Read an array from a file name.
    Returns None if the file does not exist."""
    d=None
    if os.path.exists(fitname):
        d=np.ascontiguousarray(pyfits.open(fitname)[0].data,dtype=np.double)
    return d

def linear_chi2(par,x,y) :
    """for use with some fitting techniques"""
    m = par[0]
    y0 = par[1]
    xcent = (x[0]+x[-1])/2
    ymodel = m*(x-xcent)+y0
    return(y - ymodel)

class slope_package():
    """A container for hanging on to fit products."""   
    def __init__(self):
        self.slope={}
        self.amp={}
        self.res={}
        self.range=[]
    def ingest(self,which,slope=None,amp=None,res=None):
        self.slope[which]=slope
        self.amp[which]=amp
        self.res[which]=res
    def plot(self,ax,which,norm=False,**kwargs):
        """Plots the line we found onto *ax*.
        Also update the label of the line to reflect the value
        of the slope"""
        m=self.slope[which]
        a=self.amp[which]
        ellfit=1    
        if norm:
            ellfit = np.sqrt(self.range[0]*self.range[1])
        y0 = a*(self.range[0]/ellfit)**m
        y1 = a*(self.range[1]/ellfit)**m
        label=kwargs.get('label','')
        label += r' $%0.1f$'%m
        kwargs['label']=label
        ax.plot( self.range,[y0,y1],**kwargs)



class queb_snapshot():
    """Container for projections: QUEB, T (density), H(magnetic field).
    Contains fields, their transforms, and spectra
    """
    def __init__(self,Q,U,T=None,H=None,E=None,B=None,axis='x',frame=-1, simulation=None, bin_style='dx1'):
        self.Q=Q
        self.U=U
        self.T=T
        self.H=H
        self.E=E
        self.B=B
        self.axis=axis
        self.frame=frame
        self.simulation=simulation
        self.bin_style=bin_style
    def __getitem__(self,item):
        """square bracket access method"""
        #this is a kludge
        return self.__dict__[item]
    def write(self):
        """saves E, B, and the spectra to the simulation location"""
        xd='DD'
        frb_dir = "%s/%s"%(self.simulation.directory,p49_QU2EB.frbname)
        product_dir = "%s/DD%04d.products"%(self.simulation.directory,self.frame)
        Ef= "%s/%s%04d_E%s.fits"%(frb_dir,xd,self.frame,self.axis)
        Bf= "%s/%s%04d_B%s.fits"%(frb_dir,xd,self.frame,self.axis)
        Clf= "%s/%s%04d_Cl%s.fits"%(frb_dir,xd,self.frame,self.axis)
        hdu = pyfits.PrimaryHDU(self.E)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(Ef,overwrite=True)
        hdu = pyfits.PrimaryHDU(self.B)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(Bf,overwrite=True)
        np.savetxt(Clf, list(zip(self.lbins,self.ClEE,self.ClBB, self.ClTE, self.ClEB)))

    def compute_bins_dx1(self):
        #We arrange things so that Delta=Delta x = 1.
        #Thus Delta L = 2pi/N
        #This is the best thing to use.
        self.N = np.array(self.Q.shape,dtype = np.int32)
        self.xsize = self.N[0]

        self.size2d = np.array([self.xsize]*2)
        self.Delta = self.size2d/self.N
        self.Deltal = 2*np.pi/(self.N*self.Delta) #cmbtools.Delta2l(self.Delta,self.N) #2 pi/(N Delta)

        self.lmax = self.Deltal[0]*self.N[0]/2
        self.lbins = np.arange(0,self.N[0]//2) *self.Deltal[0]
        self.lcent = self.lbins[:-1] + np.diff(self.lbins)/2.

    def compute_bins_5deg(self):
        #Earlier computations assumed a 5degree sky.
        #not robust against aliasing.  Probably shouldn't use.
        self.xsize = 5 * np.pi / 180

        self.size2d = np.array([self.xsize,self.xsize])
        self.N = np.array(np.shape(self.Q),dtype = np.int32)
        self.Delta = self.size2d/self.N
        self.Deltal = cmbtools.Delta2l(self.Delta,self.N)

        self.lmax = self.Deltal[0]*self.N[0] 
        self.lbins = np.linspace(0,self.lmax,100)
        self.lcent = self.lbins[:-1] + np.diff(self.lbins)/2.

    def compute_bins_horse_around(self):
        #Change bin size to enjoy aliasing effects.
        self.N = np.array(np.shape(self.Q),dtype = np.int32)
        self.xsize = self.N[0]

        self.size2d = np.array([self.xsize,self.xsize])
        self.Delta = self.size2d/self.N
        self.Deltal = cmbtools.Delta2l(self.Delta,self.N)

        self.lmax = self.Deltal[0]*self.N[0]
        self.lbins = np.linspace(0,self.lmax,100)
        #self.lmax = self.Deltal[0]*self.N[0]/2
        #self.lbins = np.arange(0,self.N[0]//2) *self.Deltal[0]


        self.lcent = self.lbins[:-1] + np.diff(self.lbins)/2.


    def scatter_clee(self,fname):
        fig,ax=plt.subplots(1,1,figsize=(12,8))

        N = self.Q.shape

        #build the ell 
        rxa,rya = np.mgrid[ 0:N[0]:1,
                          0:N[0]//2+1:1]

        rxa[self.N[1]//2:,:] = rxa[:self.N[1]//2,:]- N[0]/2 #self.lmax #N[0]#*self.Deltal[0]/2
        rx=rxa*self.Deltal[0] 
        ry=rya*self.Deltal[0] 
            
        ell = np.sqrt(rx**2+ry**2)

        ####
        #here is code to plot r, and ensure it's right. 
        #fig2,ax2=plt.subplots(1,1,figsize=(8,8))
        #the transposei is to make 'x' be the horizontal axis.
        #ax2.imshow(ell.transpose(),origin='lower',interpolation='nearest')
        #fig.colorbar(ppp)
        #fig.savefig(fname)
        #fig2.savefig('/home/dcollins4096/PigPen/test.png')

        #compute the spectra.
        #Take out the normalization
        spectra = (np.abs(self.Eharm)**2).flatten()
        this_ClEE = cmbtools.harm2cl(self.Eharm,self.Deltal,self.lbins)
        this_ClEE *= (2*np.pi)**2*(self.Delta[0])**2/self.Deltal[0]**2
        
        #overplot bins, to see aliasing
        for l in self.lbins:
            ax.plot( [l,l], [1e-11,100],c=[0.5]*4)
        ax.plot( [np.pi]*2, [1e-11,100],c='k')
        ax.scatter( ell.flatten(), spectra,s=0.1)
        dt.axbonk(ax,xlabel='ell',ylabel='CEE')
        ax.set_xscale('symlog',linthreshx=self.Deltal[0])
        ax.set_yscale('symlog',linthreshy=1e-9)
        ax.set_ylim( 0, spectra.max())
        ax.set_xlim( ell.min(), ell.max()*1.1)
        ax.plot( self.lcent, this_ClEE,marker='*')

        fig.savefig(fname)
        print(fname)
        plt.close(fig)
    def compute_bins(self):
        if self.bin_style=='5deg':
            self.compute_bins_5deg()
        elif self.bin_style=='horse':
            self.compute_bins_horse_around()
        else:
            self.compute_bins_dx1()

    def compute_harmonic_products(self):    
        """
        This step may prove to be expensive.  If so, insert logic to cache 
        products to and from disk within this function.  
        """
        if not self.Q.flags['C_CONTIGUOUS']:
            self.Q = np.ascontiguousarray(Q)
        if not self.U.flags['C_CONTIGUOUS']:
            self.U = np.ascontiguousarray(U)
        if not self.T.flags['C_CONTIGUOUS']:
            self.T = np.ascontiguousarray(T)
        self.N = np.array(self.Q.shape,dtype = np.int32)

        self.compute_bins()

        self.Qharm = cmbtools.map2harm(self.Q,self.Delta)
        self.Uharm = cmbtools.map2harm(self.U,self.Delta)

        self.Eharm, self.Bharm = cmbtools.QU2EB(self.Qharm,self.Uharm,self.Deltal)
        self.E = cmbtools.harm2map(self.Eharm,self.Delta)
        self.B = cmbtools.harm2map(self.Bharm,self.Delta)
        if self.T is not None:
            self.Tharm = cmbtools.map2harm(self.T,self.Delta)
        if self.H is not None:
            self.Hharm = cmbtools.map2harm(self.H,self.Delta)
        self.ClEE = cmbtools.harm2cl(self.Eharm,self.Deltal,self.lbins)
        self.ClBB = cmbtools.harm2cl(self.Bharm,self.Deltal,self.lbins)
        if self.Hharm is not None:
            self.ClHH= cmbtools.harm2cl(self.Hharm,self.Deltal,self.lbins)
        self.ClTT = cmbtools.harm2cl(self.Tharm,self.Deltal,self.lbins)
        self.ClTE = cmbtools.harm2clcross_samegrid(self.Tharm,self.Eharm,self.Deltal,self.lbins)
        self.ClTB = cmbtools.harm2clcross_samegrid(self.Tharm,self.Bharm,self.Deltal,self.lbins)
        self.ClEB = cmbtools.harm2clcross_samegrid(self.Eharm,self.Bharm,self.Deltal,self.lbins)

    def fit_eb_slopes(self,fitrange=None,slopes=None):
        """Given the fit range, determine the slope.
        The *slopes* argument a container for the slopes and amplitudes.  Developmental.
        """
        if slopes is None:
            slopes=slope_package()
        if fitrange is None:    
            fitrange = self.determine_fit_range()
        slopes.range=fitrange
        self.compute_bins()
        if self.E is None:
            self.compute_harmonic_products()
        #the star turns the output of powerlaw_fit into three arguments for 
        #slopes.ingest
        slopes.ingest("EE", *powerlaw_fit(self.lcent,self.ClEE, fitrange))
        slopes.ingest("BB", *powerlaw_fit(self.lcent,self.ClBB, fitrange))
        slopes.ingest("TT", *powerlaw_fit(self.lcent,self.ClTT, fitrange))
        return slopes
    def fit_eb_slopes_2(self,fitrange=None,slopes=None):
        """Given the fit range, determine the slope.
        The *slopes* argument a container for the slopes and amplitudes.  Developmental.
        """
        if slopes is None:
            slopes=slope_package()
        if fitrange is None:    
            fitrange = self.determine_fit_range()
        slopes.range=fitrange
        ell = self.lcent
        ellmask = (fitrange[0] < ell)*(ell < fitrange[1])
        ClE = self['ClEE']; ClB = self['ClBB']
        ellmask_finite = np.logical_and(ellmask, ClE > 0)
        ellmask_finite = np.logical_and(ellmask_finite, ClB > 0)
        x =np.log10(ell[ellmask_finite])
        y = np.log10(ClE[ellmask_finite])
        mb = np.array([0,0])
        res = leastsq(linear_chi2, mb, args=(x,y) )
        Eslope = res[0][0]
        Eamp = pow(10,res[0][1])
        slopes.ingest("EE",slope=Eslope,amp=Eamp,res=res)
        return slopes

    def determine_fit_range(self):
        """hopefully this will be made to be a more accurate representation of the fit in the future.
        """
        fitrange=np.zeros(2)
        self.compute_bins()
        if self.bin_style == 'dx1':
            fitrange[0] = 4*self.Deltal[0]
            fitrange[1]=  8*self.Deltal[0]
            fitrange[0] = self.lcent[4]
            fitrange[1]=  self.lcent[10]
        else:
            fitlmin=self.lcent[4]
            fitlmax=self.lcent[10]
            fitrange[0]=fitlmin;fitrange[1]=fitlmax
            #ellmask = (fitlmin < ell)*(ell < fitlmax)
        return fitrange
    def read_spectra(self,frame,ax='x'):
        """read 3d spectra"""
        self.vspec=dt.dpy( "%s/DD%04d.products/power_velocity.h5"%(self.simulation.directory,frame) , ['k','power'])
        self.aspec=dt.dpy( "%s/DD%04d.products/power_acceleration.h5"%(self.simulation.directory,frame) , ['k','power'])
        self.dspec=dt.dpy( "%s/DD%04d.products/power_density.h5"%(self.simulation.directory,frame) , ['k','power'])
        self.hspec=dt.dpy( "%s/DD%04d.products/power_magnetic.h5"%(self.simulation.directory,frame) , ['k','power'])


class simulation_package():
    """container for a simulation.
    Keeps track of the data location
    Produces FRBs from simulation data."""
    def __init__(self,directory=".",frames=[], prefix="RUN", 
                  plot_format='png', clobber=False):
        
        self.directory=directory
        self.frames=frames
        self.prefix=prefix
        self.plot_format=plot_format
        self.clobber=clobber #for checking before re-computing. Not used yet.

    def EBall(self):
        """compute all EB products and save them into the frb directory"""
        for frame in self.frames:
            ds = yt.load("%s/DD%04d/data%04d"%(self.directory,frame,frame))
            p49_fields.add_QU(ds)
            self.make_frbs(frame,ds=ds)
            for axis in 'xyz':
                #read and/or compute E,B, and other harmon
                this_proj=self.read_queb(frame,axis) 
                this_proj.compute_harmonic_products()
                this_proj.write()

    def make_frbs(self,frame, axes=['x','y','z'], ds=None):
        fields=[]
        for axis in axes:
          fields.append( (axis,'Q%s'%(axis))   )
          fields.append( (axis,'U%s'%(axis))   )
          fields.append( (axis,'density') )
          fields.append( (axis,'magnetic_field_strength'))

        for axis, field in fields :
            outputdir = "%s/%s/"%(self.directory,frbname)
            if not os.access(outputdir, os.F_OK):
                os.mkdir(outputdir)
            #fix names; Q and U have the name in the field, but others don't.
            if field[0] in 'QU' and field[1] in 'xyz':
                field_name = field
            else:
                field_name = field + "_"+axis
            outfile = outputdir+"/DD%.4d_%s.fits" %(frame,field_name)
            if os.access(outfile, os.F_OK) and not self.clobber:
                print("FRB exists: %s"%outfile)
            else:
                print("FRB being produced: %s"%outfile)
                res = ds.parameters['TopGridDimensions'][0] #2 + ord('x') - ord(axis)]
                proj = ds.proj(field,axis)
                frb = proj.to_frb(1,res)
                hdu = pyfits.PrimaryHDU(frb[field])
                hdulist = pyfits.HDUList([hdu])
                hdulist.writeto(outfile,clobber=True)
                print("wrote", outfile)

    def make_spectra(self,frame):
        """This makes 3d power spectra of velocity, acceleration, 
        magnetic field, and density.  Can be very slow.
        FFTs are stored in DD????.products"""
        oober = st.short_oober(self.directory, frame=frame)
        st.MakeVelocitySpectra(oober,frame)
        st.MakeAccelSpectra(oober,frame)
        st.MakeVelocitySpectra(oober,frame)
        st.MakeMagneticSpectra(oober,frame)
        st.MakeDensitySpectra(oober,frame)

    def read_queb(self,frame,ax='x',bin_style='dx1'):
        """ Read Q,U,E,B,Density,and Magnetic field FRBs.
        All values default to None if the file is not found."""
        frb_dir = "%s/%s"%(self.directory,frbname)
        product_dir = "%s/DD%04d.products"%(self.directory,frame)
        xd='DD'
        Df= "%s/%s%04d_density_%s.fits"%(frb_dir,xd,frame,ax)
        Hf= "%s/%s%04d_magnetic_field_strength_%s.fits"%(frb_dir,xd,frame,ax)
        Qf= "%s/%s%04d_Q%s.fits"%(frb_dir,xd,frame,ax)
        Uf= "%s/%s%04d_U%s.fits"%(frb_dir,xd,frame,ax)
        Ef= "%s/%s%04d_E%s.fits"%(frb_dir,xd,frame,ax)
        Bf= "%s/%s%04d_B%s.fits"%(frb_dir,xd,frame,ax)
        if not os.path.exists(Qf):
            print("Warning: no Q file exists: "+Qf)
        d=read_fits(Df)
        h=read_fits(Hf)
        q=read_fits(Qf)
        u=read_fits(Uf)
        e=read_fits(Ef)
        b=read_fits(Bf)
        ts=queb_snapshot(q,u,d,H=h, E=e,B=b,axis=ax,simulation=self,frame=frame,bin_style=bin_style)
        return ts

    def image_fields(self,frame,axis='x',ts=None):
        if ts is None:
            ts=self.read_queb(frame,axis)
            ts.compute_harmonic_products()
        for name in [ 'T','H','Q','U','E','B']:
            fig,ax=plt.subplots(1,1)
            ax.clear()
            array = ts[name]
            proj=ax.imshow(array,origin='lower',interpolation='nearest')
            fig.colorbar(proj, ax=ax)
            outname = "%s_%04d_%s_%s.png"%(self.prefix,frame,name,axis)
            fig.savefig(outname)
            print(outname)
            plt.close(fig)
        return ts

    def plot_eb(self,ts,fname='TEST.png', slopes=None):
        """plot spectra extracted from e&b.
        the function 'dostuff' treats normalization and slicing."""

        def dostuff(arr):
            return arr[1:]#/np.abs(arr[2])  #np.abs((arr/arr[2])[1:])
        def dostuff2(arr):
            return arr[1:]#/np.abs(arr[2])  #np.abs((arr/arr[2])[1:])

        fig,ax=plt.subplots(1,1)
        k = ts.vspec[0]
        ylimits=dt.extents()
        xlim=dt.extents()
        pos_k=slice(None)
        this_k=(k/k[1])[pos_k]
        this_k=(k*2*np.pi)#[1:]
        this_k = 0.5*(this_k[1:]+this_k[:-1])
        ell=ts.lcent
        #this_ell=(ell/ell[1])[pos_k]
        this_ell=ell[:-1]#[1:] 

        rEB = ts['ClEB']/(ts['ClEE']*ts['ClBB'])**0.5
        lab=r'$r_{EB}=C_{\ell}^{EB}/\sqrt{C_{\ell}^{EE}C_{\ell}^{EE}}$'
        ax.plot( this_k,dostuff2(ts['aspec'][1]),c='k',marker='*', label=r'$a$');   ylimits(ts['aspec'][1][pos_k])# print('a lim',ylimits)
        ax.plot( this_ell,dostuff2(ts['ClEE']),        marker='*', label=r'$EE$',c='g'); ylimits(rEB)# print(ylimits)
        if slopes is not None:
            slopes.plot(ax,'EE',label=r'$\alpha^{EE}$',c='r')
            
        ax.plot( this_ell,dostuff2(ts['ClBB']),        marker='*', label=r'$BB$',c='b'); ylimits(rEB)# print(ylimits)
        ax.plot( this_ell,dostuff(rEB),                marker='*', label=r'$r_{EB}$',c='m'); ylimits(rEB)# print(ylimits)
        dt.axbonk(ax,xlabel='k/k1',ylabel=lab,xscale='log',yscale='log')
        #ax.set_xscale('symlog',linthreshx=1)
        #ax.set_xlim(xlim)
        #ax.set_ylim(1e-9,1e4)
        #ax.set_ylim(ylimits)
        ax.set_yscale('symlog',linthreshy=1e-2)
        ax.set_ylim(-30,30)
        #print("ylimits",ylimits)

        title="t2 %s n%04d %s"%(self.prefix,ts['frame'],ts['axis'])
        ax.legend(loc=0)
        ax.set_title(title)
        fig.savefig(fname)
        print("saved "+fname)
        plt.close(fig)

    def plot_2eb(self,ts0,ts1,slopes0=None,slopes1=None,fname='TEST.png', slopes=None):
        """plot spectra extracted from e&b.
        the function 'dostuff' treats normalization and slicing."""

        def dostuff(arr):
            return arr[1:]#/np.abs(arr[2])  #np.abs((arr/arr[2])[1:])
        def dostuff2(arr):
            return arr[1:]#/np.abs(arr[2])  #np.abs((arr/arr[2])[1:])

        fig,axes=plt.subplots(1,2)
        ax0 = axes[0]; ax1 = axes[1]
        for ns,ts in enumerate([ts0,ts1]):
            slopes = [slopes0,slopes1][ns] #this is clunky, sorry.
            ax=axes[ns]
            ylimits=dt.extents()
            xlim=dt.extents()
            ell=ts.lcent
            this_ell=ell[:-1]#[1:] 

            lab=r'$C_\ell^{EE}$'
            ax.plot( this_ell,dostuff2(ts['ClEE']),        marker='*', label=r'$EE$',c='g')
            #the slope value is added to the label in slopes.plot
            if slopes is not None:
                slopes.plot(ax,'EE',label=r'$\alpha^{EE}$',c='r',norm=1)#ts.bin_style=='5deg')
                
            #ax.plot( this_ell,dostuff2(ts['ClBB']),        marker='*', label=r'$BB$',c='b')
            dt.axbonk(ax,xlabel='k/k1',ylabel=lab,xscale='log',yscale='log')
            ax.set_yscale('symlog',linthreshy=1e-10)
            #ax.set_yscale('log')
            ax.set_ylim(0,0.5e-3)
            ax.set_xlim(0.01,2*np.pi)
            #ax.set_ylim(0,1e-3)
            title=["new binning","old binning"][ns]
            ax.legend(loc=0)
            ax.set_title(title)
        fig.savefig(fname)
        print("saved "+fname)
        plt.close(fig)
    def plot_many_spectra(self,ts,fname='TEST.png', **kwargs):
        def dostuff(arr):
            return np.abs(arr[1:])#/np.abs(arr[2])  #np.abs((arr/arr[2])[1:])
        fig,ax=plt.subplots(1,1)
        k = ts['vspec'][0]
        ylimits=dt.extents()
        xlim=dt.extents()
        pos_k=slice(None)
        this_k=(k/k[1])[pos_k]
        this_k=(k*2*np.pi)#[1:]
        this_k = 0.5*(this_k[1:]+this_k[:-1])
        ax.plot( this_k,dostuff(ts['aspec'][1][pos_k]),marker='*', label=r'$P(a)$');   ylimits(ts['aspec'][1][pos_k])# print('a lim',ylimits)
        ax.plot( this_k,dostuff(ts['vspec'][1][pos_k]),marker='*', label=r'$P(v)$');   ylimits(ts['vspec'][1][pos_k])# print('v lim',ylimits)
        ax.plot( this_k,dostuff(ts['dspec'][1][pos_k]),marker='*', label=r'$P(\rho)$');ylimits(ts['dspec'][1][pos_k])# print('d lim',ylimits)
        ax.plot( this_k,dostuff(ts['hspec'][1][pos_k]),marker='*', label=r'$P(H)$');   ylimits(ts['hspec'][1][pos_k])# print('h lim',ylimits)
        ell=ts['lcent']
        #this_ell=(ell/ell[1])[pos_k]
        this_ell=ell[:-1]#[1:] 
        ax.plot( this_ell,dostuff(ts['ClEE'][pos_k]),marker='*', label=r'$C_{\ell}^{EE}$'); ylimits(ts['ClEE'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClBB'][pos_k]),marker='*', label=r'$C_{\ell}^{BB}$'); ylimits(ts['ClBB'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClTE'][pos_k]),marker='*', label=r'$ClTE$') ; ylimits(ts['ClTE'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClTT'][pos_k]),marker='*', label=r'$ClTT$') ; ylimits(ts['ClTT'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClHH'][pos_k]),marker='*', label=r'$ClHH$') ; ylimits(ts['ClHH'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClTB'][pos_k]),marker='*', label=r'$ClTB$') ; ylimits(ts['ClTB'])# print(ylimits)
        ax.plot( this_ell,dostuff(ts['ClEB'][pos_k]),marker='*', label=r'$ClEB$') ; ylimits(ts['ClEB'])# print(ylimits)
        dt.powerline(ax, this_ell[1]*4, this_ell[1]*10, 1, -5./3,c='k',label='-5/3')
        dt.powerline(ax, this_ell[1]*4, this_ell[1]*10, 1, -2.5,c='k',label='-2.5',linestyle='--')
        #ax.plot( this_ell, np.abs(ts['pork'])[pos_k],marker='*', label=r'pork') ; ylimits(ts['ClTE']); print(ylimits)
        #ts['ColumnDensity']= np.abs(cmbtools.harm2cl( ts['Eh'], ts['Deltal'],ts['lbins']))
        #ax.plot(this_ell,np.abs(ts['ColumnDensity'][pos_k]),c='k')
        #ax.plot( ell/ell[1], ts['ClTB'],marker='*', label=r'$ClEB$')
        xlim(this_k)
        xlim(this_ell)
        
        ax.legend(loc=1)
        ts.limits=ylimits


        #dt.axbonk(ax,xlabel='k/k1',ylabel='power',xscale='log',yscale='log')
        dt.axbonk(ax,xlabel='k/k1',ylabel='power',xscale='log',yscale='log')
        #ax.set_xscale('symlog',linthreshx=1)
        #ax.set_xlim(xlim)
        #ax.set_ylim(1e-9,1e4)
        #ax.set_yscale('symlog',linthreshy=1e-28)
        #ax.set_ylim(ylimits)
        ax.set_ylim(1e-7,50)
        print("ylimits",ylimits)

        title="t2 %s n%04d %s"%(ts['prefix'],ts['frame'],ts['axis'])
        print("POOT",title)
        ax.set_title(title)
        fig.savefig(fname)
        plt.close(fig)
