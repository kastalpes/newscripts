from GL import *
import yt
import p49_QU2EB
import p49_fields
import cmbtools
import davetools as dt
reload(dt)
frbname="frbs"
reload(p49_QU2EB)
def dpy(filename,fields):
    """Open hdf5 file *filename* and return *field*, then close the file.
    Collapses 3 function calls into 1."""
    if glob.glob(filename) == []:
        print("No such file", filename)
        return None
    fptr = h5py.File(filename,'r')
    output = []
    try:
        if hasattr(fptr,'listnames'):
            field_list = fptr.listnames()
        elif hasattr(fptr,'keys'):
            field_list = fptr.keys()

        for field in dt.ensure_list(fields):
            if field not in field_list:
                print("No field", field, "in file",filename)
                return None
            else:
                if len(fptr[field].shape) > 0:
                    output.append(fptr[field][:])
                else:
                    output.append(fptr[field].value)
    except:
        pdb.set_trace()
        raise
    finally:
        fptr.close()
    return output

def Cl2(Q,U,Density,H=None,BoxSize=1.0):
    """
    N = np.array(arr.shape,dtype = np.int32)
    arr = np.ascontiguousarray(arr,dtype=np.double)
    xsize = N.max() #5 * np.pi / 180*BoxSize
    size2d = np.array([xsize,xsize])
    Delta = size2d/N
    Deltal = cmbtools.Delta2l(Delta,N)
    harm = cmbtools.map2harm(arr,Delta)
    lmax = Deltal[0]*N[0]/2
    lbins = np.arange(N[0]//2+1)*Deltal[0] #np.linspace(0,lmax,N//2)
    lcent = lbins[:-1] + np.diff(lbins)/2.
    ClBB = cmbtools.harm2cl(harm,Deltal,lbins)
    output={'Delta':Delta,'Deltal':Deltal,'harm':harm,'lbins':lbins,'ClBB':ClBB,'lmax':lmax}
    return output
    """
    stuff=p49_QU2EB.EBfromQU(Q,U,T=Density,H=H,BoxSize=BoxSize,return_quharm=True)
    Eharm=stuff['Eh']
    Bharm=stuff['Bh']
    Tharm=stuff['Th']
    Hharm=stuff.get('Hh',None)
    N=np.array(Q.shape,dtype=np.int32)
    xsize = 64*BoxSize
    size2d = np.array([xsize,xsize])
    Delta = size2d/N
    Deltal = cmbtools.Delta2l(Delta,N)
    lmax = Deltal[0]*N[0]/2
    lbins = np.linspace(0,lmax,64*BoxSize/2)
    lcent = lbins[:-1] + np.diff(lbins)/2.
    ClEE = cmbtools.harm2cl(Eharm,Deltal,lbins)
    ClBB = cmbtools.harm2cl(Bharm,Deltal,lbins)
    if Hharm is not None:
        ClHH= cmbtools.harm2cl(Hharm,Deltal,lbins)
        stuff['ClHH']=ClHH
    ClTT = cmbtools.harm2cl(Tharm,Deltal,lbins)
    #Testing the density spectra
    #density_fft = dt.read_fft("/Users/dcollins/scratch/P49d/ca02_turb/DD0001.products/fft_density_x.float32","density")
    #density_fft_test = np.abs(density_fft)[:,:density_fft.shape[1]//2+1].astype(Eharm.dtype)
    #pork = cmbtools.harm2cl(density_fft_test,Deltal,lbins)
    #stuff['ClNN']=pork
    lbins2 = np.mgrid[0:0.5:1./128]
    Deltal2 = 1./128
    ClTE = cmbtools.harm2clcross_samegrid(Tharm,Eharm,Deltal,lbins)
    stuff.update( {'ClEE':ClEE,'ClBB':ClBB,'ClTT':ClTT,'ClTE':ClTE,'lbins':lbins,'lcent':lcent})
    return stuff
    #ClTE = cmbtools.harm2clcross_samegrid(Eharm,Nharm,Deltal,lbins)
    #return {'clee':ClEE,'clbb':ClBB,'ell':lbins,'delta':Delta,'deltal':Deltal,
    #        'Eh':Eharm,'Bh':Bharm,'Qh':Qharm,'Uh':Uharm,'Q':Q,'U':U,'xsize':xsize}
class queb_package():
    def __init__(self, clobber=False,directory=".",frames=[], prefix="RUN", fit_range=None, BoxSize=1):
        self.stuff={}
        self.clobber=clobber
        self.plot_format='pdf'
        self.clobber=clobber
        self.directory=directory
        self.frames=frames
        self.prefix=prefix
        self.fit_range=fit_range
        self.BoxSize=BoxSize
    def EBall(self):
        if 'EBcycles' not in self.stuff:
            self.stuff['EBcycles']=[]
            for ax in 'xyz':
                self.stuff['Eamp_%s'%ax]=[]
                self.stuff['Bamp_%s'%ax]=[]
                self.stuff['Eslope_%s'%ax]=[]
                self.stuff['Bslope_%s'%ax]=[]
        for frame in self.frames:
            ds = yt.load("%s/DD%04d/data%04d"%(self.directory,frame,frame))
            p49_fields.add_QU(ds)
            if frame not in self.stuff['EBcycles']:
                self.make_frbs(frame,ds=ds)
                self.QUEB(frame,BoxSize=self.BoxSize)
                self.EBslopes(frame,fit_range=self.fit_range)
#       self.dump()
    def make_frbs(self,frame, axes=['x','y','z'], ds=None):
        fields=[]
        for axis in axes:
          n0=1; p=1 #n0 in [19,39,1945] and p=0
          #fields.append( (axis,'Q%s_n0-%04d_p-%d'%(axis,n0,p))   )
          #fields.append( (axis,'U%s_n0-%04d_p-%d'%(axis,n0,p))   )
          fields.append( (axis,'Q%s'%(axis))   )
          fields.append( (axis,'U%s'%(axis))   )
          fields.append( (axis,'density') )
          fields.append( (axis,'magnetic_field_strength'))

        for axis, field in fields :
            outputdir = "%s/%s/"%(self.directory,frbname)
            if not os.access(outputdir, os.F_OK):
                os.mkdir(outputdir)
            #Hm.  Q and U have the name in the field, but others don't.
            if field[0] in 'QU' and field[1] in 'xyz':
                field_name = field
            else:
                field_name = field + "_"+axis
            outfile = outputdir+"/DD%.4d_%s.fits" %(frame,field_name)
            #move this in the 'make' conditional?
            #if ds is None:
            #    ds = self.car.load(frame)
            #    res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)] # zyx order
            if os.access(outfile, os.F_OK) and not self.clobber:
                print("FRB exists: %s"%outfile)
            else:
                #if ds is None:
                #    ds = self.car.load(frame)
                #    res = ds.parameters['TopGridDimensions'][2 + ord('x') - ord(axis)] # zyx order
                print("FRB being produced: %s"%outfile)
                res = ds.parameters['TopGridDimensions'][0] #2 + ord('x') - ord(axis)]
                proj = ds.proj(field,axis)
                frb = proj.to_frb(1,res)
                hdu = pyfits.PrimaryHDU(frb[field])
                hdulist = pyfits.HDUList([hdu])
                hdulist.writeto(outfile,clobber=True)
                print("wrote", outfile)
    def QUEB(self, frame, BoxSize=None):
        #ds = self.car.load(frame)
        frb_dir = "%s/%s/"%(self.directory,frbname)
        p49_QU2EB.QU2EB(frb_dir,frame,BoxSize=BoxSize)

    def EBslopes(self,frame, fit_range=None):
        EBSlopePower=p49_QU2EB.slopes_powers(frame,directory = self.directory + "/"+frbname, prefix=self.prefix, plot_format=self.plot_format,fit_range=fit_range )
        if 'EBcycles' not in self.stuff:
            self.stuff['EBcycles']=[]
        self.stuff['EBcycles'].append(frame)
        for ax in 'xyz':
            self.stuff['Eamp_%s'%ax  ].append(EBSlopePower['Eamp'][ax])
            self.stuff['Bamp_%s'%ax  ].append(EBSlopePower['Bamp'][ax])
            self.stuff['Eslope_%s'%ax].append(EBSlopePower['Eslope'][ax])
            self.stuff['Bslope_%s'%ax].append(EBSlopePower['Bslope'][ax])
    def pull_fft(self,frame,fit_range=None,ax='x'):
        """Pulls spectra data from several sources.
        QUEB spectra are read from fits files, and the spectra are computed.
        Fluid quantity spectra are computed elsewhere (bug collins for the code)
        and stored in DD????.products.  
        These are plotted with plot, below.
        """
        frb_dir = "%s/%s"%(self.directory,frbname)
        product_dir = "%s/DD%04d.products"%(self.directory,frame)
        xd='DD'
        Df= "%s/%s%04d_density_%s.fits"%(frb_dir,xd,frame,ax)
        Hf= "%s/%s%04d_magnetic_field_strength_%s.fits"%(frb_dir,xd,frame,ax)
        Qf= "%s/%s%04d_Q%s.fits"%(frb_dir,xd,frame,ax)
        Uf= "%s/%s%04d_U%s.fits"%(frb_dir,xd,frame,ax)
        Ef= "%s/%s%04d_E%s.fits"%(frb_dir,xd,frame,ax)
        Bf= "%s/%s%04d_B%s.fits"%(frb_dir,xd,frame,ax)
        d=np.ascontiguousarray(pyfits.open(Df)[0].data,dtype=np.double)
        h=np.ascontiguousarray(pyfits.open(Hf)[0].data,dtype=np.double)
        q=np.ascontiguousarray(pyfits.open(Qf)[0].data,dtype=np.double)
        u=np.ascontiguousarray(pyfits.open(Uf)[0].data,dtype=np.double)
        e=np.ascontiguousarray(pyfits.open(Ef)[0].data,dtype=np.double)
        b=np.ascontiguousarray(pyfits.open(Bf)[0].data,dtype=np.double)
        import cmbtools
        #x = cmbtools.map2harm(Q,np.ones(2))
        ts=Cl2(q,u,d,H=h,BoxSize=self.BoxSize)

        product_dir = "%s/DD%04d.products"%(self.directory,frame)
        #   N_fft = dt.read_fft('%s/fft_density_%s.float32'%(product_dir,ax),'density')
        #ts['N_fft']=N_fft[:,:N_fft.shape[1]//2+1]

        plt.clf()
        ppp=plt.imshow( np.log10(np.abs( ts['Th'])))
        plt.colorbar(ppp)
        plt.savefig('AbsTharm_n%04d_%s.png'%(frame,ax))
        plt.clf()
        #ppp=plt.imshow( np.log10(np.abs(ts['N_fft'])))
        #plt.colorbar(ppp)
        #plt.savefig('AbsDensity_n%04d_%s.png'%(frame,ax))
        return ts
    def make_spectra(self,frame):
        import spectra_tools as st
        reload(st)
        oober = st.short_oober(self.directory, frame=frame)
        st.MakeVelocitySpectra(oober,frame)
        st.MakeAccelSpectra(oober,frame)
        st.MakeVelocitySpectra(oober,frame)
        st.MakeMagneticSpectra(oober,frame)
        st.MakeDensitySpectra(oober,frame)



    def pull_spectra(self,frame,fit_range=None,ax='x'):
        """Pulls spectra data from several sources.
        QUEB spectra are read from fits files, and the spectra are computed.
        Fluid quantity spectra are computed elsewhere (bug collins for the code)
        and stored in DD????.products.  
        These are plotted with plot, below.
        """
        frb_dir = "%s/%s"%(self.directory,frbname)
        xd='DD'
        Df= "%s/%s%04d_density_%s.fits"%(frb_dir,xd,frame,ax)
        Hf= "%s/%s%04d_magnetic_field_strength_%s.fits"%(frb_dir,xd,frame,ax)
        Qf= "%s/%s%04d_Q%s.fits"%(frb_dir,xd,frame,ax)
        Uf= "%s/%s%04d_U%s.fits"%(frb_dir,xd,frame,ax)
        Ef= "%s/%s%04d_E%s.fits"%(frb_dir,xd,frame,ax)
        Bf= "%s/%s%04d_B%s.fits"%(frb_dir,xd,frame,ax)
        d=np.ascontiguousarray(pyfits.open(Df)[0].data,dtype=np.double)
        h=np.ascontiguousarray(pyfits.open(Hf)[0].data,dtype=np.double)
        q=np.ascontiguousarray(pyfits.open(Qf)[0].data,dtype=np.double)
        u=np.ascontiguousarray(pyfits.open(Uf)[0].data,dtype=np.double)
        e=np.ascontiguousarray(pyfits.open(Ef)[0].data,dtype=np.double)
        b=np.ascontiguousarray(pyfits.open(Bf)[0].data,dtype=np.double)
        import cmbtools
        #x = cmbtools.map2harm(Q,np.ones(2))
        ts=Cl2(q,u,d,H=h,BoxSize=self.BoxSize)
        self.make_spectra(frame)
        ts['vspec']=dpy( "%s/DD%04d.products/power_velocity.h5"%(self.directory,frame) , ['k','power'])
        ts['aspec']=dpy( "%s/DD%04d.products/power_acceleration.h5"%(self.directory,frame) , ['k','power'])
        ts['dspec']=dpy( "%s/DD%04d.products/power_density.h5"%(self.directory,frame) , ['k','power'])
        ts['hspec']=dpy( "%s/DD%04d.products/power_magnetic.h5"%(self.directory,frame) , ['k','power'])
        ts['frame']=frame
        ts['axis']=ax
        ts['prefix']=self.prefix
        product_dir = "%s/DD%04d.products"%(self.directory,frame)
        #N_fft = dt.read_fft('%s/fft_density_%s.float32'%(product_dir,ax),'density')
        #ts['N_fft']=N_fft[:,:N_fft.shape[1]//2+1]
        return ts
    def plot(self,ts,fname='TEST.png', **kwargs):
        def dostuff(arr):
            return (arr/arr[2])[1:]
        fig,ax=plt.subplots(1,1)
        k = ts['vspec'][0]
        limits=dt.extents()
        xlim=dt.extents()
        pos_k=slice(None)
        this_k=(k/k[1])[pos_k]
        this_k=(k*2*np.pi)[1:]
        ax.plot( this_k,dostuff(np.abs(ts['aspec'][1][pos_k])),marker='*', label=r'$P(a)$');   limits(ts['aspec'][1][pos_k]); print('a lim',limits)
        ax.plot( this_k,dostuff(np.abs(ts['vspec'][1][pos_k])),marker='*', label=r'$P(v)$');   limits(ts['vspec'][1][pos_k]); print('v lim',limits)
        ax.plot( this_k,dostuff(np.abs(ts['dspec'][1][pos_k])),marker='*', label=r'$P(\rho)$');limits(ts['dspec'][1][pos_k]); print('d lim',limits)
        ax.plot( this_k,dostuff(np.abs(ts['hspec'][1][pos_k])),marker='*', label=r'$P(H)$');   limits(ts['hspec'][1][pos_k]); print('h lim',limits)
        ell=ts['lcent']
        #this_ell=(ell/ell[1])[pos_k]
        this_ell=ell[1:]
        ax.plot( this_ell,dostuff(np.abs(ts['ClEE'])[pos_k]),marker='*', label=r'$C_{\ell}^{EE}$'); limits(ts['ClEE']); print(limits)
        ax.plot( this_ell,dostuff(np.abs(ts['ClBB'])[pos_k]),marker='*', label=r'$C_{\ell}^{BB}$'); limits(ts['ClBB']); print(limits)
        ax.plot( this_ell,dostuff(np.abs(ts['ClTE'])[pos_k]),marker='*', label=r'$ClTE$') ; limits(ts['ClTE']); print(limits)
        ax.plot( this_ell,dostuff(np.abs(ts['ClTT'])[pos_k]),marker='*', label=r'$ClTT$') ; limits(ts['ClTT']); print(limits)
        ax.plot( this_ell,dostuff(np.abs(ts['ClHH'])[pos_k]),marker='*', label=r'$ClHH$') ; limits(ts['ClHH']); print(limits)
        dt.powerline(ax, this_ell[1]*4, this_ell[1]*10, 1, -5./3,c='k',label='-5/3')
        dt.powerline(ax, this_ell[1]*4, this_ell[1]*10, 1, -2.5,c='k',label='-2.5',linestyle='--')
        #ax.plot( this_ell, np.abs(ts['pork'])[pos_k],marker='*', label=r'pork') ; limits(ts['ClTE']); print(limits)
        #ts['ColumnDensity']= np.abs(cmbtools.harm2cl( ts['Eh'], ts['Deltal'],ts['lbins']))
        #ax.plot(this_ell,np.abs(ts['ColumnDensity'][pos_k]),c='k')
        #ax.plot( ell/ell[1], ts['ClTB'],marker='*', label=r'$ClEB$')
        xlim(this_k)
        xlim(this_ell)
        ax.legend(loc=1)
        ts['limits']=limits


        #dt.axbonk(ax,xlabel='k/k1',ylabel='power',xscale='log',yscale='log')
        dt.axbonk(ax,xlabel='k/k1',ylabel='power',xscale='log',yscale='log')
        #ax.set_xscale('symlog',linthreshx=1)
        #ax.set_xlim(xlim)
        #ax.set_ylim(limits)
        ax.set_ylim(1e-9,1e4)
        #ax.set_yscale('symlog',linthreshy=1e-28)
        #ax.set_ylim(0,1e-4)
        ax.set_title("t2 %s n%04d %s"%(ts['prefix'],ts['frame'],ts['axis']))
        fig.savefig(fname)
        plt.close(fig)
