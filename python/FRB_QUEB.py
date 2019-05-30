from GL import *
import yt
import p49_QU2EB
import p49_fields

frbname="frbs"
class queb_package():
    def __init__(self, clobber=False):
        self.stuff={}
        self.clobber=clobber
        self.plot_format='pdf'
    def EBall(self,directory=".",frames=[], prefix="RUN"):
        self.prefix=prefix
        self.directory=directory
        if 'EBcycles' not in self.stuff:
            self.stuff['EBcycles']=[]
            for ax in 'xyz':
                self.stuff['Eamp_%s'%ax]=[]
                self.stuff['Bamp_%s'%ax]=[]
                self.stuff['Eslope_%s'%ax]=[]
                self.stuff['Bslope_%s'%ax]=[]
        for frame in frames:
            ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
            p49_fields.add_QU(ds)
            if frame not in self.stuff['EBcycles']:
                self.make_frbs(frame,ds=ds)
                self.QUEB(frame)
                self.EBslopes(frame)
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
    def QUEB(self, frame):
        #ds = self.car.load(frame)
        frb_dir = "%s/%s/"%(self.directory,frbname)
        p49_QU2EB.QU2EB(frb_dir,frame)

    def EBslopes(self,frame):
        EBSlopePower=p49_QU2EB.slopes_powers(frame,directory = self.directory + "/"+frbname, prefix=self.prefix, plot_format=self.plot_format)
        if 'EBcycles' not in self.stuff:
            self.stuff['EBcycles']=[]
        self.stuff['EBcycles'].append(frame)
        for ax in 'xyz':
            self.stuff['Eamp_%s'%ax  ].append(EBSlopePower['Eamp'][ax])
            self.stuff['Bamp_%s'%ax  ].append(EBSlopePower['Bamp'][ax])
            self.stuff['Eslope_%s'%ax].append(EBSlopePower['Eslope'][ax])
            self.stuff['Bslope_%s'%ax].append(EBSlopePower['Bslope'][ax])
