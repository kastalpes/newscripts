from go import *
import astropy.io.fits as pyfits

sims = ["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
sims = ["half_half","1_half","2_half","3_half","half_1","1_1","2_1","3_1","half_2","1_2","2_2","3_2"]
basedir="/archive2/kas14d"
outdir="./plots_to_sort"
axis='x'
density_array=[]
bx_array=[]
by_array=[]
bz_array=[]
for sim in sims:
    directory = "%s/frbs/%s"%(basedir,sim)
    #pick the last density frb
    density_frb_list = glob.glob("%s/*density_%s.fits"%(directory,axis))
    last_density = sorted(density_frb_list)[-1]
    density_array.append( pyfits.open(last_density)[0].data ) 

    bx_frb_list = glob.glob("%s/*Bx.fits"%(directory))
    last_bx = sorted(bx_frb_list)[-1]
    bx_array.append( pyfits.open(last_bx)[0].data ) 

    by_frb_list = glob.glob("%s/*By.fits"%(directory))
    last_by = sorted(by_frb_list)[-1]
    by_array.append( pyfits.open(last_by)[0].data ) 

    bz_frb_list = glob.glob("%s/*Bz.fits"%(directory))
    last_bz = sorted(bz_frb_list)[-1]
    bz_array.append( pyfits.open(last_bz)[0].data ) 

fig, axes = plt.subplots(3,4)
ax_list=axes.flatten()
fig.subplots_adjust(wspace=0, hspace=0)

for i in range(len(density_array)):
    ax_list[i].imshow( np.log10( density_array[i]),origin='lower',interpolation='nearest')
    ax_list[i].set_xticks([])
    ax_list[i].set_yticks([])
    ax_list[i].text(10,10,sims[i],c=[1.]*3)

fig.savefig("%s/multiplot.pdf"%outdir)
plt.close(fig)
