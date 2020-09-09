import pyfits
import multi_imshow
from multi_imshow import plotter2

#Script plots projections of fields together on same image
#1=B,2=H,3=rho,4=V
j=1
k=2
for axes in ['x','y','z']:
    simlist=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
    framelist=[range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(65,85),range(11,31),range(11,31),range(72,91),range(56,75),range(20,40)]
    fields=['B','H','Density','V']
    for i in range(12):
      frames=framelist[i]
      simdes=simlist[i]
      field1=fields[k]
      field2=fields[j]
      sim_dir = "/Users/Kye/512reruns/frbs/%s"%simdes
      plot_dir = "/Users/Kye/512reruns/512rerunplots/%s"%simdes
      for frame in frames:
        fname='%s/multiim_%s_%s_%04d_%s_%s'%(plot_dir,simdes,axes,frame,field1,field2)
        hdulist1=pyfits.open('%s/DD%04d_magnetic_field_strength_%s.fits'%(sim_dir,frame,axes))
        hdulist2=pyfits.open('%s/DD%04d_density_%s.fits'%(sim_dir,frame,axes))
        plotter2([hdulist2[0].data,hdulist1[0].data],fname,norm='ind',labels=[field1,field2])

