from GL import *
import h5py
import get_all_quantities as gaq
reload(gaq)
#Makes quan plots for each frame and each sim
#quan=gaq.all_quan_from_outputlog('/scratch/00369/tg456484/Paper49d_moresims/zd01_M10_MA1_512_quan')
#'/scratch/00369/tg456484/Paper49d_moresims/zd01_M10_MA1_512_quan'
# out_prefix = 'zd01'
simdic=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]
frames=[range(12,31),range(12,31),range(12,31),range(12,31),range(12,31),range(12,31),range(66,85),range(12,31),range(12,31),range(72,91),range(57,75),range(21,40)]

for i in range(12):
  direc='/Users/Kye/512reruns/frbs/%s'%simdic[i]
  plotdir='/Users/Kye/512reruns/512rerunplots/%s'%simdic[i]
#    all_files =gaq.all_quan_from_files(direc,frames[i])
#    files_to_use = all_files #all_files[slice(None,None,30)]+all_files[-1:]
#    quan = gaq.return_average_quantities(files_to_use)
  v2avg=[]
  b2avg=[]
  vaavg=[]
  maavg=[]
  v2divb2avg=[]
  timeavg=[]

  for frame in frames[i]:
    quan=h5py.File('%s/data%04d.AverageQuantities.h5'%(direc,frame))

    #These are code quantities, so you do not want a sqrt(4pi) if comparing to other code values
    v2 = np.sqrt(quan['vx_std'][:]**2+quan['vy_std'][:]**2+quan['vz_std'][:]**2)
    b2 = np.sqrt(quan['bx_std'][:]**2+quan['by_std'][:]**2+quan['bz_std'][:]**2+
                 quan['bx_avg'][:]**2+quan['by_avg'][:]**2+quan['bz_avg'][:]**2)
    va = np.sqrt(quan['alf_x_std'][:]**2+quan['alf_y_std'][:]**2+quan['alf_z_std'][:]**2)
    #Ma should be va/v2 without sqrt4pi in denominator
    ma = np.sqrt(quan['alf_x_std'][:]**2+quan['alf_y_std'][:]**2+quan['alf_z_std'][:]**2)/ np.sqrt(quan['vx_std'][:]**2+quan['vy_std'][:]**2+quan['vz_std'][:]**2)
    v2avg=np.append(v2avg,v2)
    b2avg=np.append(b2avg,b2)
    vaavg=np.append(vaavg,va)
    maavg=np.append(maavg,ma)
    v2divb2avg=np.append(v2divb2avg,v2/b2)
    timeavg=np.append(timeavg,quan['time'][:])
    plt.clf()
    plt.plot(quan['time'][:],quan['vx_avg'][:],marker="*",label='vx_avg')
    plt.plot(quan['time'][:],quan['vy_avg'][:],marker="*",label='vy_avg')
    plt.plot(quan['time'][:],quan['vz_avg'][:],marker="*",label='vz_avg')
    plt.xlabel('t');plt.ylabel(r'$\langle v_i \rangle$')
    plt.legend(loc=0)
    plt.savefig('%s/%s_time_vi_avg_%04d.png'%(plotdir,simdic[i],frame))

    plt.clf()
    plt.plot(quan['time'][:],quan['vx_std'][:],marker="*",label='vx_rms')
    plt.plot(quan['time'][:],quan['vy_std'][:],marker="*",label='vy_rms')
    plt.plot(quan['time'][:],quan['vz_std'][:],marker="*",label='vz_rms')
    plt.xlabel('t');plt.ylabel(r'$\sqrt{ \langle v_i^2\rangle - \langle v_i \rangle^2}$')
    plt.legend(loc=0)
    plt.savefig('%s/%s_time_vi_rms_%04d.png'%(plotdir,simdic[i],frame))

    plt.clf()
    plt.plot(quan['time'][:],quan['bx_avg'][:],marker="*",label='bx_avg')
    plt.plot(quan['time'][:],quan['by_avg'][:],marker="*",label='by_avg')
    plt.plot(quan['time'][:],quan['bz_avg'][:],marker="*",label='bz_avg')
    plt.xlabel('t');plt.ylabel(r'$\langle B_i \rangle$')
    plt.legend(loc=0)
    plt.savefig('%s/%s_time_bi_avg_%04d.png'%(plotdir,simdic[i],frame))

    plt.clf()
    plt.plot(quan['time'][:],quan['bx_std'][:],marker="*",label='bx_rms')
    plt.plot(quan['time'][:],quan['by_std'][:],marker="*",label='by_rms')
    plt.plot(quan['time'][:],quan['bz_std'][:],marker="*",label='bz_rms')
    plt.xlabel('t');plt.ylabel(r'$\sqrt{ \langle B_i^2\rangle - \langle B_i \rangle^2}$')
    plt.legend(loc=0)
    plt.savefig('%s/%s_time_bi_rms_%04d.png'%(plotdir,simdic[i],frame))
    plt.clf()

  plt.plot(timeavg,vaavg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$Va_{\rm{rms}}$')
  plt.savefig('%s/%s_time_va.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,v2avg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$V_{\rm{rms}}$')
  plt.savefig('%s/%s_time_mach.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,b2avg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$\sqrt{\langle B_i B_i\rangle}$')
  plt.savefig('%s/%s_time_b.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,v2avg/b2avg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$\langle v^2\rangle/\langle B^2 \rangle$')
  plt.savefig('%s/%s_time_v2divb2.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,maavg,marker="*")
  plt.xlabel('t');plt.ylabel(r'$Ma_{\rm{rms}}$')
  plt.savefig('%s/%s_time_ma.png'%(plotdir,simdic[i]))
  plt.clf()
  plt.plot(timeavg,vaavg,label=r'$V_a$')
  plt.plot(timeavg,v2avg,label=r'$V_{rms}$')
  plt.plot(timeavg,b2avg,label=r'$B_{rms}$')
  plt.plot(timeavg,maavg,label=r'$M_a$')
  plt.plot(timeavg,v2divb2avg,label=r'$V/B$')
  plt.xlabel('t')
  plt.legend(loc=0)
  plt.savefig('%s/%s_all_quan.png'%(plotdir,simdic[i]))
