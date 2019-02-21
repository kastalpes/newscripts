from GL import *
import get_all_quantities as gaq
reload(gaq)
quan=gaq.all_quan_from_outputlog('/home/dcollins/scratch/p49d_quan_test/a01b/OutputLog')
out_prefix = 'a01b'
for q in quan:
    f=quan[q][:10]
    n=len(f)
    print("%10s"%q+" %0.1f"*n%tuple(f))
if 0:
    v2 = np.sqrt(quan['vx_std']**2+quan['vy_std']**2+quan['vz_std']**2)
    b2 = np.sqrt(quan['bx_std']**2+quan['by_std']**2+quan['bz_std']**2+
                 quan['bx_avg']**2+quan['by_avg']**2+quan['bz_avg']**2)/np.sqrt(np.pi*4)
    va = np.sqrt(quan['alf_x_std']+quan['alf_y_std']+quan['alf_z_std'])//np.sqrt(np.pi*4)

    plt.clf()
    plt.plot(quan['time'],va,marker="*")
    plt.xlabel('t');plt.ylabel(r'$Va_{\rm{rms}}$')
    plt.savefig('%s_time_va.png'%out_prefix)
    plt.clf()
    plt.plot(quan['time'],v2,marker="*")
    plt.xlabel('t');plt.ylabel('$V_{\rm{rms}}$')
    plt.savefig('%s_time_mach.png'%out_prefix)
    plt.clf()
    plt.clf()
    plt.plot(quan['time'],b2,marker="*")
    plt.xlabel('t');plt.ylabel('$B_{\rm{rms}}$')
    plt.savefig('%s_time_b.png'%out_prefix)

