from GL import *
import get_all_quantities as gaq
reload(gaq)
#quan=gaq.all_quan_from_outputlog('/scratch/00369/tg456484/Paper49d_moresims/zd01_M10_MA1_512_quan')
#'/scratch/00369/tg456484/Paper49d_moresims/zd01_M10_MA1_512_quan'
# out_prefix = 'zd01'
def v2_b2_plot(out_prefix,directory):
    all_files =gaq.files_from_output(directory)
    files_to_use = all_files #all_files[slice(None,None,30)]+all_files[-1:]
    quan = gaq.return_average_quantities(files_to_use)

    #These are code quantities, so you do not want a sqrt(4pi)
    v2 = np.sqrt(quan['vx_std']**2+quan['vy_std']**2+quan['vz_std']**2)
    b2 = np.sqrt(quan['bx_std']**2+quan['by_std']**2+quan['bz_std']**2+
                 quan['bx_avg']**2+quan['by_avg']**2+quan['bz_avg']**2)
    #va = np.sqrt(quan['alf_x_std']+quan['alf_y_std']+quan['alf_z_std'])//np.sqrt(np.pi*4)

    plt.clf()
    plt.plot(quan['time'],quan['vx_avg'],marker="*",label='vx_avg')
    plt.plot(quan['time'],quan['vy_avg'],marker="*",label='vy_avg')
    plt.plot(quan['time'],quan['vz_avg'],marker="*",label='vz_avg')
    plt.xlabel('t');plt.ylabel(r'$\langle v_i \rangle$')
    plt.legend(loc=0)
    plt.savefig('%s_time_vi_avg.png'%out_prefix)

    plt.clf()
    plt.plot(quan['time'],quan['vx_std'],marker="*",label='vx_rms')
    plt.plot(quan['time'],quan['vy_std'],marker="*",label='vy_rms')
    plt.plot(quan['time'],quan['vz_std'],marker="*",label='vz_rms')
    plt.xlabel('t');plt.ylabel(r'$\sqrt{ \langle v_i^2\rangle - \langle v_i \rangle^2}$')
    plt.legend(loc=0)
    plt.savefig('%s_time_vi_rms.png'%out_prefix)

    plt.clf()
    plt.plot(quan['time'],quan['bx_avg'],marker="*",label='bx_avg')
    plt.plot(quan['time'],quan['by_avg'],marker="*",label='by_avg')
    plt.plot(quan['time'],quan['bz_avg'],marker="*",label='bz_avg')
    plt.xlabel('t');plt.ylabel(r'$\langle B_i \rangle$')
    plt.legend(loc=0)
    plt.savefig('%s_time_bi_avg.png'%out_prefix)

    plt.clf()
    plt.plot(quan['time'],quan['bx_std'],marker="*",label='bx_rms')
    plt.plot(quan['time'],quan['by_std'],marker="*",label='by_rms')
    plt.plot(quan['time'],quan['bz_std'],marker="*",label='bz_rms')
    plt.xlabel('t');plt.ylabel(r'$\sqrt{ \langle B_i^2\rangle - \langle B_i \rangle^2}$')
    plt.legend(loc=0)
    plt.savefig('%s_time_bi_rms.png'%out_prefix)

    #plt.clf()
    #plt.plot(quan['time'],va,marker="*")
    #plt.xlabel('t');plt.ylabel(r'$Va_{\rm{rms}}$')
    #plt.savefig('%s_time_va.png'%out_prefix)
    plt.clf()
    plt.plot(quan['time'],v2,marker="*")
    plt.xlabel('t');plt.ylabel(r'$V_{\rm{rms}}$')
    plt.savefig('%s_time_mach.png'%out_prefix)
    plt.clf()
    plt.plot(quan['time'],b2,marker="*")
    plt.xlabel('t');plt.ylabel(r'$\sqrt{\langle B_i B_i\rangle}$')
    plt.savefig('%s_time_b.png'%out_prefix)
    plt.clf()
    plt.plot(quan['time'],v2/b2,marker="*")
    plt.xlabel('t');plt.ylabel(r'$\langle v^2\rangle/\langle B^2 \rangle$')
    plt.savefig('%s_time_Ma1.png'%out_prefix)

