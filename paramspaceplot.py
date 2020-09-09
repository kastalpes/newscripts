from GL import *
from cycler import cycler
import h5py

plotdir='/Users/Kye/512reruns/512rerunplots/general'
simlist=["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"]

for los in ['x','y','z']:
  for runs in ['ee','bb']:
    Fptr = h5py.File("slopes_ratios_dict_%s.h5"%los,"r")
    slopes = Fptr["%s_slopes"%runs][:]
    ratios = Fptr["eb_ratios"][:]
    Fptr.close()
    slopemax=np.abs(np.min(slopes)-1) #using min because slopes are negative
    ratiomax=np.max(ratios)+1
    slopex = np.linspace(1,slopemax,10)
    ratioy = np.linspace(1,ratiomax,10)
    cm=plt.get_cmap('gist_rainbow')
    fig = plt.figure()
    axis=fig.add_subplot(111)
#axis.set_prop_cycle(color=[cm(1.*i/len(labels)) for i in range(len(labels))]) #set colormap to cycle colors i
    axis.set_prop_cycle(cycler(color=['r','g','b','r','g','b','r','g','b','r','g','b'])+cycler(marker=['.','.','.','x','x','x','*','*','*','1','1','1']))
    for i in range(len(simlist)):
        axis.plot(np.abs(slopes[i]), ratios[i],ms=5,label='%s'%simlist[i])

    axis.plot(slopex, np.zeros_like(slopex)+2.0,label='Ratio=2.0', c=[0.5]*3,marker='')
    axis.plot(np.zeros_like(ratioy)+2.43, ratioy,label='|Slope|=2.43', c=[0.5]*3,marker='')

    axes = plt.gca()
    axes.set_xlim([1,slopemax])
    axes.set_ylim([1,ratiomax])
#axis.set_title('Parameter Space of Averages')
    axis.legend(loc=2)
    axis.set_xlabel(r'Average |Slope| of $Cl%s_{%s}$($\ell$)'%(runs,los))
    axis.set_ylabel(r'Average Ratio of $Clee_{%s}$($\ell$)/$Clbb_{%s}$($\ell$)'%(los,los))
    fig.savefig("%s/ParameterSpace_%s_%s.png"%(plotdir,runs,los))
    plt.clf()
    plt.close()
