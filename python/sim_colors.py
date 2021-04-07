#color is Alfven
#line is Sonic
#marker is Sonic
from GL import *

#these need to stay in this order.
simlist=nar(["half_half","half_1","half_2","1_half","1_1","1_2","2_half","2_1","2_2","3_half","3_1","3_2"])
color_list=nar(['r','g','b','r','g','b','r','g','b','r','g','b'])
line_list=nar(['-','-','-','-.','-.','-.','--','--','--',':',':',':'])
marker_list = nar(['.','.','.','s','s','s','^','^','^','*','*','*'])

framelist=[range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(11,31),range(65,85),range(11,31),range(11,31),range(72,91),range(56,75),range(20,40)]

sim_ms = nar(['half','1','2','3'])
sim_ma = nar(['half','1','2'])
plot_order=[]
for ma in sim_ma:
    for ms in sim_ms:
        sim="%s_%s"%(ms,ma)
        plot_order.append(sim)

color=dict(zip(simlist,color_list))
linestyle = dict(zip(simlist,line_list))
marker = dict(zip(simlist,marker_list))
frames = dict(zip(simlist,framelist))

def vals_from_sim(sim):
    ms,ma = sim.split("_")
    if ms == 'half':
        ms = 0.5
    if ma == 'half':
        ma = 0.5
    ms=float(ms)
    ma=float(ma)
    return ms,ma

ms_list=[]
ma_list=[]
for sim in simlist:
    ms,ma = vals_from_sim(sim)
    ms_list.append( ms)
    ma_list.append(ma)
ms_list=nar(ms_list)
ma_list=nar(ma_list)
Ms = dict(zip(simlist,ms_list))
Ma = dict(zip(simlist,ma_list))

