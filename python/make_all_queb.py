from GL import *
import FRB_QUEB as fr
reload(fr)
pack = fr.queb_package()
#pack = fr.queb_package(directory="/Users/dcollins/scratch/P49d/ca01_cyl",frames=[0],prefix="ca01_")
#pack.EBall(directory="/Users/dcollins/scratch/P19/u05-r4-l4-128/",frames=[125],prefix="u05_",fit_range=[80,200], BoxSize=1)
pack = fr.queb_package( directory="/Users/dcollins/scratch/P49d/ca02_turb",frames=[10],prefix="ca02",fit_range=[80,200], BoxSize=2)
pack.plot_format='png'
pack.make_spectra(1)
#pack.EBall()
#for frame in list(range(10,150,10))+ list(range(1,10)):
#    tsx=pack.pull_spectra(frame=frame,ax='x')
#    tsy=pack.pull_spectra(frame=frame,ax='y')
#    pack.plot(tsx,fname = 'TESTX_%04d.png'%frame)
#    pack.plot(tsy,fname = 'TESTY_%04d.png'%frame)
#    #fftx=pack.pull_fft(frame=1,ax='x')
