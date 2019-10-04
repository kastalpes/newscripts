from GL import *
import FRB_QUEB as fr
reload(fr)
frame=10
pack = fr.simulation_package()
#pack = fr.simulation_package(directory="/Users/dcollins/scratch/P49d/ca01_cyl",frames=[0],prefix="ca01_")
#pack.EBall(directory="/Users/dcollins/scratch/P19/u05-r4-l4-128/",frames=[125],prefix="u05_",fit_range=[80,200], BoxSize=1)
<<<<<<< HEAD
pack = fr.simulation_package( directory="/Users/dcollins/scratch/P49d/ca02_turb",frames=[frame],prefix="ca02",fit_range=[80,200], BoxSize=2)
#pack = fr.simulation_package( directory="/Users/dcollins/scratch/P49d/ca01_cyl",frames=[1],prefix="ca01",fit_range=[80,200], BoxSize=2)
#pack = fr.simulation_package( directory="/Users/dcollins/repos/p49c/spiral_test_x1",frames=[0],prefix="spiral_x",fit_range=[80,200], BoxSize=2)
pack.EBall()
#pack.image_fields(frame=frame,axis='x')
#pack.make_spectra(frame)
for frame in [frame]:# list(range(10,150,10))+ list(range(1,10)):
    tsx=pack.read_spectra(frame=frame,ax='x')
    pack.plot_eb(tsx, fname='ca02_reb_x_n%04d.png'%frame)
    tsy=pack.read_spectra(frame=frame,ax='y')
    #pack.plot_many_spectra(tsx,fname = 'ca02_spectra_symlog_x_n%04d.png'%frame)
#    #pack.plot_many_spectra(tsy,fname = 'ca02_spectra_symlog_y_n%04d.png'%frame)
##    #fftx=pack.pull_fft(frame=1,ax='x')
=======
#pack = fr.queb_package( directory="/Users/dcollins/scratch/P49d/ca02_turb",frames=[10],prefix="ca02",fit_range=[80,200], BoxSize=2)
pack = fr.queb_package( directory="/scratch1/dcollins/Paper49_EBQU/ca02_turb",frames=[10],prefix="ca02",fit_range=[80,200], BoxSize=2)
pack.plot_format='png'
#pack.EBall()
for frame in range(250,500,10):
    #pack.make_spectra(frame=frame)
    pack.EBall(frame=frame)
    tsx=pack.pull_spectra(frame=frame,ax='x')
    tsy=pack.pull_spectra(frame=frame,ax='y')
    pack.plot(tsx,fname = 'TESTX_%04d.png'%frame)
    pack.plot(tsy,fname = 'TESTY_%04d.png'%frame)
    #fftx=pack.pull_fft(frame=1,ax='x')
>>>>>>> 5a7ced317edb7fa3a984d5d6c29711b2ac12cf10
