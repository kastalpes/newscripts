from GL import *
import FRB_QUEB as fr
reload(fr)
frames=list(range(0,500,10)) + list(range(0,10))
frames=[499]
#pack = fr.simulation_package(directory="/Users/dcollins/scratch/P49d/ca01_cyl",frames=[0],prefix="ca01_")
#pack.EBall(directory="/Users/dcollins/scratch/P19/u05-r4-l4-128/",frames=[125],prefix="u05_",fit_range=[80,200], BoxSize=1)
#pack = fr.simulation_package( directory="/Users/dcollins/scratch/P49d/ca02_turb",frames=[frame],prefix="ca02",fit_range=[80,200], BoxSize=2)
pack = fr.simulation_package( directory="/scratch1/dcollins/Paper49_EBQU/ca02_turb",frames=frames,prefix="ca02", BoxSize=2)
#pack = fr.simulation_package( directory="/Users/dcollins/scratch/P49d/ca01_cyl",frames=[1],prefix="ca01",fit_range=[80,200], BoxSize=2)
#pack = fr.simulation_package( directory="/Users/dcollins/repos/p49c/spiral_test_x1",frames=[0],prefix="spiral_x",fit_range=[80,200], BoxSize=2)
#pack.EBall()
for frame in frames:
    this_proj=pack.read_queb(frame=frame,ax='x')
    this_proj.read_spectra(frame)
    this_proj.compute_harmonic_products()
    fitrange=this_proj.determine_fit_range()  #something better 
    slopes=this_proj.fit_eb_slopes(fitrange=fitrange)

    #slopes=tsx.fit_slopes()
    #pack.image_fields(frame=frame,axis='x')
    pack.plot_eb(this_proj, fname='ca02_reb_x_n%04d.png'%frame)
    #pack.plot_many_spectra(tsx,fname = 'ca02_spectra_symlog_x_n%04d.png'%frame)
    #tsy=pack.read_spectra(frame=frame,ax='y')
    #pack.plot_many_spectra(tsx,fname = 'ca02_spectra_x_n%04d.png'%frame)
#    #pack.plot_many_spectra(tsy,fname = 'ca02_spectra_symlog_y_n%04d.png'%frame)
##    #fftx=pack.pull_fft(frame=1,ax='x')
