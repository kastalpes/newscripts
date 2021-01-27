#
# Import NumPy for array handling
#
import numpy as np
import math


if 0:
#
# Import plotting libraries (start Python with ipython --matplotlib)
#
#from mpl_toolkits.mplot3d import axes3d
#from matplotlib import pyplot as plt
#
# Some natural constants
#
    au  = 1.49598e13     # Astronomical Unit       [cm]
    pc  = 3.08572e18     # Parsec                  [cm]
    ms  = 1.98892e33     # Solar mass              [g]
    ts  = 5.78e3         # Solar temperature       [K]
    ls  = 3.8525e33      # Solar luminosity        [erg/s]
    rs  = 6.96e10        # Solar radius            [cm]
#
# Monte Carlo parameters
#
#
# Grid parameters
#
    nx       = 64
    ny       = 64
    nz       = 64
    sizex    = 10*au
    sizey    = 10*au
    sizez    = 10*au
#
# Model parameters
#
    radius   = 5*au
    rho0     = 1e-16
#
# Star parameters
#
    mstar    = ms
    rstar    = rs
    tstar    = ts
    pstar    = np.array([0.,0.,0.])
#
# Make the coordinates
#
    xi       = np.linspace(-sizex,sizex,nx+1)
    yi       = np.linspace(-sizey,sizey,ny+1)
    zi       = np.linspace(-sizez,sizez,nz+1)
    xc       = 0.5 * ( xi[0:nx] + xi[1:nx+1] )
    yc       = 0.5 * ( yi[0:ny] + yi[1:ny+1] )
    zc       = 0.5 * ( zi[0:nz] + zi[1:nz+1] )
#
# Make the dust density model
#
    qq       = np.meshgrid(xc,yc,zc,indexing='ij')
    xx       = qq[0]
    yy       = qq[1]
    zz       = qq[2]
    rr       = np.sqrt(xx**2+yy**2+zz**2)
    rhod     = rho0 * np.exp(-(rr**2/radius**2)/2.0)
#
# Make the wavelength_micron.inp file
#
if 1:
    lam1     = 0.1e0
    lam2     = 7.0e0
    lam3     = 25.e0
    lam4     = 1.0e4
    n12      = 20
    n23      = 100
    n34      = 30
    lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
    lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
    lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
    lam      = np.concatenate([lam12,lam23,lam34])
    nlam     = lam.size
#
# parse the dustkapscatmat file to check how many wavelengths it has
# (necessary for creating the mock alignment factor file)
#
if 1:
    with open('dustkapscatmat_pyrmg70.inp','r') as f:
        for _ in range(7): f.readline()
        dustnf = int(f.readline())
        f.readline()
        dustfreq = np.zeros(dustnf)
        f.readline()
        for inu in range(dustnf):
            s=f.readline().split()
            dustfreq[inu] = float(s[0])
#
# Now make a mock alignment factor model. This is ONLY FOR TESTING.
# 
if 1:
    nrang = 20
    muang = np.linspace(1.e0,0.e0,nrang)
    eta   = np.arccos(muang)*180./math.pi
    orth  = np.zeros(nrang) + 1.e0
    ampl  = 0.5
    para  = ( 1.e0 - ampl*np.cos(muang*math.pi) ) / ( 1.e0 + ampl)
if 0:
#
# Now the alignment vector field. This is ONLY FOR TESTING.
#
    alvec = np.zeros((nx,ny,nz,3))
#alvec[:,:,:,2] = 1.0   # Vertical field config
    rrc            = np.sqrt(xx**2+yy**2)
    alvec[:,:,:,0] = yy/rrc  # Circular field config
    alvec[:,:,:,1] = -xx/rrc # Circular field config
#
# Normalize 
#
    length = np.sqrt(alvec[:,:,:,0]*alvec[:,:,:,0]+alvec[:,:,:,1]*alvec[:,:,:,1]+alvec[:,:,:,2]*alvec[:,:,:,2])
    alvec[:,:,:,0] = np.squeeze(alvec[:,:,:,0]) / ( length + 1e-60 )
    alvec[:,:,:,1] = np.squeeze(alvec[:,:,:,1]) / ( length + 1e-60 )
    alvec[:,:,:,2] = np.squeeze(alvec[:,:,:,2]) / ( length + 1e-60 )
#
# Now include the alignment efficiency. Simple model: zero near
# midplane, then quadratically up to 1 with z, where 1 is reached
# at a distance 'radius' from the midplane. This is ONLY FOR TESTING.
#
    epsal  = (zz/radius)*(zz/radius)
    epsal[epsal>1.0] = 1.0
    alvec[:,:,:,0]  = np.squeeze(alvec[:,:,:,0]) * epsal
    alvec[:,:,:,1]  = np.squeeze(alvec[:,:,:,1]) * epsal
    alvec[:,:,:,2]  = np.squeeze(alvec[:,:,:,2]) * epsal
#
# Write the wavelength file
#
if 1:
    with open('wavelength_micron.inp','w+') as f:
        f.write('%d\n'%(nlam))
        for value in lam:
            f.write('%13.6e\n'%(value))
#
#
# Write the stars.inp file
#
if 0:
    with open('stars.inp','w+') as f:
        f.write('2\n')
        f.write('1 %d\n\n'%(nlam))
        f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
        for value in lam:
            f.write('%13.6e\n'%(value))
        f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
if 1:
    import yt
    import numpy as np
    from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter, RadMC3DSource

    ds = yt.load("ca02_units_DD0490/data0490")


if 1:
    writer = RadMC3DWriter(ds)
    writer.write_amr_grid()


# Write the density file
#
if 1:
    dust_to_gas = 0.01
    def _DustDensity(field, data):
        return dust_to_gas * data["density"]
    ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3")
    writer.write_dust_file(("gas","dust_density"), "dust_density.inp")
    def _DumbTemp(field, data):
        return ds.arr(np.ones_like(data["density"]), 'K')*10
    ds.add_field("DumbTemp", function=_DumbTemp, units="K")
    writer.write_dust_file("DumbTemp", "dumb_temp.inp")
if 1:

    #radmc3d complains if the alignment vector is >1, even by 1e-16.
    #This keeps it from complaining.
    eps_val = 1e-16
    def _bx_hat(field, data):
        eps = data.ds.quan(eps_val,'gauss')
        output = (data['magnetic_field_x']/(data['magnetic_field_strength']+eps))
        return output
    ds.add_field("bx_hat", function=_bx_hat, units="dimensionless")
    def _by_hat(field, data):
        eps = data.ds.quan(eps_val,'gauss')
        output = data['magnetic_field_y']/(data['magnetic_field_strength']+eps)
        return output
    ds.add_field("by_hat", function=_by_hat, units="dimensionless")
    def _bz_hat(field, data):
        eps = data.ds.quan(eps_val,'gauss')
        output = data['magnetic_field_z']/(data['magnetic_field_strength']+eps)
        return output
    ds.add_field("bz_hat", function=_bz_hat, units="dimensionless")
    print('write.  Kill the single 1 on the third line.')
    writer.write_dust_file(["bx_hat","by_hat","bz_hat"], "align_mine_3.inp")
#
# Dust opacity control file
#
    with open('dustopac.inp','w+') as f:
        f.write('2               Format number of this file\n')
        f.write('1               Nr of dust species\n')
        f.write('============================================================================\n')
        f.write('20              Way in which this dust species is read\n')
        f.write('0               0=Thermal grain\n')
        f.write('pyrmg70         Extension of name of dustkappa_***.inp file\n')
        f.write('----------------------------------------------------------------------------\n')

#
# Dust alignment data
#
    with open('dustkapalignfact_pyrmg70.inp','w+') as f:
        f.write('1\n')
        f.write('%d\n'%(dustnf))
        f.write('%d\n\n'%(nrang))
        for value in dustfreq:
            f.write('%13.6e\n'%(value))
        f.write('\n')
        for value in eta:
            f.write('%13.6e\n'%(value))
        f.write('\n')
        for inu in range(dustnf):
            for imu in range(nrang):
                f.write('%13.6e %13.6e\n'%(orth[imu],para[imu]))
            f.write('\n')
#
# Dust alignment direction
#
#
# Write the radmc3d.inp control file
#
    nphot    = 1000000
    with open('radmc3d.inp','w+') as f:
        f.write('nphot = %d\n'%(nphot))
        f.write('scattering_mode_max = 4\n')
        f.write('alignment_mode = -1\n')
        f.write('iranfreqmode = 1\n')

