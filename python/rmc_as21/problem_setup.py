#
# Import NumPy for array handling
#
import numpy as np
import math


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
        "We don't need any stars, we'll set the temperature by hand"
        pass

# Write the grid file
#
if 1:
    import yt
    import numpy as np
    from yt.extensions.astro_analysis.radmc3d_export.api import RadMC3DWriter, RadMC3DSource

    ds = yt.load("/data/cb1/Projects/P49_EE_BB/as21_M0.6_MA0.3_64/DD1080/data1080")


if 1:
    writer = RadMC3DWriter(ds)
    writer.write_amr_grid()


# Write the density file
#
if 1:
    dust_to_gas = 0.01
    def _DustDensity(field, data):
        return dust_to_gas * data["density"]
    ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3",sampling_type='cell')
    writer.write_dust_file(("gas","dust_density"), "dust_density.inp")
    def _DumbTemp(field, data):
        return ds.arr(np.ones_like(data["density"]), 'K')*10
    ds.add_field("DumbTemp", function=_DumbTemp, units="K", sampling_type='cell')
    writer.write_dust_file("DumbTemp", "dust_temperature.dat")

if 1:

    #radmc3d complains if the alignment vector is >1, even by 1e-16.
    #eps_val keeps it from complaining.
    eps_val = 1e-16
    def _bx_hat(field, data):
        eps = data.ds.quan(eps_val,'gauss')
        output = (data['magnetic_field_x']/(data['magnetic_field_strength']+eps))
        return output
    ds.add_field("bx_hat", function=_bx_hat, units="dimensionless",sampling_type='cell')
    def _by_hat(field, data):
        eps = data.ds.quan(eps_val,'gauss')
        output = data['magnetic_field_y']/(data['magnetic_field_strength']+eps)
        return output
    ds.add_field("by_hat", function=_by_hat, units="dimensionless",sampling_type='cell')
    def _bz_hat(field, data):
        eps = data.ds.quan(eps_val,'gauss')
        output = data['magnetic_field_z']/(data['magnetic_field_strength']+eps)
        return output
    ds.add_field("bz_hat", function=_bz_hat, units="dimensionless",sampling_type='cell')
    print('Writing the alignment file.  This is not perfect, please delete the third line from the file.')
    writer.write_dust_file(["bx_hat","by_hat","bz_hat"], "grainalign_dir.inp")
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

