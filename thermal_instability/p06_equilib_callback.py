
from yt.utilities import physical_constants as units
import numpy as np
import pdb
nar=np.array
#if 'fa01' not in dir()
#reload(taxi)
#fa03=taxi.taxi('fa03')
#fa05=taxi.taxi('fa05')
#car = fa05

mp = 1.67262171e-24  #g, from phys_constants.h
kb = 1.3806504e-16 #erg/K from phys_constants.h 
def PfromTempDen(Temp,Den):
    return Den*Temp*kb/mp
def TfromPempDen(P,Den):
    return P/(Den*kb/mp)
import equilibriumdata
class equillibrium_callback():
    def __init__(self, temps=[], mu=1.0):
        self.temps=temps
        self.MeanAtomicMass=mu
        pass
    def __call__(self,pp):
        x_field = pp.profile.x_field
        y_field = pp.profile.y_field
        verb=None
        if x_field[1] in ['density','davesNumberDensity','number_density'] and \
           y_field[1] in ['GasPressure', 'pressure']:
               verb = self.density_pressure
        elif x_field[1] in ['density','davesNumberDensity','number_density'] and\
             y_field[1] in ['Temperature']:
               verb = self.density_temperature
        else:
            print("These fields not supported: doing nothing:", x_field, y_field)

        if verb is not None:
            verb(pp)

    def density_temperature(self,pp):
        x_field = pp.profile.x_field
        y_field = pp.profile.y_field
        Equilib_NumberDensities = nar(equilibriumdata.NumberDensities)
        Equilib_GasPressure = nar(equilibriumdata.GasPressures)
        Equilib_Temperature = Equilib_GasPressure /(Equilib_NumberDensities*kb)
        Densities = None

        if 'davesNumberDensity' in x_field:
            Densities = Equilib_NumberDensities
        elif 'density' in x_field:
            Densities = Equilib_NumberDensities *self.MeanAtomicMass*units.mp

        key = list(pp.plots.keys())[0]
        this_axes = pp.plots[key].axes
        pp.save('temp.png')
        xlim = this_axes.get_xlim()
        ylim = this_axes.get_ylim()
        this_axes.plot(Densities,Equilib_Temperature,c='b')
        this_axes.set_xlim(xlim)
        this_axes.set_ylim(ylim)

    def density_pressure(self,pp):
        x_field = pp.profile.x_field
        y_field = pp.profile.y_field
        Equilib_NumberDensities = nar(equilibriumdata.NumberDensities)
        Equilib_GasPressure = nar(equilibriumdata.GasPressures)

        if 'davesNumberDensity' in x_field:
            Densities = Equilib_NumberDensities
        elif 'density' in x_field:
            Densities = Equilib_NumberDensities *self.MeanAtomicMass*units.mp

        if  y_field[1] in ['GasPressure', 'pressure']:
            Pressures = Equilib_GasPressure
        elif 'PoverK' in y_field:
            Pressures = Equilib_GasPressure/units.kboltz
        key = list(pp.plots.keys())[0]
        this_axes = pp.plots[key].axes
        pp.save('temp.png')
        this_axes.plot(Densities,Pressures,c='b')
        for T in self.temps:
            this_axes.plot(Densities,PfromTempDen(T,Densities),c='k')

