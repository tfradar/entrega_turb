import isa
from air_model import RealGas
from isentropic_gas import IsentropicGas
import numpy as np

air = RealGas(cp_option='naca',gamma_option='standard')
gas = IsentropicGas(selected_cp_air_model='naca', selected_gamma_air_model= 'standard')

altitud = np.array([500, 2000, 4000, 5000, 7000, 9000, 11000])
mach = np.array([0.3, .45, .6, .8, .85, .9, .95])

T0 = np.zeros(altitud.size * mach.size)
p0 = np.zeros_like(T0)

area = np.pi * 0.4**2
ii = 0
for h in altitud:

    for m in mach:
        T0[ii] = isa.temp_isa(h)
        p0[ii] = isa.pres_isa(h)
        rho0 = p0[ii] / (T0[ii] * air.R_air)
        v0 = m * gas.sound_speed(T0[ii])
        G0 = rho0 * area * v0
        ii += 1
        print(G0)
