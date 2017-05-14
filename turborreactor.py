"""
Explicación de qué hace este .py aquí.
"""

from air_model import RealGas
from isentropic_gas import IsentropicGas
import isa
air = RealGas(cp_option='naca', gamma_option='standard')
gas = IsentropicGas(selected_cp_air_model='naca', selected_gamma_air_model='standard')


# Cálculo del difusor
def difusor(mach, p0, T0):
    # Difusor:
    T2t = gas.stagnation_temperature_from_mach(mach, T0)
    p2t = gas.stagnation_pressure_from_mach(mach, p0, T0)
    return T2t, p2t


# Función rendimiento tobera
def rend_tobera(x):
    return -1/3 * ((x-8000)/11000)**2 + 1


# Calculo de la tobera
def tobera(h, p5t, T5t, T2t):
    # Tobera
    p9 = isa.pres_isa(h)
    Ttob = (T5t + T2t) * 0.5
    T9 = (((p9/p5t)**((air.gamma_air(Ttob) - 1)/air.gamma_air(Ttob)) - 1) * rend_tobera(h) + 1) * T5t

    return T9, p9


def empuje(G, c, v9, v0):
    return (G + c)*v9 - G*v0

