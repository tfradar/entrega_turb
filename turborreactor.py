"""
Explicación de qué hace este .py aquí.
"""

from air_model import RealGas
from isentropic_gas import IsentropicGas
import isa
air = RealGas(cp_option='naca', gamma_option='standard')
gas = IsentropicGas(selected_cp_air_model='naca', selected_gamma_air_model='standard')
import generador_gas


# Cálculo del difusor
def difusor(mach, p0, T0):
    # Difusor:
    T2t = gas.stagnation_temperature_from_mach(mach, T0)
    p2t = gas.stagnation_pressure_from_mach(mach, p0, T0)
    return T2t, p2t


# Función rendimiento tobera
def rend_tobera(x):

    return - 1/3 * ((x - 8000) / 11000) ** 2 + 1


# Calculo de la tobera
def tobera(h, p5t, T5t, T2t):
    # Tobera
    p9 = isa.pres_isa(h)
    Ttob = (T5t + T2t) * 0.5
    T9 = (((p9 / p5t) ** ((air.gamma_air(Ttob) - 1) / air.gamma_air(Ttob)) - 1) * rend_tobera(h) + 1) * T5t

    return T9, p9


def actuaciones(G0, c, v9, v0):
    #Empuje
    E = (G0 + c) * v9 - G0 * v0
    #Empuje neto
    Eneto = G0 * (v9 - v0)
    # Impulso específico total (sin aproximar la presión de salida)
    Ie = E / G0
    # Consumo específico total (sin aproximar la presión de salida)
    Ce = c / E

    return E, Ie, Ce, Eneto,

# Rendimiento motor, propulsor y motopropulsor:
def rendimiento_TB(Eneto, G, c, v9, v0, ):
    # Rendimiento motor:
    eta_m = (Eneto * v0 + 0.5 * (G + c) * (v9 - v0) ** 2 - 0.5 * c * v0 ** 2) / (c * generador_gas.heating_value)

    # Rendimiento propulsor:
    eta_p = Eneto * v0 / (Eneto * v0 + 0.5 * (G + c) * (v9 - v0) ** 2 - 0.5 * c * v0 ** 2)

    # Rendimiento motopropulsor:
    eta_mp = eta_m * eta_p

    return eta_m, eta_p, eta_mp