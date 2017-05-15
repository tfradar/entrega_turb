"""
Archivo con las partes del turbohélice

Difusor, 
"""

from air_model import RealGas
from isentropic_gas import IsentropicGas
import isa
import numpy as np
import generador_gas

air = RealGas(cp_option='naca', gamma_option='standard')
gas = IsentropicGas(selected_cp_air_model='naca', selected_gamma_air_model='standard')

# Función rendimiento de la turbina
def rend_turb(x):
    if x < 1000:
        return 0.88
    elif x > 2000:
        return 0.95
    else:
        return ((x - 1000) * 0.1 / 1000) + 0.88

# Cálculo del difusor
def difusor(mach, p0, T0):
    # Difusor:
    T2t = gas.stagnation_temperature_from_mach(mach, T0)
    p2t = gas.stagnation_pressure_from_mach(mach, p0, T0)
    return T2t, p2t

# Función rendimiento de la hélice
def rendimiento_helice(mach_i, alt_i):
    altitud = []
    rendimientos = np.zeros([12,16])
    mach = [0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    for ii, alt in enumerate(np.linspace(0, 11000, 12)):
        x = np.linspace(0, 1, 16)
        altitud.append(alt)
        rendimientos[ii] = (- (x - 0.2) **2 ) * (10-ii) / 10
        rendimientos[ii] = rendimientos[ii] + .95 ** 2 - 0.975 *ii / 11
        rendimientos[ii] = rendimientos[ii] * ((ii + 12) / 11)
    rendimientos[11] = rendimientos[11] * 0
    altitud = np.asarray(altitud)
    indice = np.where(altitud <= alt_i)
    return np.interp(mach_i, mach, rendimientos[indice[-1][-1]])

# Función Hélice
def helice(v0, G0, v6, W56, m, h):
    Eneto = G0 * (v6 - v0)
    Phelice = Eneto * v0
    rend_mec = 0,95
    rend_helice = rendimiento_helice(m, h)
    Wprop = Phelice/ rend_helice /rend_mec

    return Wprop

# Función turbina
def turbina(T5t, G0, p5t, p0, Wprop):
    # Turbina
    p45t = p5t
    T45t = T5t
    T5_t = T45t - Wprop / G0 / air.cp_air(T5t)
    Tturb = 0.5 * (T5_t + T45t)
    W56 = G0 * air.cp_air(Tturb) * (T5_t - T45t)
    p5_t = ((T5_t / T5t - 1) / rend_turb(T5t) + 1) ** (air.gamma_air(Tturb) /
                                                     (air.gamma_air(Tturb) - 1)) * p45t
    p6_ = p0
    T5_ = (p6_ / p45t) ** ((air.gamma_air(Tturb) - 1) / air.gamma_air(Tturb)) * T5_t
    v8 = (np.sqrt(2 * air.cp_air(Tturb) * (T5t - T5_)))

    return T5_, T5_t, p5_t, W56, v8

def actuaciones(G0, c, v8, v0):
    #Empuje
    E = (G0 + c) * v8 - G0 * v0
    #Empuje neto
    Eneto = G0 * (v8 - v0)
    # Impulso específico total (sin aproximar la presión de salida)
    Ie = E / G0
    # Consumo específico total (sin aproximar la presión de salida)
    Ce = c / E

    return E, Eneto, Ie, Ce,

# Rendimiento motor, propulsor y motopropulsor:
def rendimiento_TB(Eneto, c, v8, v0, G0):
    # Rendimiento motor:
    eta_m = (Eneto * v0 + 0.5 * (G0 + c) * (v8 - v0) ** 2 - 0.5 * c * v0 ** 2) / (c * generador_gas.heating_value)

    # Rendimiento propulsor:
    eta_p = Eneto * v0 / (Eneto * v0 + 0.5 * (G0 + c) * (v8 - v0) ** 2 - 0.5 * c * v0 ** 2)

    # Rendimiento motopropulsor:
    eta_mp = eta_m * eta_p

    return eta_m, eta_p, eta_mp
