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
    rendimientos = np.zeros([12, 16])
    mach = [0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    for ii, alt in enumerate(np.linspace(0, 11000, 12)):
        x = np.linspace(0, 1, 16)
        altitud.append(alt)
        rendimientos[ii] = (- (x - 0.2) ** 2) * (10 - ii) / 10
        rendimientos[ii] = rendimientos[ii] + .95 ** 2 - 0.975 * ii / 11
        rendimientos[ii] = rendimientos[ii] * ((ii + 12) / 11)
    rendimientos[11] = rendimientos[11] * 0
    altitud = np.asarray(altitud)
    indice = np.where(altitud <= alt_i)
    return np.interp(mach_i, mach, rendimientos[indice[-1][-1]])

def tubofan(G0, v0, rho0, T0, P0, T2t, T5t, p2t, p4t):
    G_secundario = 5 * G0
    Gtotal = G_secundario + G0
    Afan = G0 / (rho0 * v0)
    dfan = np.sqrt(Afan/np.pi*4)
    T45t = T5t
    T5_t = (T2t + T45t)*0.5
    Tturb2 = 0.5*(T45t + T5t)
    p5_t = p4t * (1 + 1 / rend_turb(Tturb2) * (T5_t / T45t - 1)) ** (air.gamma_air(Tturb2) / (air.gamma_air(Tturb2) - 1))
    T13t = 1/5 * air.cp_air(Tturb2) * (T45t - T5_t) /air.cp_air(T2t) + T2t
    Tfan = 0.5*(T2t+T13t)
    P13t = p2t * (0.8 * (T13t / T2t - 1) + 1) ** (air.gamma_air(Tfan) / (air.gamma_air(Tfan) - 1))
    Ttob = (T5t + T0) / 2
    Ttob2 = (T13t + T0) / 2
    v19 = np.sqrt(2 * air.cp_air(Ttob2)*T13t*(1-(P0/P13t) ** ((air.gamma_air(Ttob2) - 1) / (air.gamma_air(Ttob2)))))
    v9 = np.sqrt(2*air.cp_air(Ttob)*T5_t*(1-(P0/p5_t)**((air.gamma_air(Ttob)-1) / (air.gamma_air(Ttob)))))

    return v9, v19, G_secundario

def actuaciones(G0, v9, v19, v0, G_secundario, c):
    E_primario = G0 * (v9 - v0)
    E_secundario = G_secundario * (v19 - v0)
    Eneto = E_primario + E_secundario
    Ie= Eneto / G0
    Ce = c / Eneto

    return Eneto, Ie, Ce, E_primario, E_secundario

def rendimiento_TB(c, v19, v9, v0, G0, E_primario, E_secundario, G_secundario):
    eta_m = (G0 * 0.5 * (v9 ** 2 - v0 ** 2) + G_secundario * 0.5 * (v19 ** 2 - v0 ** 2)) / c / generador_gas.heating_value
    eta_p = (E_primario + E_secundario)*v0 / (G0 * 0.5 * (v9 ** 2 - v0 ** 2) + G_secundario * 0.5 * (v19 ** 2 - v0 ** 2))
    eta_mp = eta_m * eta_p

    return eta_m, eta_p, eta_mp

