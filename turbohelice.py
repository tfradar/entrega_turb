"""
Explicación de qué hace este .py aquí.
"""

from air_model import RealGas
from isentropic_gas import IsentropicGas
import isa
import numpy as np

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

# Función turbina
def turbina(T5t, G0, p5t, p0, T2t):
    # Turbina
    p45t = p5t
    T45t = T5t
    T5_t = (T2t + T45t) / 2
    Tturb = 0.5 * (T5_t + T45t)
    W56 = G0 * air.cp_air(Tturb) * (T5_t - T45t)
    p5_t = ((T5_t / T5t - 1) / rend_turb(T5t) + 1) ** (air.gamma_air(Tturb) /
                                                     (air.gamma_air(Tturb) - 1)) * p45t
    p6_ = p0
    T5_ = (p6_ / p45t) ** ((air.gamma_air(Tturb) - 1) / air.gamma_air(Tturb)) * T5_t
    v6 = (np.sqrt(2 * air.cp_air(Tturb) * (T5t - T5_)))

    return T5_t, p45t, W56, p5_t, T5_, v6

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
    Potencia_TJ = Eneto * v0
    Phelice = Potencia_TJ
    rend_mec = Phelice / W56
    rend_helice = rendimiento_helice(m, h)
    Wprop = Phelice/ rend_helice /rend_mec

    return Wprop
