"""
Explicación de qué hace este .py aquí.
"""

import isa
import numpy as np
from air_model import RealGas
from isentropic_gas import IsentropicGas
import generador_gas
import turbohelice
import turborreactor

# Carga de funciones:
air = RealGas(cp_option='naca', gamma_option='standard')
gas = IsentropicGas(selected_cp_air_model='naca', selected_gamma_air_model='standard')

# Creación de los vectores de altitud y mach:
height = np.array([500, 2000, 4000, 5000, 7000, 9000, 11000])
mach = np.array([0.3, .45, .6, .8, .85, .9, .95])

# Creación de las variables resultado
T0 = np.zeros(height.size * mach.size)
p0 = np.zeros_like(T0)

v0 = np.zeros_like(T0)
G0 = np.zeros_like(T0)

E = np.zeros_like(T0)
Eneto = np.zeros_like(T0)
Ie = np.zeros_like(T0)
Ce = np.zeros_like(T0)
eta_m = np.zeros_like(T0)
eta_p = np.zeros_like(T0)
eta_mp = np.zeros_like(T0)


area = np.pi * 0.4 ** 2

ii = 0
for h in height:
    for m in mach:
        # Condiciones ambiente:
        T0[ii] = isa.temp_isa(h)
        p0[ii] = isa.pres_isa(h)
        rho0 = p0[ii] / (T0[ii] * air.R_air)
        v0[ii] = m * gas.sound_speed(T0[ii])
        G0[ii] = rho0 * area * v0[ii]

        # Difusor:
        T2t, p2t = turbohelice.difusor(m, p0[ii], T0[ii])

        # Generador de Gas:
        T45t, p45t, c, p4t = generador_gas.generador_gas(T2t, p2t, G0[ii])

        # Turbina turbohélice
        T5_, T5_t, p5_t, W56, v8 = turbohelice.turbina(T45t, G0[ii], p45t, p0[ii], T2t)

        # Actuaciones:
        # Empuje, Impulso específico, Consumo específico
        E, Ie, Ce, Eneto = turbohelice.actuaciones(G0[ii], c, v8, v0[ii])
        # Rendimiento
        eta_m, eta_p, eta_mp = turbohelice.rendimiento_TB(Eneto, c, v8, v0[ii], G0[ii])

        print(ii + 1, E)

        ii += 1


