"""
# Se importa un generador de gas (compresor, cámara de combustión, turbina)
# Se le añade otra turbina acoplada a una hélice propulsiva de turbohelice.py
# Saca los rendimientos, impulsos y empujes para distintas condiciones de vuelo
"""

import isa
import numpy as np
from air_model import RealGas
from isentropic_gas import IsentropicGas
import generador_gas
import turbohelice
from matplotlib import pyplot as plt

# Carga de funciones:
air = RealGas(cp_option='naca', gamma_option='standard')
gas = IsentropicGas(selected_cp_air_model='naca', selected_gamma_air_model='standard')

# Funciones matemáticas
area = np.pi * 0.4 ** 2

# Creación de los vectores de altitud y mach:
height = np.array([500, 2000, 4000, 5000, 7000, 9000, 11000])
mach = np.array([0.3, .45, .6, .8, .85, .9, .95])

# Creación de las variables resultado
T0 = np.zeros(height.size * mach.size)
p0 = np.zeros_like(T0)
v0 = np.zeros_like(T0)
G0 = np.zeros_like(T0)

# Creación de las listas resultado
Tp_lista = []
Ce_lista = []
eta_m_lista = []
eta_p_lista = []
eta_mp_lista = []

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
        T5_t, p5_t, W56 = turbohelice.turbina(T45t, G0[ii], p45t, p0[ii], T2t)

        # Hélice
        Phelice, T5t, Pot_H_u = turbohelice.helice(T0[ii], v0[ii], G0[ii], m, h, c, T45t)

        # Actuaciones:
        # Empuje, Impulso específico, Consumo específico
        Tp, Ce = turbohelice.actuaciones(c, v0[ii], Phelice)
        # Rendimiento
        eta_m, eta_p, eta_mp = turbohelice.rendimiento_TB(c, Pot_H_u, Tp)

        # Reunimos los resultados en listas
        Tp_lista.append(Tp)
        Ce_lista.append(Ce)
        eta_m_lista.append(eta_m)
        eta_p_lista.append(eta_p)
        eta_mp_lista.append(eta_mp)

        ii += 1
print(eta_mp_lista)
#TFD: Me falta ordenar y plotear las listas
