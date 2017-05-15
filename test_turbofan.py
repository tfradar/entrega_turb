"""
# Se importa un generador de gas (compresor, cámara de combustión, turbina)
# Se le añade un una turbina y un fan de turbofan.py
# Saca los rendimientos, impulsos y empujes para distintas condiciones de vuelo
"""

import isa
import numpy as np
from air_model import RealGas
from isentropic_gas import IsentropicGas
import generador_gas
import turbofan
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
E_primario_lista = []
E_secundario_lista = []
Eneto_lista = []
Ie_lista = []
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

        # Difusor
        T2t, p2t = turbofan.difusor(m, p0[ii], T0[ii])

        # Generador de Gas
        T5t, p5t, c, p4t = generador_gas.generador_gas(T2t, p2t, G0[ii])

        # Turbofan
        v9, v19, G_secundario = turbofan.tubofan(G0[ii], v0[ii], rho0, T0[ii], p0[ii], T2t, T5t, p2t, p4t)

        # Actuaciones
        Eneto, Ie, Ce, E_primario, E_secundario = turbofan.actuaciones(G0[ii], v9, v19, v0[ii], G_secundario, c)

        # Rendimientos
        eta_m, eta_p, eta_mp = turbofan.rendimiento_TB(c, v19, v9, v0[ii], G0[ii], E_primario, E_secundario, G_secundario)

        # Reunimos los resultados en listas
        E_primario_lista.append(E_primario)
        E_secundario_lista.append(E_secundario)
        Eneto_lista.append(Eneto)
        Ie_lista.append(Ie)
        Ce_lista.append(Ce)
        eta_m_lista.append(eta_m)
        eta_p_lista.append(eta_p)
        eta_mp_lista.append(eta_mp)

        ii += 1

#TFD: Me falta ordenar y plotear las listas
