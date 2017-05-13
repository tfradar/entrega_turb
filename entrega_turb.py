import isa
from air_model import RealGas
from isentropic_gas import IsentropicGas
import numpy as np

air = RealGas(cp_option='naca',gamma_option='standard')
gas = IsentropicGas(selected_cp_air_model='naca', selected_gamma_air_model= 'standard')


## Datos Cámara de combustión
rend_comb = 0.95 #rendimiento de la combustión
L=42e6 #J. poder calorífico del combustible
f = 0.02 #dosado = c/G
pi_comb = 0.95 #relación de presiones tras la cámara de combustión

## Función rendimiento de la turbina
def rend_turb(x):
    if x < 1000:
        return 0.88
    elif x > 2000:
        return 0.95
    else:
        return ((x - 1000) * 0.1 / 1000) + 0.88

## Función rendimiento tobera
def rend_tobera(x):
    return -1/3 * ((x-8000)/11000)**2 + 1


## Función rendimiento compresor
def rend_compresor(G0):
    if G0 < 25:
        rend_c = 0.6 + (G0 - 0) * 0.1/25
    elif (25<=G0<35):
        rend_c = 0.7 + (G0 - 25) * 0.1 /10
    elif (35<=G0<45):
        rend_c = 0.8 + (G0 - 35) * 0.05 / 10
    elif (45<=G0<55):
        rend_c = 0.85 + (G0 - 45) * 0.05 / 10
    elif (55<=G0<65):
        rend_c = 0.9
    elif (65<=G0<80):
        rend_c = 0.9 - (G0 - 65) * 0.05 / 15
    elif (80<=G0<95):
        rend_c = 0.85 - (G0 - 80) * 0.05 / 15
    elif (95<=G0):
        rend_c = 0.8 - (G0 - 95) * 0.1 / 15
    else:
        print('Error en el rendimiento del compresor')
        return 0
    return rend_c

##Función relación de compresión normal
def relacion_compresor_normal(G0):
    if G0 < 25:
        rel_c = 1 + (G0 - 0) * 9/25
    elif (25<=G0<35):
        rel_c = 10 + (G0 - 25) * 5 /10
    elif (35<=G0<45):
        rel_c = 15 + (G0 - 35) * 10 / 10
    elif (45<=G0<55):
        rel_c = 25 + (G0 - 45) * 15 / 10
    elif (55<=G0<65):
        rel_c = 40
    elif (65<=G0<80):
        rel_c = 40 + (G0 - 65) * 5 / 15
    elif (80<=G0<95):
        rel_c = 45 + (G0 - 80) * 5 / 15
    elif (95<=G0):
        rel_c = 50
    else:
        print('Error en la relación del compresor')
        return 0
    return rel_c

##Función rendimiento de la hélice
def rendimiento_helice(mach_i, alt_i):
    altitud = []
    rendimientos = np.zeros([12, 16])
    mach = [0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    for ii, alt in enumerate(np.linspace(0, 11000, 12)):
        x = np.linspace(0, 1, 16)
        altitud.append(alt)
        rendimientos[ii] = (-(x - 0.2) ** 2) * (10 - ii) / 10
        rendimientos[ii] = rendimientos[ii] + .95 ** 2 - 0.975 * ii / 11
        rendimientos[ii] = rendimientos[ii] * ((ii + 12) / 11)
    rendimientos[11] = rendimientos[11] * 0
    altitud = np.asarray(altitud)
    indice = np.where(altitud <= alt_i)
    return np.interp(mach_i, mach, rendimientos[indice[-1][-1]])

##Función generador de gas
def generador_gas(P0, T0, G0):
    # Compresor:
    P3t = P0 * relacion_compresor_normal(G0)
    T3t = T0 * (1 + 1 / rend_comb * (-1 + pi_comb ** ((air.gamma_air(T2t) - 1) / air.gamma_air(T2t))))
    # Cámara de combustión
    T4t = f * L * rend_comb / (1 + f) / air.cp_air(T3t) + T3t / (1 + f)
    c = f * G0
    Tcomp = 0.5 * (T3t + T2t)
    P4t = pi_comb * P3t

    # Turbina
    T5t = T4t - air.cp_air(Tcomp) / air.cp_air((T3t + T4t) / 2) * (T3t - T2t)
    Tturb = 0.5 * (T5t + T4t)
    W45 = G0 * air.cp_air(Tturb) * (T5t - T4t)
    W23 = G0 * air.cp_air(Tcomp) * (T3t - T2t)
    P5t = ((T5t / T4t - 1) / rend_turb(T4t) + 1) ** (air.gamma_air(Tturb) / (air.gamma_air(Tturb) - 1)) * P4t

##Función Turborreactor
def turborreactor(m, T0, G0, p0):
    # Difusor:
    T2t = gas.stagnation_temperature_from_mach(m, T0)
    P2t = gas.stagnation_pressure_from_mach(m, p0, T0)

    # Compresor:
    P3t = P2t * relacion_compresor_normal(G0)
    T3t = T2t * (1 + 1 / rend_comb * (-1 + pi_comb ** ((air.gamma_air(T2t) - 1) / air.gamma_air(T2t))))
    # Cámara de combustión
    T4t = f * L * rend_comb / (1 + f) / air.cp_air(T3t) + T3t / (1 + f)
    c = f * G0
    Tcomp = 0.5 * (T3t + T2t)
    P4t = pi_comb * P3t

    # Turbina
    T5t = T4t - air.cp_air(Tcomp) / air.cp_air((T3t + T4t) / 2) * (T3t - T2t)
    Tturb = 0.5 * (T5t + T4t)
    W45 = G0 * air.cp_air(Tturb) * (T5t - T4t)
    W23 = G0 * air.cp_air(Tcomp) * (T3t - T2t)
    P5t = ((T5t / T4t - 1) / rend_turb(T4t) + 1) ** (air.gamma_air(Tturb) / (air.gamma_air(Tturb) - 1)) * P4t

    # Tobera
    P9 = isa.pres_isa(h)
    T9 = T5t * (
    rend_tobera(h) * ((P9 / P5t) ** (air.gamma_air((T5t + T2t) / 2) - 1) / air.gamma_air((T5t + T2t) / 2) - 1) + 1)

##Función Turbohélice
def turbohelice(m, G0, T0, p0):
    # Difusor:
    T2t = gas.stagnation_temperature_from_mach(m, T0)
    P2t = gas.stagnation_pressure_from_mach(m, p0, T0)

    # Compresor:
    P3t = P2t * relacion_compresor_normal(G0)
    T3t = T2t * (1 + 1 / rend_comb * (-1 + pi_comb ** ((air.gamma_air(T2t) - 1) / air.gamma_air(T2t))))

    # Cámara de combustión
    T4t = f * L * rend_comb / (1 + f) / air.cp_air(T3t) + T3t / (1 + f)
    c = f * G0
    Tcomp = 0.5 * (T3t + T2t)
    P4t = pi_comb * P3t

    # Turbina 1
    T5t = T4t - air.cp_air(Tcomp) / air.cp_air((T3t + T4t) / 2) * (T3t - T2t)
    Tturb = 0.5 * (T5t + T4t)
    W45 = G0 * air.cp_air(Tturb) * (T5t - T4t)
    W23 = G0 * air.cp_air(Tcomp) * (T3t - T2t)
    P5t = ((T5t / T4t - 1) / rend_turb(T4t) + 1) ** (air.gamma_air(Tturb) / (air.gamma_air(Tturb) - 1)) * P4t
    
    # Turbina 2
    T6t = T5t - air.cp_air(Tcomp) / air.cp_air((T4t + T5t) / 2) * (T4t - T3t)
    Tturb2 = 0.5 * (T6t + T5t)
    W45 = G0 * air.cp_air(Tturb2) * (T6t - T5t)
    W23 = G0 * air.cp_air(Tcomp) * (T4t - T3t)
    P6t = ((T6t / T5t - 1) / 0.95 + 1) ** (air.gamma_air(Tturb2) / (air.gamma_air(Tturb2) - 1)) * P5t

    
    
height = np.array([500, 2000, 4000, 5000, 7000, 9000, 11000])
mach = np.array([0.3, .45, .6, .8, .85, .9, .95])

T0 = np.zeros(height.size * mach.size)
p0 = np.zeros_like(T0)


area = np.pi * 0.4**2
ii = 0
for h in height:

    for m in mach:
        T0[ii] = isa.temp_isa(h)
        p0[ii] = isa.pres_isa(h)
        rho0 = p0[ii] / (T0[ii] * air.R_air)
        v0 = m * gas.sound_speed(T0[ii])
        G0 = rho0 * area * v0

        ii += 1
