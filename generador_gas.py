from air_model import RealGas
from isentropic_gas import IsentropicGas

air = RealGas(cp_option='naca', gamma_option='standard')
gas = IsentropicGas(selected_cp_air_model='naca', selected_gamma_air_model='standard')


# Datos Cámara de combustión
rend_comb = 0.95  # rendimiento de la combustión
heating_value = 42e6  # J. poder calorífico del combustible
f = 0.02  # dosado = c/G
pi_comb = 0.95  # relación de presiones tras la cámara de combustión


# Función rendimiento de la turbina
def rend_turb(x):
    if x < 1000:
        return 0.88
    elif x > 2000:
        return 0.95
    else:
        return ((x - 1000) * 0.1 / 1000) + 0.88


# Función rendimiento compresor
def rend_compresor(G0):
    if G0 < 25:
        rend_c = 0.6 + (G0 - 0) * 0.1/25
    elif 25 <= G0 < 35:
        rend_c = 0.7 + (G0 - 25) * 0.1 /10
    elif 35 <= G0 < 45:
        rend_c = 0.8 + (G0 - 35) * 0.05 / 10
    elif 45 <= G0 < 55:
        rend_c = 0.85 + (G0 - 45) * 0.05 / 10
    elif 55 <= G0 < 65:
        rend_c = 0.9
    elif 65 <= G0 < 80:
        rend_c = 0.9 - (G0 - 65) * 0.05 / 15
    elif 80 <= G0 < 95:
        rend_c = 0.85 - (G0 - 80) * 0.05 / 15
    elif 95 <= G0:
        rend_c = 0.8 - (G0 - 95) * 0.1 / 15
    else:
        print('Error en el rendimiento del compresor')
        return 0
    return rend_c


# Función relación de compresión normal
def relacion_compresor_normal(G0):
    if G0 < 25:
        rel_c = 1 + (G0 - 0) * 9/25
    elif 25 <= G0 < 35:
        rel_c = 10 + (G0 - 25) * 5 /10
    elif 35 <= G0 < 45:
        rel_c = 15 + (G0 - 35) * 10 / 10
    elif 45 <= G0 < 55:
        rel_c = 25 + (G0 - 45) * 15 / 10
    elif 55 <= G0 < 65:
        rel_c = 40
    elif 65 <= G0 < 80:
        rel_c = 40 + (G0 - 65) * 5 / 15
    elif 80 <= G0 < 95:
        rel_c = 45 + (G0 - 80) * 5 / 15
    elif 95 <= G0:
        rel_c = 50
    else:
        print('Error en la relación del compresor')
        return 0
    return rel_c


# Función generador de gas
def generador_gas(T2t, P2t, G0):
    # Compresor:
    P3t = P2t * relacion_compresor_normal(G0)
    T3t = T2t * (1 + 1 / rend_comb * (-1 + pi_comb ** ((air.gamma_air(T2t) - 1) / air.gamma_air(T2t))))
    # Cámara de combustión
    T4t = f * heating_value * rend_comb / (1 + f) / air.cp_air(T3t) + T3t / (1 + f)
    Tcomp = 0.5 * (T3t + T2t)
    P4t = pi_comb * P3t
    c = f * G0

    # Turbina
    T5t = T4t - air.cp_air(Tcomp) / air.cp_air((T3t + T4t) / 2) * (T3t - T2t)
    Tturb = 0.5 * (T5t + T4t)
    W45 = G0 * air.cp_air(Tturb) * (T5t - T4t)
    W23 = G0 * air.cp_air(Tcomp) * (T3t - T2t)
    P5t = ((T5t / T4t - 1) / rend_turb(T4t) + 1) ** (air.gamma_air(Tturb) /
                                                     (air.gamma_air(Tturb) - 1)) * P4t
    return T5t, P5t, W45, W23, c