# Urrios 2016: multicellular memory

import numpy as np

def not_cell(state, params):
    L_X, x, y, N_X, N_Y = state
    delta_L, gamma_X, n_y, theta_X, eta_x, omega_x, m_x, delta_x, rho_x = params

    f = gamma_X * (y ** n_y)/(1 + (theta_X*y)**n_y )
    dL_X_dt = f - delta_L * L_X

    dx_dt = N_X * (eta_x * (1/(1+ (omega_x*L_X)**m_x))) - N_Y * (delta_x * x) - rho_x * x

    return dL_X_dt, dx_dt


def not_cell_a(state, params):
    return not_cell(state, params)


"""
    L_A, a, b, N_A, N_B = state
    delta_L, gamma_A, n_b, theta_A, eta_a, omega_a, m_a, delta_a, rho_a= params

    f_b = gamma_A * (b ** n_b)/(1 + (theta_A*b)**n_b )
    dL_A_dt = f_b - delta_L * L_A

    da_dt = N_A * (eta_a * (1/(1+ (omega_a*L_A)**m_b))) - N_B * (delta_a * a) - rho_a * a
"""

def not_cell_b(state, params):
    return not_cell(state, params)

def population(state, params):
    N = state
    r = params

    dN = r * N * (1 - N)

    return dN


def toggle_model(state, T, params):
    L_A, L_B, a, b, N_A, N_B = state

    state_A = L_A, a, b, N_A, N_B
    state_B = L_B, b, a, N_B, N_A
    
    delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, r_A, r_B = params

    params_A = delta_L, gamma_A, n_b, theta_A, eta_a, omega_a, m_a, delta_a, rho_a
    params_B = delta_L, gamma_B, n_a, theta_B, eta_b, omega_b, m_b, delta_b, rho_b

    dL_A_dt, da_dt = not_cell(state_A, params_A)
    dL_B_dt, db_dt = not_cell(state_B, params_B)

    dN_A_dt = population(N_A, r_A)
    dN_B_dt = population(N_B, r_B)
        
    return np.array([dL_A_dt, dL_B_dt, da_dt, db_dt, dN_A_dt, dN_B_dt])

def toggle_model_ODE(T, state, params):
    return toggle_model(state, T, params)