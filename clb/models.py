# Urrios 2016: multicellular memory

import numpy as np

def not_cell(state, params):
    L_X, x, y, N_X, N_Y = state
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x = params

    f = gamma_L_X * (y ** n_y)/(1 + (theta_L_X*y)**n_y )
    dL_X_dt = f - delta_L * L_X

    dx_dt = N_X * (eta_x * (1/(1+ (omega_x*L_X)**m_x))) - N_Y * (delta_x * x) - rho_x * x

    return dL_X_dt, dx_dt

    
def yes_cell(state, params):
    x, y, N_X, N_Y = state
    gamma_x, n_y, theta_x, delta_x, rho_x = params

    dx_dt = N_X * gamma_x * (y ** n_y)/(1 + (theta_x*y)**n_y ) - N_Y * (delta_x * x) - rho_x * x
    
    return dx_dt

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


# L_A ... intermediate
# a ... out
# b ... in
# N_A ... number of cells
def not_cell_wrapper(state, params):
    L_A, a, b, N_A = state
    
    state_A = L_A, a, b, N_A, N_A
    params_A = params

    return not_cell(state_A, params_A)

# a ... out
# b ... in
# N_A ... number of cells
def yes_cell_wrapper(state, params):
    a, b, N_A = state

    state_A = a, b, N_A, N_A
    params_A = params

    return yes_cell(state_A, params_A)

def not_model(state, T, params):
    L_A, a, b, N_A = state

    delta_L, gamma_L_A, n_b, theta_L_A, eta_a, omega_a, m_a, delta_a, delta_b, rho_a, rho_b, r_A = params

    state_not = L_A, a, b, N_A
    params_not = delta_L, gamma_L_A, n_b, theta_L_A, eta_a, omega_a, m_a, delta_a, rho_a
    dL_A_dt, da_dt = not_cell_wrapper(state_not, params_not)
    
    db_dt = 0#- N_A * delta_b * b - rho_b * b

    dN_A_dt = population(N_A, r_A)

    return np.array([dL_A_dt, da_dt, db_dt, dN_A_dt])

def yes_model(state, T, params):
    a, b, N_A = state
    
    gamma_a, n_b, theta_a, delta_a, delta_b, rho_a, rho_b, r_A = params

    state_yes = a, b, N_A
    params_yes = gamma_a, n_b, theta_a, delta_a, rho_a
    da_dt = yes_cell_wrapper(state_yes, params_yes)
    
    db_dt = 0 #- N_A * delta_b * b - rho_b * b

    dN_A_dt = population(N_A, r_A)

    return np.array([da_dt, db_dt, dN_A_dt])


def toggle_model_ODE(T, state, params):
    return toggle_model(state, T, params)

def not_model_ODE(T, state, params):
    return not_model(state, T, params)


def yes_model_ODE(T, state, params):
    return yes_model(state, T, params)

def MUX_4_1_model(state, T, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x


    S0, S1, I0, I1, I2, I3 = state[:6]
    I0_out, I1_out, I2_out, I3_out = state[6:10]
    L_I0_I0, L_I1_S0, L_I1_I1, L_I2_S1, L_I2_I2, L_I3_S0, L_I3_S1, L_I3_I3, L_I0, L_I1, L_I2, L_I3 = state[10:22]
    N_I0_S0, N_I0_S1, N_I0_I0, N_I1_S0, N_I1_S1, N_I1_I1, N_I2_S0, N_I2_S1, N_I2_I2, N_I3_S0, N_I3_S1, N_I3_I3, N_I0, N_I1, N_I2, N_I3 = state[22:38]
    out = state[38]
    
        
    """
     I0
    """
    dI0_out = 0
    
    # yes S0: I0_S0
    state_yes_I0_S0 = I0_out, S0, N_I0_S0
    dI0_out += yes_cell_wrapper(state_yes_I0_S0, params_yes)
    dN_I0_S0 = population(N_I0_S0, r_X)

    # yes S1: I0_S1
    state_yes_I0_S1 = I0_out, S1, N_I0_S1
    dI0_out += yes_cell_wrapper(state_yes_I0_S1, params_yes)
    dN_I0_S1 = population(N_I0_S1, r_X)

    # not I0: I0_I0
    state_not_I0_I0 = L_I0_I0, I0_out, I0, N_I0_I0
    dL_I0_I0, dd = not_cell_wrapper(state_not_I0_I0, params_not)    
    dI0_out += dd
    dN_I0_I0 = population(N_I0_I0, r_X)


    """
     I1
    """
    dI1_out = 0

    # not S0: I1_S0
    state_not_I1_S0 = L_I1_S0, I1_out, S0, N_I1_S0
    dL_I1_S0, dd = not_cell_wrapper(state_not_I1_S0, params_not)    
    dI1_out += dd
    dN_I1_S0 = population(N_I1_S0, r_X)

    # yes S1: I1_S1
    state_yes_I1_S1 = I1_out, S1, N_I1_S1
    dI1_out += yes_cell_wrapper(state_yes_I1_S1, params_yes)
    dN_I1_S1 = population(N_I1_S1, r_X)

    # not I1: I1_I1
    state_not_I1_I1 = L_I1_I1, I1_out, I1, N_I1_I1
    dL_I1_I1, dd = not_cell_wrapper(state_not_I1_I1, params_not)    
    dI1_out += dd
    dN_I1_I1 = population(N_I1_I1, r_X)


    """
    I2
    """
    dI2_out = 0

    # yes S0: I2_S0
    state_yes_I2_S0 = I2_out, S0, N_I2_S0
    dI2_out += yes_cell_wrapper(state_yes_I2_S0, params_yes)
    dN_I2_S0 = population(N_I2_S0, r_X)

    # not S1: I2_S1
    state_not_I2_S1 = L_I2_S1, I2_out, S1, N_I2_S1
    dL_I2_S1, dd = not_cell_wrapper(state_not_I2_S1, params_not)    
    dI2_out += dd
    dN_I2_S1 = population(N_I2_S1, r_X)
    
    # not I2: I2_I2
    state_not_I2_I2 = L_I2_I2, I2_out, I2, N_I2_I2
    dL_I2_I2, dd = not_cell_wrapper(state_not_I2_I2, params_not)    
    dI2_out += dd
    dN_I2_I2 = population(N_I2_I2, r_X)

    """
    I3
    """
    dI3_out = 0

    # not S0: I3_S0
    state_not_I3_S0 = L_I3_S0, I3_out, S0, N_I3_S0
    dL_I3_S0, dd = not_cell_wrapper(state_not_I3_S0, params_not)    
    dI3_out += dd
    dN_I3_S0 = population(N_I3_S0, r_X)
    
    # not S1: I3_S1
    state_not_I3_S1 = L_I3_S1, I3_out, S1, N_I3_S1
    dL_I3_S1, dd = not_cell_wrapper(state_not_I3_S1, params_not)    
    dI3_out += dd
    dN_I3_S1 = population(N_I3_S1, r_X)

    # not I3: I3_I3
    state_not_I3_I3 = L_I3_I3, I3_out, I3, N_I3_I3
    dL_I3_I3, dd = not_cell_wrapper(state_not_I3_I3, params_not)    
    dI3_out += dd
    dN_I3_I3 = population(N_I3_I3, r_X)

    """
    out
    """
    dout = 0

    # not I0: I0
    state_not_I0 = L_I0, out, I0_out, N_I0
    dL_I0, dd = not_cell_wrapper(state_not_I0, params_not)    
    dout += dd
    dN_I0 = population(N_I0, r_X)

    # not I1: I1
    state_not_I1 = L_I1, out, I1_out, N_I1
    dL_I1, dd = not_cell_wrapper(state_not_I1, params_not)    
    dout += dd
    dN_I1 = population(N_I1, r_X)

    # not I2: I2
    state_not_I2 = L_I2, out, I2_out, N_I2
    dL_I2, dd = not_cell_wrapper(state_not_I2, params_not)    
    dout += dd
    dN_I2 = population(N_I2, r_X)

    # not I3: I3
    state_not_I3 = L_I3, out, I3_out, N_I3
    dL_I3, dd = not_cell_wrapper(state_not_I3, params_not)    
    dout += dd
    dN_I3 = population(N_I3, r_X)

    dS0, dS1, dI0, dI1, dI2, dI3 = 0, 0, 0, 0, 0, 0

    dstate = np.array([dS0, dS1, dI0, dI1, dI2, dI3,
              dI0_out, dI1_out, dI2_out, dI3_out,
              dL_I0_I0, dL_I1_S0, dL_I1_I1, dL_I2_S1, dL_I2_I2, dL_I3_S0, dL_I3_S1, dL_I3_I3, dL_I0, dL_I1, dL_I2, dL_I3,
              dN_I0_S0, dN_I0_S1, dN_I0_I0, dN_I1_S0, dN_I1_S1, dN_I1_I1, dN_I2_S0, dN_I2_S1, dN_I2_I2, dN_I3_S0, dN_I3_S1, dN_I3_I3, dN_I0, dN_I1, dN_I2, dN_I3,
              dout])

    return dstate
    
def MUX_4_1_model_ODE(T, state, params):
    return MUX_4_1_model(state, T, params)