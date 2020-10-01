# Urrios 2016: multicellular memory + Macia 2016

import numpy as np

def not_cell(state, params):
    L_X, x, y = state
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, N_X, N_Y = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X

    f = gamma_L_X * (y ** n_y)/(1 + (theta_L_X*y)**n_y )
    dL_X_dt = N_X * (f - delta_L * L_X)

    dx_dt = N_X * (eta_x * (1/(1+ (omega_x*L_X)**m_x))) - N_Y * (delta_x * x) - rho_x * x

    return dL_X_dt, dx_dt

def yes_cell(state, params):
    x, y = state
    gamma_x, n_y, theta_x, delta_x, rho_x, N_X, N_Y = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X

    dx_dt = N_X * gamma_x * (y ** n_y)/(1 + (theta_x*y)**n_y ) - N_Y * (delta_x * x) - rho_x * x
    
    return dx_dt


def toggle_model(state, T, params):
    L_A, L_B, a, b = state

    state_A = L_A, a, b
    state_B = L_B, b, a
    
    delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, N_A, N_B = params

    params_A = delta_L, gamma_A, n_b, theta_A, eta_a, omega_a, m_a, delta_a, rho_a, N_A, N_B
    params_B = delta_L, gamma_B, n_a, theta_B, eta_b, omega_b, m_b, delta_b, rho_b, N_B, N_A

    dL_A_dt, da_dt = not_cell(state_A, params_A)
    dL_B_dt, db_dt = not_cell(state_B, params_B)
           
    return np.array([dL_A_dt, dL_B_dt, da_dt, db_dt])


# L_A ... intermediate
# a ... out
# b ... in
def not_cell_wrapper(state, params):
    L_A, a, b = state
    
    state_A = L_A, a, b
    
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, N_X = params
    params_A = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, N_X, N_X
    

    return not_cell(state_A, params_A)

# a ... out
# b ... in
def yes_cell_wrapper(state, params):
    a, b = state

    state_A = a, b
    gamma_x, n_y, theta_x, delta_x, rho_x, N_X = params
    params_A = gamma_x, n_y, theta_x, delta_x, rho_x, N_X, N_X
    

    return yes_cell(state_A, params_A)

def not_model(state, T, params):
    L_A, a, b = state

    delta_L, gamma_L_A, n_b, theta_L_A, eta_a, omega_a, m_a, delta_a, delta_b, rho_a, rho_b, N_A = params

    state_not = L_A, a, b
    params_not = delta_L, gamma_L_A, n_b, theta_L_A, eta_a, omega_a, m_a, delta_a, rho_a, N_A
    dL_A_dt, da_dt = not_cell_wrapper(state_not, params_not)
    
    db_dt = 0#- N_A * delta_b * b - rho_b * b

    return np.array([dL_A_dt, da_dt, db_dt])

def yes_model(state, T, params):
    a, b = state
    
    gamma_a, n_b, theta_a, delta_a, delta_b, rho_a, rho_b, N_A = params

    state_yes = a, b
    params_yes = gamma_a, n_b, theta_a, delta_a, rho_a, N_A
    da_dt = yes_cell_wrapper(state_yes, params_yes)
    
    db_dt = 0 #- N_A * delta_b * b - rho_b * b

    return np.array([da_dt, db_dt])

def MUX_4_1_model(state, T, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x = params[:11]
    N_I0_S0, N_I0_S1, N_I0_I0, N_I1_S0, N_I1_S1, N_I1_I1, N_I2_S0, N_I2_S1, N_I2_I2, N_I3_S0, N_I3_S1, N_I3_I3, N_I0, N_I1, N_I2, N_I3 = params[11:]


    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x


    I0, I1, I2, I3, S0, S1 = state[:6]
    I0_out, I1_out, I2_out, I3_out = state[6:10]
    L_I0_I0, L_I1_S0, L_I1_I1, L_I2_S1, L_I2_I2, L_I3_S0, L_I3_S1, L_I3_I3, L_I0, L_I1, L_I2, L_I3 = state[10:22]
    
    out = state[-1]
    
        
    """
     I0
    """
    dI0_out = 0
    
    # yes S0: I0_S0
    state_yes_I0_S0 = I0_out, S0
    dI0_out += yes_cell_wrapper(state_yes_I0_S0, params_yes + (N_I0_S0,))
    
    # yes S1: I0_S1
    state_yes_I0_S1 = I0_out, S1
    dI0_out += yes_cell_wrapper(state_yes_I0_S1, params_yes + (N_I0_S1,))    

    # not I0: I0_I0
    state_not_I0_I0 = L_I0_I0, I0_out, I0
    dL_I0_I0, dd = not_cell_wrapper(state_not_I0_I0, params_not + (N_I0_I0,))    
    dI0_out += dd
    

    """
     I1
    """
    dI1_out = 0

    # not S0: I1_S0
    state_not_I1_S0 = L_I1_S0, I1_out, S0
    dL_I1_S0, dd = not_cell_wrapper(state_not_I1_S0, params_not + (N_I1_S0,))    
    dI1_out += dd
    
    # yes S1: I1_S1
    state_yes_I1_S1 = I1_out, S1
    dI1_out += yes_cell_wrapper(state_yes_I1_S1, params_yes + (N_I1_S1,))
    
    # not I1: I1_I1
    state_not_I1_I1 = L_I1_I1, I1_out, I1
    dL_I1_I1, dd = not_cell_wrapper(state_not_I1_I1, params_not + (N_I1_I1,))    
    dI1_out += dd
    
    """
    I2
    """
    dI2_out = 0

    # yes S0: I2_S0
    state_yes_I2_S0 = I2_out, S0
    dI2_out += yes_cell_wrapper(state_yes_I2_S0, params_yes + (N_I2_S0,))
    
    # not S1: I2_S1
    state_not_I2_S1 = L_I2_S1, I2_out, S1
    dL_I2_S1, dd = not_cell_wrapper(state_not_I2_S1, params_not + (N_I2_S1,))    
    dI2_out += dd
        
    # not I2: I2_I2
    state_not_I2_I2 = L_I2_I2, I2_out, I2
    dL_I2_I2, dd = not_cell_wrapper(state_not_I2_I2, params_not + (N_I2_I2,))    
    dI2_out += dd
    

    """
    I3
    """
    dI3_out = 0

    # not S0: I3_S0
    state_not_I3_S0 = L_I3_S0, I3_out, S0
    dL_I3_S0, dd = not_cell_wrapper(state_not_I3_S0, params_not + (N_I3_S0,))    
    dI3_out += dd
        
    # not S1: I3_S1
    state_not_I3_S1 = L_I3_S1, I3_out, S1
    dL_I3_S1, dd = not_cell_wrapper(state_not_I3_S1, params_not + (N_I3_S1,))    
    dI3_out += dd
    
    # not I3: I3_I3
    state_not_I3_I3 = L_I3_I3, I3_out, I3
    dL_I3_I3, dd = not_cell_wrapper(state_not_I3_I3, params_not + (N_I3_I3,))    
    dI3_out += dd
    
    """
    out
    """
    dout = 0

    # not I0: I0
    state_not_I0 = L_I0, out, I0_out
    dL_I0, dd = not_cell_wrapper(state_not_I0, params_not + (N_I0,))    
    dout += dd
    
    # not I1: I1
    state_not_I1 = L_I1, out, I1_out
    dL_I1, dd = not_cell_wrapper(state_not_I1, params_not + (N_I1,))    
    dout += dd
    
    # not I2: I2
    state_not_I2 = L_I2, out, I2_out
    dL_I2, dd = not_cell_wrapper(state_not_I2, params_not + (N_I2,))     
    dout += dd
    

    # not I3: I3
    state_not_I3 = L_I3, out, I3_out
    dL_I3, dd = not_cell_wrapper(state_not_I3, params_not + (N_I3,))    
    dout += dd
    
    dI0, dI1, dI2, dI3, dS0, dS1 = 0, 0, 0, 0, 0, 0

    dstate = np.array([dI0, dI1, dI2, dI3, dS0, dS1,
              dI0_out, dI1_out, dI2_out, dI3_out,
              dL_I0_I0, dL_I1_S0, dL_I1_I1, dL_I2_S1, dL_I2_I2, dL_I3_S0, dL_I3_S1, dL_I3_I3, dL_I0, dL_I1, dL_I2, dL_I3,              
              dout])

    return dstate



def CLB_model(state, T, params):
    
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, gamma_x, theta_x, rho_I0_a, rho_I0_b, rho_I1_a, rho_I1_b, rho_I2_a, rho_I2_b, rho_I3_a, rho_I3_b = params[:21]
    # toggle switches
    N_TS0, N_not_TS0, N_TS1, N_not_TS1, N_TS2, N_not_TS2, N_TS3, N_not_TS3 = params[21:29]
    # mux   
    N_I0_S0, N_I0_S1, N_I0_I0, N_I1_S0, N_I1_S1, N_I1_I1, N_I2_S0, N_I2_S1, N_I2_I2, N_I3_S0, N_I3_S1, N_I3_I3, N_I0, N_I1, N_I2, N_I3 = params[29:]
 
  
    """
    latches
    """
    #########
    # params

    # set params for symmetric toggle switch topology
    gamma_L_Y, theta_L_Y = gamma_L_X, theta_L_X
    n_x, m_y = n_y, m_x
    eta_y, omega_y = eta_x, omega_x
 
    params_toggle =  [delta_L, gamma_L_X, gamma_L_Y, n_x, n_y, theta_L_X, theta_L_Y, eta_x, eta_y, omega_x, omega_y, m_x, m_y, delta_x, delta_y, rho_x, rho_y]
    
    # degradation rates for induction of switches are specific for each toggle switch    
    params_toggle_I0 =  params_toggle.copy()
    params_toggle_I0[-2:] = rho_I0_a, rho_I0_b
    params_toggle_I0 += [N_TS0, N_not_TS0]

    params_toggle_I1 =  params_toggle.copy()
    params_toggle_I1[-2:] = rho_I1_a, rho_I1_b
    params_toggle_I1 += [N_TS1, N_not_TS1]

    params_toggle_I2 =  params_toggle.copy()
    params_toggle_I2[-2:] = rho_I2_a, rho_I2_b
    params_toggle_I2 += [N_TS2, N_not_TS2]

    params_toggle_I3 =  params_toggle.copy()
    params_toggle_I3[-2:] = rho_I3_a, rho_I3_b   
    params_toggle_I3 += [N_TS3, N_not_TS3]

    #########
    # states
    
    # latch I0
    I0_L_A, I0_L_B, I0_a, I0_b = state[:4]
    state_toggle_IO = I0_L_A, I0_L_B, I0_a, I0_b

    # latch I1
    I1_L_A, I1_L_B, I1_a, I1_b = state[4:8]
    state_toggle_I1 = I1_L_A, I1_L_B, I1_a, I1_b

    # latch I2
    I2_L_A, I2_L_B, I2_a, I2_b = state[8:12]
    state_toggle_I2 = I2_L_A, I2_L_B, I2_a, I2_b

    # latch I3
    I3_L_A, I3_L_B, I3_a, I3_b = state[12:16]
    state_toggle_I3 = I3_L_A, I3_L_B, I3_a, I3_b

    #########
    # models
    dstate_toggle_IO = toggle_model(state_toggle_IO, T, params_toggle_I0)
    dstate_toggle_I1 = toggle_model(state_toggle_I1, T, params_toggle_I1)
    dstate_toggle_I2 = toggle_model(state_toggle_I2, T, params_toggle_I2)
    dstate_toggle_I3 = toggle_model(state_toggle_I3, T, params_toggle_I3)

    dstate_toggles = np.append(np.append(np.append(dstate_toggle_IO, dstate_toggle_I1, axis=0), dstate_toggle_I2, axis = 0), dstate_toggle_I3, axis = 0)

    """
    mux
    """
    #########
    # params
    params_mux = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x
    params_mux += N_I0_S0, N_I0_S1, N_I0_I0, N_I1_S0, N_I1_S1, N_I1_I1, N_I2_S0, N_I2_S1, N_I2_I2, N_I3_S0, N_I3_S1, N_I3_I3, N_I0, N_I1, N_I2, N_I3 

    #########
    # state
    I0, I1, I2, I3 = I0_a, I1_a, I2_a, I3_a
    state_mux = np.append([I0, I1, I2, I3], state[16:], axis=0)

    ########
    # model
    dstate_mux = MUX_4_1_model(state_mux, T, params_mux)
    dstate_mux = dstate_mux[4:] # ignore dI0, dI1, dI2, dI3

    """
    return
    """
    dstate = np.append(dstate_toggles, dstate_mux, axis = 0)
    return dstate




"""
wrappers for scipy.integrate.ode
"""

def toggle_model_ODE(T, state, params):
    return toggle_model(state, T, params)


def not_model_ODE(T, state, params):
    return not_model(state, T, params)
    
def yes_model_ODE(T, state, params):
    return yes_model(state, T, params)

def MUX_4_1_model_ODE(T, state, params):
    return MUX_4_1_model(state, T, params)

def CLB_model_ODE(T, state, params):
    return CLB_model(state, T, params)

