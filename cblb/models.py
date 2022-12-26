# Urrios 2016: multicellular memory + Macia 2016

import numpy as np

def merge_N(x,y):
    x1 = np.append(x, np.zeros([x.shape[0], y.shape[1]]), axis=1)
    y1 = np.append(np.zeros([y.shape[0], x.shape[1]]), y, axis=1)
    return np.append(x1,y1,axis=0)

# a model of inverter
def not_cell(state, params):
    L_X, x, y, N_X, N_Y = state
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X


    f = gamma_L_X * (y ** n_y)/(1 + (theta_L_X*y)**n_y )
    dL_X_dt = N_X * (f - delta_L * L_X)

    dx_dt = N_X * (eta_x * (1/(1+ (omega_x*L_X)**m_x))) - N_Y * (delta_x * x) - rho_x * x

    return dL_X_dt, dx_dt

def not_cell_stochastic(state, params, Omega):
    L_X, x, y, N_X, N_Y = state
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X

    #Omega *= N_X # reaction space volume is proportional to the number of cells

    gamma_L_X *= Omega
    eta_x *= Omega
    theta_L_X /= Omega
    omega_x /= Omega
    

    p = [0]*5
    
    p[0] = N_X * gamma_L_X * (y ** n_y)/(1 + (theta_L_X*y)**n_y ) / Omega
    #p[0] = gamma_L_X * (y ** n_y)/(1 + (theta_L_X*y)**n_y ) / Omega # N_x already included in reaction space volume (Omega)
    
    p[1] = N_X * delta_L * L_X 
    #p[1] = delta_L * L_X # N_x already included in reaction space volume (Omega)

    p[2] = N_X * (eta_x * (1/(1+ (omega_x*L_X)**m_x)))
    #p[2] = (eta_x * (1/(1+ (omega_x*L_X)**m_x))) # N_x already included in reaction space volume (Omega)

    p[3] = N_Y * (delta_x * x)
    #p[3] = (delta_x * x) # N_y already included in reaction space volume (Omega)
    
    p[4] = rho_x * x

    return p

# a model of driver
def yes_cell(state, params):
    x, y, N_X, N_Y = state
    gamma_x, n_y, theta_x, delta_x, rho_x = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X


    dx_dt = N_X * gamma_x * (y ** n_y)/(1 + (theta_x*y)**n_y ) - N_Y * (delta_x * x) - rho_x * x
    
    return dx_dt

def yes_cell_stochastic(state, params, Omega):
    x, y, N_X, N_Y = state
    gamma_x, n_y, theta_x, delta_x, rho_x = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X

    #Omega *= N_X # reaction space volume is proportional to the number of cells

    gamma_x *= Omega
    theta_x /= Omega
        

    p = [0]*3

    p[0] = N_X * gamma_x * (y ** n_y)/(1 + (theta_x*y)**n_y )
    #p[0] = gamma_x * (y ** n_y)/(1 + (theta_x*y)**n_y ) # N_x already included in reaction space volume (Omega)
    
    p[1] = N_Y * (delta_x * x) 
    #p[1] = delta_x * x # N_y already included in reaction space volume (Omega)

    p[2] = rho_x * x
    
    return p


def population(state, params):
    N = state
    r = params

    dN = r * N * (1 - N)    

    return dN

def population_stochastic(state, params, Omega):
    N = state
    r = params
    
    p = [0]*2

    p[0] = r * N
    p[1] = r * Omega * N**2

    return p


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

def toggle_model_stochastic(state, params, Omega):
    L_A, L_B, a, b, N_A, N_B = state

    state_A = L_A, a, b, N_A, N_B
    state_B = L_B, b, a, N_B, N_A
    
    delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, r_A, r_B = params

    params_A = delta_L, gamma_A, n_b, theta_A, eta_a, omega_a, m_a, delta_a, rho_a
    params_B = delta_L, gamma_B, n_a, theta_B, eta_b, omega_b, m_b, delta_b, rho_b

    
    p1 = not_cell_stochastic(state_A, params_A, Omega)
    p2 = not_cell_stochastic(state_B, params_B, Omega)

    #p3 = population_stochastic(N_A, r_A, Omega)
    #p4 = population_stochastic(N_B, r_B, Omega)
        
    return p1 + p2

def toggle_generate_stoichiometry():
    #
        # x axis ... species # Y = L_A, L_B, a, b, N_A, N_B
        # y axis ... reactions
        #
        idx_L_A, idx_L_B, idx_a, idx_b, idx_N_A, idx_N_B = 0,1,2,3,4,5

        N = np.zeros((6, 10))

        # reaction 0
        r = 0
        # 0 --> L_A
        N[idx_L_A, r] = 1

        # reaction 1
        r = 1
        # L_A --> 0
        N[idx_L_A, r] = -1

        # reaction 2
        r = 2
        # 0 --> a
        N[idx_a, r] = 1

        # reaction 3
        r = 3
        # a --> 0
        N[idx_a, r] = -1

        # reaction 4
        r = 4
        # a --> 0
        N[idx_a, r] = -1

        # reaction 5
        r = 5
        # 0 --> L_B
        N[idx_L_B, r] = 1

        # reaction 6
        r = 6
        # L_B --> 0
        N[idx_L_B, r] = -1

        # reaction 7
        r = 7
        # 0 --> b
        N[idx_b, r] = 1

        # reaction 8
        r = 8
        # b --> 0
        N[idx_b, r] = -1

        # reaction 9
        r = 9
        # b --> 0
        N[idx_b, r] = -1

        return N


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

def MUX_4_1_model(state, T, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x


    I0, I1, I2, I3, S0, S1 = state[:6]
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

    dI0, dI1, dI2, dI3, dS0, dS1 = 0, 0, 0, 0, 0, 0

    dstate = np.array([dI0, dI1, dI2, dI3, dS0, dS1,
              dI0_out, dI1_out, dI2_out, dI3_out,
              dL_I0_I0, dL_I1_S0, dL_I1_I1, dL_I2_S1, dL_I2_I2, dL_I3_S0, dL_I3_S1, dL_I3_I3, dL_I0, dL_I1, dL_I2, dL_I3,
              dN_I0_S0, dN_I0_S1, dN_I0_I0, dN_I1_S0, dN_I1_S1, dN_I1_I1, dN_I2_S0, dN_I2_S1, dN_I2_I2, dN_I3_S0, dN_I3_S1, dN_I3_I3, dN_I0, dN_I1, dN_I2, dN_I3,
              dout])

    return dstate

def ENCODER_4_2_model(state, T, params):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x


    I0, I1, I2, I3 = state[:4]
    M0_out, M1_out, M2_out, M3_out = state[4:8]
    #L...
    L_M0_I1, L_O0_M0 = state[8:10]
    #N...
    N_M0_I1, N_M0_I2, N_M1_I3, N_M2_I2, N_M2_I3, N_M3_I0, N_M3_I1, N_M3_I2, N_M3_I3, N_O0_M0, N_O0_M1, N_O1, N_V  = state[10:23]
    O0_out, O1_out, V_out = state[23:26]
    """
    M0
    """
    dM0_out = 0

    # not I1: M0_I1
    state_not_M0_I1 = L_M0_I1, M0_out, I1, N_M0_I1
    dL_M0_I1, dd = not_cell_wrapper(state_not_M0_I1, params_not)
    dM0_out += dd
    dN_M0_I1 = population(N_M0_I1, r_X)

    # yes I2: M0_I2
    state_yes_M0_I2 = M0_out, I2, N_M0_I2
    dM0_out += yes_cell_wrapper(state_yes_M0_I2, params_yes)
    dN_M0_I2 = population(N_M0_I2, r_X)


    """
    M1
    """
    dM1_out = 0

    # yes I3: M1_I3
    state_yes_M1_I3 = M1_out, I3, N_M1_I3
    dM1_out += yes_cell_wrapper(state_yes_M1_I3, params_yes)
    dN_M1_I3 = population(N_M1_I3, r_X)


    """
    M2
    """
    dM2_out = 0

    # yes I2: M2_I2
    state_yes_M2_I2 = M2_out, I2, N_M2_I2
    dM2_out += yes_cell_wrapper(state_yes_M2_I2, params_yes)
    dN_M2_I2 = population(N_M2_I2, r_X)

    # yes I3: M2_I3
    state_yes_M2_I3 = M2_out, I3, N_M2_I3
    dM2_out += yes_cell_wrapper(state_yes_M2_I3, params_yes)
    dN_M2_I3 = population(N_M2_I3, r_X)

    """
    M3
    """
    dM3_out = 0

    # yes I0: M3_I0
    state_yes_M3_I0 = M3_out, I0, N_M3_I0
    dM3_out += yes_cell_wrapper(state_yes_M3_I0, params_yes)
    dN_M3_I0 = population(N_M3_I0, r_X)

    # yes I1: M3_I1
    state_yes_M3_I1 = M3_out, I1, N_M3_I1
    dM3_out += yes_cell_wrapper(state_yes_M3_I1, params_yes)
    dN_M3_I1 = population(N_M3_I1, r_X)

    # yes I2: M3_I2
    state_yes_M3_I2 = M3_out, I2, N_M3_I2
    dM3_out += yes_cell_wrapper(state_yes_M3_I2, params_yes)
    dN_M3_I2 = population(N_M3_I2, r_X)

    # yes I3: M3_I3
    state_yes_M3_I3 = M3_out, I3, N_M3_I3
    dM3_out += yes_cell_wrapper(state_yes_M3_I3, params_yes)
    dN_M3_I3 = population(N_M3_I3, r_X)


    """
    O0
    """
    dO0_out = 0

    # not M0: O0_M0
    state_not_O0_M0 = L_O0_M0, O0_out, M0_out, N_O0_M0
    dL_O0_M0, dd = not_cell_wrapper(state_not_O0_M0, params_not)
    dO0_out += dd
    dN_O0_M0 = population(N_O0_M0, r_X)

    # yes M1: O0_M1
    state_yes_O0_M1 = O0_out, M1_out, N_O0_M1
    dO0_out += yes_cell_wrapper(state_yes_O0_M1, params_yes)
    dN_O0_M1 = population(N_O0_M1, r_X)


    """
    O1
    """
    dO1_out = 0

    #yes M2: O1_M2
    state_yes_O1_M2 = O1_out, M2_out, N_O1
    dO1_out += yes_cell_wrapper(state_yes_O1_M2, params_yes)
    dN_O1 = population(N_O1, r_X)


    """
    V
    """
    dV_out = 0

    #yes M3: V_M3
    state_yes_V_M3 = V_out, M3_out, N_V
    dV_out += yes_cell_wrapper(state_yes_V_M3, params_yes)
    dN_V = population(N_V, r_X)






    dI0, dI1, dI2, dI3 = 0, 0, 0, 0

    dstate = np.array([
        dI0, dI1, dI2, dI3,
        dM0_out, dM1_out, dM2_out, dM3_out,
        dL_M0_I1, dL_O0_M0,
        dN_M0_I1, dN_M0_I2, dN_M1_I3, dN_M2_I2, dN_M2_I3, dN_M3_I0, dN_M3_I1, dN_M3_I2, dN_M3_I3, dN_O0_M0, dN_O0_M1, dN_O1, dN_V,
        dO0_out, dO1_out, dV_out
    ])

    return dstate

def MUX_4_1_generate_stoichiometry():

    """
    I0, I1, I2, I3, S0, S1 = state[:6]
    I0_out, I1_out, I2_out, I3_out = state[6:10]
    L_I0_I0, L_I1_S0, L_I1_I1, L_I2_S1, L_I2_I2, L_I3_S0, L_I3_S1, L_I3_I3, L_I0, L_I1, L_I2, L_I3 = state[10:22]
    N_I0_S0, N_I0_S1, N_I0_I0, N_I1_S0, N_I1_S1, N_I1_I1, N_I2_S0, N_I2_S1, N_I2_I2, N_I3_S0, N_I3_S1, N_I3_I3, N_I0, N_I1, N_I2, N_I3 = state[22:38]
    out = state[38]
    """
    #I0, I1, I2, I3, S0, S1 = range(6)
    I0_out, I1_out, I2_out, I3_out = range(6,10)
    L_I0_I0, L_I1_S0, L_I1_I1, L_I2_S1, L_I2_I2, L_I3_S0, L_I3_S1, L_I3_I3, L_I0, L_I1, L_I2, L_I3 = range(10,22)
    #N_I0_S0, N_I0_S1, N_I0_I0, N_I1_S0, N_I1_S1, N_I1_I1, N_I2_S0, N_I2_S1, N_I2_I2, N_I3_S0, N_I3_S1, N_I3_I3, N_I0, N_I1, N_I2, N_I3 = range(22,38)
    out = 38

    #
    # x axis ... species 
    # y axis ... reactions
    #
    N = np.zeros((39, 72))

    """
    # yes S0: I0_S0
    """

    r = 0
    # reaction 0
    # 0 --> I0_out     
    N[I0_out, r] = 1
    
    
    r += 1        
    # reaction 1
    # I0_out --> 0
    N[I0_out, r] = -1
    
    r += 1
    # reaction 2
    # I0_out --> 0
    N[I0_out, r] = -1
    
    """
    # yes S1: I0_S1 
    """

    r += 1
    # reaction 3
    # 0 --> I0_out
    N[I0_out, r] = 1

    r += 1
    # reaction 4
    # I0_out --> 0
    N[I0_out, r] = -1

    r += 1
    # reaction 5
    # I0_out --> 0
    N[I0_out, r] = -1

    """
    # not I0: I0_I0 
    """

    r += 1
    # reaction 6
    # 0 --> L_I0_I0
    N[L_I0_I0, r] = 1

    r += 1
    # reaction 7
    # L_I0_I0 --> 0
    N[L_I0_I0, r] = -1
    
    r += 1
    # reaction 8
    # 0 --> I0_out     
    N[I0_out, r] = 1

    r += 1
    # reaction 9
    # I0_out --> 0
    N[I0_out, r] = -1

    r += 1
    # reaction 10
    # I0_out --> 0
    N[I0_out, r] = -1

    """
    # not S0: I1_S0
    """

    r += 1
    # reaction 11
    # 0 --> L_I1_S0
    N[L_I1_S0, r] = 1

    r += 1
    # reaction 12
    # L_I1_S0 --> 0
    N[L_I1_S0, r] = -1

    r += 1
    # reaction 13
    # 0 --> I1_out
    N[I1_out, r] = 1
    
    r += 1
    # reaction 14
    # I1_out --> 0
    N[I1_out, r] = -1

    r += 1
    # reaction 15
    # I1_out --> 0
    N[I1_out, r] = -1

    """
    # yes S1: I1_S1
    """

    r += 1
    # reaction 16
    # 0 --> I1_out
    N[I1_out, r] = 1

    r += 1
    # reaction 17
    # I1_out --> 0
    N[I1_out, r] = -1

    r += 1
    # reaction 18
    # I1_out --> 0
    N[I1_out, r] = -1

    """
    # not I1: I1_I1
    """

    r += 1
    # reaction 19
    # 0 --> L_I1_I1
    N[L_I1_I1, r] = 1

    r += 1
    # reaction 20
    # L_I1_I1 --> 0
    N[L_I1_I1, r] = -1

    r += 1
    # reaction 21
    # 0 --> I1_out
    N[I1_out, r] = 1

    r += 1
    # reaction 22
    # I1_out --> 0
    N[I1_out, r] = -1

    r += 1
    # reaction 23
    # I1_out --> 0
    N[I1_out, r] = -1

    """
    # yes S0: I2_S0
    """

    r += 1
    # reaction 24
    # 0 --> I2_out
    N[I2_out, r] = 1

    r += 1
    # reaction 25
    # I2_out --> 0
    N[I2_out, r] = -1

    r += 1
    # reaction 26
    # I2_out --> 0
    N[I2_out, r] = -1

    """
    # not S1: I2_S1
    """

    r += 1
    # reaction 27
    # 0 --> L_I2_S1
    N[L_I2_S1, r] = 1

    r += 1
    # reaction 28
    # L_I2_S1 --> 0
    N[L_I2_S1, r] = -1

    r += 1
    # reaction 29
    # 0 --> I2_out
    N[I2_out, r] = 1
    
    r += 1
    # reaction 30
    # I2_out --> 0
    N[I2_out, r] = -1

    r += 1
    # reaction 31
    # I2_out --> 0
    N[I2_out, r] = -1

    """
    # not I2: I2_I2
    """

    r += 1
    # reaction 32
    # 0 --> L_I2_I2
    N[L_I2_I2, r] = 1

    r += 1
    # reaction 33
    # L_I2_I2 --> 0
    N[L_I2_I2, r] = -1

    r += 1
    # reaction 34
    # 0 --> I2_out
    N[I2_out, r] = 1

    r += 1
    # reaction 35
    # I2_out --> 0
    N[I2_out, r] = -1

    r += 1
    # reaction 36
    # I2_out --> 0
    N[I2_out, r] = -1

    """ 
    # not S0: I3_S0
    """

    r += 1
    # reaction 37
    # 0 --> L_I3_S0
    N[L_I3_S0, r] = 1

    r += 1
    # reaction 38
    # 0 --> L_I3_S0
    N[L_I3_S0, r] = -1

    r += 1
    # reaction 39
    # 0 --> I3_out
    N[I3_out, r] = 1


    r += 1
    # reaction 40
    # I3_out --> 0
    N[I3_out, r] = -1

    r += 1
    # reaction 41
    # I3_out --> 0
    N[I3_out, r] = -1

    """
    # not S1: I3_S1
    """

    r += 1
    # reaction 42
    # 0 --> L_I3_S1
    N[L_I3_S1, r] = 1

    r += 1
    # reaction 43
    # L_I3_S1 --> 0
    N[L_I3_S1, r] = -1

    r += 1
    # reaction 44
    # 0 --> I3_out
    N[L_I3_S1, r] = 1

    r += 1
    # reaction 45
    # I3_out --> 0
    N[I3_out, r] = -1

    r += 1
    # reaction 46
    # I3_out --> 0
    N[I3_out, r] = -1

    """
    # not I3: I3_I3
    """

    r += 1
    # reaction 47
    # 0 --> L_I3_I3
    N[L_I3_I3, r] = 1

    r += 1
    # reaction 48
    # L_I3_I3 --> 0
    N[L_I3_I3, r] = -1

    r += 1
    # reaction 49
    # 0 --> I3_out
    N[I3_out, r] = 1

    r += 1
    # reaction 50
    # I3_out --> 0
    N[I3_out, r] = -1

    r += 1
    # reaction 51
    # I3_out --> 0
    N[I3_out, r] = -1

    """ 
    # not I0: I0
    """

    r += 1
    # reaction 52
    # 0 --> L_I0
    N[L_I0, r] = 1

    r += 1
    # reaction 53
    # L_I0 --> 0
    N[L_I0, r] = -1

    r += 1
    # reaction 54
    # 0 --> out
    N[out, r] = 1

    r += 1
    # reaction 55
    # out --> 0
    N[out, r] = -1

    r += 1
    # reaction 56
    # out --> 0
    N[out, r] = -1

    """
    # not I1: I1
    """

    r += 1
    # reaction 57
    # 0 --> L_I1
    N[L_I1, r] = 1

    r += 1
    # reaction 58
    # L_I1 --> 0
    N[L_I1, r] = -1

    r += 1
    # reaction 59
    # 0 --> out
    N[out, r] = 1

    r += 1
    # reaction 60
    # out --> 0
    N[out, r] = -1

    r += 1
    # reaction 61
    # out --> 0
    N[out, r] = -1

    """
    # not I2: I2
    """

    r += 1
    # reaction 62
    # 0 --> L_I2
    N[L_I2, r] = 1

    r += 1
    # reaction 63
    # L_I2 --> 0
    N[L_I2, r] = -1

    r += 1
    # reaction 64
    # 0 --> out
    N[out, r] = 1

    r += 1
    # reaction 65
    # out --> 0
    N[out, r] = -1

    r += 1
    # reaction 66
    # out --> 0
    N[out, r] = -1

    """
    # not I3: I3
    """

    r += 1
    # reaction 67
    # 0 --> L_I3
    N[L_I3, r] = 1

    r += 1
    # reaction 68
    # L_I3 --> 0
    N[L_I3, r] = -1

    r += 1
    # reaction 69
    # 0 --> out
    N[out, r] = 1

    r += 1
    # reaction 70
    # out --> 0
    N[out, r] = -1

    r += 1
    # reaction 71
    # out --> 0
    N[out, r] = -1

    return N

def ENCODER_4_2_generate_stoichiometry():
    raise NotImplementedError

def MUX_4_1_model_stochastic(state, params, Omega):
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x


    I0, I1, I2, I3, S0, S1 = state[:6]
    I0_out, I1_out, I2_out, I3_out = state[6:10]
    L_I0_I0, L_I1_S0, L_I1_I1, L_I2_S1, L_I2_I2, L_I3_S0, L_I3_S1, L_I3_I3, L_I0, L_I1, L_I2, L_I3 = state[10:22]
    N_I0_S0, N_I0_S1, N_I0_I0, N_I1_S0, N_I1_S1, N_I1_I1, N_I2_S0, N_I2_S1, N_I2_I2, N_I3_S0, N_I3_S1, N_I3_I3, N_I0, N_I1, N_I2, N_I3 = state[22:38]
    out = state[38]
    
        
    """
     I0
    """
        
    # yes S0: I0_S0
    state_yes_I0_S0 = I0_out, S0, N_I0_S0, N_I0_S0
    p_I0_S0 = yes_cell_stochastic(state_yes_I0_S0, params_yes, Omega)
    

    # yes S1: I0_S1
    state_yes_I0_S1 = I0_out, S1, N_I0_S1, N_I0_S1
    p_I0_S1 = yes_cell_stochastic(state_yes_I0_S1, params_yes, Omega)
    

    # not I0: I0_I0
    state_not_I0_I0 = L_I0_I0, I0_out, I0, N_I0_I0, N_I0_I0
    p_I0_I0 = not_cell_stochastic(state_not_I0_I0, params_not, Omega)    
    
    
    """
     I1
    """
    
    # not S0: I1_S0
    state_not_I1_S0 = L_I1_S0, I1_out, S0, N_I1_S0, N_I1_S0
    p_I1_S0 = not_cell_stochastic(state_not_I1_S0, params_not, Omega)    
    
    
    # yes S1: I1_S1
    state_yes_I1_S1 = I1_out, S1, N_I1_S1, N_I1_S1
    p_I1_S1 = yes_cell_stochastic(state_yes_I1_S1, params_yes, Omega)
    
    # not I1: I1_I1
    state_not_I1_I1 = L_I1_I1, I1_out, I1, N_I1_I1, N_I1_I1
    p_I1_I1 = not_cell_stochastic(state_not_I1_I1, params_not, Omega)    
    

    """
    I2
    """
    
    # yes S0: I2_S0
    state_yes_I2_S0 = I2_out, S0, N_I2_S0, N_I2_S0
    p_I2_S0 = yes_cell_stochastic(state_yes_I2_S0, params_yes, Omega)
    

    # not S1: I2_S1
    state_not_I2_S1 = L_I2_S1, I2_out, S1, N_I2_S1, N_I2_S1
    p_I2_S1= not_cell_stochastic(state_not_I2_S1, params_not, Omega)    
    
    
    # not I2: I2_I2
    state_not_I2_I2 = L_I2_I2, I2_out, I2, N_I2_I2, N_I2_I2
    p_I2_I2 = not_cell_stochastic(state_not_I2_I2, params_not, Omega)    
       

    """
    I3
    """
    # not S0: I3_S0
    state_not_I3_S0 = L_I3_S0, I3_out, S0, N_I3_S0, N_I3_S0
    p_I3_S0 = not_cell_stochastic(state_not_I3_S0, params_not, Omega)    
        
    
    # not S1: I3_S1
    state_not_I3_S1 = L_I3_S1, I3_out, S1, N_I3_S1, N_I3_S1
    p_I3_S1 = not_cell_stochastic(state_not_I3_S1, params_not, Omega)    
       

    # not I3: I3_I3
    state_not_I3_I3 = L_I3_I3, I3_out, I3, N_I3_I3, N_I3_I3
    p_I3_I3 = not_cell_stochastic(state_not_I3_I3, params_not, Omega)    
       

    """
    out
    """
    # not I0: I0
    state_not_I0 = L_I0, out, I0_out, N_I0, N_I0
    p_I0 = not_cell_stochastic(state_not_I0, params_not, Omega)    
        

    # not I1: I1
    state_not_I1 = L_I1, out, I1_out, N_I1, N_I1
    p_I1 = not_cell_stochastic(state_not_I1, params_not, Omega)    
    

    # not I2: I2
    state_not_I2 = L_I2, out, I2_out, N_I2, N_I2
    p_I2 = not_cell_stochastic(state_not_I2, params_not, Omega)    
        

    # not I3: I3
    state_not_I3 = L_I3, out, I3_out, N_I3, N_I3
    p_I3 = not_cell_stochastic(state_not_I3, params_not, Omega)    
       
    
    return (p_I0_S0 + p_I0_S1 + p_I0_I0 + 
           p_I1_S0 + p_I1_S1 + p_I1_I1 +
           p_I2_S0 + p_I2_S1 + p_I2_I2 +
           p_I3_S0 + p_I3_S1 + p_I3_I3 +
           p_I0 + p_I1 + p_I2 + p_I3)

def ENCODER_4_2_model_stochastic(state, params, Omega):
    raise NotImplementedError

def CLB_generate_stoichiometry():
    N_toggle_IO = toggle_generate_stoichiometry()
    N_toggle_I1 = toggle_generate_stoichiometry()
    N_toggle_I2 = toggle_generate_stoichiometry()
    N_toggle_I3 = toggle_generate_stoichiometry()


    N_mux = MUX_4_1_generate_stoichiometry()
    # skip first four rows (I0, I1, I2, I3)
    N_mux = N_mux[4:,:]

    return merge_N(merge_N(merge_N(merge_N(N_toggle_IO, N_toggle_I1), N_toggle_I2), N_toggle_I3), N_mux)

def CLB_4_2_ENCODER_generate_stoichiometry():
    raise NotImplementedError

def CLB_model(state, T, params):
    
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, gamma_x, theta_x, r_X, r_Y, rho_I0_a, rho_I0_b, rho_I1_a, rho_I1_b, rho_I2_a, rho_I2_b, rho_I3_a, rho_I3_b = params
 
  
    """
    latches
    """
    #########
    # params

    # set params for symmetric toggle switch topology
    gamma_L_Y, theta_L_Y = gamma_L_X, theta_L_X
    n_x, m_y = n_y, m_x
    eta_y, omega_y = eta_x, omega_x
 
    params_toggle =  [delta_L, gamma_L_X, gamma_L_Y, n_x, n_y, theta_L_X, theta_L_Y, eta_x, eta_y, omega_x, omega_y, m_x, m_y, delta_x, delta_y, rho_x, rho_y, r_X, r_Y]

    # degradation rates for induction of switches are specific for each toggle switch    
    params_toggle_I0 =  params_toggle.copy()
    params_toggle_I0[-4:-2] = rho_I0_a, rho_I0_b
    params_toggle_I1 =  params_toggle.copy()
    params_toggle_I1[-4:-2] = rho_I1_a, rho_I1_b
    params_toggle_I2 =  params_toggle.copy()
    params_toggle_I2[-4:-2] = rho_I2_a, rho_I2_b
    params_toggle_I3 =  params_toggle.copy()
    params_toggle_I3[-4:-2] = rho_I3_a, rho_I3_b   

    #########
    # states
    
    # latch I0
    I0_L_A, I0_L_B, I0_a, I0_b, I0_N_a, I0_N_b = state[:6]
    state_toggle_IO = I0_L_A, I0_L_B, I0_a, I0_b, I0_N_a, I0_N_b

    # latch I1
    I1_L_A, I1_L_B, I1_a, I1_b, I1_N_a, I1_N_b = state[6:12]
    state_toggle_I1 = I1_L_A, I1_L_B, I1_a, I1_b, I1_N_a, I1_N_b

    # latch I2
    I2_L_A, I2_L_B, I2_a, I2_b, I2_N_a, I2_N_b = state[12:18]
    state_toggle_I2 = I2_L_A, I2_L_B, I2_a, I2_b, I2_N_a, I2_N_b

    # latch I3
    I3_L_A, I3_L_B, I3_a, I3_b, I3_N_a, I3_N_b = state[18:24]
    state_toggle_I3 = I3_L_A, I3_L_B, I3_a, I3_b, I3_N_a, I3_N_b

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
    params_mux = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X

    #########
    # state
    I0, I1, I2, I3 = I0_a, I1_a, I2_a, I3_a
    state_mux = np.append([I0, I1, I2, I3], state[24:], axis=0)

    ########
    # model
    dstate_mux = MUX_4_1_model(state_mux, T, params_mux)
    dstate_mux = dstate_mux[4:] # ignore dI0, dI1, dI2, dI3

    """
    return
    """
    dstate = np.append(dstate_toggles, dstate_mux, axis = 0)
    return dstate

def CLB_4_2_ENCODER_model(state, T, params):
    
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, gamma_x, theta_x, r_X, r_Y, rho_I0_a, rho_I0_b, rho_I1_a, rho_I1_b, rho_I2_a, rho_I2_b, rho_I3_a, rho_I3_b = params
 
  
    """
    latches
    """
    #########
    # params

    # set params for symmetric toggle switch topology
    gamma_L_Y, theta_L_Y = gamma_L_X, theta_L_X
    n_x, m_y = n_y, m_x
    eta_y, omega_y = eta_x, omega_x
 
    params_toggle =  [delta_L, gamma_L_X, gamma_L_Y, n_x, n_y, theta_L_X, theta_L_Y, eta_x, eta_y, omega_x, omega_y, m_x, m_y, delta_x, delta_y, rho_x, rho_y, r_X, r_Y]

    # degradation rates for induction of switches are specific for each toggle switch    
    params_toggle_I0 =  params_toggle.copy()
    params_toggle_I0[-4:-2] = rho_I0_a, rho_I0_b
    params_toggle_I1 =  params_toggle.copy()
    params_toggle_I1[-4:-2] = rho_I1_a, rho_I1_b
    params_toggle_I2 =  params_toggle.copy()
    params_toggle_I2[-4:-2] = rho_I2_a, rho_I2_b
    params_toggle_I3 =  params_toggle.copy()
    params_toggle_I3[-4:-2] = rho_I3_a, rho_I3_b   

    #########
    # states
    
    # latch I0
    I0_L_A, I0_L_B, I0_a, I0_b, I0_N_a, I0_N_b = state[:6]
    state_toggle_IO = I0_L_A, I0_L_B, I0_a, I0_b, I0_N_a, I0_N_b

    # latch I1
    I1_L_A, I1_L_B, I1_a, I1_b, I1_N_a, I1_N_b = state[6:12]
    state_toggle_I1 = I1_L_A, I1_L_B, I1_a, I1_b, I1_N_a, I1_N_b

    # latch I2
    I2_L_A, I2_L_B, I2_a, I2_b, I2_N_a, I2_N_b = state[12:18]
    state_toggle_I2 = I2_L_A, I2_L_B, I2_a, I2_b, I2_N_a, I2_N_b

    # latch I3
    I3_L_A, I3_L_B, I3_a, I3_b, I3_N_a, I3_N_b = state[18:24]
    state_toggle_I3 = I3_L_A, I3_L_B, I3_a, I3_b, I3_N_a, I3_N_b

    #########
    # models
    dstate_toggle_IO = toggle_model(state_toggle_IO, T, params_toggle_I0)
    dstate_toggle_I1 = toggle_model(state_toggle_I1, T, params_toggle_I1)
    dstate_toggle_I2 = toggle_model(state_toggle_I2, T, params_toggle_I2)
    dstate_toggle_I3 = toggle_model(state_toggle_I3, T, params_toggle_I3)

    dstate_toggles = np.append(np.append(np.append(dstate_toggle_IO, dstate_toggle_I1, axis=0), dstate_toggle_I2, axis = 0), dstate_toggle_I3, axis = 0)

    """
    decoder
    """
    #########
    # params
    params_dec = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X

    #########
    # state
    I0, I1, I2, I3 = I0_a, I1_a, I2_a, I3_a
    state_dec = np.append([I0, I1, I2, I3], state[24:], axis=0)

    ########
    # model
    dstate_dec = ENCODER_4_2_model(state_dec, T, params_dec)
    dstate_dec = dstate_dec[4:] # ignore dI0, dI1, dI2, dI3

    """
    return
    """
    dstate = np.append(dstate_toggles, dstate_dec, axis = 0)
    return dstate

def CLB_model_stochastic(state, params, Omega):
    
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, gamma_x, theta_x, r_X, r_Y, rho_I0_a, rho_I0_b, rho_I1_a, rho_I1_b, rho_I2_a, rho_I2_b, rho_I3_a, rho_I3_b = params
 
  
    """
    latches
    """
    #########
    # params

    # set params for symmetric toggle switch topology
    gamma_L_Y, theta_L_Y = gamma_L_X, theta_L_X
    n_x, m_y = n_y, m_x
    eta_y, omega_y = eta_x, omega_x
 
    params_toggle =  [delta_L, gamma_L_X, gamma_L_Y, n_x, n_y, theta_L_X, theta_L_Y, eta_x, eta_y, omega_x, omega_y, m_x, m_y, delta_x, delta_y, rho_x, rho_y, r_X, r_Y]

    # degradation rates for induction of switches are specific for each toggle switch    
    params_toggle_I0 =  params_toggle.copy()
    params_toggle_I0[-4:-2] = rho_I0_a, rho_I0_b
    params_toggle_I1 =  params_toggle.copy()
    params_toggle_I1[-4:-2] = rho_I1_a, rho_I1_b
    params_toggle_I2 =  params_toggle.copy()
    params_toggle_I2[-4:-2] = rho_I2_a, rho_I2_b
    params_toggle_I3 =  params_toggle.copy()
    params_toggle_I3[-4:-2] = rho_I3_a, rho_I3_b   

    #########
    # states
    
    # latch I0
    I0_L_A, I0_L_B, I0_a, I0_b, I0_N_a, I0_N_b = state[:6]
    state_toggle_IO = I0_L_A, I0_L_B, I0_a, I0_b, I0_N_a, I0_N_b

    # latch I1
    I1_L_A, I1_L_B, I1_a, I1_b, I1_N_a, I1_N_b = state[6:12]
    state_toggle_I1 = I1_L_A, I1_L_B, I1_a, I1_b, I1_N_a, I1_N_b

    # latch I2
    I2_L_A, I2_L_B, I2_a, I2_b, I2_N_a, I2_N_b = state[12:18]
    state_toggle_I2 = I2_L_A, I2_L_B, I2_a, I2_b, I2_N_a, I2_N_b

    # latch I3
    I3_L_A, I3_L_B, I3_a, I3_b, I3_N_a, I3_N_b = state[18:24]
    state_toggle_I3 = I3_L_A, I3_L_B, I3_a, I3_b, I3_N_a, I3_N_b

    #########
    # models
    p_toggle_IO = toggle_model_stochastic(state_toggle_IO, params_toggle_I0, Omega)
    p_toggle_I1 = toggle_model_stochastic(state_toggle_I1, params_toggle_I1, Omega)
    p_toggle_I2 = toggle_model_stochastic(state_toggle_I2, params_toggle_I2, Omega)
    p_toggle_I3 = toggle_model_stochastic(state_toggle_I3, params_toggle_I3, Omega)

    """
    mux
    """
    #########
    # params
    params_mux = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X

    #########
    # state
    I0, I1, I2, I3 = I0_a, I1_a, I2_a, I3_a
    state_mux = np.append([I0, I1, I2, I3], state[24:], axis=0)

    ########
    # model
    p_mux = MUX_4_1_model_stochastic(state_mux, params_mux, Omega)

    """
    return
    """

    return p_toggle_IO + p_toggle_I1 + p_toggle_I2 + p_toggle_I3 + p_mux

def CLB_4_2_ENCODER_model_stochastic(state, params, Omega):
    raise NotImplementedError

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

def ENCODER_4_2_model_ODE(T, state, params):
    return ENCODER_4_2_model(state, T, params)

def CLB_4_2_ENCODER_model_ODE(T, state, params):
    return CLB_4_2_ENCODER_model(T, state, params)
