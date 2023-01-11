import numpy as np
import matplotlib.pyplot as plt

from parameters import *
from models import *
from binarize import binarize

def simulate_stochastic_encoder(params, Y0, Omega, T_end, dt = 1): 
            
        state = np.array(Y0)
        
        Y_total = np.zeros([1+T_end//dt, len(state)])
        T = np.zeros(1+T_end//dt)
        t = 0 

        Y_total[0, :] = state
        T[0] = t
   
        """N = MUX_4_1_generate_stoichiometry()"""
        N = ENCODER_4_2_generate_stoichiometry()

        i = 1
        last_time = t

        while t < T_end: 
            

            #choose two random numbers 
            r = np.random.uniform(size=2)
            r1 = r[0]
            r2 = r[1]           

            """a = MUX_4_1_model_stochastic(state, params, Omega)"""
            a = ENCODER_4_2_model_stochastic(state, params, Omega)
                    
            asum = np.cumsum(a)
            a0 = np.sum(a)  
            #get tau
            tau = (1.0/a0)*np.log(1.0/r1)    
        
            #print(t)       
            #select reaction 
            reaction_number = np.argwhere(asum > r2*a0)[0,0] #get first element     
        
            #update concentrations
            state = state + N[:,reaction_number]      
            
            #update time
            t = t + tau  
   

            if (t - last_time >= dt) or (t >= T_end):
                last_time = t
                Y_total[i, :] = state
                T[i] = t                
                i += 1
            


        return T[:i], Y_total[:i,:]     

rho_x = 0 
params = [delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X]
# reaction space volume for the whole cell population
# N_cells should be set to 1
Omega = 10


# I0, I1, I2, I3
inp = [0, 1, 0, 0]
I = np.array(inp)*100
"""
# S0, S1
S = np.array([1, 0])
"""


"""Y0 = np.zeros(39)"""
Y0 = np.zeros(26)
"""Y0[22:38] = 1 # number of cells"""
Y0[10:23] = 1 # number of cells

Y0[:4] = I
"""Y0[4:6] = S"""


T, Y = simulate_stochastic_encoder(params, Y0, Omega, 100)

"""
out = Y[:,-1]
#plt.plot(T,out)
#plt.show()

I0_a, I0_b = Y[:,2], Y[:,3]
I1_a, I1_b = Y[:,8], Y[:,9]
I2_a, I2_b = Y[:,14], Y[:,15]
I3_a, I3_b = Y[:,20], Y[:,21]

out = Y[:,-1]

# plot

I0_out = Y[:,6]
I1_out = Y[:,7]
I2_out = Y[:,8]
I3_out = Y[:,9]

out = Y[:,-1]"""

M0_out = Y[:,4]
M1_out = Y[:,5]
M2_out = Y[:,6]
M3_out = Y[:,7]

O0_out = Y[:,-3]
O1_out = Y[:,-2]
V_out = Y[:,-1]

# plot
"""
ax1 = plt.subplot(241)
ax1.plot(T, I0_out)
ax1.legend(["I0_out"])

ax2 = plt.subplot(242)
ax2.plot(T, I1_out)
ax2.legend(["I1_out"])


ax3 = plt.subplot(243)
ax3.plot(T, I2_out)
ax3.legend(["I2_out"])

ax4 = plt.subplot(244)
ax4.plot(T, I3_out)
ax4.legend(["I3_out"])

ax5 = plt.subplot(212)
ax5.plot(T,out)
ax5.set_title('out')
"""
# plot M outputs
"""
ax1 = plt.subplot(411)
ax1.plot(T, M0_out)
ax1.legend(["M0_out"])

ax2 = plt.subplot(412)
ax2.plot(T, M1_out)
ax2.legend(["M1_out"])

ax3 = plt.subplot(413)
ax3.plot(T, M2_out)
ax3.legend(["M2_out"])

ax4 = plt.subplot(414)
ax4.plot(T, M3_out)
ax4.legend(["M3_out"])
plt.show()
"""

ax5 = plt.subplot(311)
ax5.plot(T, binarize(O0_out))
ax5.legend(["O0"])
ax5.set_xlabel("Time [min]")
ax5.set_ylabel("Molecules")

ax6 = plt.subplot(312)
ax6.plot(T, binarize(O1_out))
ax6.legend(["O1"])
ax6.set_xlabel("Time [min]")
ax6.set_ylabel("Molecules")

ax7 = plt.subplot(313)
ax7.plot(T, binarize(V_out))
ax7.legend(["V"])
ax7.set_xlabel("Time [min]")
ax7.set_ylabel("Molecules")

plt.suptitle("Input: "+str(inp))

"""plt.suptitle(f"S = [{S[1]},{S[0]}]")"""
plt.show()
