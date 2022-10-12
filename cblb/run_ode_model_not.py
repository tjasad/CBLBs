from scipy.integrate import ode
import matplotlib.pyplot as plt

# import models and parameter values
from models import *
from parameters import *

# increased degradation rates are set to 0
rho_x = 0
rho_y = 0

# prepare a tuple with parameter values
params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, r_X 


# set simulation parameters
t_end = 1500 # simulation duration
N = t_end # number of samples

# set initial conditions
# Y = L_A, a, b, N_A
Y0 = np.zeros(4)
Y0[2] = 0 # b
Y0[3] = 1 # N_A

# prepare a vector of timepoints (samples) for plotting
#T = np.linspace(0, t_end, N)

###
# simulations
t1 = t_end # simulate until this time is reached
dt = t_end/N # timestep for simulation
T = np.arange(0,t1+dt,dt) # prepare a vector of timepoints (samples) for plotting
Y = np.zeros([1+N,4]) # prepare a vector in which observed variables will be stored in each timepoint
Y[0,:] = Y0 # set initial conditions

r = ode(not_model_ODE).set_integrator('zvode', method='bdf') # prepare a solver using scipy.integrate.ode; you can change a solver in case convergence fails (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html)
r.set_initial_value(Y0, T[0]).set_f_params(params) # give initial conditions to the solver

i = 1
while r.successful() and r.t < t1:  # run simulations until...
    Y[i,:] = r.integrate(r.t+dt)
    i += 1


####
# analysis of results

# Y = L_A, a, b, N_A
L_A = Y[:,0] 
a = Y[:,1]
b = Y[:,2]
N_A = Y[:,3]

ax1 = plt.subplot(311)
ax1.plot(T,L_A)
ax1.legend(["L_A"])

ax2 = plt.subplot(312)
ax2.plot(T,a)
ax2.plot(T,b)
ax2.legend(["a", "b"])

ax3 = plt.subplot(313)
ax3.plot(T,N_A)
ax3.legend(["N_A"])


plt.show()
