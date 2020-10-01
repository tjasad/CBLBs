from scipy.integrate import ode
import matplotlib.pyplot as plt

from models_population import *
from parameters import *

rho_x = 0
rho_y = 0

N_A = 1
params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, N_A


# simulation parameters
t_end = 1500
N = t_end

# Y = L_A, a, b
Y0 = np.zeros(3)
Y0[2] = 0 # b


T = np.linspace(0, t_end, N)



t1 = t_end
dt = t_end/N
T = np.arange(0,t1+dt,dt)
Y = np.zeros([1+N,len(Y0)])
Y[0,:] = Y0

r = ode(not_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T[0]).set_f_params(params)

i = 1
while r.successful() and r.t < t1:
    Y[i,:] = r.integrate(r.t+dt)
    i += 1

# Y = L_A, a, b
L_A = Y[:,0]
a = Y[:,1]
b = Y[:,2]


ax1 = plt.subplot(211)
ax1.plot(T,L_A)
ax1.legend(["L_A"])

ax2 = plt.subplot(212)
ax2.plot(T,a)
ax2.plot(T,b)
ax2.legend(["a", "b"])


plt.show()
