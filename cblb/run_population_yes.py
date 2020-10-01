from scipy.integrate import ode
import matplotlib.pyplot as plt

from models_population import *
from parameters import *

rho_x = 0
rho_y = 0

N_A = 0.5
params = gamma_x, n_y, theta_x, delta_x, delta_y, rho_x, rho_y, N_A

# simulation parameters
t_end = 1500
N = t_end

# Y = a, b
Y0 = np.zeros(2)
Y0[1] = 10


T = np.linspace(0, t_end, N)

t1 = t_end
dt = t_end/N
T = np.arange(0,t1+dt,dt)
Y = np.zeros([1+N,len(Y0)])
Y[0,:] = Y0

r = ode(yes_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T[0]).set_f_params(params)

i = 1
while r.successful() and r.t < t1:
    Y[i,:] = r.integrate(r.t+dt)
    i += 1


a = Y[:,0]
b = Y[:,1]


ax1 = plt.subplot(111)
ax1.plot(T,a)
ax1.plot(T,b)
ax1.legend(["a", "b"])


plt.show()
