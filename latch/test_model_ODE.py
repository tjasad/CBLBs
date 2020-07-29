from scipy.integrate import ode
import matplotlib.pyplot as plt

from models import *
from parameters import *

#params = [delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, r_A, r_B]
#params = [delta_L, gamma_A, gamma_A, n_a, n_a, theta_A, theta_A, eta_a, eta_a, omega_a, omega_a, m_a, m_a, delta_a, delta_a, rho_a, rho_a, r_A, r_A]
params = [delta_L, gamma_B, gamma_B, n_b, n_b, theta_B, theta_B, eta_b, eta_b, omega_b, omega_b, m_b, m_b, delta_b, delta_b, rho_b, rho_b, r_B, r_B]

# simulation parameters
t_end = 500
N = t_end*2

# Y = L_A, L_B, a, b, N_A, N_B
Y0 = np.array([0.0]*6)
Y0[-2] = 0.1
Y0[-1] = 0.1
Y0[2] = 0
T1 = np.linspace(0, t_end, N)
T2 = np.linspace(0, t_end, N)
T3 = np.linspace(0, t_end, N)
T4 = np.linspace(0, t_end, N)
T5 = np.linspace(0, t_end, N)

# 1

params[-4] = rho_a
params[-3] = 0

t1 = t_end
dt = t_end/N
T1 = np.arange(0,t1+dt,dt)
Y1 = np.zeros([1+N,6])
Y1[0,:] = Y0

r = ode(toggle_model_ODE).set_integrator('lsoda', method='adams')
r.set_initial_value(Y0, T1[0]).set_f_params(params)

i = 1
while r.successful() and r.t < t1:
    Y1[i,:] = r.integrate(r.t+dt)
    i += 1

# 2 

params[-4] = 0
params[-3] = 0

t1 = t_end
dt = t_end/N
T2 = np.arange(0,t1+dt,dt)
Y0 = Y1[-1,:]
Y2 = np.zeros([1+N,6])
Y2[0,:] = Y0

r = ode(toggle_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T2[0]).set_f_params(params)
i = 1
while r.successful() and r.t < t1:
    Y2[i,:] = r.integrate(r.t+dt)
    i += 1

T2 += t_end

# 3

params[-4] = 0
params[-3] = rho_b

t1 = t_end
dt = t_end/N
T3 = np.arange(0,t1+dt,dt)
Y0 = Y2[-1,:]
Y3 = np.zeros([1+N,6])
Y3[0,:] = Y0

r = ode(toggle_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T3[0]).set_f_params(params)
i = 1
while r.successful() and r.t < t1:
    Y3[i,:] = r.integrate(r.t+dt)
    i += 1

T3 += 2*t_end

# 4

params[-4] = 0
params[-3] = 0

t1 = t_end
dt = t_end/N
T4 = np.arange(0,t1+dt,dt)
Y0 = Y3[-1,:]
Y4 = np.zeros([1+N,6])
Y4[0,:] = Y0

r = ode(toggle_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T4[0]).set_f_params(params)
i = 1
while r.successful() and r.t < t1:
    Y4[i,:] = r.integrate(r.t+dt)
    i += 1

T4 += 3*t_end

# 5

params[-4] = rho_a
params[-3] = 0

t1 = t_end
dt = t_end/N
T5 = np.arange(0,t1+dt,dt)
Y0 = Y4[-1,:]
Y5 = np.zeros([1+N,6])
Y5[0,:] = Y0

r = ode(toggle_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T5[0]).set_f_params(params)
i = 1
while r.successful() and r.t < t1:
    Y5[i,:] = r.integrate(r.t+dt)
    i += 1

T5 += 4*t_end

# 6

params[-4] = 0
params[-3] = 0

t1 = t_end
dt = t_end/N
T6 = np.arange(0,t1+dt,dt)
Y0 = Y5[-1,:]
Y6 = np.zeros([1+N,6])
Y6[0,:] = Y0

r = ode(toggle_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T6[0]).set_f_params(params)
i = 1
while r.successful() and r.t < t1:
    Y6[i,:] = r.integrate(r.t+dt)
    i += 1

T6 += 5*t_end


"""
Y1 = odeint(toggle_model, Y0, T1, args=((tuple(params),)))

params[-4] = 0
Y2 = odeint(toggle_model, Y1[-1], T2, args=((tuple(params),)))

params[-3] = rho_b
Y3 = odeint(toggle_model, Y2[-2], T3, args=((tuple(params),)))

params[-3] = 0
Y4 = odeint(toggle_model, Y3[-2], T4, args=((tuple(params),)))

params[-3] = 0
Y5 = odeint(toggle_model, Y4[-2], T5, args=((tuple(params),)))

T2 += T1[-1]
T3 += T2[-1]
T4 += T3[-1]
T5 += T4[-1]

Y = np.append(np.append(np.append(np.append(Y1, Y2, axis=0), Y3, axis=0), Y4, axis=0), Y5, axis=0)
T = np.append(np.append(np.append(np.append(T1, T2, axis=0), T3, axis=0), T4, axis=0), T5, axis=0)
"""

Y = np.append(np.append(np.append(np.append(np.append(Y1, Y2, axis=0), Y3, axis=0), Y4, axis=0), Y5, axis=0), Y6, axis=0)
T = np.append(np.append(np.append(np.append(np.append(T1, T2, axis=0), T3, axis=0), T4, axis=0), T5, axis=0), T6, axis=0)

L_A = Y[:,0]
L_B = Y[:,1]
a = Y[:,2]
b = Y[:,3]
N_A = Y[:,4]
N_B = Y[:,5]

ax1 = plt.subplot(311)
ax1.plot(T,L_A)
ax1.plot(T,L_B)
ax1.legend(["L_A", "L_B"])

ax2 = plt.subplot(312)
ax2.plot(T,a)
ax2.plot(T,b)
ax2.legend(["a", "b"])

ax3 = plt.subplot(313)
ax3.plot(T,N_A)
ax3.plot(T,N_B)
ax3.legend(["N_A", "N_B"])

#plt.plot(Y)
#plt.legend(["L_A", "L_B", "a", "b", "N_A", "N_B"])
#plt.xticks([0, len(T)-1], [0, T[-1]])

plt.show()
