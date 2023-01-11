from scipy.integrate import ode
import matplotlib.pyplot as plt

from models import *
from parameters import *
from binarize import *


rho_x = 0
rho_y = 0

params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X

O0_list = []
O1_list = []
V_list = []

# I0, I1, I2, I3
I_list = [
    np.array([0,0,0,0]),
    np.array([1,0,0,0]),
    np.array([0,1,0,0]),
    np.array([0,0,1,0]),
    np.array([0,0,0,1]),
]
#I = np.array([0, 1, 0, 0])
"""# S0, S1
S = np.array([1, 1])"""

for I in I_list:


    # simulation parameters
    t_end = 1500
    N = t_end


    """Y0 = np.zeros(39)"""
    Y0 = np.zeros(26)
    """Y0[22:38] = 1 # number of cells"""
    Y0[10:23] = 1 # number of cells

    Y0[:4] = I
    """Y0[4:6] = S"""


    T = np.linspace(0, t_end, N)

    t1 = t_end
    dt = t_end/N
    T = np.arange(0,t1+dt,dt)
    """Y = np.zeros([1+N,39])"""
    Y = np.zeros([1+N,26])
    Y[0,:] = Y0

    r = ode(ENCODER_4_2_model_ODE).set_integrator('zvode', method='bdf')
    r.set_initial_value(Y0, T[0]).set_f_params(params)

    i = 1
    while r.successful() and r.t < t1:
        Y[i,:] = r.integrate(r.t+dt)
        i += 1

    O0_list.append(Y[:,-3])
    O1_list.append(Y[:,-2])
    V_list.append(Y[:,-1])

l = len(I_list)

if l > 1:
    fig, axs = plt.subplots(l)
    fig.tight_layout(rect=(0.025,0.025,1,1))
    fig.supylabel('Value')
    fig.supxlabel('Time [min]')

    for i, (O0, O1, V) in enumerate(zip(O0_list, O1_list, V_list)):
        axs[i].plot(T,binarize(O0), label="O0")
        axs[i].plot(T,binarize(O1), label="O1")
        axs[i].plot(T,binarize(V), label="V")
        axs[i].set_title("Input: "+str(I_list[i]))
        
else:
    plt.plot(T,binarize(O0_list[0]), label="O0")    
    plt.plot(T,binarize(O1_list[0]), label="O1")
    plt.plot(T,binarize(V_list[0]), label="V")
    plt.xlabel('Concentrations [nM]')
    plt.ylabel('Time [min]')
    plt.suptitle("Input: "+str(I_list[0]))





plt.legend()
plt.show()


"""
# Y = a, b, N_A
a = Y[:,0]
b = Y[:,1]
N_A = Y[:,2]

ax1 = plt.subplot(211)
ax1.plot(T,a)
ax1.plot(T,b)
ax1.legend(["a", "b"])

ax2 = plt.subplot(212)
ax2.plot(T,N_A)
ax2.legend(["N_A"])


plt.show()
"""