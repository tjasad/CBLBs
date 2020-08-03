from parameters import *
from model_clb import model_clb

clb = model_clb()

rho = 5
#candidate = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, gamma_x, theta_x, rho
candidate = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, gamma_x, theta_x, rho

out = clb.simulate(candidate, plot_on=True)
print(clb.getFitness(out))

#print(clb.eval(candidate))