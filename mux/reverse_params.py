import numpy as np
import matplotlib.pyplot as plt
import pickle


from scipy.optimize import curve_fit

from mux import *
from transfer import *
from data import *

bounds_ID = ([min_gamma, min_alpha1, min_alpha2, min_omega1, min_omega2, min_n],
          [max_gamma, max_alpha1, max_alpha2, max_omega1, max_omega2, max_n])

bounds_NOT = ([min_gamma, min_alpha1, min_alpha2, min_omega1, min_omega2, min_n],
          [max_gamma, max_alpha1, max_alpha2+4, max_omega1, max_omega2, max_n])

direct_params = {}

plot_on = True

for name, max_val in zip(names[:], max_vals[:]):
    x = np.logspace(0, max_val, 1000)
    x = x.astype(np.double)

    """
    ID
    """

    y = T_f(x, *params[cells["ID_"+name]])
    GFP_x = iT_f(y,*params[cells["NOT_GFP"]])

    idx = np.isfinite(GFP_x)
    x1 = x[idx]
    y1 = y[idx]
    GFP_x = GFP_x[idx]
    
    try:      
        popt1, pcov = curve_fit(T_f_fit, x1, GFP_x, bounds=bounds_ID)
        #popt11, pcov = curve_fit(T_f_fit, x1, GFP_x)
        #print(popt1)
        gamma, alpha1, alpha2, omega1, omega2, n = popt1
        alpha = alpha1 * 10**(-alpha2)
        omega = omega1 * 10**(-omega2)

        y_fit1 = T_f(x1, gamma, alpha, omega, n)


        #y_fit11 = T_f_fit(x1, *popt11)

        if plot_on:
            plt.plot(x1,y1, label = "ID_"+name)
            plt.plot(x1,GFP_x, label = "GFP_x")
            plt.plot(x1, y_fit1, label = "fit")
            #plt.plot(x1, y_fit11, label = "fit no bounds")
            plt.legend()

            ax = plt.gca()
            ax.set_xscale('log')    
            
            plt.show()

        direct_params["ID_"+name] = (gamma, alpha, omega, n)
        
    except:
        print("ID_"+name)

    """
    NOT
    """
    
    not_y = T_f(x, *params[cells["NOT_"+name]])
    not_GFP_x = iT_f(not_y,*params[cells["NOT_GFP"]])

    idx = np.isfinite(not_GFP_x)
    x2 = x[idx]
    y2 = not_y[idx]
    not_GFP_x = not_GFP_x[idx]
       
    try:
        popt2, pcov = curve_fit(T_f_fit, x2, not_GFP_x, bounds=bounds_NOT)
        #popt22, pcov = curve_fit(T_f_fit, x2, not_GFP_x)
        gamma, alpha1, alpha2, omega1, omega2, n = popt2
        alpha = alpha1 * 10**(-alpha2)
        omega = omega1 * 10**(-omega2)

        y_fit2 = T_f(x2, gamma, alpha, omega, n)
        
        #y_fit2 = T_f_fit(x2, *popt2)
        #y_fit22 = T_f_fit(x2, *popt22)

        if plot_on:
            plt.plot(x2,y2, label = "NOT_"+name)
            plt.plot(x2,not_GFP_x, label = "NOT_GFP_x")
            plt.plot(x2, y_fit2, label = "fit")
            #plt.plot(x2, y_fit22, label = "fit no bounds")
            plt.legend()

            ax = plt.gca()
            ax.set_xscale('log')    
            
            plt.show()

        direct_params["NOT_"+name] = (gamma, alpha, omega, n)
    except:
        print("NOT_"+name)

with open('direct_params.pickle', 'wb') as handle:
    pickle.dump(direct_params, handle)

