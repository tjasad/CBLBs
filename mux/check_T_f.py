from mux import *
from transfer import *
import numpy as np
import matplotlib.pyplot as plt

from parameters import *


for name, max_val in zip(names[:], max_vals[:]):
    x = np.logspace(0, max_val, 1000)
    x = x.astype(np.double)

    print("ID_"+name, end="...")
    print(*params[cells["ID_"+name]])
    print("NOT_ID_"+name, end="...")
    print(*params[cells["NOT_"+name]])
    print("######")    

    y = T_f(x, *params[cells["ID_"+name]])
    inv_y = iT_f(y,*params[cells["ID_"+name]])
    #print(x-inv_y)


    not_y = T_f(x, *params[cells["NOT_"+name]])
    not_inv_y = iT_f(not_y,*params[cells["NOT_"+name]])
    
    #print(x-not_inv_y)

    gamma, alpha, omega, n = params[cells["NOT_"+name]]
    l = alpha-(not_y/gamma)    

    plt.plot(l,x, '.', label='x')
    plt.plot(l, not_inv_y, 'x',label='x_approx')
    plt.legend()
    plt.show()


    plt.plot(x,y, label = "ID_"+name)
    plt.plot(inv_y,y,'o', label = "ID_"+name+" approx.")
    plt.plot(x,not_y, label = "NOT_"+name)
    plt.plot(not_inv_y,not_y,'x', label = "NOT_"+name+" approx.")
    plt.legend()

    ax = plt.gca()
    ax.set_xscale('log')    
    
    plt.show()

    #plt.plot(x,not_inv_y)
    #plt.show()

    #break
