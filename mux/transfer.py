# transfer function
def T_f(I, gamma, alpha, omega, n):
    return gamma * (1 + alpha * omega * (I**n)) / (1 + omega * (I**n))

def T_f_fit(I, gamma, alpha1, alpha2, omega1, omega2, n):
    #return gamma * (1 + alpha * omega * (I**n)) / (1 + omega * (I**n))
    alpha = alpha1 * 10**(-alpha2)
    omega = omega1 * 10**(-omega2)

    return gamma * (1 + alpha * omega * (I**n)) / (1 + omega * (I**n))


# inverse transfer function
def iT_f(O, gamma, alpha, omega, n):
    #return (((O/gamma) - 1)/(alpha*omega - (O/gamma) * omega))**(1/n)
    #if alpha-(O/gamma) == 0: # pri izpeljavi (invertiranju) smo predpostavili, da je alpha-(O/gamma) != 0 (to ni nujno res)
    #    return 1000
    #else:
    return ((O-gamma)/(alpha*gamma*omega-omega*O))**(1/n)
