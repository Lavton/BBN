from constants import *

def p_n__g_d(T, units="eV"):
    """Experimental computation and obsevational analysis of primordial nucleosynthesis, 
    sj 1993, M.S.Smith and L.H. Wawano 227"""
    if units=="eV":
        T = T/k_b
    T /= 10**9
    coef = 4.72*10**4
    in_breckets = 1. - \
        0.850*T**(1/2.) + \
        0.490*T - \
        0.0962*T**(3/2.) + \
        8.47*10**(-3)*T**2 - \
        2.80*10**(-4)*T**(5/2.)

    return coef*in_breckets


if __name__ == '__main__':
    p_n__g_d(1.)
    import matplotlib as mpl
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    plt.xscale('log')
    plt.yscale('log')

    # Ts = np.logspace(math.log10(10**-4), math.log10(10**2), num=200)
    # Ts *= 10**6
    Ts = np.logspace(math.log10(10**8), math.log10(10**11), num=200)*k_b

    # Ts = Ts*k_b/(10**6)
    plt.plot(Ts/10**6, [p_n__g_d(T)*10**4 for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{n\to p+e+\nu}$')
    plt.show()

