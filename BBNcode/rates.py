from constants import *
import constants
import math
from math import pi

def __photon_density__(T):
    # T in eV
    zeta_3 = 1.202
    # return (2 * (8*pi)/c**3 * zeta_3) / (h/(2*pi*constants.to_norm_tempreture(T)))**3
    return 0.244*(constants.to_norm_tempreture(T)/(h*c))**3

def __proton_mass_density__(T):
    return (m_p/5.60958835719e+32) * __photon_density__(T) * nu_n

def rat_scale(T):
    return __proton_mass_density__(T) # * constants.N_a

def p_n__g_d(T):
    """Experimental computation and obsevational analysis of primordial nucleosynthesis, 
    sj 1993, M.S.Smith and L.H. Wawano 227"""
    sigma_v = 6 * 10**-20
    return constants.less_time(sigma_v*__photon_density__(T)/10**20) #*N_a
    # if units=="eV":
    #     T = T/k_b
    # T /= 10**9
    # coef = 4.72*10**4
    # in_breckets = 1. - \
    #     0.850*T**(1/2.) + \
    #     0.490*T - \
    #     0.0962*T**(3/2.) + \
    #     8.47*10**(-3)*T**2 - \
    #     2.80*10**(-4)*T**(5/2.)

    # print(coef*in_breckets)

    return coef*in_breckets

def d_g__n_p(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    zn = 1.440e-5 * T9**(3./2) * nu_n * math.exp(25.815/T9)
    return 2 * p_n__g_d(T) / zn


def bind_en(Z, A, alpha):
    alpha0=7.297352533e-3
    b=-14.04*A+89* (A/2-Z)**2 /A+15 * A**(2./3.)
    b+=0.675 * (Z**2) * A**(-1./3.) *alpha/alpha0

    if A%2:
        if Z%2 == 0:
            b+=33.5 * A**-0.75
        else:
            b-=33.5 * A**-0.75

    return b

def pngd(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    T912 = T9 ** 0.5
    T932 = T9 ** (3./2)
    T92 = T9 ** 2
    T952 = T9 ** (5./2)
    alpha0=7.297352533e-3
    forw=1.-0.8504*T912+0.4895*T9-0.09623*T932;
    forw=47420*(forw+0.008471*T92-0.00028*T952);

    forw *= ro_b*alpha/alpha0
    alpha=alpha0*(1+dylaton)

def dgpn(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    T912 = T9 ** 0.5
    T932 = T9 ** (3./2)
    T92 = T9 ** 2
    T952 = T9 ** (5./2)
    alpha0=7.297352533e-3
    forw = pngd(T)
    m_d_amu=2.01355321271
    del_E = E_d_t9=constants.E_d_t0*Bind_en(1,m_d_amu,alpha)/Bind_en(1,m_d_amu,alpha0);
    if (T9 >= abs(del_E/250)):
        rev=4.71e9*T932*math.exp(del_E/T9)*forw/ro_b
    else:
        rev=0



if __name__ == '__main__':
    print ('photon density today ~{}'.format(__photon_density__(constants.less_tempreture(2.7255, units="K"))))
    # p_n__g_d(1.)
    # import matplotlib as mpl
    # import numpy as np
    # import math
    # import matplotlib.pyplot as plt
    # plt.xscale('log')
    # plt.yscale('log')

    # # Ts = np.logspace(math.log10(10**-4), math.log10(10**2), num=200)
    # # Ts *= 10**6
    # Ts = np.logspace(math.log10(10**8), math.log10(10**11), num=20)*k_b

    # # Ts = Ts*k_b/(10**6)
    # plt.plot(Ts/10**6, [p_n__g_d(T)*10**4 for T in Ts], 
    #     linewidth=2.0, label=r'$\lambda_{n\to p+e+\nu}$')
    # plt.show()

