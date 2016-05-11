from constants import *


def p_n__g_d(T, *, units="eV"):
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