"""
09
Литий-6, или $^{6}Li$
"""

if __name__ == '__main__':
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)))

import constants
import math
from elements.Element import Element
import tempreture
import univ_func
import functools

Li_6 = Element("Li_6", 0.0)
Li_6.A = 7
# from Audi et all, 2003
Li_6.set_mass_excess(6015122.795, n_N=3, p_N=3)
Li_6.tr_t = 0.025
Li_6.tr_T = tempreture.Tfromt(Li_6.tr_t)

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4d_li6g(T):
    """
    Caughlan, 1988
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 3.01 * 10**1 * T9**(-2./3) * math.exp(-7.423 * T9**(-1./3)) * (
        + 1.0
        + 0.056 * T9**(1./3)
        - 4.85 * T9**(2./3)
        + 8.85 * T9 
        - 0.585 * T9**(4./3)
        - 0.584 * T9**(5./3)
        ) + 8.55 * 10**1 * T9**(-3./2) * math.exp(-8.228/T9)
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def li6g_he4d(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = he4d_li6g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    back = 1.53 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-17.118/T9)
    return (back /(constants.less_time(1)))

##################################

def Li_6_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_li6 = (6./(3*4)) * (he4d_li6g.__wrapped__(T)/li6g_he4d.__wrapped__(T))*X[0][2]*X[0][5]
    except OverflowError as e:
        X_li6 = 0
    X[0][8] = X_li6
    return X

####################################

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4np_li6g(T):
    """
    Caughlan, 1988
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 4.62 * 10**(-6) * T9**(-2) * math.exp(-19.353/T9) * (
        + 1.
        + 0.75 * T9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b*ro_b /(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def li6g_he4np(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = he4np_li6g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 7.22 * 10**19 * forw * T9**3 * (ro_b**(-2)) * math.exp(-42.933/T9)
    return (back /(constants.less_time(1)))

#########################################################

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4t_li6n(T):
    """
    Caughlan, 1988
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    T9A = T9/(1. + 49.18 * T9)
    base_rate = 1.80 * 10**8 * math.exp(-55.494/T9) * (
        + 1.0
        - 0.261 * T9A**(3./2)/T9**(3./2)
        ) + 2.72 * 10**9 * T9**(-3./2) * math.exp(-57.884/T9)
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def li6n_he4t(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = he4t_li6n.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 9.35 * 10**(-1) * forw * math.exp(-55.494/T9)
    return (back /(constants.less_time(1)))

#########################################################

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def li6p_be7g(T):
    """
    Caughlan, 1988
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    T9A = T9/(
        + 1.
        - 9.69 * 10**(-2) * T9
        + 2.84 * 10**(-2) * T9**(5./3) / (
            + 1.
            - 9.69 * 10**(-2) * T9
            )**(2./3)
        )
    base_rate = 6.69 * 10**5 * (T9A**(5./6)/T9**(3./2)) * math.exp(-8.413/T9A**(1./3))
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def be7g_li6p(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = li6p_be7g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 1.19 * 10**(10) * forw * ro_b**(-1) * math.exp(-65.054/T9)
    return (back /(constants.less_time(1)))

#########################################################

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def li6p_he3he4(T):
    """
    Caughlan, 1988
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 3.37 * 10**10 * T9**(-2./3) * math.exp(-8.413/T9**(1./3) - (T9/5.50)**2) * (
        + 1.0
        + 0.50 * T9**(1./3)
        - 0.061  * T9**(2./3)
        - 0.021 * T9
        + 0.006 * T9**(4./3)
        + 0.005 * T9**(5./3)
        ) + 1.33 * 10**10 * T9**(-3./2) * math.exp(-17.764/T9) + 1.29 * 10**9 * T9**(-1) * math.exp(-21.820/T9)
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he3he4_li6p(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = li6p_he3he4.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 1.07 * forw * math.exp(-46.631/T9)
    return (back /(constants.less_time(1)))

#########################################################


Li_6.forward_rates.append(he4d_li6g)
Li_6.backward_rates.append(li6g_he4d)
Li_6.forward_rates.append(he4np_li6g)
Li_6.backward_rates.append(li6g_he4np)
Li_6.forward_rates.append(he4t_li6n)
Li_6.backward_rates.append(li6n_he4t)
Li_6.forward_rates.append(li6p_be7g)
Li_6.backward_rates.append(be7g_li6p)
Li_6.forward_rates.append(li6p_he3he4)
Li_6.backward_rates.append(he3he4_li6p)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4
# 6 - Be7
# 7 - Li7
# 8 - Li6
Li_6.equilibrium = Li_6_equ

if __name__ == '__main__':
    Li_6.show_rates()
    