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
# Li_6.mass_excess = constants.less_tempreture(2161062.7, units="eV")
Li_6.set_mass_excess(6015122.795, n_N=3, p_N=3)
Li_6.tr_t = 0.020
Li_6.tr_T = tempreture.Tfromt(Li_6.tr_t)
# Li_6.tr_T = constants.less_tempreture(2*10**10, units="K")

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def the4_li7g(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 5.28 * 10**5 * T9**(-2./3) * math.exp(-8.08 * T9**(-1./3)) * (
        + 1.0
        + 0.0516 * T9**(1./3)
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def li7g_the4(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = the4_li7g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    back = 1.12 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-28.63/T9)
    return (back /(constants.less_time(1)))

##################################

def Li_6_equ(X, T):
    T9 = constants.to_norm_tempreture(-T, units="T9")
    try:
        X_be7 = (7./(3*4)) * (the4_li7g.__wrapped__(T)/li7g_the4.__wrapped__(T))*X[0][4]*X[0][5]
    except OverflowError as e:
        X_be7 = 0
    X[0][6] = X_be7
    return X

####################################

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def nbe7_pli7(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = (
        6.74 * 10**9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def pli7_nbe7(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = nbe7_pli7.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    # back = forw * math.exp(-Li_6.mass_excess/T)
    ro_b = univ_func.rat_scale(T)

    back = forw * math.exp(-19.07/T9)
    return (back /(constants.less_time(1)))

#########################################################

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def pli7_he4he4(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = (
        1.2 * 10**7 * T9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4he4_pli7(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = pli7_he4he4.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    # back = forw * math.exp(-Li_6.mass_excess/T)
    ro_b = univ_func.rat_scale(T)

    back = 4.64 * forw * math.exp(-220.4/T9)
    return (back /(constants.less_time(1)))

#########################################################


Li_6.forward_rates.append(the4_li7g)
Li_6.backward_rates.append(li7g_the4)
Li_6.forward_rates.append(nbe7_pli7)
Li_6.backward_rates.append(pli7_nbe7)
Li_6.forward_rates.append(pli7_he4he4)
Li_6.backward_rates.append(he4he4_pli7)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4
# 6 - Be7
# 7 - Li7
# 8 - Li6
Li_6.ode_elem = {
}

Li_6.jacob = { 
}

Li_6.equilibrium = Li_6_equ


if __name__ == '__main__':
    Li_6.show_rates()
    