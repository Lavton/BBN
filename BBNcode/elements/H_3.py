"""
05
тритий, или $^{3}H$
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

H_3 = Element("H_3", 0.0)
H_3.A = 3
# from Audi et all, 2003
H_3.set_mass_excess(3016049.2777, n_N=2, p_N=1)
H_3.tr_t =  0.0008 * 6./4
H_3.tr_T = tempreture.Tfromt(H_3.tr_t)

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def nd_tg(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = (
        75.5
        + 1250 * T9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def tg_nd(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = nd_tg.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 1.63 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-72.62/T9)
    return (back /(constants.less_time(1)))

#############################################

def H_3_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_h3 = (3./2) * (nd_tg.__wrapped__(T)/tg_nd.__wrapped__(T))*X[0][0]*X[0][2]
    except OverflowError as e:
        X_h3 = 0
    X[0][4] = X_h3
    return X

#############################################

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def nhe3_pt(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = (
        7.06 * (10**8)
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def pt_nhe3(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = nhe3_pt.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = forw * math.exp(-8.864/T9)
    return (back /(constants.less_time(1)))


#############################################

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def dd_pt(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 3.9 * (10**8) * (T9**(-2./3)) * math.exp(-4.26*(T9**(-1./3))) * (
        + 1
        + 0.0979 * T9**(1./3)
        + 0.642 * T9**(2./3)
        + 0.440 * T9
        )
    return forw * ro_b * (1./(constants.less_time(1)))

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def pt_dd(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = dd_pt.__wrapped__(T) / constants.to_norm_time(1)
    back = 1.73 * forw * math.exp(-46.80 * (T9**(-1)))
    return back * (1./(constants.less_time(1)))

H_3.forward_rates.append(nd_tg)
H_3.backward_rates.append(tg_nd)
H_3.forward_rates.append(nhe3_pt)
H_3.backward_rates.append(pt_nhe3)
H_3.forward_rates.append(dd_pt)
H_3.backward_rates.append(pt_dd)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3


H_3.equilibrium = H_3_equ


if __name__ == '__main__':
    H_3.show_rates()