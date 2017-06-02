"""
04
гелий-3, или $^{3}He$
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

He_3 = Element("He_3", 0.0)
He_3.A = 3
# from Audi et all, 2003
He_3.set_mass_excess(3016029.3191, n_N=1, p_N=2)
He_3.tr_t =  0.0011 * (constants.nu_0/constants.nu_n)
He_3.tr_T = tempreture.Tfromt(He_3.tr_t)

@He_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def d_pg_he3(T):
    """
    Vagoner
    table 1 reac 2
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 2.23 * 10**3 * T9**(-2./3) * math.exp(-3.72 * T9**(-1./3)) * (
        + 1.
        + 0.112 * T9**(1./3)
        + 3.38 ** T9**(2./3)
        + 2.65 * T9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@He_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he3_gp_d(T):
    """
    Vagoner
    t2, r2
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = d_pg_he3.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 1.63*10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-63.75/T9)
    return (back /(constants.less_time(1)))


def He_3_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_he = (3./2) * (d_pg_he3.__wrapped__(T)/he3_gp_d.__wrapped__(T))*X[0][1]*X[0][2]
    except OverflowError as e:
        X_he = 0
    X[0][3] = X_he
    return X

@He_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def dd_nhe3(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 3.9*(10**8)*ro_b*(T9**(-2./3)) * math.exp(-4.26*(T9**(-1./3))) * (
        + 1
        + 0.0979 * T9**(1./3)
        + 0.642 * T9**(2./3)
        + 0.440 * T9
        )
    return forw * (1./(constants.less_time(1)))

@He_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def nhe3_dd(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = dd_nhe3.__wrapped__(T) / constants.to_norm_time(1)
    back = 1.73 * forw * math.exp(-37.94 * (T9**(-1)))
    return back * (1./(constants.less_time(1)))


He_3.forward_rates.append(d_pg_he3)
He_3.backward_rates.append(he3_gp_d)
He_3.forward_rates.append(dd_nhe3)
He_3.backward_rates.append(nhe3_dd)


# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3

He_3.equilibrium = He_3_equ


if __name__ == '__main__':
    He_3.show_rates()