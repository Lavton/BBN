"""
дейтерий, или $^{2}H$
"""

if __name__ == '__main__':
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)))

import constants
import math
import tempreture
from elements.Element import Element
import functools
import univ_func

H_2 = Element("H_2", 0.0)
H_2.A = 2
# from Audi et all, 2003
H_2.mass_excess = constants.less_tempreture(2161062.7, units="eV")


# change from equilibrium state
H_2.tr_t = 0.0015
H_2.tr_T = tempreture.Tfromt(H_2.tr_t)

@H_2.equilib_zeroize
@functools.lru_cache(maxsize=8)
def H_2_forw_rate(T):
    """
    Smith et all
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 4.742 * 10**4 * (
        + 1. 
        - 0.8540 * T9**(1./2)
        + 0.4895 * T9
        - 0.09623 * T9**(3./2)
        + 8.471*1e-3 * T9**2
        - 2.80*1e-4 * T9**(5./2)
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@H_2.equilib_zeroize
@functools.lru_cache(maxsize=8)
def H_2_backward_rate(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H_2_forw_rate.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    back = forw * math.exp(-H_2.mass_excess/T)
    ro_b = univ_func.rat_scale(T)

    back = 4.68*10**9 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-H_2.mass_excess/T)
    return (back /(constants.less_time(1)))


def H_2_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_n = (2./1) * (H_2_forw_rate.__wrapped__(T)/H_2_backward_rate.__wrapped__(T))*X[0][0]*X[0][1]
    except OverflowError as e:
        X_n = 0
    X[0][2] = X_n
    return X


H_2.forward_rates.append(H_2_forw_rate)
H_2.backward_rates.append(H_2_backward_rate)


H_2.equilibrium = H_2_equ


if __name__ == '__main__':
    H_2.show_rates()
