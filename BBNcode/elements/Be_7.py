"""
07
берилий-7, или $^{7}Be$
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

Be_7 = Element("Be_7", 0.0)
Be_7.A = 7
# from Audi et all, 2003
Be_7.set_mass_excess(7016929.83, n_N=3, p_N=4)
Be_7.tr_t = 0.02
Be_7.tr_T = tempreture.Tfromt(Be_7.tr_t)

@Be_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he3he4_be7g(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 4.8 * 10**6 * T9**(-2./3) * math.exp(-12.8 * T9**(-1./3)) * (
        + 1.0
        + 0.0326 * T9**(1./3)
        - 0.219 * T9**(2./3)
        - 0.0499 * T9
        + 0.0258 * T9**(4./3)
        + 0.0150 * T9**(5./3)
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Be_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def be7g_he3he4(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = he3he4_be7g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    back = 1.12 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-18.42/T9)
    return (back /(constants.less_time(1)))

##################################

def Be_7_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_be7 = (7./(3*4)) * (he3he4_be7g.__wrapped__(T)/be7g_he3he4.__wrapped__(T))*X[0][3]*X[0][5]
    except OverflowError as e:
        X_be7 = 0
    X[0][6] = X_be7
    return X

####################################

@Be_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def nBe7_He4He4(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = (
        1.2 * 10**7 * T9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Be_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4He4_nBe7(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = nBe7_He4He4.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 4.64 * forw * math.exp(-220.4/T9)
    return (back /(constants.less_time(1)))

#########################################################

Be_7.forward_rates.append(he3he4_be7g)
Be_7.backward_rates.append(be7g_he3he4)
Be_7.forward_rates.append(nBe7_He4He4)
Be_7.backward_rates.append(He4He4_nBe7)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4
# 6 - Be7

Be_7.equilibrium = Be_7_equ


if __name__ == '__main__':
    Be_7.show_rates()
    