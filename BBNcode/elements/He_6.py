"""
10
Гелий-6, или $^{6}He$
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

He_6 = Element("He_6", 0.0)
He_6.A = 7
# from Audi et all, 2003
He_6.set_mass_excess(6018889.1, n_N=3, p_N=3)
He_6.tr_t = 0.025
He_6.tr_T = tempreture.Tfromt(He_6.tr_t)

@He_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4nn_he6_g(T):
    """
    Caughlan, 1988
    """
    k = 1000.0
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = k * 4.04 * 10**(-11) * T9**(-2.) * math.exp(-9.585 / T9) * (
        + 1.0
        + 0.138 * T9 
        )
    ro_b = univ_func.rat_scale(T)
    # return 0
    return base_rate * ro_b*ro_b/(constants.less_time(1))

@He_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he6g_he4nn(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = he4nn_he6_g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    back = 1.08 * 10**20 * forw * (ro_b**(-2)) * T9**(3) * math.exp(-11.319/T9)
    return (back /(constants.less_time(1)))

##################################

def He_6_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_he6 = (6./4) * (he4nn_he6_g.__wrapped__(T)/he6g_he4nn.__wrapped__(T))*X[0][2]*X[0][5]
    except OverflowError as e:
        X_he6 = 0
    X[0][9] = X_he6
    return X

####################################

@He_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he6_li6(T):
    """
    http://periodictable.com/Isotopes/002.6/index3.full.html
    """
    return (1/0.001)/(constants.less_time(1))

@He_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def li6_he6(T):
    return (0 /(constants.less_time(1)))

#########################################################

He_6.forward_rates.append(he4nn_he6_g)
He_6.backward_rates.append(he6g_he4nn)
He_6.forward_rates.append(he6_li6)
He_6.backward_rates.append(li6_he6)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4
# 6 - Be7
# 7 - Li7
# 8 - Li6
# 9 - He6
He_6.equilibrium = He_6_equ

if __name__ == '__main__':
    He_6.show_rates()
    