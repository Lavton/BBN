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
H_3.tr_t =  0.0017 * (constants.nu_0/constants.nu_n)
H_3.tr_T = tempreture.Tfromt(H_3.tr_t)

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def H2n_H3g(T):
    """Wagoner 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 66.2*(
        + 1.
        + 18.9 * T9)
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))    

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def H3g_H2n(T):
    """Wagoner, 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H2n_H3g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 1.63 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-72.62/T9)
    return (back /(constants.less_time(1)))

H_3.reactions.append((
    ("n", "D"), 
    ("T",), 
    H2n_H3g, 
    H3g_H2n
    ))

#############################################

def H_3_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_h3 = (3./2) * (H2n_H3g.__wrapped__(T)/H3g_H2n.__wrapped__(T))*X[0][0]*X[0][2]
    except OverflowError as e:
        X_h3 = 0
    X[0][4] = X_h3
    X[0][1] -= X_h3
    return X

#############################################


@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He3n_H3p(T):
    """Wagoner, 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 7.06 * 10**8 * (
            + 1.0 
            - 0.597 * T9**(1./2)
            + 0.183 * T9
            ) + T9**(-3./2) * 1.29 * 10**11 * math.exp(-20.61/T9)

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def H2n_H3p(T):
    """Wagoner, 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He3n_H3p.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = forw * math.exp(-8.864/T9)
    return (back /(constants.less_time(1)))

H_3.reactions.append((
    ("n", "He_3"), 
    ("p" + "T"), 
    He3n_H3p, 
    H2n_H3p
    ))

#############################################

H_3.add_interpo("H2H2_H3p","""0.001 1.45E−08 1.33E−08 1.57E−08 0.14 5.57E+05 5.26E+05 5.91E+05
0.002 6.02E−05 5.52E−05 6.54E−05 0.15 6.47E+05 6.12E+05 6.87E+05
0.003 3.27E−03 3.00E−03 3.55E−03 0.16 7.42E+05 7.02E+05 7.87E+05
0.004 4.00E−02 3.67E−02 4.34E−02 0.18 9.44E+05 8.95E+05 1.00E+06
0.005 2.36E−01 2.17E−01 2.57E−01 0.2 1.16E+06 1.10E+06 1.23E+06
0.006 9.12E−01 8.38E−01 9.89E−01 0.25 1.75E+06 1.66E+06 1.84E+06
0.007 2.67E+00 2.45E+00 2.90E+00 0.3 2.38E+06 2.26E+06 2.51E+06
0.008 6.46E+00 5.94E+00 7.00E+00 0.35 3.04E+06 2.89E+06 3.20E+06
0.009 1.36E+01 1.25E+01 1.47E+01 0.4 3.72E+06 3.54E+06 3.91E+06
0.01 2.57E+01 2.37E+01 2.79E+01 0.45 4.41E+06 4.20E+06 4.63E+06
0.011 4.49E+01 4.13E+01 4.87E+01 0.5 5.11E+06 4.87E+06 5.36E+06
0.012 7.34E+01 6.76E+01 7.96E+01 0.6 6.52E+06 6.23E+06 6.82E+06
0.013 1.14E+02 1.05E+02 1.23E+02 0.7 7.93E+06 7.61E+06 8.27E+06
0.014 1.69E+02 1.56E+02 1.83E+02 0.8 9.34E+06 8.99E+06 9.71E+06
0.015 2.42E+02 2.23E+02 2.62E+02 0.9 1.07E+07 1.04E+07 1.11E+07
0.016 3.35E+02 3.09E+02 3.63E+02 1. 1.21E+07 1.17E+07 1.25E+07
0.018 5.96E+02 5.49E+02 6.45E+02 1.25 1.54E+07 1.50E+07 1.59E+07
0.02 9.77E+02 9.01E+02 1.06E+03 1.5 1.86E+07 1.81E+07 1.90E+07
0.025 2.62E+03 2.41E+03 2.83E+03 1.75 2.15E+07 2.11E+07 2.19E+07
0.03 5.51E+03 5.10E+03 5.95E+03 2. 2.43E+07 2.38E+07 2.47E+07
0.04 1.62E+04 1.50E+04 1.74E+04 2.5 2.94E+07 2.89E+07 2.98E+07
0.05 3.46E+04 3.21E+04 3.72E+04 3. 3.40E+07 3.35E+07 3.44E+07
0.06 6.14E+04 5.71E+04 6.59E+04 3.5 3.82E+07 3.77E+07 3.86E+07
0.07 9.68E+04 9.02E+04 1.04E+05 4. 4.21E+07 4.17E+07 4.25E+07
0.08 1.41E+05 1.31E+05 1.50E+05 5. 4.91E+07 4.88E+07 4.95E+07
0.09 1.93E+05 1.80E+05 2.06E+05 6. 5.54E+07 5.51E+07 5.58E+07
0.1 2.52E+05 2.37E+05 2.69E+05 7. 6.10E+07 6.07E+07 6.14E+07
0.11 3.19E+05 3.00E+05 3.40E+05 8. 6.59E+07 6.55E+07 6.63E+07
0.12 3.92E+05 3.69E+05 4.17E+05 9. 7.00E+07 6.96E+07 7.03E+07
0.13 4.72E+05 4.45E+05 5.02E+05 10. 7.33E+07 7.29E+07 7.36E+07
    """)


@H_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
@H_3.nacreII
def H2H2_H3p(T):
    """NACRE I and II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    ro_b = univ_func.rat_scale(T)
    base_rate = 4.66 * 10**8 * T9**(-2./3) * math.exp(-4.259/T9**(1./3)) * (
        + 1.0
        + 0.759 * T9
        - 0.061 * T9**2
        + 2.78 * 10**(-3) * T9**3
        )
    return base_rate * ro_b * (1./(constants.less_time(1)))

def H3p_H2H2(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H2H2_H3p.__wrapped__(T) / constants.to_norm_time(1)
    back = 1.73 * forw * math.exp(-46.799 * (T9**(-1)))
    return back * (1./(constants.less_time(1)))

H_3.reactions.append((
    ("D", "D"),
    ("p", "T"), 
    H2H2_H3p, 
    H3p_H2H2
    ))

# H_3.forward_rates.append(nd_tg)
# H_3.backward_rates.append(tg_nd)
# H_3.forward_rates.append(nhe3_pt)
# H_3.backward_rates.append(pt_nhe3)
# H_3.forward_rates.append(dd_pt)
# H_3.backward_rates.append(pt_dd)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3


H_3.equilibrium = H_3_equ
H_3.names = ["H_3", "t", "^3H", "T"]

if __name__ == '__main__':
    # H_3.show_rates()
    import numpy as np
    from tempreture import tfromT
    grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
    grid2 = grid
    Ts = constants.less_tempreture(grid2, units="K")
    ts = np.array([tfromT(T) for T in Ts])
    for T in Ts:
        pass
        # print(dd_pt.__wrapped__(T)/H2H2_H3p.__wrapped__(T))
        # print(H2H2_H3p.__wrapped__(T)/ta.__wrapped__(T))
    # ta()

# @H_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def nd_tg(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = (
#         75.5
#         + 1250 * T9
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @H_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def tg_nd(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = nd_tg.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = 1.63 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-72.62/T9)
#     return (back /(constants.less_time(1)))

# @H_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def nhe3_pt(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = (
#         7.06 * (10**8)
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @H_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def pt_nhe3(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = nhe3_pt.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = forw * math.exp(-8.864/T9)
#     return (back /(constants.less_time(1)))

# @H_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def dd_pt(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     forw = 3.9 * (10**8) * (T9**(-2./3)) * math.exp(-4.26*(T9**(-1./3))) * (
#         + 1
#         + 0.0979 * T9**(1./3)
#         + 0.642 * T9**(2./3)
#         + 0.440 * T9
#         )
#     return forw * ro_b * (1./(constants.less_time(1)))

# @H_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def pt_dd(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = dd_pt.__wrapped__(T) / constants.to_norm_time(1)
#     back = 1.73 * forw * math.exp(-46.80 * (T9**(-1)))
#     return back * (1./(constants.less_time(1)))