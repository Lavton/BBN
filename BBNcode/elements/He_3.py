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

# @He_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def d_pg_he3(T):
#     """
#     Vagoner
#     table 1 reac 2
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = 2.23 * 10**3 * T9**(-2./3) * math.exp(-3.72 * T9**(-1./3)) * (
#         + 1.
#         + 0.112 * T9**(1./3)
#         + 3.38 ** T9**(2./3)
#         + 2.65 * T9
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @He_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he3_gp_d(T):
#     """
#     Vagoner
#     t2, r2
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = d_pg_he3.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = 1.63*10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-63.75/T9)
#     return (back /(constants.less_time(1)))

@He_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
@He_3.nacreII
def H2p_He3g(T):
    """NACRE"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    if T9 <= 0.11:
        base_rate = 1.81 * 10**3 * T9**(-2./3) * math.exp(-3.721/T9**(1./3)) * (
            + 1.
            + 14.3 * T9
            - 90.5 * T9**2
            + 395 * T9**3
            )
    else:
        base_rate = 2.58 * 10**3 * T9**(-2./3) * math.exp(-3.721/T9**(1./3)) * (
            + 1.
            + 3.96 * T9
            + 0.116 * T9**2
            )

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

He_3.add_interpo("H2p_He3g",
"""
0.001 1.35E−11 1.09E−11 1.62E−11 0.14 1.15E+01 1.03E+01 1.27E+01
0.002 1.87E−08 1.51E−08 2.25E−08 0.15 1.33E+01 1.20E+01 1.47E+01
0.003 6.07E−07 4.93E−07 7.27E−07 0.16 1.52E+01 1.37E+01 1.67E+01
0.004 5.38E−06 4.38E−06 6.43E−06 0.18 1.93E+01 1.75E+01 2.12E+01
0.005 2.52E−05 2.06E−05 3.01E−05 0.2 2.37E+01 2.16E+01 2.60E+01
0.006 8.15E−05 6.67E−05 9.71E−05 0.25 3.62E+01 3.31E+01 3.94E+01
0.007 2.07E−04 1.70E−04 2.47E−04 0.3 5.01E+01 4.59E+01 5.44E+01
0.008 4.47E−04 3.68E−04 5.31E−04 0.35 6.52E+01 5.98E+01 7.08E+01
0.009 8.55E−04 7.04E−04 1.01E−03 0.4 8.13E+01 7.45E+01 8.84E+01
0.01 1.49E−03 1.23E−03 1.77E−03 0.45 9.84E+01 9.00E+01 1.07E+02
0.011 2.43E−03 2.01E−03 2.87E−03 0.5 1.16E+02 1.06E+02 1.27E+02
0.012 3.73E−03 3.09E−03 4.40E−03 0.6 1.54E+02 1.40E+02 1.69E+02
0.013 5.47E−03 4.54E−03 6.45E−03 0.7 1.95E+02 1.77E+02 2.14E+02
0.014 7.73E−03 6.42E−03 9.10E−03 0.8 2.38E+02 2.16E+02 2.61E+02
0.015 1.06E−02 8.80E−03 1.24E−02 0.9 2.84E+02 2.57E+02 3.12E+02
0.016 1.41E−02 1.17E−02 1.65E−02 1. 3.32E+02 2.99E+02 3.64E+02
0.018 2.33E−02 1.95E−02 2.74E−02 1.25 4.57E+02 4.12E+02 5.03E+02
0.02 3.60E−02 3.02E−02 4.22E−02 1.5 5.91E+02 5.31E+02 6.52E+02
0.025 8.59E−02 7.25E−02 1.00E−01 1.75 7.30E+02 6.54E+02 8.07E+02
0.03 1.66E−01 1.41E−01 1.93E−01 2. 8.73E+02 7.79E+02 9.68E+02
0.04 4.35E−01 3.73E−01 5.01E−01 2.5 1.17E+03 1.04E+03 1.30E+03
0.05 8.62E−01 7.44E−01 9.86E−01 3. 1.46E+03 1.30E+03 1.63E+03
0.06 1.45E+00 1.26E+00 1.65E+00 3.5 1.76E+03 1.56E+03 1.97E+03
0.07 2.20E+00 1.93E+00 2.50E+00 4. 2.07E+03 1.84E+03 2.30E+03
0.08 3.11E+00 2.74E+00 3.51E+00 5. 2.67E+03 2.39E+03 2.95E+03
0.09 4.17E+00 3.69E+00 4.69E+00 6. 3.27E+03 2.95E+03 3.60E+03
0.1 5.37E+00 4.77E+00 6.02E+00 7. 3.86E+03 3.50E+03 4.23E+03
0.11 6.71E+00 5.98E+00 7.49E+00 8. 4.45E+03 4.05E+03 4.86E+03
0.12 8.18E+00 7.31E+00 9.10E+00 9. 5.02E+03 4.58E+03 5.48E+03
0.13 9.77E+00 8.76E+00 1.08E+01 10. 5.59E+03 5.10E+03 6.08E+03
""")

@He_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He3g_H2p(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H2p_He3g.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)

    back = 1.63 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-63.752/T9)
    return (back /(constants.less_time(1)))


def He_3_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_he = (3./2) * (H2p_He3g.__wrapped__(T)/He3g_H2p.__wrapped__(T))*X[0][1]*X[0][2]
    except OverflowError as e:
        X_he = 0
    X[0][3] = X_he
    X[0][1] -= X_he
    return X

# @He_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def dd_nhe3(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     forw = 3.9*(10**8)*ro_b*(T9**(-2./3)) * math.exp(-4.26*(T9**(-1./3))) * (
#         + 1
#         + 0.0979 * T9**(1./3)
#         + 0.642 * T9**(2./3)
#         + 0.440 * T9
#         )
#     return forw * (1./(constants.less_time(1)))

He_3.add_interpo("H2H2_He3n", 
"""
0.001 1.43E−08 1.28E−08 1.58E−08 0.14 5.48E+05 5.01E+05 5.92E+05
0.002 5.95E−05 5.31E−05 6.57E−05 0.15 6.38E+05 5.83E+05 6.88E+05
0.003 3.23E−03 2.89E−03 3.57E−03 0.16 7.32E+05 6.70E+05 7.89E+05
0.004 3.95E−02 3.53E−02 4.36E−02 0.18 9.36E+05 8.59E+05 1.01E+06
0.005 2.34E−01 2.09E−01 2.58E−01 0.2 1.16E+06 1.06E+06 1.24E+06
0.006 9.01E−01 8.05E−01 9.95E−01 0.25 1.77E+06 1.64E+06 1.90E+06
0.007 2.64E+00 2.36E+00 2.91E+00 0.3 2.46E+06 2.28E+06 2.63E+06
0.008 6.38E+00 5.70E+00 7.04E+00 0.35 3.20E+06 2.97E+06 3.41E+06
0.009 1.34E+01 1.20E+01 1.48E+01 0.4 3.98E+06 3.70E+06 4.23E+06
0.01 2.54E+01 2.27E+01 2.80E+01 0.45 4.78E+06 4.46E+06 5.08E+06
0.011 4.44E+01 3.97E+01 4.89E+01 0.5 5.59E+06 5.23E+06 5.93E+06
0.012 7.25E+01 6.49E+01 7.99E+01 0.6 7.25E+06 6.80E+06 7.67E+06
0.013 1.12E+02 1.01E+02 1.24E+02 0.7 8.90E+06 8.38E+06 9.40E+06
0.014 1.67E+02 1.49E+02 1.84E+02 0.8 1.05E+07 9.95E+06 1.11E+07
0.015 2.39E+02 2.14E+02 2.63E+02 0.9 1.22E+07 1.15E+07 1.28E+07
0.016 3.31E+02 2.96E+02 3.64E+02 1. 1.37E+07 1.30E+07 1.45E+07
0.018 5.88E+02 5.27E+02 6.48E+02 1.25 1.76E+07 1.66E+07 1.85E+07
0.02 9.64E+02 8.64E+02 1.06E+03 1.5 2.12E+07 2.01E+07 2.23E+07
0.025 2.58E+03 2.31E+03 2.84E+03 1.75 2.46E+07 2.34E+07 2.59E+07
0.03 5.44E+03 4.88E+03 5.97E+03 2. 2.79E+07 2.65E+07 2.93E+07
0.04 1.60E+04 1.44E+04 1.75E+04 2.5 3.39E+07 3.23E+07 3.55E+07
0.05 3.40E+04 3.07E+04 3.73E+04 3. 3.93E+07 3.76E+07 4.10E+07
0.06 6.04E+04 5.45E+04 6.60E+04 3.5 4.43E+07 4.25E+07 4.60E+07
0.07 9.52E+04 8.61E+04 1.04E+05 4. 4.87E+07 4.70E+07 5.05E+07
0.08 1.38E+05 1.25E+05 1.51E+05 5. 5.65E+07 5.48E+07 5.81E+07
0.09 1.89E+05 1.72E+05 2.06E+05 6. 6.26E+07 6.10E+07 6.42E+07
0.1 2.48E+05 2.25E+05 2.69E+05 7. 6.73E+07 6.58E+07 6.88E+07
0.11 3.13E+05 2.85E+05 3.40E+05 8. 7.27E+07 7.11E+07 7.44E+07
0.12 3.86E+05 3.51E+05 4.17E+05 9. 7.72E+07 7.54E+07 7.89E+07
0.13 4.64E+05 4.23E+05 5.02E+05 10. 8.13E+07 8.03E+07 8.40E+07
"""
    )

@He_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
@He_3.nacreII
def H2H2_He3n(T):
    """NACRE"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 4.67 * 10**8 * T9**(-2./3) * math.exp(-4.259/T9**(1./3)) * (
            + 1.
            + 1.079 * T9
            - 0.1124 * T9**2
            + 5.68 * 10**(-3) * T9**3
        )

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

# @He_3.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def nhe3_dd(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = dd_nhe3.__wrapped__(T) / constants.to_norm_time(1)
#     back = 1.73 * forw * math.exp(-37.94 * (T9**(-1)))
#     return back * (1./(constants.less_time(1)))

@He_3.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He3n_H2H2(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H2H2_He3n.__wrapped__(T) / constants.to_norm_time(1)
    back = 1.73 * forw * math.exp(-37.936 * (T9**(-1)))
    return back * (1./(constants.less_time(1)))


# He_3.forward_rates.append(d_pg_he3)
# He_3.backward_rates.append(he3_gp_d)
# He_3.forward_rates.append(dd_nhe3)
# He_3.backward_rates.append(nhe3_dd)


# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3

He_3.equilibrium = He_3_equ


if __name__ == '__main__':
    # He_3.show_rates()
    import numpy as np
    from tempreture import tfromT
    grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
    grid2 = grid
    Ts = constants.less_tempreture(grid2, units="K")
    ts = np.array([tfromT(T) for T in Ts])
    for T in Ts:
        pass
        # print(H2p_He3g.__wrapped__(T)/H2p_He3g1.__wrapped__(T))
