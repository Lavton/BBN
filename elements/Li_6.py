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
from elements.He_6 import He_6

Li_6 = Element("Li_6", 0.0)
Li_6.A = 6
# from Audi et all, 2003
Li_6.set_mass_excess(6015122.795, n_N=3, p_N=3)
Li_6.tr_t = 0.0077
Li_6.tr_T = tempreture.Tfromt(Li_6.tr_t)

Li_6.add_interpo("He4H2_Li6g", """
0.002 2.98E−23 1.35E−23 4.26E−23 0.15 1.13E−04 6.29E−05 1.50E−04
0.003 4.01E−20 1.83E−20 5.72E−20 0.16 1.49E−04 8.40E−05 1.99E−04
0.004 3.73E−18 1.71E−18 5.32E−18 0.18 2.46E−04 1.41E−04 3.26E−04
0.005 9.34E−17 4.29E−17 1.33E−16 0.2 3.80E−04 2.20E−04 4.99E−04
0.006 1.09E−15 5.00E−16 1.54E−15 0.25 9.10E−04 5.43E−04 1.18E−03
0.007 7.67E−15 3.55E−15 1.09E−14 0.3 1.79E−03 1.09E−03 2.30E−03
0.008 3.85E−14 1.78E−14 5.46E−14 0.35 3.09E−03 1.93E−03 3.94E−03
0.009 1.50E−13 6.98E−14 2.13E−13 0.4 4.89E−03 3.12E−03 6.16E−03
0.01 4.84E−13 2.26E−13 6.86E−13 0.45 7.24E−03 4.70E−03 9.05E−03
0.011 1.35E−12 6.30E−13 1.91E−12 0.5 1.02E−02 6.75E−03 1.27E−02
0.012 3.33E−12 1.56E−12 4.71E−12 0.6 1.84E−02 1.26E−02 2.24E−02
0.013 7.48E−12 3.52E−12 1.06E−11 0.7 3.04E−02 2.15E−02 3.65E−02
0.014 1.55E−11 7.31E−12 2.19E−11 0.8 4.78E−02 3.51E−02 5.66E−02
0.015 3.01E−11 1.42E−11 4.25E−11 0.9 7.26E−02 5.50E−02 8.48E−02
0.016 5.51E−11 2.61E−11 7.77E−11 1. 1.06E−01 8.29E−02 1.23E−01
0.018 1.61E−10 7.65E−11 2.26E−10 1.25 2.36E−01 1.93E−01 2.69E−01
0.02 4.04E−10 1.93E−10 5.68E−10 1.5 4.28E−01 3.59E−01 4.84E−01
0.025 2.55E−09 1.23E−09 3.58E−09 1.75 6.67E−01 5.64E−01 7.53E−01
0.03 1.04E−08 5.07E−09 1.46E−08 2. 9.32E−01 7.92E−01 1.05E+00
0.04 8.06E−08 3.99E−08 1.12E−07 2.5 1.49E+00 1.26E+00 1.70E+00
0.05 3.45E−07 1.74E−07 4.77E−07 3. 2.03E+00 1.71E+00 2.34E+00
0.06 1.05E−06 5.34E−07 1.44E−06 3.5 2.57E+00 2.13E+00 2.99E+00
0.07 2.54E−06 1.31E−06 3.49E−06 4. 3.10E+00 2.54E+00 3.64E+00
0.08 5.29E−06 2.77E−06 7.23E−06 5. 4.22E+00 3.39E+00 5.04E+00
0.09 9.85E−06 5.21E−06 1.34E−05 6. 5.50E+00 4.36E+00 6.63E+00
0.1 1.68E−05 8.98E−06 2.28E−05 7. 6.98E+00 5.49E+00 8.46E+00
0.11 2.69E−05 1.45E−05 3.63E−05 8. 8.67E+00 6.79E+00 1.05E+01
0.12 4.08E−05 2.22E−05 5.49E−05 9. 1.05E+01 8.23E+00 1.28E+01
0.13 5.92E−05 3.25E−05 7.94E−05 10. 1.25E+01 9.79E+00 1.53E+01
0.14 8.29E−05 4.59E−05 1.11E−04
    """)

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
@Li_6.nacreII
def He4H2_Li6g(T):
    """NACRE I and II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 14.82 * T9**(-2./3) * math.exp(-7.435 * T9**(-1./3)) * (
        + 1.0
        + 6.572 * T9 
        + 0.076 * T9**2
        + 0.0248 * T9**3
        ) + 82.8 * T9**(-3./2) * math.exp(-7.904/T9)

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def Li6g_He4H2(T):
    """NACRE"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He4H2_Li6g.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 1.53 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-17.104/T9) / (
        + 1.0
        + 2.333 * math.exp(-25.369/T9)
        )
    return (back /(constants.less_time(1)))

Li_6.reactions.append((
    ("He_4", "D"),
    ("Li_6",),
    He4H2_Li6g,
    Li6g_He4H2
    ))

##################################

def Li_6_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_li6 = (6./(3*4)) * (He4H2_Li6g.__wrapped__(T)/Li6g_He4H2.__wrapped__(T))*X[0][2]*X[0][5]
    except OverflowError as e:
        X_li6 = 0
    X[0][8] = X_li6
    X[0][1] -= X_li6
    return X

####################################

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4np_Li6g(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 4.62 * 10**(-6) * T9**(-2) * math.exp(-19.353/T9) * (
        + 1.
        + 0.75 * T9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b*ro_b /(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def Li6g_He4np(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He4np_Li6g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 7.22 * 10**19 * forw * T9**3 * (ro_b**(-2)) * math.exp(-42.933/T9)
    return (back /(constants.less_time(1)))

Li_6.reactions.append((
    ("He_4", "n", "p"),
    ("Li_6",), 
    He4np_Li6g, 
    Li6g_He4np
    ))
#########################################################

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4H3_Li6n(T):
    """Caughlan, 1988"""
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
def Li6n_He4H3(T):
    """Caughlan, 1988"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He4H3_Li6n.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 9.35 * 10**(-1) * forw * math.exp(-55.494/T9)
    return (back /(constants.less_time(1)))

Li_6.reactions.append((
    ("He_4", "T"), 
    ("Li_6", "n"),
    He4H3_Li6n,
    Li6n_He4H3
    ))
#########################################################

Li_6.add_interpo("Li6p_Be7g", """
0.001 2.23E−29 1.89E−29 3.92E−29 0.14 3.07E−01 2.67E−01 4.94E−01
0.002 5.38E−22 4.55E−22 9.45E−22 0.15 4.25E−01 3.72E−01 6.79E−01
0.003 1.91E−18 1.62E−18 3.36E−18 0.16 5.74E−01 5.03E−01 9.08E−01
0.004 3.28E−16 2.77E−16 5.75E−16 0.18 9.74E−01 8.59E−01 1.51E+00
0.005 1.26E−14 1.07E−14 2.21E−14 0.2 1.54E+00 1.36E+00 2.35E+00
0.006 2.04E−13 1.72E−13 3.57E−13 0.25 3.85E+00 3.42E+00 5.58E+00
0.007 1.87E−12 1.58E−12 3.27E−12 0.3 7.73E+00 6.82E+00 1.07E+01
0.008 1.16E−11 9.82E−12 2.03E−11 0.35 1.34E+01 1.17E+01 1.78E+01
0.009 5.44E−11 4.59E−11 9.49E−11 0.4 2.11E+01 1.82E+01 2.70E+01
0.01 2.05E−10 1.73E−10 3.57E−10 0.45 3.07E+01 2.60E+01 3.83E+01
0.011 6.52E−10 5.50E−10 1.14E−09 0.5 4.21E+01 3.53E+01 5.15E+01
0.012 1.82E−09 1.53E−09 3.16E−09 0.6 7.00E+01 5.72E+01 8.31E+01
0.013 4.53E−09 3.82E−09 7.89E−09 0.7 1.03E+02 8.30E+01 1.21E+02
0.014 1.03E−08 8.71E−09 1.80E−08 0.8 1.41E+02 1.12E+02 1.63E+02
0.015 2.18E−08 1.84E−08 3.80E−08 0.9 1.83E+02 1.43E+02 2.10E+02
0.016 4.32E−08 3.64E−08 7.51E−08 1. 2.26E+02 1.76E+02 2.60E+02
0.018 1.45E−07 1.22E−07 2.51E−07 1.25 3.43E+02 2.65E+02 3.95E+02
0.02 4.09E−07 3.44E−07 7.08E−07 1.5 4.64E+02 3.59E+02 5.40E+02
0.025 3.25E−06 2.74E−06 5.62E−06 1.75 5.85E+02 4.56E+02 6.92E+02
0.03 1.57E−05 1.32E−05 2.71E−05 2. 7.04E+02 5.53E+02 8.50E+02
0.04 1.55E−04 1.30E−04 2.65E−04 2.5 9.32E+02 7.45E+02 1.19E+03
0.05 7.80E−04 6.54E−04 1.33E−03 3. 1.14E+03 9.28E+02 1.54E+03
0.06 2.66E−03 2.24E−03 4.51E−03 3.5 1.34E+03 1.10E+03 1.90E+03
0.07 7.06E−03 5.95E−03 1.19E−02 4. 1.52E+03 1.26E+03 2.25E+03
0.08 1.58E−02 1.33E−02 2.64E−02 5. 1.84E+03 1.55E+03 2.91E+03
0.09 3.10E−02 2.63E−02 5.15E−02 6. 2.11E+03 1.80E+03 3.46E+03
0.1 5.53E−02 4.73E−02 9.15E−02 7. 2.34E+03 2.02E+03 3.91E+03
0.11 9.18E−02 7.88E−02 1.51E−01 8. 2.54E+03 2.20E+03 4.27E+03
0.12 1.44E−01 1.24E−01 2.34E−01 9. 2.71E+03 2.36E+03 4.56E+03
0.13 2.14E−01 1.85E−01 3.47E−01 10. 2.86E+03 2.50E+03 4.79E+03

    """)

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
@Li_6.nacreII
def Li6p_Be7g(T):
    """ NACRE I and II """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 1.25 * 10**6 * T9**(-2./3) * math.exp(-8.415*T9**(-1./3)) * (
        + 1.0
        - 0.252 * T9 
        + 5.19 * 10**(-2) * T9**2
        - 2.92 * 10**(-3) * T9**3
        )

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def Be7g_Li6p(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = Li6p_Be7g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 1.19 * 10**10 * T9**(3./2) * forw * ro_b**(-1) * math.exp(-65.054/T9) * (
        + 1.0 
        + 2.333 * math.exp(-25.369 / T9)
        ) / (
        + 1.0 
        + 0.5 * math.exp(-4.979 / T9)
        )
    return (back /(constants.less_time(1)))

Li_6.reactions.append((
    ("Li_6", "p"),
    ("Be_7",),
    Li6p_Be7g,
    Be7g_Li6p
    ))
#########################################################


Li_6.add_interpo("Li6p_He3He4", """
0.001 1.00E−24 8.72E−25 1.14E−24 0.14 1.26E+04 1.10E+04 1.39E+04
0.002 2.24E−17 1.95E−17 2.54E−17 0.15 1.74E+04 1.52E+04 1.92E+04
0.003 8.03E−14 6.99E−14 9.09E−14 0.16 2.33E+04 2.04E+04 2.57E+04
0.004 1.37E−11 1.19E−11 1.55E−11 0.18 3.90E+04 3.41E+04 4.29E+04
0.005 5.28E−10 4.60E−10 5.98E−10 0.2 6.05E+04 5.30E+04 6.67E+04
0.006 8.53E−09 7.42E−09 9.65E−09 0.25 1.45E+05 1.27E+05 1.60E+05
0.007 7.83E−08 6.82E−08 8.86E−08 0.3 2.78E+05 2.43E+05 3.07E+05
0.008 4.87E−07 4.24E−07 5.50E−07 0.35 4.65E+05 4.06E+05 5.15E+05
0.009 2.28E−06 1.98E−06 2.57E−06 0.4 7.07E+05 6.18E+05 7.85E+05
0.01 8.58E−06 7.47E−06 9.69E−06 0.45 1.00E+06 8.76E+05 1.12E+06
0.011 2.73E−05 2.38E−05 3.09E−05 0.5 1.35E+06 1.18E+06 1.50E+06
0.012 7.61E−05 6.63E−05 8.59E−05 0.6 2.20E+06 1.91E+06 2.44E+06
0.013 1.90E−04 1.65E−04 2.14E−04 0.7 3.22E+06 2.80E+06 3.58E+06
0.014 4.33E−04 3.77E−04 4.89E−04 0.8 4.38E+06 3.81E+06 4.87E+06
0.015 9.15E−04 7.97E−04 1.03E−03 0.9 5.67E+06 4.92E+06 6.29E+06
0.016 1.81E−03 1.58E−03 2.04E−03 1. 7.05E+06 6.11E+06 7.83E+06
0.018 6.07E−03 5.29E−03 6.84E−03 1.25 1.08E+07 9.32E+06 1.20E+07
0.02 1.71E−02 1.49E−02 1.93E−02 1.5 1.48E+07 1.28E+07 1.65E+07
0.025 1.37E−01 1.19E−01 1.54E−01 1.75 1.90E+07 1.63E+07 2.11E+07
0.03 6.61E−01 5.76E−01 7.43E−01 2. 2.33E+07 1.99E+07 2.59E+07
0.04 6.50E+00 5.68E+00 7.29E+00 2.5 3.21E+07 2.74E+07 3.59E+07
0.05 3.28E+01 2.86E+01 3.66E+01 3. 4.13E+07 3.51E+07 4.63E+07
0.06 1.12E+02 9.75E+01 1.25E+02 3.5 5.07E+07 4.31E+07 5.70E+07
0.07 2.96E+02 2.59E+02 3.30E+02 4. 6.01E+07 5.11E+07 6.78E+07
0.08 6.59E+02 5.76E+02 7.33E+02 5. 7.83E+07 6.66E+07 8.86E+07
0.09 1.29E+03 1.13E+03 1.44E+03 6. 9.48E+07 8.07E+07 1.07E+08
0.1 2.31E+03 2.02E+03 2.56E+03 7. 1.09E+08 9.30E+07 1.24E+08
0.11 3.82E+03 3.34E+03 4.22E+03 8. 1.21E+08 1.03E+08 1.37E+08
0.12 5.95E+03 5.20E+03 6.57E+03 9. 1.31E+08 1.12E+08 1.49E+08
0.13 8.83E+03 7.73E+03 9.75E+03 10. 1.39E+08 1.19E+08 1.58E+08

    """)

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
@Li_6.nacreII
def Li6p_He3He4(T):
    """NACRE I and II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 3.54 * 10**10 * T9**(-2./3) * math.exp(-8.415*T9**(-1./3)) * (
        + 1.0
        - 0.137 * T9 
        + 2.41 * 10**(-2) * T9**2
        - 1.28 * 10**(-3) * T9**3
        )

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_6.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He3He4_Li6p(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = Li6p_He3He4.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)

    back = 1.07 * forw * math.exp(-46.648/T9) * (
        + 1.0 
        + 2.333 * math.exp(-25.369/T9)
        )
    return (back /(constants.less_time(1)))

Li_6.reactions.append((
    ("Li_6", "p"),
    ("He_4", "He_3"),
    Li6p_He3He4,
    He3He4_Li6p
    ))
#########################################################


# Li_6.forward_rates.append(he4d_li6g)
# Li_6.backward_rates.append(li6g_he4d)
# Li_6.forward_rates.append(he4np_li6g)
# Li_6.backward_rates.append(li6g_he4np)
# Li_6.forward_rates.append(he4t_li6n)
# Li_6.backward_rates.append(li6n_he4t)
# Li_6.forward_rates.append(li6p_be7g)
# Li_6.backward_rates.append(be7g_li6p)
# Li_6.forward_rates.append(li6p_he3he4)
# Li_6.backward_rates.append(he3he4_li6p)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4
# 6 - Be7
# 7 - Li7
# 8 - Li6
# Li_6.equilibrium = Li_6_equ
Li_6.names = ["Li_6", "^6Li"]

if __name__ == '__main__':
    # Li_6.show_rates()
    import numpy as np
    from tempreture import tfromT
    grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
    grid2 = grid
    Ts = constants.less_tempreture(grid2, units="K")
    ts = np.array([tfromT(T) for T in Ts])
    for T in np.copy(Ts):
        pass
        print(li6p_he3he4.__wrapped__(T)/Li6p_He3He4.__wrapped__(T))
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from tempreture import tfromT
    plt.cla()
    plt.clf()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{tempreture} (MeV)')
    plt.ylabel(r'\textbf{\lambda}')
    plt.plot(ts, [li6p_he3he4(T) for T in np.copy(Ts)], label='forw1')
    plt.plot(ts, [Li6p_He3He4(T) for T in np.copy(Ts)], label='forw2')
    plt.legend()
    plt.show()
    plt.cla()
    plt.clf()
    

# @Li_6.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4d_li6g(T):
#     """
#     Caughlan, 1988
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = 3.01 * 10**1 * T9**(-2./3) * math.exp(-7.423 * T9**(-1./3)) * (
#         + 1.0
#         + 0.056 * T9**(1./3)
#         - 4.85 * T9**(2./3)
#         + 8.85 * T9 
#         - 0.585 * T9**(4./3)
#         - 0.584 * T9**(5./3)
#         ) + 8.55 * 10**1 * T9**(-3./2) * math.exp(-8.228/T9)
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @Li_6.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def li6g_he4d(T):
#     """Caughlan, 1988"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = he4d_li6g.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     back = 1.53 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-17.118/T9)
#     return (back /(constants.less_time(1)))

# @Li_6.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def li6p_be7g(T):
#     """
#     Caughlan, 1988
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     T9A = T9/(
#         + 1.
#         - 9.69 * 10**(-2) * T9
#         + 2.84 * 10**(-2) * T9**(5./3) / (
#             + 1.
#             - 9.69 * 10**(-2) * T9
#             )**(2./3)
#         )
#     base_rate = 6.69 * 10**5 * (T9A**(5./6)/T9**(3./2)) * math.exp(-8.413/T9A**(1./3))
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @Li_6.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def be7g_li6p(T):
#     """Caughlan, 1988"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = li6p_be7g.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = 1.19 * 10**(10) * forw * ro_b**(-1) * math.exp(-65.054/T9)
#     return (back /(constants.less_time(1)))

# @Li_6.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def li6p_he3he4(T):
#     """
#     Caughlan, 1988
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = 3.37 * 10**10 * T9**(-2./3) * math.exp(-8.413/T9**(1./3) - (T9/5.50)**2) * (
#         + 1.0
#         + 0.50 * T9**(1./3)
#         - 0.061  * T9**(2./3)
#         - 0.021 * T9
#         + 0.006 * T9**(4./3)
#         + 0.005 * T9**(5./3)
#         ) + 1.33 * 10**10 * T9**(-3./2) * math.exp(-17.764/T9) + 1.29 * 10**9 * T9**(-1) * math.exp(-21.820/T9)
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @Li_6.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he3he4_li6p(T):
#     """Caughlan, 1988"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = li6p_he3he4.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = 1.07 * forw * math.exp(-46.631/T9)
#     return (back /(constants.less_time(1)))
