"""
08
Литий-7, или $^{7}Li$
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

Li_7 = Element("Li_7", 0.0)
Li_7.A = 7
# from Audi et all, 2003
Li_7.set_mass_excess(7016004.55, n_N=4, p_N=3)
Li_7.tr_t = 0.030
Li_7.tr_T = tempreture.Tfromt(Li_7.tr_t)

Li_7.add_interpo("H3He4_Li7g", """
0.002 7.06E−21 6.52E−21 7.88E−21 0.15 6.64E−01 6.18E−01 7.38E−01
0.003 1.79E−17 1.65E−17 2.00E−17 0.16 8.76E−01 8.15E−01 9.73E−01
0.004 2.48E−15 2.29E−15 2.77E−15 0.18 1.43E+00 1.33E+00 1.58E+00
0.005 8.21E−14 7.59E−14 9.17E−14 0.2 2.16E+00 2.01E+00 2.40E+00
0.006 1.18E−12 1.09E−12 1.32E−12 0.25 4.93E+00 4.60E+00 5.46E+00
0.007 9.87E−12 9.12E−12 1.10E−11 0.3 9.14E+00 8.54E+00 1.01E+01
0.008 5.68E−11 5.25E−11 6.34E−11 0.35 1.49E+01 1.39E+01 1.65E+01
0.009 2.49E−10 2.30E−10 2.77E−10 0.4 2.21E+01 2.06E+01 2.44E+01
0.01 8.85E−10 8.18E−10 9.88E−10 0.45 3.07E+01 2.86E+01 3.40E+01
0.011 2.68E−09 2.48E−09 2.99E−09 0.5 4.06E+01 3.79E+01 4.51E+01
0.012 7.15E−09 6.61E−09 7.98E−09 0.6 6.39E+01 5.93E+01 7.13E+01
0.013 1.72E−08 1.59E−08 1.92E−08 0.7 9.12E+01 8.41E+01 1.02E+02
0.014 3.78E−08 3.49E−08 4.22E−08 0.8 1.22E+02 1.11E+02 1.37E+02
0.015 7.73E−08 7.15E−08 8.62E−08 0.9 1.55E+02 1.41E+02 1.75E+02
0.016 1.49E−07 1.37E−07 1.66E−07 1. 1.90E+02 1.72E+02 2.16E+02
0.018 4.72E−07 4.36E−07 5.26E−07 1.25 2.85E+02 2.54E+02 3.26E+02
0.02 1.27E−06 1.18E−06 1.42E−06 1.5 3.87E+02 3.42E+02 4.44E+02
0.025 9.25E−06 8.56E−06 1.03E−05 1.75 4.92E+02 4.32E+02 5.65E+02
0.03 4.17E−05 3.86E−05 4.65E−05 2. 5.99E+02 5.24E+02 6.86E+02
0.04 3.70E−04 3.42E−04 4.12E−04 2.5 8.15E+02 7.13E+02 9.28E+02
0.05 1.73E−03 1.60E−03 1.92E−03 3. 1.03E+03 9.04E+02 1.17E+03
0.06 5.55E−03 5.14E−03 6.18E−03 3.5 1.25E+03 1.11E+03 1.40E+03
0.07 1.40E−02 1.30E−02 1.56E−02 4. 1.46E+03 1.28E+03 1.62E+03
0.08 3.00E−02 2.78E−02 3.34E−02 5. 1.86E+03 1.63E+03 2.04E+03
0.09 5.68E−02 5.27E−02 6.32E−02 6. 2.27E+03 1.99E+03 2.42E+03
0.1 9.83E−02 9.12E−02 1.09E−01 7. 2.51E+03 2.19E+03 2.73E+03
0.11 1.58E−01 1.47E−01 1.76E−01 8. 2.74E+03 2.38E+03 2.96E+03
0.12 2.41E−01 2.24E−01 2.68E−01 9. 2.91E+03 2.51E+03 3.13E+03
0.13 3.50E−01 3.26E−01 3.90E−01 10. 3.02E+03 2.61E+03 3.25E+03
0.14 4.90E−01 4.56E−01 5.45E−01

    """)

@Li_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
@Li_7.nacreII
def H3He4_Li7g(T):
    """NACRE I and II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 8.20 * 10**5 * T9**(-2./3) * math.exp(-8.081/T9**(1./3)) * (
        + 1.0
        - 0.389 * T9 
        + 0.134 * T9**2
        - 1.81 * 10**(-2) * T9**3
        + 9.23 * 10**(-4) * T9**4
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def Li7g_H3He4(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H3He4_Li7g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    back = 1.11 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-28.625/T9) / (
        + 1.0
        + 0.5 * math.exp(-5.543/T9)
        )
    return (back /(constants.less_time(1)))

Li_7.reactions.append((
    ("T", "He_4"),
    ("Li_7",),
    H3He4_Li7g,
    Li7g_H3He4
    ))

##################################

def Li_7_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_li7 = (7./(3*4)) * (H3He4_Li7g.__wrapped__(T)/Li7g_H3He4.__wrapped__(T))*X[0][4]*X[0][5]
    except OverflowError as e:
        X_li7 = 0
    X[0][7] = X_li7
    X[0][1] -= X_li7
    return X

####################################

@Li_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def Be7n_Li7p(T):
    """Shepiro"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = (
        + 6.8423 * 10**9 
        - 1.4988 * 10**10 * T9**(1./2)
        + 1.76749 * 10**10 * T9 
        - 1.05769 * 10**10 * T9**(3./2)
        + 2.6622 * 10**9 * T9**2
        + 2.74476 * 10**8 * T9**(5./2)
        - 3.35616 * 10**8 * T9**3
        + 7.64252 * 10**7 * T9**(7./2)
        - 5.93091 * 10**6 * T9**4
        - 2.28294 * 10**7 * math.exp(-0.0503518/T9) * T9**(-3./2)
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def Li7p_Be7n(T):
    """Caughlan"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = Be7n_Li7p.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 9.98 * 10**(-1) * forw * math.exp(-19.81/T9)
    return (back /(constants.less_time(1)))

Li_7.reactions.append((
    ("n", "Be_7"),
    ("p", "Li_7"),
    Be7n_Li7p,
    Li7p_Be7n
    ))

#########################################################


Li_7.add_interpo("Li7p_He4He4", """
0.001 9.17E−27 7.79E−27 1.12E−26 0.14 2.60E+02 2.32E+02 2.87E+02
0.002 2.38E−19 2.02E−19 2.90E−19 0.15 3.62E+02 3.25E+02 4.00E+02
0.003 9.06E−16 7.71E−16 1.10E−15 0.16 4.92E+02 4.42E+02 5.42E+02
0.004 1.61E−13 1.37E−13 1.96E−13 0.18 8.42E+02 7.60E+02 9.25E+02
0.005 6.41E−12 5.45E−12 7.78E−12 0.2 1.34E+03 1.21E+03 1.46E+03
0.006 1.06E−10 9.02E−11 1.29E−10 0.25 3.35E+03 3.05E+03 3.66E+03
0.007 9.93E−10 8.46E−10 1.20E−09 0.3 6.73E+03 6.14E+03 7.31E+03
0.008 6.28E−09 5.36E−09 7.61E−09 0.35 1.17E+04 1.07E+04 1.27E+04
0.009 2.98E−08 2.55E−08 3.61E−08 0.4 1.84E+04 1.68E+04 1.99E+04
0.01 1.14E−07 9.74E−08 1.38E−07 0.45 2.68E+04 2.46E+04 2.91E+04
0.011 3.68E−07 3.15E−07 4.45E−07 0.5 3.71E+04 3.40E+04 4.02E+04
0.012 1.04E−06 8.89E−07 1.25E−06 0.6 6.31E+04 5.78E+04 6.83E+04
0.013 2.62E−06 2.25E−06 3.16E−06 0.7 9.58E+04 8.79E+04 1.04E+05
0.014 6.04E−06 5.18E−06 7.27E−06 0.8 1.35E+05 1.24E+05 1.46E+05
0.015 1.29E−05 1.11E−05 1.55E−05 0.9 1.79E+05 1.65E+05 1.94E+05
0.016 2.58E−05 2.22E−05 3.10E−05 1. 2.29E+05 2.10E+05 2.48E+05
0.018 8.80E−05 7.56E−05 1.05E−04 1.25 3.72E+05 3.41E+05 4.03E+05
0.02 2.53E−04 2.17E−04 3.01E−04 1.5 5.39E+05 4.93E+05 5.84E+05
0.025 2.09E−03 1.80E−03 2.47E−03 1.75 7.25E+05 6.63E+05 7.88E+05
0.03 1.04E−02 8.98E−03 1.23E−02 2. 9.31E+05 8.49E+05 1.01E+06
0.04 1.08E−01 9.30E−02 1.26E−01 2.5 1.40E+06 1.28E+06 1.53E+06
0.05 5.63E−01 4.88E−01 6.52E−01 3. 1.98E+06 1.79E+06 2.17E+06
0.06 1.98E+00 1.72E+00 2.28E+00 3.5 2.68E+06 2.42E+06 2.94E+06
0.07 5.40E+00 4.71E+00 6.17E+00 4. 3.51E+06 3.16E+06 3.85E+06
0.08 1.23E+01 1.08E+01 1.40E+01 5. 5.52E+06 4.97E+06 6.08E+06
0.09 2.47E+01 2.16E+01 2.78E+01 6. 7.86E+06 7.08E+06 8.65E+06
0.1 4.48E+01 3.94E+01 5.03E+01 7. 1.03E+07 9.29E+06 1.13E+07
0.11 7.53E+01 6.66E+01 8.42E+01 8. 1.27E+07 1.14E+07 1.40E+07
0.12 1.19E+02 1.06E+02 1.33E+02 9. 1.49E+07 1.35E+07 1.64E+07
0.13 1.79E+02 1.60E+02 1.99E+02 10. 1.70E+07 1.53E+07 1.86E+07

    """)

@Li_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
@Li_7.nacreII
def Li7p_He4He4(T):
    """NACRE I and II """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 7.20 * 10**8 * T9**(-2./3) * math.exp(-8.473*T9**(-1./3) - (T9/6.5)**2) * (
        + 1.0
        + 1.05 * T9 
        - 0.653 * T9**2
        + 0.185 * T9**3
        - 2.12 * 10**(-2) * T9**4
        + 9.30 * 10**(-4) * T9**5
        ) + 9.85 * 10**6 * T9**0.576 * math.exp(-10.415/T9)

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Li_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4He4_Li7p(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = Li7p_He4He4.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)

    back = 4.69 * forw * math.exp(-201.32/T9) / (
        + 1.0 
        + 0.5 * math.exp(-5.543/T9)
        )
    return (back /(constants.less_time(1)))

Li_7.reactions.append((
    ("p", "Li_7"),
    ("He_4", "He_4"),
    Li7p_He4He4,
    He4He4_Li7p
    ))
#########################################################


# Li_7.forward_rates.append(the4_li7g)
# Li_7.backward_rates.append(li7g_the4)
# Li_7.forward_rates.append(nbe7_pli7)
# Li_7.backward_rates.append(pli7_nbe7)
# Li_7.forward_rates.append(pli7_he4he4)
# Li_7.backward_rates.append(he4he4_pli7)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4
# 6 - Be7
# 7 - Li7
# Li_7.equilibrium = Li_7_equ
Li_7.names = ["Li_7", "^7Li"]

if __name__ == '__main__':
    # Li_7.show_rates()
    import numpy as np
    from tempreture import tfromT
    grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
    grid2 = grid
    Ts = constants.less_tempreture(grid2, units="K")
    ts = np.array([tfromT(T) for T in Ts])
    for T in np.copy(Ts):
        pass
        print(Be7n_Li7p.__wrapped__(T)/Be7n_Li7p1.__wrapped__(T))
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
    plt.plot(ts, [Be7n_Li7p(T) for T in np.copy(Ts)], label='forw1')
    plt.plot(ts, [Be7n_Li7p1(T) for T in np.copy(Ts)], label='forw2')
    plt.legend()
    plt.show()
    plt.cla()
    plt.clf()
    
    

# @Li_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def the4_li7g(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = 5.28 * 10**5 * T9**(-2./3) * math.exp(-8.08 * T9**(-1./3)) * (
#         + 1.0
#         + 0.0516 * T9**(1./3)
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @Li_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def li7g_the4(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = the4_li7g.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     back = 1.12 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-28.63/T9)
#     return (back /(constants.less_time(1)))

# @Li_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def nbe7_pli7(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = (
#         6.74 * 10**9
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @Li_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def pli7_nbe7(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = nbe7_pli7.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = forw * math.exp(-19.07/T9)
#     return (back /(constants.less_time(1)))

# @Li_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def pli7_he4he4(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = 1.42 * 10**9 * T9**(-2./3) * math.exp(-8.47*T9**(-1./3)) * (
#         + 1.
#         + 0.0493 * T9**(1./3)
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @Li_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4he4_pli7(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = pli7_he4he4.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = 4.64 * forw * math.exp(-201.3/T9)
#     return (back /(constants.less_time(1)))
