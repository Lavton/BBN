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
Be_7.tr_t = 0.03
Be_7.tr_T = tempreture.Tfromt(Be_7.tr_t)

# @Be_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he3he4_be7g(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = 4.8 * 10**6 * T9**(-2./3) * math.exp(-12.8 * T9**(-1./3)) * (
#         + 1.0
#         + 0.0326 * T9**(1./3)
#         - 0.219 * T9**(2./3)
#         - 0.0499 * T9
#         + 0.0258 * T9**(4./3)
#         + 0.0150 * T9**(5./3)
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @Be_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def be7g_he3he4(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = he3he4_be7g.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     back = 1.12 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-18.42/T9)
#     return (back /(constants.less_time(1)))

Be_7.add_interpo("He3He4_Be7g", """
0.005 5.14E−25 4.57E−25 5.62E−25 0.16 9.98E−04 8.93E−04 1.08E−03
0.006 3.79E−23 3.37E−23 4.14E−23 0.18 2.28E−03 2.04E−03 2.47E−03
0.007 1.17E−21 1.04E−21 1.28E−21 0.2 4.62E−03 4.14E−03 5.02E−03
0.008 1.99E−20 1.77E−20 2.17E−20 0.25 1.89E−02 1.70E−02 2.05E−02
0.009 2.17E−19 1.93E−19 2.37E−19 0.3 5.49E−02 4.93E−02 5.94E−02
0.01 1.70E−18 1.51E−18 1.86E−18 0.35 1.28E−01 1.15E−01 1.38E−01
0.011 1.03E−17 9.13E−18 1.12E−17 0.4 2.54E−01 2.29E−01 2.74E−01
0.012 5.03E−17 4.48E−17 5.50E−17 0.45 4.54E−01 4.09E−01 4.89E−01
0.013 2.09E−16 1.86E−16 2.28E−16 0.5 7.44E−01 6.72E−01 8.02E−01
0.014 7.51E−16 6.68E−16 8.20E−16 0.6 1.67E+00 1.51E+00 1.80E+00
0.015 2.40E−15 2.14E−15 2.63E−15 0.7 3.16E+00 2.86E+00 3.39E+00
0.016 6.96E−15 6.20E−15 7.60E−15 0.8 5.30E+00 4.79E+00 5.70E+00
0.018 4.56E−14 4.06E−14 4.98E−14 0.9 8.17E+00 7.38E+00 8.79E+00
0.02 2.30E−13 2.05E−13 2.51E−13 1. 1.18E+01 1.07E+01 1.27E+01
0.025 5.86E−12 5.21E−12 6.39E−12 1.25 2.44E+01 2.19E+01 2.64E+01
0.03 6.88E−11 6.12E−11 7.50E−11 1.5 4.19E+01 3.73E+01 4.57E+01
0.04 2.46E−09 2.19E−09 2.69E−09 1.75 6.40E+01 5.66E+01 7.02E+01
0.05 3.11E−08 2.77E−08 3.39E−08 2. 9.04E+01 7.93E+01 9.98E+01
0.06 2.14E−07 1.91E−07 2.33E−07 2.5 1.54E+02 1.34E+02 1.72E+02
0.07 9.92E−07 8.85E−07 1.08E−06 3. 2.32E+02 1.98E+02 2.61E+02
0.08 3.50E−06 3.13E−06 3.82E−06 3.5 3.21E+02 2.73E+02 3.64E+02
0.09 1.02E−05 9.06E−06 1.11E−05 4. 4.21E+02 3.56E+02 4.79E+02
0.1 2.53E−05 2.26E−05 2.76E−05 5. 6.53E+02 5.49E+02 7.44E+02
0.11 5.62E−05 5.02E−05 6.11E−05 6. 9.15E+02 7.71E+02 1.04E+03
0.12 1.14E−04 1.02E−04 1.24E−04 7. 1.19E+03 1.01E+03 1.36E+03
0.13 2.13E−04 1.90E−04 2.31E−04 8. 1.47E+03 1.26E+03 1.67E+03
0.14 3.75E−04 3.35E−04 4.07E−04 9. 1.74E+03 1.49E+03 1.96E+03
0.15 6.25E−04 5.59E−04 6.80E−04 10. 1.99E+03 1.71E+03 2.24E+03

    """)

@Be_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
@Be_7.nacreII
def He3He4_Be7g(T):
    """Nacre I and II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 5.46 * 10**6 * T9**(-2./3) * math.exp(-12.827*T9**(-1./3)) * (
        + 1.0
        - 0.307 * T9 
        + 8.81 * 10**(-2) * T9**2
        - 1.06 * 10**(-2) * T9**3
        + 4.46 * 10**(-4) * T9**4
        )

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Be_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def Be7g_He3He4(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He3He4_Be7g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    back_base = 1.11 * 10**10 * math.exp(-18.407/T9) / (
        + 1.0 
        + 0.5 * math.exp(-4.979/T9)
        )
    back = back_base * forw * (ro_b**(-1)) * T9**(3./2)
    return (back /(constants.less_time(1)))


##################################

def Be_7_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_be7 = (7./(3*4)) * (He3He4_Be7g.__wrapped__(T)/Be7g_He3He4.__wrapped__(T))*X[0][3]*X[0][5]
    except OverflowError as e:
        X_be7 = 0
    X[0][6] = X_be7
    X[0][1] -= X_be7
    return X

####################################

# @Be_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def nbe7_He4He4(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = (
#         1.2 * 10**7 * T9
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @Be_7.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def He4he4_nBe7(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = nbe7_He4He4.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = 4.64 * forw * math.exp(-220.4/T9)
#     return (back /(constants.less_time(1)))

@Be_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def nBe7_He4He4(T):
    """ Wagoner, 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 2.05 * 10**4 * (
        + 1.0 
        + 3760 * T9
        )

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@Be_7.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4He4_nBe7(T):
    """Wagoner, 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = nBe7_He4He4.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)

    back = 4.70 * forw * math.exp(-220.39/T9)
    return (back /(constants.less_time(1)))


#########################################################

# Be_7.forward_rates.append(he3he4_be7g)
# Be_7.backward_rates.append(be7g_he3he4)
# Be_7.forward_rates.append(nbe7_He4He4)
# Be_7.backward_rates.append(He4He4_nBe7)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4
# 6 - Be7

Be_7.equilibrium = Be_7_equ


if __name__ == '__main__':
    # Be_7.show_rates()
    import numpy as np
    from tempreture import tfromT
    grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
    grid2 = grid
    Ts = constants.less_tempreture(grid2, units="K")
    ts = np.array([tfromT(T) for T in Ts])
    for T in np.copy(Ts):
        pass
        print(nbe7_He4He4.__wrapped__(T)/nBe7_He4He4.__wrapped__(T))
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
    plt.plot(ts, [nbe7_He4He4(T) for T in np.copy(Ts)], label='forw1')
    plt.plot(ts, [nBe7_He4He4(T) for T in np.copy(Ts)], label='forw2')
    plt.legend()
    plt.show()
    plt.cla()
    plt.clf()
    
    