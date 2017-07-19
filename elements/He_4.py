"""
06
гелий-4, или $^{4}He$
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

He_4 = Element("He_4", 0.0)
He_4.A = 4
# from Audi et all, 2003
He_4.set_mass_excess(4002603.25415, n_N=2, p_N=2)
He_4.tr_t =  0.0016 * max(constants.nu_0/constants.nu_n, 1.0)
He_4.tr_T = tempreture.Tfromt(He_4.tr_t)

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def H3p_He4g(T):
    """Shapiro, 2004"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 2.2 * 10**4 * T9**(-2./3) * math.exp(-3.869/T9**(1./3)) * (
        + 1.0
        + 0.108 * T9**(1./3)
        + 1.68 * T9**(2./3)
        + 1.26 * T9
        + 0.551 * T9**(4./3)
        + 1.06 * T9**(5./3)
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))


@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4g_H3p(T):
    """Wagoner, 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H3p_He4g.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)

    back = 2.61 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-229.94/T9)
    return (back /(constants.less_time(1)))

He_4.reactions.append((
    ("p", "T"),
    ("He_4",), 
    H3p_He4g, 
    He4g_H3p
    ))

##################################

def He_4_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        X_he4 =  (4./3) * (H3p_He4g.__wrapped__(T)/He4g_H3p.__wrapped__(T))*X[0][1]*X[0][4]
    except OverflowError as e:
        X_he4 = 0
    X[0][5] = X_he4
    X[0][1] -= X_he4
    return X

####################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He3n_He4g(T):
    """Wagoner 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 6.62 * (
        + 1.0 
        + 905 * T9
        )

    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4g_He3n(T):
    """Wagoner, 1969"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He3n_He4g.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)

    back = 2.61 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-238.81/T9)
    return (back /(constants.less_time(1)))

He_4.reactions.append((
    ("n", "He_3"),
    ("He_4",), 
    He3n_He4g, 
    He4g_He3n
    ))

#####################################

He_4.add_interpo("H2H2_He4g", """
0.001 1.46E−15 1.08E−15 1.72E−15 0.14 4.57E−02 3.41E−02 5.45E−02
0.002 6.02E−12 4.47E−12 7.12E−12 0.15 5.26E−02 3.93E−02 6.27E−02
0.003 3.30E−10 2.45E−10 3.90E−10 0.16 5.98E−02 4.47E−02 7.13E−02
0.004 4.05E−09 3.01E−09 4.79E−09 0.18 7.48E−02 5.60E−02 8.93E−02
0.005 2.39E−08 1.77E−08 2.83E−08 0.2 9.06E−02 6.79E−02 1.08E−01
0.006 9.19E−08 6.82E−08 1.09E−07 0.25 1.32E−01 9.96E−02 1.58E−01
0.007 2.68E−07 1.99E−07 3.18E−07 0.3 1.76E−01 1.33E−01 2.11E−01
0.008 6.47E−07 4.80E−07 7.66E−07 0.35 2.21E−01 1.67E−01 2.66E−01
0.009 1.36E−06 1.01E−06 1.61E−06 0.4 2.67E−01 2.02E−01 3.23E−01
0.01 2.57E−06 1.90E−06 3.04E−06 0.45 3.15E−01 2.37E−01 3.82E−01
0.011 4.47E−06 3.31E−06 5.29E−06 0.5 3.63E−01 2.73E−01 4.44E−01
0.012 7.28E−06 5.41E−06 8.62E−06 0.6 4.65E−01 3.47E−01 5.78E−01
0.013 1.13E−05 8.36E−06 1.33E−05 0.7 5.74E−01 4.25E−01 7.27E−01
0.014 1.67E−05 1.24E−05 1.98E−05 0.8 6.92E−01 5.07E−01 8.92E−01
0.015 2.38E−05 1.77E−05 2.82E−05 0.9 8.18E−01 5.94E−01 1.07E+00
0.016 3.29E−05 2.44E−05 3.90E−05 1. 9.54E−01 6.86E−01 1.27E+00
0.018 5.83E−05 4.33E−05 6.91E−05 1.25 1.33E+00 9.38E−01 1.85E+00
0.02 9.52E−05 7.07E−05 1.13E−04 1.5 1.76E+00 1.22E+00 2.51E+00
0.025 2.52E−04 1.87E−04 2.99E−04 1.75 2.23E+00 1.53E+00 3.25E+00
0.03 5.27E−04 3.91E−04 6.25E−04 2. 2.75E+00 1.88E+00 4.05E+00
0.04 1.52E−03 1.13E−03 1.80E−03 2.5 3.92E+00 2.67E+00 5.81E+00
0.05 3.19E−03 2.37E−03 3.79E−03 3. 5.27E+00 3.62E+00 7.74E+00
0.06 5.57E−03 4.14E−03 6.62E−03 3.5 6.80E+00 4.76E+00 9.81E+00
0.07 8.65E−03 6.44E−03 1.03E−02 4. 8.50E+00 6.07E+00 1.20E+01
0.08 1.24E−02 9.22E−03 1.47E−02 5. 1.23E+01 9.17E+00 1.67E+01
0.09 1.68E−02 1.25E−02 1.99E−02 6. 1.66E+01 1.28E+01 2.17E+01
0.1 2.16E−02 1.61E−02 2.58E−02 7. 2.12E+01 1.67E+01 2.69E+01
0.11 2.71E−02 2.02E−02 3.22E−02 8. 2.58E+01 2.08E+01 3.20E+01
0.12 3.29E−02 2.45E−02 3.92E−02 9. 3.04E+01 2.49E+01 3.70E+01
0.13 3.91E−02 2.92E−02 4.66E−02 10. 3.48E+01 2.88E+01 4.18E+01
    """)

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
@He_4.nacreII
def H2H2_He4g(T):
    """NACRE I and II"""
    #WARNING: сильно ниже скорость! 0.0389719772608
    T9 = constants.to_norm_tempreture(T, units="T9")
    ro_b = univ_func.rat_scale(T)
    base_rate = 42.1 * T9**(-2./3) * math.exp(-4.259/T9**(1./3)) * (
        + 1. 
        + 0.514 * T9
        + 0.399 * T9**2
        - 1.18 * 10**(-2) * T9**3
        )

    return base_rate * ro_b * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4g_H2H2(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H2H2_He4g.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 4.53 * 10**10 * forw * ro_b**(-1) * T9**(3./2) * math.exp(-276.74 * (T9**(-1)))
    return back * (1./(constants.less_time(1)))

He_4.reactions.append((
    ("D", "D"), 
    ("He_4",), 
    H2H2_He4g, 
    He4g_H2H2
    ))

#########################################

He_4.add_interpo("He3H2_He4p", """
0.001 3.54E−19 3.21E−19 3.87E−19 0.14 2.88E+05 2.64E+05 3.15E+05
0.002 6.14E−13 5.58E−13 6.72E−13 0.15 3.85E+05 3.53E+05 4.22E+05
0.003 6.36E−10 5.78E−10 6.96E−10 0.16 5.03E+05 4.61E+05 5.51E+05
0.004 5.02E−08 4.56E−08 5.49E−08 0.18 8.07E+05 7.40E+05 8.84E+05
0.005 1.11E−06 1.01E−06 1.22E−06 0.2 1.22E+06 1.11E+06 1.33E+06
0.006 1.18E−05 1.07E−05 1.29E−05 0.25 2.77E+06 2.53E+06 3.04E+06
0.007 7.74E−05 7.04E−05 8.47E−05 0.3 5.19E+06 4.71E+06 5.69E+06
0.008 3.64E−04 3.31E−04 3.99E−04 0.35 8.51E+06 7.68E+06 9.34E+06
0.009 1.35E−03 1.23E−03 1.48E−03 0.4 1.27E+07 1.14E+07 1.39E+07
0.01 4.15E−03 3.78E−03 4.54E−03 0.45 1.76E+07 1.57E+07 1.93E+07
0.011 1.11E−02 1.01E−02 1.21E−02 0.5 2.31E+07 2.06E+07 2.54E+07
0.012 2.64E−02 2.40E−02 2.89E−02 0.6 3.53E+07 3.13E+07 3.89E+07
0.013 5.74E−02 5.22E−02 6.29E−02 0.7 4.83E+07 4.26E+07 5.33E+07
0.014 1.16E−01 1.05E−01 1.26E−01 0.8 6.12E+07 5.37E+07 6.77E+07
0.015 2.18E−01 1.98E−01 2.39E−01 0.9 7.35E+07 6.42E+07 8.15E+07
0.016 3.89E−01 3.54E−01 4.26E−01 1. 8.49E+07 7.39E+07 9.43E+07
0.018 1.08E+00 9.87E−01 1.19E+00 1.25 1.09E+08 9.40E+07 1.21E+08
0.02 2.62E+00 2.38E+00 2.87E+00 1.5 1.27E+08 1.08E+08 1.42E+08
0.025 1.52E+01 1.39E+01 1.67E+01 1.75 1.39E+08 1.18E+08 1.56E+08
0.03 5.82E+01 5.30E+01 6.37E+01 2. 1.48E+08 1.25E+08 1.67E+08
0.04 4.07E+02 3.71E+02 4.46E+02 2.5 1.57E+08 1.31E+08 1.78E+08
0.05 1.62E+03 1.48E+03 1.77E+03 3. 1.61E+08 1.33E+08 1.83E+08
0.06 4.63E+03 4.23E+03 5.07E+03 3.5 1.61E+08 1.32E+08 1.83E+08
0.07 1.07E+04 9.80E+03 1.17E+04 4. 1.59E+08 1.30E+08 1.82E+08
0.08 2.14E+04 1.96E+04 2.35E+04 5. 1.53E+08 1.24E+08 1.76E+08
0.09 3.85E+04 3.52E+04 4.22E+04 6. 1.46E+08 1.18E+08 1.68E+08
0.1 6.38E+04 5.85E+04 6.99E+04 7. 1.38E+08 1.11E+08 1.59E+08
0.11 9.94E+04 9.11E+04 1.09E+05 8. 1.31E+08 1.05E+08 1.51E+08
0.12 1.47E+05 1.35E+05 1.61E+05 9. 1.24E+08 9.85E+07 1.42E+08
0.13 2.09E+05 1.92E+05 2.29E+05 10. 1.17E+08 9.27E+07 1.34E+08
    """)

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
@He_4.nacreII
def He3H2_He4p(T):
    """Caughlan, 1988; NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    base_rate = 5.86 * 10**10 / T9**(2./3) * math.exp(-7.181/T9**(1./3) - (T9/0.315)**2) * (
        + 1.
        + 0.58 * T9**(1./3)
        + 0.142 * T9**(2./3)
        + 5.78 * 10**(-2) * T9 
        + 2.25 * T9**(4./3)
        + 2.32 * T9**(5./3)
        ) + 4.36 * 10**8 / T9**(1./2) * math.exp(-1.720/T9)

    return base_rate * ro_b * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4p_He3H2(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He3H2_He4p.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 5.54 * forw * math.exp(-212.99/T9)
    return back * (1./(constants.less_time(1)))

He_4.reactions.append((
    ("D", "He_3"),
    ("He_4", "p"), 
    He3H2_He4p, 
    He4p_He3H2
    ))

#############################################

He_4.add_interpo("H3H2_He4n", """
0.001 1.87E−07 1.67E−07 2.09E−07 0.14 1.05E+08 9.66E+07 1.13E+08
0.002 1.37E−03 1.22E−03 1.53E−03 0.15 1.21E+08 1.12E+08 1.31E+08
0.003 9.98E−02 8.90E−02 1.11E−01 0.16 1.38E+08 1.28E+08 1.49E+08
0.004 1.47E+00 1.31E+00 1.64E+00 0.18 1.73E+08 1.60E+08 1.85E+08
0.005 9.93E+00 8.86E+00 1.10E+01 0.2 2.06E+08 1.91E+08 2.21E+08
0.006 4.24E+01 3.78E+01 4.71E+01 0.25 2.83E+08 2.63E+08 3.04E+08
0.007 1.35E+02 1.20E+02 1.50E+02 0.3 3.47E+08 3.23E+08 3.72E+08
0.008 3.49E+02 3.12E+02 3.87E+02 0.35 3.98E+08 3.70E+08 4.25E+08
0.009 7.80E+02 7.01E+02 8.62E+02 0.4 4.36E+08 4.07E+08 4.66E+08
0.01 1.56E+03 1.41E+03 1.71E+03 0.45 4.65E+08 4.34E+08 4.96E+08
0.011 2.85E+03 2.59E+03 3.12E+03 0.5 4.85E+08 4.53E+08 5.17E+08
0.012 4.87E+03 4.44E+03 5.31E+03 0.6 5.10E+08 4.77E+08 5.42E+08
0.013 7.85E+03 7.19E+03 8.53E+03 0.7 5.19E+08 4.86E+08 5.51E+08
0.014 1.21E+04 1.11E+04 1.31E+04 0.8 5.19E+08 4.87E+08 5.50E+08
0.015 1.79E+04 1.65E+04 1.93E+04 0.9 5.13E+08 4.82E+08 5.44E+08
0.016 2.56E+04 2.37E+04 2.75E+04 1. 5.04E+08 4.74E+08 5.34E+08
0.018 4.81E+04 4.48E+04 5.16E+04 1.25 4.76E+08 4.48E+08 5.04E+08
0.02 8.30E+04 7.74E+04 8.87E+04 1.5 4.47E+08 4.21E+08 4.73E+08
0.025 2.47E+05 2.32E+05 2.63E+05 1.75 4.19E+08 3.95E+08 4.43E+08
0.03 5.70E+05 5.34E+05 6.06E+05 2. 3.94E+08 3.71E+08 4.17E+08
0.04 1.93E+06 1.81E+06 2.06E+06 2.5 3.52E+08 3.32E+08 3.72E+08
0.05 4.62E+06 4.30E+06 4.93E+06 3. 3.19E+08 3.00E+08 3.37E+08
0.06 8.98E+06 8.34E+06 9.63E+06 3.5 2.92E+08 2.75E+08 3.09E+08
0.07 1.52E+07 1.41E+07 1.64E+07 4. 2.70E+08 2.54E+08 2.86E+08
0.08 2.34E+07 2.16E+07 2.52E+07 5. 2.37E+08 2.23E+08 2.52E+08
0.09 3.35E+07 3.09E+07 3.61E+07 6. 2.14E+08 2.01E+08 2.27E+08
0.1 4.52E+07 4.17E+07 4.88E+07 7. 1.97E+08 1.84E+08 2.09E+08
0.11 5.85E+07 5.40E+07 6.30E+07 8. 1.83E+08 1.72E+08 1.95E+08
0.12 7.30E+07 6.73E+07 7.86E+07 9. 1.73E+08 1.62E+08 1.84E+08
0.13 8.84E+07 8.16E+07 9.52E+07 10. 1.64E+08 1.54E+08 1.75E+08
    """)

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
@He_4.nacreII
def H3H2_He4n(T):
    """nacre I and II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    ro_b = univ_func.rat_scale(T)
    base_rate = 8.29 * 10**10 * T9**(-2./3) * math.exp(-4.524*T9**(-1./3) - (T9/0.08)**2) * (
        + 1.0
        + 17.2 * T9 
        + 175 * T9**2
        ) + 8.12 * 10**8 * T9**(-0.712) * math.exp(-0.506/T9)

    return base_rate * ro_b * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4n_H3H2(T):
    """NACRE II"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H3H2_He4n.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 5.54 * forw * math.exp(-204.12/T9)
    return back * (1./(constants.less_time(1)))

He_4.reactions.append((
    ("D", "T"),
    ("He_4", "n"), 
    H3H2_He4n, 
    He4n_H3H2
    ))

#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he3he3_he4pp(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 1.19 * 10**10 * T9**(-2./3) * math.exp(-12.25 * T9**(-1./3)) * (
        + 1. 
        + 0.0340 * T9**(1./3)
        )
    return forw * ro_b * (1./(constants.less_time(1)))

He_4.add_interpo("He3He3_He4pp", """
0.003 2.80E−25 2.55E−25 3.04E−25 0.15 1.74E+01 1.61E+01 1.92E+01
0.004 5.55E−22 5.07E−22 6.04E−22 0.16 2.71E+01 2.51E+01 3.00E+01
0.005 1.22E−19 1.12E−19 1.33E−19 0.18 5.96E+01 5.52E+01 6.59E+01
0.006 7.48E−18 6.83E−18 8.14E−18 0.2 1.17E+02 1.08E+02 1.30E+02
0.007 1.99E−16 1.82E−16 2.17E−16 0.25 4.50E+02 4.16E+02 4.98E+02
0.008 2.98E−15 2.72E−15 3.24E−15 0.3 1.24E+03 1.15E+03 1.38E+03
0.009 2.93E−14 2.68E−14 3.19E−14 0.35 2.79E+03 2.58E+03 3.09E+03
0.01 2.09E−13 1.91E−13 2.28E−13 0.4 5.40E+03 4.99E+03 5.98E+03
0.011 1.17E−12 1.07E−12 1.27E−12 0.45 9.42E+03 8.70E+03 1.04E+04
0.012 5.34E−12 4.89E−12 5.82E−12 0.5 1.52E+04 1.40E+04 1.68E+04
0.013 2.08E−11 1.90E−11 2.27E−11 0.6 3.30E+04 3.04E+04 3.64E+04
0.014 7.08E−11 6.47E−11 7.71E−11 0.7 6.10E+04 5.61E+04 6.73E+04
0.015 2.15E−10 1.97E−10 2.35E−10 0.8 1.01E+05 9.24E+04 1.11E+05
0.016 5.95E−10 5.44E−10 6.48E−10 0.9 1.53E+05 1.40E+05 1.69E+05
0.018 3.59E−09 3.28E−09 3.91E−09 1. 2.19E+05 2.01E+05 2.42E+05
0.02 1.68E−08 1.54E−08 1.84E−08 1.25 4.46E+05 4.06E+05 4.93E+05
0.025 3.71E−07 3.40E−07 4.05E−07 1.5 7.60E+05 6.90E+05 8.40E+05
0.03 3.91E−06 3.58E−06 4.27E−06 1.75 1.16E+06 1.05E+06 1.28E+06
0.04 1.19E−04 1.10E−04 1.31E−04 2. 1.64E+06 1.49E+06 1.81E+06
0.05 1.35E−03 1.24E−03 1.48E−03 2.5 2.83E+06 2.56E+06 3.13E+06
0.06 8.49E−03 7.82E−03 9.32E−03 3. 4.30E+06 3.89E+06 4.74E+06
0.07 3.68E−02 3.39E−02 4.04E−02 3.5 6.01E+06 5.43E+06 6.60E+06
0.08 1.23E−01 1.13E−01 1.35E−01 4. 7.92E+06 7.16E+06 8.69E+06
0.09 3.40E−01 3.14E−01 3.74E−01 5. 1.22E+07 1.11E+07 1.34E+07
0.1 8.13E−01 7.52E−01 8.95E−01 6. 1.70E+07 1.54E+07 1.85E+07
0.11 1.74E+00 1.61E+00 1.92E+00 7. 2.21E+07 2.01E+07 2.40E+07
0.12 3.41E+00 3.15E+00 3.76E+00 8. 2.73E+07 2.50E+07 2.96E+07
0.13 6.21E+00 5.75E+00 6.85E+00 9. 3.27E+07 2.99E+07 3.53E+07
0.14 1.07E+01 9.86E+00 1.18E+01 10. 3.80E+07 3.49E+07 4.09E+07

    """)

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
@He_4.nacreII
def He3He3_He4pp(T):
    """NACRE"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    base_rate = 5.59 * 10**10 * T9**(-2./3) * math.exp(-12.277 * T9**(-1./3)) * (
        + 1.0
        - 0.135 * T9 
        + 2.54 * 10**(-2) * T9**2
        - 1.29 * 10**(-3) * T9**3
        )
    return base_rate * ro_b * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4pp_He3He3(T):
    """NACRE """
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He3He3_He4pp.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 3.92 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-149.23/T9)
    return back * (1./(constants.less_time(1)))

He_4.reactions.append((
    ("He_3", "He_3"),
    ("He_4", "p", "p"),
    He3He3_He4pp,
    He4pp_He3He3
    ))

#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def H3H3_He4nn(T):
    """Caughlan-Fowler"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    ro_b = univ_func.rat_scale(T)
    base_rate = 1.67 * 10**9 * T9**(-2./3) * math.exp(-4.872/T9**(1./3)) * (
        + 1.0
        + 0.086 * T9**(1./3)
        - 0.455 * T9**(2./3)
        - 0.272 * T9
        + 0.148 * T9**(4./3)
        + 0.225 * T9**(5./3)
        )
    return base_rate * ro_b * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4nn_H3H3(T):
    """Caughlan-Fowler"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H3H3_He4nn.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 3.38 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-131.504/T9)
    return back * (1./(constants.less_time(1)))

He_4.reactions.append((
    ("T", "T"),
    ("He_4", "n", "n"),
    H3H3_He4nn,
    He4nn_H3H3
    ))

#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He3H3_He4pn(T):
    """Caughlan-Fowler"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    T9A = T9/(1.0 + 0.115 * T9)
    base_rate = 7.71 * 10**9 * T9A**(5./6)/T9**(3./2) * math.exp(-7.733/T9A**(1./3))
    return base_rate * ro_b * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4pn_He3H3(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He3H3_He4pn.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 3.39 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-140.367/T9)
    return back * (1./(constants.less_time(1)))

He_4.reactions.append((
    ("He_3", "T"),
    ("He_4", "p", "n"),
    He3H3_He4pn,
    He4pn_He3H3
    ))
#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He3H3_He4H2(T):
    """Caughlan-Fowler"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    T9A = T9/(1.0 + 0.128 * T9)
    ro_b = univ_func.rat_scale(T)
    base_rate = 5.46 * 10**9 * T9A**(5./6)/T9**(3./2) * math.exp(-7.733/T9A**(1./3))

    return base_rate * ro_b * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def He4H2_He3H2(T):
    """Caughlan-Fowler"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = He3H3_He4H2.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 1.60 * forw * math.exp(-166.182/T9)
    return back * (1./(constants.less_time(1)))

He_4.reactions.append((
    ("He_3", "T"),
    ("He_4", "D"),
    He3H3_He4H2,
    He4H2_He3H2
    ))

#############################################


# He_4.forward_rates.append(pt_he4g)
# He_4.backward_rates.append(he4g_pt)
# He_4.forward_rates.append(nhe3_he4g)
# He_4.backward_rates.append(he4g_nhe3)
# He_4.forward_rates.append(dd_he4g)
# He_4.backward_rates.append(he4g_dd)
# He_4.forward_rates.append(dhe3_he4p)
# He_4.backward_rates.append(he4p_dhe3)
# He_4.forward_rates.append(dt_he4n)
# He_4.backward_rates.append(he4n_dt)
# He_4.forward_rates.append(he3he3_he4pp)
# He_4.backward_rates.append(he4pp_he3he3)
# He_4.forward_rates.append(tt_he4nn)
# He_4.backward_rates.append(he4nn_tt)
# He_4.forward_rates.append(he3t_he4pn)
# He_4.backward_rates.append(he4pn_he3t)
# He_4.forward_rates.append(he3t_he4d)
# He_4.backward_rates.append(he4d_he3t)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4

He_4.equilibrium = He_4_equ
He_4.names = ["He_4", "a", "^4He"]


if __name__ == '__main__':
    # He_4.show_rates()
    import numpy as np
    from tempreture import tfromT
    grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=320)
    grid2 = grid
    Ts = constants.less_tempreture(grid2, units="K")
    ts = np.array([tfromT(T) for T in Ts])
    for T in np.copy(Ts):
        pass
        print(he3he3_he4pp.__wrapped__(T)/He3He3_He4pp.__wrapped__(T))
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
    plt.plot(ts, [he3he3_he4pp(T) for T in np.copy(Ts)], label='forw1')
    plt.plot(ts, [He3He3_He4pp(T) for T in np.copy(Ts)], label='forw2')
    plt.legend()
    plt.show()
    plt.cla()
    plt.clf()
    

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def pt_he4g(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = 2.87 * 10**4 * T9**(-2./3) *math.exp(-3.87*T9**(-1./3)) * (
#         + 1.
#         + 0.108 * T9**(1./3)
#         + 0.466 * T9**(2./3)
#         + 0.352 * T9 
#         + 0.300 * T9**(4./3)
#         + 0.576 * T9**(5./3)
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))


# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4g_pt(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = pt_he4g.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = 2.59 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-229.9/T9)
#     return (back /(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def nhe3_he4g(T):
#     """
#     Wagoner
#     """
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     base_rate = (
#         6.0 * 10**3 * T9
#         )
#     ro_b = univ_func.rat_scale(T)
#     return base_rate * ro_b/(constants.less_time(1))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4g_nhe3(T):
#     """Wagoner, 1966"""
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = nhe3_he4g.__wrapped__(T) / constants.to_norm_time(1)
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)

#     back = 2.60 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-238.8/T9)
#     return (back /(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def dd_he4g(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     forw = 24.1 * (T9**(-2./3)) * math.exp(-4.26*(T9**(-1./3))) * (
#         + T9**(2./3)
#         + 0.685 * T9
#         + 0.152 * T9**(4./3)
#         + 0.265 * T9**(5./3)
#         )
#     return forw * ro_b * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4g_dd(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = dd_he4g.__wrapped__(T) / constants.to_norm_time(1)
#     ro_b = univ_func.rat_scale(T)
#     back = 4.50 * 10**10 * forw * ro_b**(-1) * T9**(3./2) * math.exp(-276.7 * (T9**(-1)))
#     return back * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def dhe3_he4p(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     forw = 2.60 * 10**9 * T9**(-3./2) * math.exp(-2.99/T9)
#     return forw * ro_b * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4p_dhe3(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = dhe3_he4p.__wrapped__(T) / constants.to_norm_time(1)
#     ro_b = univ_func.rat_scale(T)
#     back = 5.50 * forw * math.exp(-213.0/T9)
#     return back * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def dt_he4n(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     forw = 1.38 * 10**9 * T9**(-3./2) * math.exp(-0.745/T9)
#     return forw * ro_b * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4n_dt(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = dt_he4n.__wrapped__(T) / constants.to_norm_time(1)
#     ro_b = univ_func.rat_scale(T)
#     back = 5.50 * forw * math.exp(-204.1/T9)
#     return back * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4pp_he3he3(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = he3he3_he4pp.__wrapped__(T) / constants.to_norm_time(1)
#     ro_b = univ_func.rat_scale(T)
#     back = 3.37 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-149.2/T9)
#     return back * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def tt_he4nn(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     forw = 1.10 * 10**9 * T9**(-2./3) * math.exp(-4.87*T9**(-1./3)) * (
#         + 1. 
#         + 0.0857 * T9**(1./3)
#         )
#     return forw * ro_b * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4nn_tt(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = tt_he4nn.__wrapped__(T) / constants.to_norm_time(1)
#     ro_b = univ_func.rat_scale(T)
#     back = 3.37 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-131.5/T9)
#     return back * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he3t_he4pn(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     forw = 5.60 * 10**9 * T9**(-2./3) * math.exp(-7.72*T9**(-1./3)) * (
#         + 1.
#         + 0.0540 * T9**(1./3)
#         )
#     return forw * ro_b * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4pn_he3t(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = he3t_he4pn.__wrapped__(T) / constants.to_norm_time(1)
#     ro_b = univ_func.rat_scale(T)
#     back = 3.37 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-140.4/T9)
#     return back * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he3t_he4d(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     E = constants.to_norm_tempreture(T, units="MeV")
#     ro_b = univ_func.rat_scale(T)
#     forw = 3.88 * 10**9 * T9**(-2./3) * math.exp(-7.72*T9**(-1./3)) * (
#         +1.
#         + 0.0540 * T9**(1./3)
#         )
#     return forw * ro_b * (1./(constants.less_time(1)))

# @He_4.equilib_zeroize
# @functools.lru_cache(maxsize=8)
# def he4d_he3t(T):
#     T9 = constants.to_norm_tempreture(T, units="T9")
#     forw = he3t_he4d.__wrapped__(T) / constants.to_norm_time(1)
#     ro_b = univ_func.rat_scale(T)
#     back = 1.59 * forw * math.exp(-166.2/T9)
#     return back * (1./(constants.less_time(1)))
