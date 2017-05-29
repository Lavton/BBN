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
# He_4.mass_excess = constants.less_tempreture(2161062.7, units="eV")
He_4.set_mass_excess(4002603.25415, n_N=2, p_N=2)
He_4.tr_t =  0.0017
He_4.tr_T = tempreture.Tfromt(He_4.tr_t)
# He_4.tr_T = constants.less_tempreture(2*10**10, units="K")

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def pt_he4g(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 2.87 * 10**4 * T9**(-2./3) *math.exp(-3.87*T9**(-1./3)) * (
        + 1
        + 0.108 * T9**(1./3)
        + 0.466 * T9**(2./3)
        + 0.352 * T9 
        + 0.300 * T9**(4./3)
        +0.576 * T9**(5./3)
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4g_pt(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = pt_he4g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    # back = forw * math.exp(-He_4.mass_excess/T)
    ro_b = univ_func.rat_scale(T)

    back = 2.59 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-229.9/T9)
    return (back /(constants.less_time(1)))

##################################

def He_4_equ(X, T):
    T9 = constants.to_norm_tempreture(-T, units="T9")
    try:
        X_he4 = (4./3) * (pt_he4g.__wrapped__(T)/he4g_pt.__wrapped__(T))*X[0][0]*X[0][4]
    except OverflowError as e:
        X_he4 = 0
    X[0][5] = X_he4
    return X

####################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def nhe3_he4g(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = (
        6.0 * 10**3 * T9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4g_nhe3(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = nhe3_he4g.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    # back = forw * math.exp(-He_4.mass_excess/T)
    ro_b = univ_func.rat_scale(T)

    back = 2.60 * 10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-238.8/T9)
    return (back /(constants.less_time(1)))

#####################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def dd_he4g(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 24.1 * ro_b * (T9**(-2./3)) * math.exp(-4.26*(T9**(-1./3))) * (
        + T9**(2./3)
        + 0.685 * T9
        + 0.152 * T9**(4./3)
        + 0.265 * T9**(5./3)
        )
    return forw * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4g_dd(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = dd_he4g.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 4.50 * 10**10 * forw * ro_b**(-1) * T9**(3./2) * math.exp(-276.7 * (T9**(-1)))
    return back * (1./(constants.less_time(1)))

#########################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def dhe3_he4p(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 2.60 * 10**9 * T9**(-3./2) * math.exp(-2.99/T9)
    return forw * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4p_dhe3(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = dhe3_he4p.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 5.50 * forw * math.exp(-213.0/T9)
    return back * (1./(constants.less_time(1)))

#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def dt_he4n(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 1.38 * 10**9 * ro_b * T9**(-3./2) * math.exp(-0.745/T9)
    return forw * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4n_dt(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = dt_he4n.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 5.50 * forw * math.exp(-204.1/T9)
    return back * (1./(constants.less_time(1)))

#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he3he3_he4pp(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 1.19 * 10**10 * ro_b * T9**(-2./3) * math.exp(-12.25*T9**(-1./3)) * (
        1 + 0.0340 * T9**(1./3))
    return forw * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4pp_he3he3(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = he3he3_he4pp.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 3.37 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-149.2/T9)
    return back * (1./(constants.less_time(1)))

#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def tt_he4nn(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 1.10 * 10**9 * ro_b * T9**(-2./3) * math.exp(-4.87*T9**(-1./3)) * (
        1 + 0.0857 * T9**(1./3))
    return forw * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4nn_tt(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = tt_he4nn.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 3.37 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-131.5/T9)
    return back * (1./(constants.less_time(1)))

#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he3t_he4pn(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 5.60 * 10**9 * ro_b * T9**(-2./3) * math.exp(-7.72*T9**(-1./3)) * (
        1 + 0.0540 * T9**(1./3))
    return forw * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4pn_he3t(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = he3t_he4pn.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 3.37 * 10**(-10) * forw * ro_b * T9**(-3./2) * math.exp(-140.4/T9)
    return back * (1./(constants.less_time(1)))

#############################################

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he3t_he4d(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 3.88 * 10**9 * ro_b * T9**(-2./3) * math.exp(-7.72*T9**(-1./3)) * (
        1 + 0.0540 * T9**(1./3))
    return forw * (1./(constants.less_time(1)))

@He_4.equilib_zeroize
@functools.lru_cache(maxsize=8)
def he4d_he3t(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = he3t_he4d.__wrapped__(T) / constants.to_norm_time(1)
    ro_b = univ_func.rat_scale(T)
    back = 1.59 * forw * math.exp(-166.2/T9)
    return back * (1./(constants.less_time(1)))

#############################################


He_4.forward_rates.append(pt_he4g)
He_4.backward_rates.append(he4g_pt)
He_4.forward_rates.append(nhe3_he4g)
He_4.backward_rates.append(he4g_nhe3)
He_4.forward_rates.append(dd_he4g)
He_4.backward_rates.append(he4g_dd)
He_4.forward_rates.append(dhe3_he4p)
He_4.backward_rates.append(he4p_dhe3)
He_4.forward_rates.append(dt_he4n)
He_4.backward_rates.append(he4n_dt)
He_4.forward_rates.append(he3he3_he4pp)
He_4.backward_rates.append(he4pp_he3he3)
He_4.forward_rates.append(tt_he4nn)
He_4.backward_rates.append(he4nn_tt)
He_4.forward_rates.append(he3t_he4pn)
He_4.backward_rates.append(he4pn_he3t)
He_4.forward_rates.append(he3t_he4d)
He_4.backward_rates.append(he4d_he3t)

# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
# 5 - He4
He_4.ode_elem = {
    # "n": (lambda X, T: 
    #     -X[0]*(X[2]/2)*nd_tg(T) + (X[4]/3)*tg_nd(T)
    #     -X[0]*(X[3]/3)*nhe3_pt(T) + X[1]*(X[4]/3)*pt_nhe3(T)
    #     ),
    # "H_1": (lambda X, T:
    #     +X[0]*(X[3]/3)*nhe3_pt(T) - X[1]*(X[4]/3)*pt_nhe3(T)
    #     ),
    # "H_2": (lambda X, T:  -X[0]*(X[2]/2)*nd_tg(T) + (X[4]/3)*tg_nd(T)),
    # "He_3": (lambda X, T:
    #     -X[0]*(X[3]/3)*nhe3_pt(T) + X[1]*(X[4]/3)*pt_nhe3(T)
    #     ),
    # "He_4": (lambda X, T: 
    #     X[0]*(X[2]/2)*nd_tg(T) - (X[4]/3)*tg_nd(T)
    #     +X[0]*(X[3]/3)*nhe3_pt(T) - X[1]*(X[4]/3)*pt_nhe3(T)
    #     )
}

He_4.jacob = { 
    # "n": { #берём первую строку и дифференцируем по каждому
    #     "n": (lambda X, T: 
    #         -(X[2]/2)*nd_tg(T)
    #         -(X[3]/3)*nhe3_pt(T)
    #         ),
    #     "H_1": (lambda X, T:
    #         +(X[4]/3)*pt_nhe3(T)
    #         ),
    #     "H_2": (lambda X, T: -X[0]*(1./2)*nd_tg(T)),
    #     "He_3": (lambda X, T: 
    #         -X[0]*(1./3)*nhe3_pt(T)
    #         ),
    #     "He_4": (lambda X, T: 
    #         (1./3)*tg_nd(T)
    #         +X[1]*(1./3)*pt_nhe3(T)
    #         )
    # },
    # "H_1": {
    #     "n": (lambda X, T:
    #         +(X[3]/3)*nhe3_pt(T)
    #         ),
    #     "H_1": (lambda X, T:
    #         -(X[4]/3)*pt_nhe3(T)
    #         ), 
    #     "He_3": (lambda X, T: 
    #         +X[0]*(1./3)*nhe3_pt(T)
    #         ),
    #     "He_4": (lambda X, T: 
    #         -X[1]*(1./3)*pt_nhe3(T)
    #         )
    # },
    # "H_2": {
    #     "n": (lambda X, T:  -(X[2]/2)*nd_tg(T)),
    #     "H_2": (lambda X, T: -X[0]*(1./2)*nd_tg(T)),
    #     "He_4": (lambda X, T: (1./3)*tg_nd(T))
    # },
    # "He_3": {
    #     "n": (lambda X, T:
    #         -(X[3]/3)*nhe3_pt(T)
    #         ),
    #     "H_1": (lambda X, T:
    #         +(X[4]/3)*pt_nhe3(T)
    #         ), 
    #     "He_3": (lambda X, T: 
    #         -X[0]*(1./3)*nhe3_pt(T)
    #         ),
    #     "He_4": (lambda X, T: 
    #         +X[1]*(1./3)*pt_nhe3(T)
    #         )
    # },

    # "He_4": {
    #     "n": (lambda X, T: 
    #         +(X[2]/2)*nd_tg(T)
    #         +(X[3]/3)*nhe3_pt(T)
    #         ),
    #     "H_1": (lambda X, T:
    #         -(X[4]/3)*pt_nhe3(T)
    #         ), 
    #     "He_3": (lambda X, T: 
    #         +X[0]*(1./3)*nhe3_pt(T)
    #         ),

    #     "H_2": (lambda X, T: +X[0]*(1./2)*nd_tg(T)),
    #     "He_4": (lambda X, T: 
    #         -(1./3)*tg_nd(T)
    #         -X[1]*(1./3)*pt_nhe3(T)
    #         )
    # }
}

He_4.equilibrium = He_4_equ


if __name__ == '__main__':
    He_4.show_rates()
    