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

He_3 = Element("He_3", 0.0)
He_3.A = 3
# from Audi et all, 2003
He_3.set_mass_excess(3016029.3191, n_N=1, p_N=2)
# He_3.tr_T = constants.less_tempreture(2*10**7, units="K")
He_3.tr_t =  0.009
He_3.tr_T = tempreture.Tfromt(He_3.tr_t)

@He_3.equilib_zeroize
def d_pg_he3(T):
    """
    Smith et all
    table 1 reac 2
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 2.65*(10**3)*(T9**(-2./3))*\
        math.exp(-3.720/(T9**(1./3))) * (
            +1.
            + 0.112 * T9**(1./3)
            + 1.99 * T9**(2./3)
            + 1.56 * T9
            + 0.162 * T9**(4./3)
            + 0.324 * T9**(5./3)
            )

    ro_b = univ_func.rat_scale(T)
    # base_rate = 0 # отладка
    return base_rate * ro_b/(constants.less_time(1))

@He_3.equilib_zeroize
def he3_gp_d(T):
    """
    Vagoner
    t2, r2
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = d_pg_he3.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)

    back = 1.63*10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-He_3.mass_excess/T)
    return (back /(constants.less_time(1)))


def He_3_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    try:
        # print(T, he3_gp_d.__wrapped__(T))
        X_he = (3./2) * (d_pg_he3.__wrapped__(T)/he3_gp_d.__wrapped__(T))*X[0][1]*X[0][2]
        # X_he = 4.463*(10**-11)*(T9**3.)*(constants.nu_n**2)*math.exp(89.564/T9)*X[0][0]*(X[0][1]**2)
    except OverflowError as e:
        X_he = 0
    X[0][3] = X_he
    return X

@He_3.equilib_zeroize
def dd_nhe3(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = d_pg_he3.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    ro_b = univ_func.rat_scale(T)
    forw = 3.9*(10**8)*ro_b*(T9**(-2./3)) * math.exp(-4.26*(T9**(-1./3))) * (
        + 1
        + 0.0979 * T9**(1./3)
        + 0.642 * T9**(2./3)
        + 0.440 * T9
        )
    return forw * (1./(constants.less_time(1)))

@He_3.equilib_zeroize
def nhe3_dd(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = dd_nhe3.__wrapped__(T) / constants.to_norm_time(1)
    back = 1.73 * forw * math.exp(-37.94 * (T9**(-1)))
    return back * (1./(constants.less_time(1)))


He_3.forward_rates.append(d_pg_he3)
He_3.backward_rates.append(he3_gp_d)
He_3.forward_rates.append(dd_nhe3)
He_3.backward_rates.append(nhe3_dd)


# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
He_3.ode_elem = {
    "n": (lambda X, T: (X[2]/2)*(X[2]/2)*dd_nhe3(T) - X[0]*(X[3]/3)*nhe3_dd(T)),
    "H_1": (lambda X, T: -X[1]*(X[2]/2)*d_pg_he3(T) + (X[3]/3)*he3_gp_d(T)),
    "H_2": (
        lambda X, T: -X[1]*(X[2]/2)*d_pg_he3(T) + (X[3]/3)*he3_gp_d(T)
        -2*(X[2]/2)*(X[2]/2)*dd_nhe3(T) + 2*X[0]*(X[3]/3)*nhe3_dd(T)
        ),
    "He_3": (
        lambda X, T: +X[1]*(X[2]/2)*d_pg_he3(T) - (X[3]/3)*he3_gp_d(T)
        + (X[2]/2)*(X[2]/2)*dd_nhe3(T) - X[0]*(X[3]/3)*nhe3_dd(T)
        )
}

He_3.jacob = { 
    "n": {
        "n": (lambda X, T: -(X[3]/3)*nhe3_dd(T)),
        "H_2": (lambda X, T: 2*(1./2)*(X[2]/2)*dd_nhe3(T)),
        "He_3": (lambda X, T: -X[0]*(1./3)*nhe3_dd(T)),
    },
    "H_1": {
        "H_1": (lambda X, T: -(X[2]/2)*d_pg_he3(T)),
        "H_2": (lambda X, T: -X[1]*(1./2)*d_pg_he3(T)),
        "He_3": (lambda X, T: +(1./3)*he3_gp_d(T))
    },
    "H_2": { #берём первую строку и дифференцируем по каждому
        "n": (lambda X, T: 2*(X[3]/3)*nhe3_dd(T)),
        "H_1": (lambda X, T: -(X[2]/2)*d_pg_he3(T)),
        "H_2": (
            lambda X, T: -X[1]*(1./2)*d_pg_he3(T)
            -2*2*(X[2]/2)*(1./2)*dd_nhe3(T)
            ),
        "He_3": (
            lambda X, T: +(1./3)*he3_gp_d(T)
            +2*X[0]*(1./3)*nhe3_dd(T)
            )
    },
    "He_3": {
        "n": (lambda X, T: -(X[3]/3)*nhe3_dd(T)),
        "H_1": (lambda X, T: (X[2]/2)*d_pg_he3(T)),
        "H_2": (
            lambda X, T: X[1]*(1./2)*d_pg_he3(T)
            +2*(1./2)*(X[2]/2)*dd_nhe3(T)
            ),
        "He_3": (
            lambda X, T: -(1./3)*he3_gp_d(T)
            -X[0]*(1./3)*nhe3_dd(T)
            )
    }
}

He_3.equilibrium = He_3_equ


if __name__ == '__main__':
    He_3.show_rates()
    exit()
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    Ts = constants.less_tempreture(np.logspace(math.log10(1e-2*1e6), math.log10(1*1e6), num=200), units="eV")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{tempreture} (MeV)')
    plt.ylabel(r'\textbf{\univ_func}')

    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [He_3_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [He_3_backward_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{back}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()
    import nTOp
    plt.cla()
    plt.clf()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{tempreture} (MeV)')
    plt.ylabel(r'\textbf{\lambda}')

    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [nTOp.lambda_n__p(T) for T in Ts], 
        linewidth=1.0, label=r'$\lambda_{n\to p}$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [nTOp.lambda_p__n(T) for T in Ts],
        linewidth=1.0, label=r'$\lambda_{p\to n}$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [He_3_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [He_3_backward_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{back}$')

    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()


    ########## смотрим Wagoner, 1966 для сравнения ###########
    def pn(T):
        coef = 2.5*10**4
        rho_b = univ_func.__proton_mass_density__(T)
        return rho_b*coef / (constants.less_time(1))

    def lambda_d(T):
        T9 = constants.to_norm_tempreture(T, "T9")
        _pn_ = pn(T) / constants.to_norm_time(1)
        rho_b = univ_func.__proton_mass_density__(T)
        l = 4.68*10**9 * _pn_ * (rho_b**(-1)) * T9**(3./2) * math.exp(-25.82/T9)
        #l = math.exp(-25.82/T9) # otladka
        return l/constants.less_time(1)

    plt.cla()
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [pn(T) for T in Ts],
        linewidth=1.0, label=r'$H^2_{forw}(Wag)$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [He_3_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()



    plt.cla()
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [lambda_d(T) for T in Ts],
        linewidth=1.0, label=r'$H^2_{back}(Wag)$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [He_3_backward_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{bacl}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()