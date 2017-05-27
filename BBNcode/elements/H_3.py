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

H_3 = Element("H_3", 0.0)
H_3.A = 3
# from Audi et all, 2003
# H_3.mass_excess = constants.less_tempreture(2161062.7, units="eV")
H_3.set_mass_excess(3016049.2777, n_N=2, p_N=1)
H_3.tr_t =  0.012
H_3.tr_T = tempreture.Tfromt(H_3.tr_t)
# H_3.tr_T = constants.less_tempreture(2*10**10, units="K")

@H_3.equilib_zeroize
def nd_tg(T):
    """
    Wagoner
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = (
        75.5
        + 1250 * T9
        )
    ro_b = univ_func.rat_scale(T)
    return base_rate * ro_b/(constants.less_time(1))

@H_3.equilib_zeroize
def tg_nd(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = nd_tg.__wrapped__(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    # back = forw * math.exp(-H_3.mass_excess/T)
    ro_b = univ_func.rat_scale(T)

    back = 1.63*10**10 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-72.62/T9)
    return (back /(constants.less_time(1)))


def H_3_equ(X, T):
    T9 = constants.to_norm_tempreture(-T, units="T9")
    try:
        # print(T, he3_gp_d.__wrapped__(T))
        X_h3 = (3./2) * (nd_tg.__wrapped__(T)/tg_nd.__wrapped__(T))*X[0][0]*X[0][2]
    except OverflowError as e:
        X_h3 = 0
    X[0][4] = X_h3
    return X

H_3.forward_rates.append(nd_tg)
H_3.backward_rates.append(tg_nd)


# 0 - n
# 1 - H1
# 2 - H2
# 3 - He3
# 4 - H3
H_3.ode_elem = {
    "n": (lambda X, T: -X[0]*(X[2]/2)*nd_tg(T) + (X[4]/3)*tg_nd(T)),
    "H_2": (lambda X, T:  -X[0]*(X[2]/2)*nd_tg(T) + (X[4]/3)*tg_nd(T)),
    "H_3": (lambda X, T: X[0]*(X[2]/2)*nd_tg(T) - (X[4]/3)*tg_nd(T))
}

H_3.jacob = { 
    "n": { #берём первую строку и дифференцируем по каждому
        "n": (lambda X, T: -(X[2]/2)*nd_tg(T)),
        "H_2": (lambda X, T: -X[0]*(1./2)*nd_tg(T)),
        "H_3": (lambda X, T: (1./3)*tg_nd(T))
    },
    "H_2": {
        "n": (lambda X, T:  -(X[2]/2)*nd_tg(T)),
        "H_2": (lambda X, T: -X[0]*(1./2)*nd_tg(T)),
        "H_3": (lambda X, T: (1./3)*tg_nd(T))
    },
    "H_3": {
        "n": (lambda X, T: +(X[2]/2)*nd_tg(T)),
        "H_2": (lambda X, T: +X[0]*(1./2)*nd_tg(T)),
        "H_3": (lambda X, T: -(1./3)*tg_nd(T))
    }
}

H_3.equilibrium = H_3_equ


if __name__ == '__main__':
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

    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_3_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_3_backward_rate(T) for T in Ts], 
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
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_3_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_3_backward_rate(T) for T in Ts], 
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
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_3_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()



    plt.cla()
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [lambda_d(T) for T in Ts],
        linewidth=1.0, label=r'$H^2_{back}(Wag)$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_3_backward_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{bacl}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()