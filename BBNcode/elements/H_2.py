"""
дейтерий, или $^{2}H$
"""

if __name__ == '__main__':
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)))

import constants
import math
import tempreture
from elements.Element import Element

import univ_func

H_2 = Element("H_2", 0.0)
H_2.A = 2
# from Audi et all, 2003
H_2.mass_excess = constants.less_tempreture(2161062.7, units="eV")
H_2.tr_T = constants.less_tempreture(2.5*10**9, units="K")
H_2.tr_t = tempreture.tfromT(H_2.tr_T)
print(H_2.tr_t)
# exit()


def H_2_forw_rate(T):
    """
    Smith et all
    """
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 4.742 * 10**4 * (
        + 1. 
        - 0.8540 * T9**(1./2)
        + 0.4895 * T9
        - 0.09623 * T9**(3./2)
        + 8.471*1e-3 * T9**2
        - 2.80*1e-4 * T9**(5./2)
        )
    ro_b = univ_func.rat_scale(T)
    # base_rate = 0
    return base_rate * ro_b/(constants.less_time(1)) if T < H_2.tr_T else 0

def H_2_backward_rate(T):
    """Wagoner, 1966"""
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H_2_forw_rate(T) / constants.to_norm_time(1)
    E = constants.to_norm_tempreture(T, units="MeV")
    back = forw * math.exp(-H_2.mass_excess/T)
    ro_b = univ_func.rat_scale(T)

    back = 4.68*10**9 * forw * (ro_b**(-1)) * T9**(3./2) * math.exp(-H_2.mass_excess/T)
    return (back /(constants.less_time(1))) if T < H_2.tr_T else 0


def H_2_equ(X, T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    # print('IN EQ', T9)
    try:
        X_n = 1.440*(10**-5)*(T9**(3./2))*constants.nu_n*math.exp(25.815/T9)*X[0][0]*X[0][1]
    except OverflowError as e:
        X_n = 0
    X[0][2] = X_n
    return X


H_2.forward_rates.append(H_2_forw_rate)
H_2.backward_rates.append(H_2_backward_rate)


H_2.ode_elem = {
    "n": (lambda X, T: -X[0]*X[1]*H_2_forw_rate(T) + X[2]*H_2_backward_rate(T)/2),
    "H_1": (lambda X, T: -X[0]*X[1]*H_2_forw_rate(T) + X[2]*H_2_backward_rate(T)/2),
    "H_2": (lambda X, T: +X[0]*X[1]*H_2_forw_rate(T) - X[2]*H_2_backward_rate(T)/2)
}

H_2.jacob = { 
    "n": { #берём первую строку и дифференцируем по каждому
        "n": (lambda X, T: -X[1]*H_2_forw_rate(T)),
        "H_1": (lambda X, T: -X[0]*H_2_forw_rate(T)),
        "H_2": (lambda X, T: +H_2_backward_rate(T)/2)
    },
    "H_1": {
        "n": (lambda X, T:  -X[1]*H_2_forw_rate(T)),
        "H_1": (lambda X, T: -X[0]*H_2_forw_rate(T)),
        "H_2": (lambda X, T: +H_2_backward_rate(T)/2)
    },
    "H_2": {
        "n": (lambda X, T: X[1]*H_2_forw_rate(T)),
        "H_1": (lambda X, T: X[0]*H_2_forw_rate(T)),
        "H_2": (lambda X, T: - H_2_backward_rate(T)/2)
    }
}

H_2.equilibrium = H_2_equ


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

    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_2_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_2_backward_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{back}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()
    H_2.show_rates()
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
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_2_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_2_backward_rate(T) for T in Ts], 
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
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [pn(T) for T in Ts],
        linewidth=1.0, label=r'$H^2_{forw}(Wag)$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_2_forw_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{forw}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()



    plt.cla()
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [lambda_d(T) for T in Ts],
        linewidth=1.0, label=r'$H^2_{back}(Wag)$')
    plt.plot(constants.to_norm_tempreture(Ts, units="eV")*1e-6, [H_2_backward_rate(T) for T in Ts], 
        linewidth=2.0, label=r'$H^2_{bacl}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()