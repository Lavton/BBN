
import constants
import math
from elements.Element import Element

import rates
def H_2_forw_rate(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    base_rate = 4.742 * 10**4 * (
        + 1. 
        - 0.8540 * T9**(1./2)
        + 0.4895 * T9
        - 0.09623 * T9**(3./2)
        + 0.008471 * T9**2
        - 0.00028 * T9**(5./2)
        )
    ro_b = rates.__proton_mass_density__(T)
    return base_rate * ro_b if T9 < 10 else 0

def H_2_backward_rate(T):
    T9 = constants.to_norm_tempreture(T, units="T9")
    forw = H_2_forw_rate(T)
    back = forw / (1.440*(10**-5)*(T9**(3./2))*constants.nu_n*math.exp(25.815/T9))
    return (back if T9 > 25.8/250 else 0) if T9 < 10 else 0


def H_2_equ(X, T):
    T9 = constants.to_norm_tempreture(-T, units="T9")
    X_n = 1.440*(10**-5)*(T9**(3./2))*constants.nu_n*math.exp(25.815/T9)*X[0][0]*X[0][1]
    X[0][2] = X_n
    return X

H_2 = Element("H_2", 0.0)
H_2.A = 2


H_2.tr_T = constants.less_tempreture(10**10, units="K")


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

# from Audi et all, 2003
H_2.mass_excess = constants.less_tempreture(13135721.6)