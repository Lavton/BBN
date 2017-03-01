# import elements.n as n
import constants
import math
import nTOp

class Element():
    """elements"""

    def __init__(self, str_view, X_0):
        self.X_0 = X_0
        self.str_view = str_view
        self.X_c = X_0
        self.ode_elem = {}
        self.jacob = {}    
        self.equilibrium = None

n = Element("n", 0.5)
H_1 = Element("H_1", 0.5)
H_2 = Element("H_2", 0.0)

n.ode_elem = {
    n.str_view : (lambda X, T: -nTOp.lambda_n__p(T) * X[0])
}


H_1.ode_elem = {
    n.str_view: (lambda X, T: nTOp.lambda_p__n(T) * X[1]),
    H_1.str_view: (lambda X, T: nTOp.lambda_n__p(T) * X[0] - nTOp.lambda_p__n(T) * X[1])
}


n.jacob = {
    n.str_view: {
        n.str_view: (lambda X, T: -nTOp.lambda_n__p(T))
    }
}


H_1.jacob = {
    n.str_view: {
        H_1.str_view: (lambda X, T: +nTOp.lambda_p__n(T))
    },
    H_1.str_view: {
        n.str_view: (lambda X, T: +nTOp.lambda_n__p(T)),
        H_1.str_view: (lambda X, T: -nTOp.lambda_p__n(T))
    }
}

H_2.ode_elem = {
    
}

H_2.jacob = {
    
}
def H_2_equ(X, T):
    T9 = constants.to_norm_tempreture(-T, units="T9")
    X_n = 1.440*(10**-5)*(T9**(3./2))*constants.nu_n*math.exp(25.815/T9)*X[0][0]*X[0][1]
    X[0][2] = X_n
    return X

H_2.equilibrium = H_2_equ

if __name__ == '__main__':
    h_1 = H_1()
    print(h_1.X_c)
    print(H_1.str_view)