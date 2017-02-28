# import elements.n as n
import nTOp

class Element():
    """elements"""

    def __init__(self, str_view, X_0):
        self.X_0 = X_0
        self.str_view = str_view
        self.X_c = X_0
        self.ode_elem = {}
        self.jacob = {}    

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

if __name__ == '__main__':
    h_1 = H_1()
    print(h_1.X_c)
    print(H_1.str_view)