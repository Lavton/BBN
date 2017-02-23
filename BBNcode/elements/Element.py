# import elements.n as n
import nTOp

class Element():
    """elements"""

    def __init__(self, str_view, X_0):
        self.X_0 = X_0
        self.str_view = str_view
        self.X_c = X_0
        
    def ode_elem(T):
        return {self.str_view: -nTOp.lambda_n__p(T)*self.X_c}
        # self.X_c 


n = Element("n", 0.5)
H_1 = Element("H_1", 0.5)

n.ode_elem = lambda X, T: {
    n.str_view: -nTOp.lambda_n__p(T) * X[n.str_view]
}

H_1.ode_elem = lambda X, T: {
    n.str_view: lambda_p__n(T) * X[H_1.str_view],
    H_1.str_view: lambda_n__p(T) * X[n.str_view] - lambda_p__n(T) * X[H_1.str_view]
}

if __name__ == '__main__':
    h_1 = H_1()
    print(h_1.X_c)
    print(H_1.str_view)