import elements.n as n
import nTOp

class H_1():
    """protons"""

    X_c = 0.5
    str_view = "p"

    def __init__(self):
        self.X_0 = 0.5
        
    def ode_elem(T):
        return {self.str_view: -nTOp.lambda_n__p(T)*self.X_c}
        # self.X_c 



if __name__ == '__main__':
    h_1 = H_1()
    print(h_1.X_c)
    print(H_1.str_view)