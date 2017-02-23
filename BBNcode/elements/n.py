import nTOp

class n_class():
    """neitrons"""

    X_c = 0.5
    str_view = "n"

    def __init__(self):
        self.X_0 = 0.5
        
    def ode_elem(T):
        return {self.str_view: -nTOp.lambda_n__p(T)*self.X_c}
        # self.X_c 



if __name__ == '__main__':
    n = n_class()
    print(n.X_c)