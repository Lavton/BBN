# import elements.n as n
import constants
import math
import nTOp
# from elements.H_2 import H_2_forw_rate, H_2_equ, H_2_backward_rate

class Element():
    """elements"""

    def __init__(self, str_view, X_0):
        self.X_0 = X_0
        self.str_view = str_view
        self.X_c = X_0
        self.ode_elem = {}
        self.A = 1
        self.jacob = {}    
        self.equilibrium = None



if __name__ == '__main__':
    h_1 = H_1()
    print(h_1.X_c)
    print(H_1.str_view)