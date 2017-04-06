import constants
import math
import nTOp

class Element():
    """
    Общий класс для всех элементов.
    у каждого есть начальное значение, строковое представление, массовое число и функции самого подсчёта
    """

    def __init__(self, str_view, X_0):
        self.X_0 = X_0
        self.str_view = str_view
        self.X_c = X_0
        self.ode_elem = {}
        self.A = 1
        self.jacob = {}    
        self.equilibrium = None
        self.is_ode_state = False # = in equilibrium


if __name__ == '__main__':
    h_1 = H_1()
    print(h_1.X_c)
    print(H_1.str_view)