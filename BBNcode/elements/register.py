import sys
import os
from collections import defaultdict
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)))


# import elements.n as n
# import elements.H_1 as H_1
import elements.Element as el

# n = n.n_class()
# H_1 = H_1.H_1()

# X_0 = [n.X_0, H_1.X_0]


class Registrator():
    """docstring for Registrator"""
    def __init__(self):
        self.X_0 = []
        self.elements = []
        self.rev_element_list = {}
        self.ode_form = []


    def registrate(self, element):
        self.X_0.append(element.X_0)
        # for i in range(len(self.elements)):
            # self.ode_form[i] 
        self.elements.append(element)
        self.rev_element_list[element.str_view] = len(self.elements) - 1

    def finish_registration(self):
        self.ode_funcs = []
        for i in range(len(self.elements)):
            element = self.elements[i]
            # ode_func = []
            self.ode_funcs.append([])
            for key, value in element.ode_elem.items():
                self.ode_funcs[self.rev_element_list[key]].append(value)
            # self.X_0.append(element.X_0)
        for i in range(len(self.ode_funcs)):
            self.ode_funcs[i] = lambda X, T: sum(map(lambda f: f(X, T), self.ode_funcs[i]))

        self.jacob_funcs = [[[lambda X, T: 0] for _ in self.elements] for _ in self.elements]
        for i in range(len(self.elements)):
            element = self.elements[i]
            for key, value in element.jacob.items():
                for k, v in value.items():
                    self.jacob_funcs[self.rev_element_list[key]][self.rev_element_list[k]].append(v)
        for i in range(len(self.jacob_funcs)):
            for j in range(len(self.jacob_funcs[i])):
                self.jacob_funcs[i][j] = lambda X, T: sum(map(lambda f: f(X, T), self.jacob_funcs[i][j]))

    def ode_int(X, T):

        # dX = defaultdict(lambda: 0)
        # X_new = defaultdict(lambda: 0)
        dX = []
        X_new = []
        # for i in range(len(self.elements)):
            # X_new.a = X[i]
        for element in self.elements:
            for equ in element.ode_elem(X, T)
                dX[element.str_view] += equ
        f_dX = []
        for i in range(len(self.elements)):
            f_dX.append(dX[self.elements[i]])
        return f_dX

registrator = Registrator()
registrator.registrate(el.n)
registrator.registrate(el.H_1)

registrator.finish_registration()

X_0 = registrator.X_0