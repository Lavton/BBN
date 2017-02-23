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
        self.ode_form = []


    def registrate(self, element):
        self.X_0.append(element.X_0)
        # for i in range(len(self.elements)):
            # self.ode_form[i] 
        self.elements.append(element)

    # def finish_registration(self):
        # for i in range(len(self.elements)):
            # element = self.elements[i]
            # self.X_0.append(element.X_0)

            # def o_form(X, T):
                # dX = defaultdict(lambda : 0)
                # for j in range(i):
                    # dX[element.str_view] += element.ode_elem(X, T).get(element.str_view, 0)
                # return dX
            # self.ode_form.append(o_form)

    def ode_int(X, T):
        dX = defaultdict(lambda: 0)
        X_new = defaultdict(lambda: 0)
        for i in range(len(self.elements)):
            X_new[self.elements[i].str_view] = X[i]
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

# registrator.finish_registration()

X_0 = registrator.X_0