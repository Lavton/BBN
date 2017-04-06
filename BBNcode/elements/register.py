"""
объединение всех элементов
"""

import sys
import os
from collections import defaultdict
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)))
import elements.Element as el
from elements.n import n 
from elements.H_1 import H_1
from elements.H_2 import H_2



class Registrator():
    def __init__(self):
        self.X_0 = []
        self.elements = []
        self.rev_element_list = {}
        self.ode_form = []


    def registrate(self, element):
        self.X_0.append(element.X_0)
        self.elements.append(element)
        self.rev_element_list[element.str_view] = len(self.elements) - 1

    def finish_registration(self):
        self.ode_funcs = []
        for i in range(len(self.elements)):
            element = self.elements[i]
            self.ode_funcs.append([])
            for key, value in element.ode_elem.items():
                self.ode_funcs[self.rev_element_list[key]].append(value)

        self.jacob_funcs = [[[lambda X, T: 0] for _ in self.elements] for _ in self.elements]
        for i in range(len(self.elements)):
            element = self.elements[i]
            for key, value in element.jacob.items():
                for k, v in value.items():
                    self.jacob_funcs[self.rev_element_list[key]][self.rev_element_list[k]].append(v)

    def sode_int(self, X, T):
        """
        к этому моменту 
        self.ode_funcs = [
        [lambda X, T: ... ,lambda X, T: ... ,lambda X, T: ... ,],
        [lambda X, T: ... ,lambda X, T: ... ,lambda X, T: ... ,],
        [lambda X, T: ... ,lambda X, T: ... ,lambda X, T: ... ,]
        ]
        """
        dX = []
        for i in range(len(self.elements)):
            ode_f = self.ode_funcs[i]
            dX.append(sum(map(lambda f: f(X, T), ode_f)) * self.elements[i].A)
        return dX

    def jacob(self, X, T):
        """
        к этому моменту 
        self.jacob_funcs = [
        [lambda X, T: [lambda X,T:],[] , lambda X, T: [lambda X,T:],[] ,],
        [lambda X, T: [lambda X,T:],[] , lambda X, T: [lambda X,T:],[] ,],
        [lambda X, T: [lambda X,T:],[] , lambda X, T: [lambda X,T:],[] ,],
        ]
        """
        res_jacob = []
        for i in range(len(self.elements)):
            j = self.jacob_funcs[i]
            jacob_row = []
            for J in j:
                jacob_row.append(sum(map(lambda f: f(X, T), J)) * self.elements[i].A)
            res_jacob.append(jacob_row)
        return res_jacob

    def calc_plot(self, plt, ts, X_ans, num_of_el=0):
        if not num_of_el:
            num_of_el = len(self.elements)

        for i in range(num_of_el):
            plt.plot(ts, X_ans[:,i], 
            linewidth=2.0, label="$"+self.elements[i].str_view+"$")


registrator = Registrator()
registrator.registrate(n)
registrator.registrate(H_1)
registrator.registrate(H_2)

registrator.finish_registration()

X_0 = registrator.X_0