"""
объединение всех элементов
"""

import sys
import os
from collections import defaultdict
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)))
import elements.Element as el
import nTOp
import constants
from elements.n import n 
from elements.H_1 import H_1
from elements.H_2 import H_2
from elements.He_3 import He_3
import tempreture
import logging
import functools
import Cacher

def check_jacob_online(X, T, res_jacob):
    def approx_jacobian(x,func,epsilon,*args):
        from numpy import asfarray, zeros
        """Approximate the Jacobian matrix of callable function func

           * Parameters
             x       - The state vector at which the Jacobian matrix is
    desired
             func    - A vector-valued function of the form f(x,*args)
             epsilon - The peturbation used to determine the partial derivatives
             *args   - Additional arguments passed to func

           * Returns
             An array of dimensions (lenf, lenx) where lenf is the length
             of the outputs of func, and lenx is the number of

           * Notes
             The approximation is done using forward differences

        """
        x0 = asfarray(x)
        f0 = func(*((x0,)+args))
        jac = zeros([len(x0),len(f0)])
        dx = zeros(len(x0))
        for i in range(len(x0)):
           dx[i] = epsilon
           jac[i] = (func(*((x0+dx,)+args)) - f0)/epsilon
           dx[i] = 0.0
        return jac.transpose()
    import numpy as np
    import tempreture
    apj = approx_jacobian(X, lambda X: np.array(registrator.sode_int(X,T)), 1e-6)
    # print(ap_j)
    total_s = 0
    for i in range(len(res_jacob)):
        total_s += sum([abs(res_jacob[i][j] - apj[i][j])/(1+abs(apj[i][j])) for j in range(len(res_jacob[i]))])
    if total_s >= 1e-2:
        logging.error(("BUG JACOB:", total_s, tempreture.tfromT(T), X))


class Registrator():
    def __init__(self):
        self.X_0 = []
        self.elements = []
        self.rev_element_list = {} # string -> order number
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

    # @Cacher.cacher.np_array_to_list_decor
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
        logging.debug("t = {}, X = {}".format(tempreture.tfromT(T), X))
        return dX

    # @Cacher.cacher.np_array_to_list_decor
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
        if logging.root.level==logging.DEBUG:
            check_jacob_online(X, T, res_jacob)
        logging.debug("t = {}, X = {}".format(tempreture.tfromT(T), X))
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
registrator.registrate(He_3)
registrator.finish_registration()

X_0 = registrator.X_0

if __name__ == '__main__':
    def check_jacob():
        from scipy.misc import derivative
        from tempreture import Tfromt
        st = "[[  1.33289087e-01   8.66677827e-01   1.45546654e-04]] i = 159/320 0.201898823095 raz 0.000112460505732"
        st = st.replace("]","").split()
        X = [  1.67386032e-01,   8.32614019e-01,  -5.11260695e-08,-2.21515705e-12]
         # [float(st[1]), float(st[2]), float(st[3])]
        # X = [1.88997594e-01, 8.11002406e-01, 1.70829884e-12]
        # X = [0.3, 0.3, 0.3]
        X[0] = 1-X[1]-X[2]-X[3]
        # t = float(st[7])
        # 
        t = 0.011797924341382172
        T = Tfromt(t)
        import numpy as np
        jaaaa = np.array(registrator.jacob(X, T))
        odeeee = registrator.sode_int(X, T)
        print("analitic")
        print(jaaaa)
        print("func")
        print(odeeee)

        def approx_jacobian(x,func,epsilon,*args):
            from numpy import asfarray, zeros
            """Approximate the Jacobian matrix of callable function func

               * Parameters
                 x       - The state vector at which the Jacobian matrix is
        desired
                 func    - A vector-valued function of the form f(x,*args)
                 epsilon - The peturbation used to determine the partial derivatives
                 *args   - Additional arguments passed to func

               * Returns
                 An array of dimensions (lenf, lenx) where lenf is the length
                 of the outputs of func, and lenx is the number of

               * Notes
                 The approximation is done using forward differences

            """
            x0 = asfarray(x)
            f0 = func(*((x0,)+args))
            jac = zeros([len(x0),len(f0)])
            dx = zeros(len(x0))
            for i in range(len(x0)):
               dx[i] = epsilon
               jac[i] = (func(*((x0+dx,)+args)) - f0)/epsilon
               dx[i] = 0.0
            return jac.transpose()

        # for i in range(len(odeeee)):
            # for j in range(len(X)):
                # print(derivative())
        X = np.array(X)
        ap_j = approx_jacobian(X, lambda X: np.array(registrator.sode_int(X,T)), 1e-10)
        print("approx")
        print(ap_j)

        for i in range(len(jaaaa)):
            print([(jaaaa[i][j] - ap_j[i][j])/(1+abs(ap_j[i][j])) for j in range(len(jaaaa[i]))])

    check_jacob()
    exit()
    import numpy as np
    import math
    import constants
    import nTOp
    from tempreture import tfromT
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=160)
    grid2 = grid
    # grid2 = np.array(sorted(list(set(list(np.logspace(math.log10(grid[100]), math.log10(10**7), num=50))+list(grid))), reverse=True))
    Ts = constants.less_tempreture(grid2, units="K")
    # переводим в отрицательную шкалу, чтобы Ts[i] > Ts[i-1]
    ts = np.array([tfromT(T) for T in Ts])
    plt.cla()
    plt.clf()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1e-5,1e10])
    plt.xlabel(r'\textbf{tempreture}')
    plt.ylabel(r'\textbf{\lambda}')

    plt.plot(ts, [nTOp.lambda_n__p(T) for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{n\to p}$')
    plt.plot(ts, [nTOp.lambda_p__n(T) for T in Ts],
        linewidth=2.0, label=r'$\lambda_{p\to n}$')
    import elements.H_2 as eH2
    plt.plot(ts, [eH2.H_2_forw_rate(T) for T in Ts],
        linewidth=1.0, label=r'$p,n \to d,\gamma$')
    plt.plot(ts, [eH2.H_2_backward_rate(T) for T in Ts],
        linewidth=1.0, label=r'$\lambda_{d\to p,n}$')
    plt.legend()
    # plt.gca().invert_xaxis()
    plt.show()
