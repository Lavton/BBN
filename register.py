"""
объединение всех элементов
"""

import sys
import os
from collections import defaultdict, Counter
import elements.Element as el
import nTOp
import constants
import elements.n
from elements.n import n 
import elements.H_1
from elements.H_1 import H_1
import elements.H_2
from elements.H_2 import H_2
import elements.He_3
from elements.He_3 import He_3
import elements.H_3
from elements.H_3 import H_3
import elements.He_4
from elements.He_4 import He_4
import elements.Be_7
from elements.Be_7 import Be_7
import elements.Li_7
from elements.Li_7 import Li_7
import elements.Li_6
from elements.Li_6 import Li_6
import elements.He_6
from elements.He_6 import He_6
import tempreture
import logging
import functools
import Cacher

def check_jacob_online(X, T, res_jacob):
    def approx_jacobian(x,func,epsilon,*args):
        # from StackOverflow
        from numpy import asfarray, zeros
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
    apj = approx_jacobian(X, lambda X: np.array(registrator.sode_int(X,T)), 1e-8)
    total_s = 0
    for i in range(len(res_jacob)):
        total_s += sum([abs(res_jacob[i][j] - apj[i][j])/(max(abs(res_jacob[i][j]), abs(apj[i][j]))+1e-13) for j in range(len(res_jacob[i]))])
    total_s /= len(res_jacob)*len(res_jacob)
    if total_s >= 10e-1:
        logging.error(("BUG JACOB:", total_s, tempreture.tfromT(T), X))


class Registrator():
    def __init__(self):
        self.X_0 = []
        self.elements = []
        self.string_reactions = []
        self.reactions = []
        self.rev_element_list = {} # string -> order number
        self.ode_form = []


    def registrate(self, element):
        self.X_0.append(element.X_0)
        self.elements.append(element)
        self.rev_element_list[element.str_view] = len(self.elements) - 1
        self.reactions += element.reactions

    def _factorial(self, n):
        if n == 0: 
            return 1
        else: 
            return self._factorial(n-1)

    def finish_registration(self):
        # элементы, участвующие в реакции, добавляются с помощью регистратора
        # если какой-либо элемент так не добавлен
        # но присутствует в реакции - мы эту реакцию исключаем
        i_for_del = set()
        i = 0
        for r in self.reactions:
            left, right, forw, backw = r
            self.string_reactions.append(" + ".join(left) + " <--> " + " + ".join(right))
            for l in left:
                if self.get_num_on_name(l) == -1:
                    i_for_del.add(i)
            for r in right:
                if self.get_num_on_name(r) == -1:
                    i_for_del.add(i)
            i += 1
        i_for_del = sorted(list(i_for_del))
        for i in i_for_del:
            self.reactions.pop(i)
        
        dXss = defaultdict(list)
        for r in self.reactions:
            left, right, forw, backw = r
            c_left = Counter(left)
            c_right = Counter(right)

            # хотим немного магии: записываем нужные нам вещи в строковом представлении
            forw2 = forw.__module__+"."+forw.__name__
            backw2 = backw.__module__+"."+backw.__name__
            lefts = left
            rights = right
            # объединяются элементы, ответственные за распад при прямой реакции
            forw_part = []
            for l in lefts:
                forw_part.append(l)
            factor = self._factorial(max(c_left.values()))
            forw_part = (forw_part, forw2 + "(T)/"+str(factor))
            # объединяются элементы, ответственные за синтез при прямой реакции
            back_part = []
            for l in rights:
                back_part.append(l)
            factor = self._factorial(max(c_right.values()))
            back_part = (back_part, backw2 + "(T)/"+str(factor))
            
            # для прямой реакции распад идёт со знаком "-"
            # синтез со знаком "+"
            for left in c_left:
                dXss[self.get_num_on_name(left)].append((" - ", self._factorial(c_left[left]), forw_part))
                dXss[self.get_num_on_name(left)].append((" + ", self._factorial(c_left[left]), back_part))
            # для обратной реакции - распад с "+", синтез с "-"
            for right in c_right:
                dXss[self.get_num_on_name(right)].append((" - ", self._factorial(c_right[right]), back_part))
                dXss[self.get_num_on_name(right)].append((" + ", self._factorial(c_right[right]), forw_part))
        dX = [[] for _ in range(max(dXss.keys()) + 1)]
        for k in dXss:
            dX[k] = dXss[k]
        equation = []
        for i in range(len(dX)):
            lines = "0"
            for j in range(len(dX[i])):
                lines += "{sign}{coef}*{seq}*{rate}".format(
                    sign=dX[i][j][0], 
                    coef=dX[i][j][1],
                    seq=self._to_ode(dX[i][j][2][0]),
                    rate=dX[i][j][2][1]
                    )
            equation.append("my_dx[{}] = ({})*{}.A".format(i, lines, self.elements[i].names[0]))
            logging.debug("add '{}' to the equation formula".format(equation[-1]))

        # получаем функцию для диффура
        ode_line = (
"""
def _magic_ode(X, T):
    my_dx = [0 for _ in range(len(X))]
    
    {}
    
    return my_dx
    
self._magic_ode = _magic_ode
""".format("\n    ".join(equation))
    )

        exec(ode_line)
        logging.info(ode_line)

        ############
        # аналогично для jacob
        equations = []
        for i in range(len(dX)):
            equ_jac = ["0" for _ in range(len(dX))]
            for j in range(len(dX[i])):
                pre_jac = self._to_jacob(dX[i][j][2][0])
                for prj in pre_jac:
                    equ_jac[self.get_num_on_name(prj[0])] += "{sign}{coef}*{seq}*{rate}*{el}.A".format(
                    sign=dX[i][j][0], 
                    coef=dX[i][j][1],
                    seq=prj[1],
                    rate=dX[i][j][2][1],
                    el=self.elements[i].names[0]
                )
            equations.append(equ_jac)
        jac_equations = []
        for i in range(len(equations)):
            for j in range(len(equations[i])):
                jac_equations.append("my_j[{}][{}] = {}".format(i, j, equations[i][j]))

        jacob_line = (
"""
def _magic_jacob(X, T):
    my_j = [[0 for _ in range(len(X))] for _ in range(len(X))]
    {}
    
    return my_j

self._magic_jacob = _magic_jacob
""".format("\n    ".join(jac_equations))
    )

        exec(jacob_line)
        logging.info(jacob_line)

        logging.info("elements: \n{}".format(
            "\n".join([e.names[0] for e in self.elements])
            ))
        logging.info("reactions: \n{}".format("\n".join(self.string_reactions)))


    def _to_ode(self, seq):
        return "*".join(map(self._get_seq, seq))

    def _to_jacob(self, seq):
        c_seq = Counter(seq)
        res = []
        for c_s in c_seq:
            part = []
            for el in c_seq:
                coef = self._factorial(c_seq[c_s]) if el==c_s else 1
                power = coef - 1 if el==c_s else 1
                n_line = "({}*{}**{})".format(coef, self._get_seq(el), power)
                if el==c_s:
                    n_line += "*(1.0/{}.A)".format(
                        self.elements[self.get_num_on_name(el)].names[0]
                        )
                part.append(n_line)
            res.append((c_s, "*".join(part)))
        return res

    def _get_seq(self, el):
        num = self.get_num_on_name(el)
        return "(X[{}]/{}.A)".format(num, self.elements[num].names[0])

    def get_num_on_name(self, name):
        for j in range(len(self.elements)):
            if name in self.elements[j].names:
                return j
        logging.error("no such element name")
        return -1


    def sode_int(self, X, T):
        dX = self._magic_ode(X, T)
        return dX

    def jacob(self, X, T):
        res_jacob = self._magic_jacob(X, T)
        if logging.root.level==logging.DEBUG:
            check_jacob_online(X, T, res_jacob)
        logging.debug("t = {}, X = {}".format(tempreture.tfromT(T), X))
        return res_jacob

    def calc_plot(self, plt, ts, X_ans, num_of_el=0):
        if not num_of_el:
            num_of_el = len(self.elements)

        for i in range(num_of_el):
            plt.plot(ts, X_ans[:,i], 
            linewidth=1.8, label="$"+self.elements[i].str_view+"$")


registrator = Registrator()
registrator.registrate(n)
registrator.registrate(H_1)
registrator.registrate(H_2)
registrator.registrate(He_3)
registrator.registrate(H_3)
registrator.registrate(He_4)
registrator.registrate(Be_7)
registrator.registrate(Li_7)
registrator.registrate(Li_6)
# registrator.registrate(He_6)
registrator.finish_registration()
X_0 = registrator.X_0

if __name__ == '__main__':
    def check_jacob():
        from scipy.misc import derivative
        from tempreture import Tfromt
        st = "[[  1.33289087e-01   8.66677827e-01   1.45546654e-04]] i = 159/320 0.201898823095 raz 0.000112460505732"
        st = st.replace("]","").split()
        X = [ 2.12487653e-01,   7.87512347e-01,   8.68469663e-13, -4.06999347e-19,   2.97767158e-16]
         # [float(st[1]), float(st[2]), float(st[3])]
        # X = [1.88997594e-01, 8.11002406e-01, 1.70829884e-12]
        # X = [0.3, 0.3, 0.3]
        X[0] = 1-sum(X[1:])
        # t = float(st[7])
        # 
        t = 0.0017
        T = Tfromt(t)
        import numpy as np
        jaaaa = np.array(registrator.jacob(X, T))
        odeeee = registrator.sode_int(X, T)
        print(X, t, T)
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
            print([(jaaaa[i][j] - ap_j[i][j])/(max(abs(jaaaa[i][j]), abs(ap_j[i][j]))+1e-13) for j in range(len(jaaaa[i]))])

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
