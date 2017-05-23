r"""
модуль для общёта скоростей реакций перехода
протонов в нейтроны и обратно
"""

import Cacher
import constants
import numpy as np
import math
from constants import *
from scipy import integrate
from math import pi, sqrt, exp
from tempreture import tfromT, TnuFromT


__lambda0__ = 1.63616
# __a0__ = (4*g**2)/(pi**3*(c*h*2*pi)**6*(2*pi*h))

def __x__(T):
    return m_e/constants.to_norm_tempreture(T)

def __x_nu__(T):
    return m_e/constants.to_norm_tempreture(TnuFromT(T))

def __q__():
    return Q/m_e

def __eps__(E_e):
    return E_e/m_e

__consta__ = (1/(constants.less_time(t_n)*__lambda0__))

def __under_int_lamda_n__p_e_nu(eps, x, x_nu, q):
    r"""
    подынтыгральное выражение из
    "Кинетика первичного НС".
    3.4a
    скорость процесса n->p+e+\nu
    """
    try:
        chis = eps*((eps-q)**2)*sqrt(eps**2-1)
        znam = ((1+exp(-eps*x))*(1+exp((eps-q)*x_nu)))
        res = chis/znam
    except OverflowError:
        res = 0.0
    return res


def __lamda_n__p_e_nu(T):
    r"""
    "Кинетика первичного НС".
    3.4a
    скорость процесса n->p+e+\nu
    """
    return __consta__*integrate.quad(lambda eps: __under_int_lamda_n__p_e_nu(eps, __x__(T), __x_nu__(T), __q__()), 1.0, __q__())[0]


def __under_int_lamda_nu_n__p_e(eps, x, x_nu, q):
    r"""
    подынтыгральное выражение из
    "Кинетика первичного НС".
    3.4б
    скорость процесса \nu+n->p+e
    """
    try:
        chis = eps*((eps-q)**2)*sqrt(eps**2-1)
        znam = ((1+exp(-eps*x))*(1+exp((eps-q)*x_nu)))
        res = chis/znam
    except OverflowError:
        res = 0.0
    return res


def __lamda_nu_n__p_e(T, max_y=np.inf):
    r"""
    "Кинетика первичного НС".
    3.4б
    скорость процесса \nu+n->p+e
    """
    return __consta__*integrate.quad(lambda eps: __under_int_lamda_nu_n__p_e(eps, __x__(T), __x_nu__(T), __q__()), __q__(), max_y)[0]

def __under_int_lamda_e_n__p_nu(eps, x, x_nu, q):
    r"""
    подынтыгральное выражение из
    "Кинетика первичного НС".
    3.4в
    скорость процесса e+n->p+\nu
    """
    try:
        chis = eps*((eps+q)**2)*sqrt(eps**2-1)
        znam = ((1+exp(eps*x))*(1+exp(-(eps+q)*x_nu)))
        res = chis/znam
    except OverflowError:
        res = 0.0
    return res


def __lamda_e_n__p_nu(T, max_y=np.inf):
    r"""
    "Кинетика первичного НС".
    3.4в
    скорость процесса e+n->p+\nu
    """
    return __consta__*integrate.quad(lambda eps: __under_int_lamda_e_n__p_nu(eps, __x__(T), __x_nu__(T), __q__()), 1.0, max_y)[0]


@Cacher.cacher.sql_base_cache
def lambda_n__p(T):
    r"""
    суммарная скорость перехода нейтронов в протоны, см 3.6а
    """
    return __lamda_n__p_e_nu(T)+__lamda_nu_n__p_e(T)+__lamda_e_n__p_nu(T)


###################################################################
def __under_int_lamda_p_e_nu__n(eps, x, x_nu, q):
    r"""
    подынтыгральное выражение из
    "Кинетика первичного НС".
    3.4г
    скорость процесса p+e+\nu->n
    """
    try:
        chis = eps*((eps-q)**2)*sqrt(eps**2-1)
        znam = ((1+exp(eps*x))*(1+exp((q-eps)*x_nu)))
        res = chis/znam
    except OverflowError:
        res = 0.0
    return res

def __lamda_p_e_nu__n(T):
    r"""
    "Кинетика первичного НС".
    3.4г
    скорость процесса p+e+\nu->n
    """
    return __consta__*integrate.quad(lambda eps: __under_int_lamda_p_e_nu__n(eps, __x__(T), __x_nu__(T), __q__()), 1.0, __q__())[0]

def __under_int_lamda_p_e__n_nu(eps, x, x_nu, q):
    r"""
    подынтыгральное выражение из
    "Кинетика первичного НС".
    3.4д
    скорость процесса p+e -> \nu+n
    """
    try:
        chis = eps*((eps-q)**2)*sqrt(eps**2-1)
        znam = ((1+exp(eps*x))*(1+exp((q-eps)*x_nu)))
        res = chis/znam
    except OverflowError:
        res = 0.0
    return res


def __lamda_p_e__n_nu(T, max_y=np.inf):
    r"""
    "Кинетика первичного НС".
    3.4д
    скорость процесса p+e -> \nu+n
    """
    return __consta__*integrate.quad(lambda eps: __under_int_lamda_p_e__n_nu(eps, __x__(T), __x_nu__(T), __q__()), __q__(), max_y)[0]

def __under_int_lamda_p_nu__e_n(eps, x, x_nu, q):
    r"""
    подынтыгральное выражение из
    "Кинетика первичного НС".
    3.4е
    скорость процесса p+\nu -> e+n
    """
    try:
        chis = eps*((eps+q)**2)*sqrt(eps**2-1)
        znam = ((1+exp(-eps*x))*(1+exp((eps+q)*x_nu)))
        res = chis/znam
    except OverflowError:
        res = 0.0
    return res


def __lamda_p_nu__e_n(T, max_y=np.inf):
    r"""
    "Кинетика первичного НС".
    3.4е
    скорость процесса p+\nu -> e+n
    """
    return __consta__*integrate.quad(lambda eps: __under_int_lamda_p_nu__e_n(eps, __x__(T), __x_nu__(T), __q__()), 1.0, max_y)[0]


@Cacher.cacher.sql_base_cache
def lambda_p__n(T):
    r"""
    суммарная скорость перехода протонов в нейтроны, см 3.6б
    """
    return __lamda_p_e_nu__n(T)+__lamda_p_e__n_nu(T)+__lamda_p_nu__e_n(T)


if __name__ == '__main__':
    import sys
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    if len(sys.argv)-1:
        T = sys.argv[1]
        T = eval(T)
        # exit()
    grid = np.logspace(math.log10(9.8*10**10), math.log10(10**7), num=160)
    grid2 = grid
    # grid2 = np.array(sorted(list(set(list(np.logspace(math.log10(grid[100]), math.log10(10**7), num=50))+list(grid))), reverse=True))
    Ts = constants.less_tempreture(grid2, units="K")
    # переводим в отрицательную шкалу, чтобы Ts[i] > Ts[i-1]
    ts = np.array([tfromT(T) for T in Ts])

    # # переход из нейтронов
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel(r'\textbf{tempreture} (K)')
    # plt.ylabel(r'\textbf{\lambda}')

    # plt.plot(Ts, [__lamda_n__p_e_nu(T) for T in Ts], 
    #     linewidth=2.0, label=r'$\lambda_{n\to p+e+\nu}$')
    # plt.plot(Ts, [__lamda_nu_n__p_e(T) for T in Ts], 
    #     linewidth=2.0, label=r'$\lambda_{n+\nu \to p+e}$')
    # plt.plot(Ts, [__lamda_e_n__p_nu(T) for T in Ts], 
    #     linewidth=2.0, label=r'$\lambda_{n+e \to p+\nu }$')
    # plt.plot(Ts, [lambda_n__p(T) for T in Ts],
    #     'r--', label=r'$\lambda_{n\to p}$')
    # plt.legend()
    # plt.gca().invert_xaxis()
    # plt.show()


    # # # переход из протонов
    # plt.cla()
    # plt.clf()
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel(r'\textbf{tempreture} (K)')
    # plt.ylabel(r'\textbf{\lambda}')

    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [__lamda_p_e_nu__n(T) for T in Ts], 
    #     linewidth=2.0, label=r'$\lambda_{p+e+\nu\to n}$')
    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [__lamda_p_e__n_nu(T) for T in Ts], 
    #     linewidth=2.0, label=r'$\lambda_{p+e\to n+\nu}$')
    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [__lamda_p_nu__e_n(T) for T in Ts], 
    #     linewidth=2.0, label=r'$\lambda_{p+\nu \to e+n}$')
    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [lambda_p__n(T) for T in Ts],
    #     'r--', label=r'$\lambda_{p\to n}$')
    # plt.legend()
    # plt.gca().invert_xaxis()
    # plt.show()

    # print ([(lambda_n__p(T), T) for T in Ts])
    # # сравнение скоростей
    # plt.cla()
    # plt.clf()
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel(r'\textbf{tempreture} (K)')
    # plt.ylabel(r'\textbf{\lambda}')

    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [lambda_n__p(T) for T in Ts], 
    #     linewidth=2.0, label=r'$\lambda_{n\to p}$')
    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [lambda_p__n(T) for T in Ts],
    #     linewidth=2.0, label=r'$\lambda_{p\to n}$')
    # plt.legend()
    # plt.gca().invert_xaxis()
    # plt.show()


    # асимптотика
    plt.cla()
    plt.clf()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    # plt.ylim([-3, 20])
    # plt.yscale('log')
    plt.xlabel(r'\textbf{tempreture} (K)')
    # plt.ylabel(r'\textbf{\lambda}')

    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [(lambda_p__n(T)/lambda_n__p(T))/exp(-constants.less_tempreture(Q)/T) for T in Ts], 
    #     linewidth=2.0)
    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [(lambda_p__n(T)/lambda_n__p(T)) for T in Ts], 
    #     linewidth=2.0)
    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [-constants.less_tempreture(Q, units="eV")/T for T in Ts], 
    #     linewidth=2.0)
    # plt.plot(constants.to_norm_tempreture(Ts, units="K"), [(lambda_p__n(T)/lambda_n__p(T))/exp(-constants.less_tempreture(Q, units="eV")/T) for T in Ts], 
    #     linewidth=2.0)
    ls = [(lambda_p__n(T)/lambda_n__p(T)) for T in Ts]
    l1s = [(lambda_p__n(T)) for T in Ts]
    l2s = [(lambda_n__p(T)) for T in Ts]
    eqs = [exp(-constants.less_tempreture(Q, units="eV")/T) for T in Ts]
    Is = [i for i in range(len(Ts))]
    # dis = [ls[i]/eqs[i] for i in range(len(Ts))]
    tms = [__x_nu__(T) for T in Ts]
    # import pickle 
    # with open("/tmp/new_result.pickle", "wb") as f:
        # pickle.dump((Is, constants.to_norm_tempreture(Ts, units="K"), l1s, l2s, ls, 
            # eqs, dis, tms),f)
    # for i in range(len(ls)):
        # print("i {} la {} exp {}: {}".format(i, ls[i], eqs[i], ls[i]/ eqs[i]))
    # plt.xlim([1e-5,1e3])
    plt.yscale("log")
    plt.plot(ts, l1s, 
        linewidth=2.0, label="lamdasPN")
    plt.plot(ts, l2s, 
        linewidth=2.0, label="lamdasNP")
    # plt.plot(ts, ls, 
        # linewidth=2.0, label="lamdas")
    # plt.plot(ts, eqs, 
        # linewidth=2.0, label="Q/k")
    # plt.plot(ts, dis, label="div")
    plt.legend()

    # plt.gca().invert_xaxis()
    plt.show()

    # print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    # print([(lambda_p__n(T)/lambda_n__p(T))/exp(-constants.less_tempreture(Q)/T) for T in Ts])