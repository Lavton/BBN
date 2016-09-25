r"""
модуль для общёта скоростей реакций перехода
протонов в нейтроны и обратно
"""

import numpy as np
import math
from constants import *
from scipy import integrate
from math import pi, sqrt, exp
from tempreture import tfromT, TnuFromT
import functools

if sql_enabled: # если кеширование в БД есть
    import sqlite3
    db = sqlite3.connect('cache.db')
    cur = db.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Lambdas(tempeture REAL, speed REAL NULL, title VARCHAR(32) NULL);")


def sql_lambda_cache(func):
    """
    кеширование скоростей протекания реакций
    """
    @functools.wraps(func)
    def inner(*args, **kwargs):
        coef = 1.0
        if kwargs.get("units", "") == "K":
            coef = k_b # если температура была дана в Кельвинах, переводин в eV
        cached = cur.execute("SELECT speed FROM Lambdas WHERE tempeture={temp} AND title=\"{title}\";".format(temp=args[0]*coef, title=func.__name__))
        res = cur.fetchone() # ищем результат в кеше
        if res: # если нашли температуру
            return res[0]

        if not res: # не нашли - вставим
            res = func(*args, **kwargs)
            cur.execute("INSERT INTO Lambdas (tempeture, speed, title) VALUES ({temp}, {speed_v}, \"{title_v}\");".
                format(title_v=func.__name__, temp=args[0]*coef, speed_v=res))
            db.commit()
            return res

    return inner if sql_enabled else func


__lambda0__ = 1.63616
# __a0__ = (4*g**2)/(pi**3*(c*h*2*pi)**6*(2*pi*h))

def __x__(T):
    return m_e/T

def __x_nu__(T):
    return m_e/TnuFromT(T)

def __q__():
    return Q/m_e

def __eps__(E_e):
    return E_e/m_e

__consta__ = (1/(t_n*__lambda0__))

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


@sql_lambda_cache
def lambda_n__p(T, *, units="eV"):
    r"""
    суммарная скорость перехода нейтронов в протоны, см 3.6а
    """
    if units=="K":
        T = T*k_b
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


@sql_lambda_cache
def lambda_p__n(T, *, units="eV"):
    r"""
    суммарная скорость перехода протонов в нейтроны, см 3.6б
    """
    if units=="K":
        T = T*k_b
    return __lamda_p_e_nu__n(T)+__lamda_p_e__n_nu(T)+__lamda_p_nu__e_n(T)





# def __Kp_under_int__(T, T_nu, eps):

# def __Kp__(T, T_nu, max_y=np.inf):
    # return integrate.quad(lambda eps: __Kp_under_int__(T, T_nu, eps), 1, max_y)
if __name__ == '__main__':
    import sys
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    if len(sys.argv)-1:
        T = sys.argv[1]
        T = eval(T)
        # print("{:.1E}: {:.4E}\n".format(T, tfromT(T, units="K")))
        # exit()
    Ts = np.logspace(math.log10(10**8), math.log10(10**11), num=200)

    # переход из нейтронов
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{tempreture} (K)')
    plt.ylabel(r'\textbf{\lambda} (c^{-1})')

    plt.plot(Ts, [__lamda_n__p_e_nu(T*k_b) for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{n\to p+e+\nu}$')
    plt.plot(Ts, [__lamda_nu_n__p_e(T*k_b) for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{n+\nu \to p+e}$')
    plt.plot(Ts, [__lamda_e_n__p_nu(T*k_b) for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{n+e \to p+\nu }$')
    plt.plot(Ts, [lambda_n__p(T, units="K") for T in Ts],
        'r--', label=r'$\lambda_{n\to p}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()


    # переход из протонов
    plt.cla()
    plt.clf()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{tempreture} (K)')
    plt.ylabel(r'\textbf{\lambda} (c^{-1})')

    plt.plot(Ts, [__lamda_p_e_nu__n(T*k_b) for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{p+e+\nu\to n}$')
    plt.plot(Ts, [__lamda_p_e__n_nu(T*k_b) for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{p+e\to n+\nu}$')
    plt.plot(Ts, [__lamda_p_nu__e_n(T*k_b) for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{p+\nu \to e+n}$')
    plt.plot(Ts, [lambda_p__n(T, units="K") for T in Ts],
        'r--', label=r'$\lambda_{p\to n}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()

    print ([(lambda_n__p(T, units="K"), T) for T in Ts])
    # сравнение скоростей
    plt.cla()
    plt.clf()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{tempreture} (K)')
    plt.ylabel(r'\textbf{\lambda} (c^{-1})')

    plt.plot(Ts, [lambda_n__p(T, units="K") for T in Ts], 
        linewidth=2.0, label=r'$\lambda_{n\to p}$')
    plt.plot(Ts, [lambda_p__n(T, units="K") for T in Ts],
        linewidth=2.0, label=r'$\lambda_{p\to n}$')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()


    # асимптотика
    plt.cla()
    plt.clf()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{tempreture} (K)')
    # plt.ylabel(r'\textbf{\lambda} (c^{-1})')

    plt.plot(Ts, [(lambda_p__n(T, units="K")/lambda_n__p(T, units="K"))/exp(-Q/(k_b*T)) for T in Ts], 
        linewidth=2.0)
    plt.gca().invert_xaxis()
    plt.show()