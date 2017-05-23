r"""
Высчитывание времени с момента БВ по температуре.
Формулы взяты из книги "Космология" Стивена Вейнберга,
Переводное издание 2013 года.
Глава 3.1
"""
import constants
import math
import numpy as np
import scipy
from constants import *
from scipy import integrate
from scipy.misc import derivative
import Cacher
import logging

def __s_beaut_integrand_f__(x, y):
    r"""
    Подынтыгральное выражение для плотности энтропии фотонов, электронов и позитронов
    из уравнения 3.1.18 Вайнберга
    $$
    y^2\left(\sqrt{y^2+x^2}+\frac{y^2}{3\sqrt{y^2+x^2}}\right)\frac{1}{\exp\sqrt{y^2+x^2}+1}
    $$
    """
    try:
        res = (y**2)*(math.sqrt(y**2+x**2)+(y**2)/(3*math.sqrt(y**2+x**2)))*1/(math.exp(math.sqrt(y**2+x**2))+1)
    except OverflowError:
        res = 0.0
    return res

def __S_beaut__(x, max_y=np.inf):
    r"""
    плотность энтропии фотонов, электронов и позитронов,
    уравнение 3.1.18
    $$S(x)=1+\frac{45}{2\pi^4}\int\limits_0^\infty dy ...$$
    """
    integ = scipy.integrate.quad(lambda y: __s_beaut_integrand_f__(x, y), 0, max_y)
    return 1 + (45/(2*math.pi**4))*integ[0]


def __eps_beaut_integrand_f__(x, y):
    r"""
    подынтыгральное выражение для полной плотности энергии
    из уравнения 3.1.23
    $$
    \frac{y^2\sqrt{y^2+x^2}}{\exp\sqrt{y^2+x^2}+1}
    $$
    """
    try:
        res = y**2*math.sqrt(y**2+x**2)/(math.exp(math.sqrt(y**2+x**2))+1)
    except OverflowError:
        res = 0.0
    return res

def __epsilon__(x, max_y=np.inf):
    r"""
    выражение для полной плотности энергии
    из уравнения 3.1.23
    $$
    1+\frac{21}{8}\left(\frac{4S(x)}{11}\right)^{4/3}+
    \frac{30}{pi^4}\int\limits_0^\infty dy ...
    $$
    """
    integ = scipy.integrate.quad(lambda y: __eps_beaut_integrand_f__(x, y), 0, max_y)
    return 1+(21/8.)*(4.*__S_beaut__(x)/11)**(4./3)+(30/(math.pi**4))*integ[0]

def __under_int__(x):
    r"""
    подынтыгральное уравнение для ненормированной температуры фотонов
    уравнение 3.1.24

    $$
    \left(3-\frac{xS'(x)}{S(x)}\right)\epsilon^{-1/2}(x)x
    $$
    """
    return (3-x*derivative(__S_beaut__, x)/__S_beaut__(x))*(__epsilon__(x)**(-1/2.))*x

def __tfromT__(T):
    r"""
    не нормированная температура фотонов
    уравнение 3.1.24
    $$
    t=t_e \int\limits_0^\infty dx ...
    $$
    """
    t_e = 4.3694
    T_ = constants.to_norm_tempreture(T) / constants.m_e
    integ = t_e*scipy.integrate.quad(__under_int__, 0, 1/T_)[0]

    return constants.less_time(integ)

__t0__ = __tfromT__(constants.less_tempreture(10**11, units="K"))
@Cacher.cacher.sql_base_cache
def tfromT(T):
    r"""
    зависимость времени от температуры. По умолчанию в эВ, 
       доступные единицы
       eV, K
    нормировка на 10**11 K
    """
    return __tfromT__(T) - __t0__

@Cacher.cacher.sql_base_cache
def Tfromt(t):
    T0 = constants.less_tempreture(10**11, units="K")
    Tn = constants.less_tempreture(10**6, units="K")
    c_Tm = lambda Tn, T0: (math.sqrt(2./(Tn**(-2)+T0**(-2))))
    # Tm = (T0+Tn)/2
    Tm = c_Tm(Tn, T0)
    t0 = tfromT(T0)
    tn = tfromT(Tn)
    tm = tfromT(Tm)
    # for _ in range(10):
    i = 0
    if t <= t0:
        return T0
    while abs(t-tm)/t > 0.001:
        if tm >= t >= t0:
            Tn = Tm
            tn = tfromT(Tn)
        elif tm < t <= tn:
            T0 = Tm 
            t0 = tfromT(T0)
        # Tm = (T0+Tn)/2
        Tm = c_Tm(Tn, T0)
        tm = tfromT(Tm)
    logging.debug(t)

    return Tm

@Cacher.cacher.sql_base_cache
def TnuFromT(T):
    r"""
    температура нейтрино. При T>>m_e равна температуре фотонов, 
    при T<<m_e T/T_{\nu}=(11/4)^(1/3)
    в общем случае
    T_\nu=(4/11)^{1/3}TS^{1/3}(m_e/T)
    """
    Tu = T
    T_ = constants.to_norm_tempreture(T) / constants.m_e
    return (4*__S_beaut__(1/T_)/11.0)**(1/3)*T

# @Cacher.cacher.sql_base_cache
def derriviate_T_from_t(T):
    return 1.
    T_ = T
    if type(T) is np.ndarray:
        T_ = np.copy(T)
        for i in range(len(T_)):
            T_[i] = 1.0/derivative(tfromT, T_[i])
    else:
        T_ = 1.0/derivative(tfromT, T_)
    return T_

if __name__ == '__main__':
    import sys
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    if len(sys.argv)-1:
        T = sys.argv[1]
        T = eval(T)
        print("{:.1E}: {:.4E}\n".format(T, tfromT(T, units="K")))
        exit()
    Ts = constants.less_tempreture(np.logspace(math.log10(10**7), math.log10(10**11), num=200), units="K")
    import datetime
    a = datetime.datetime.now()

    ts = np.array([tfromT(T) for T in Ts])

    b = datetime.datetime.now()
    print(b-a)
    print("DONE")
    Is = [i for i in range(len(Ts))]
    tms = [(TnuFromT(T)) for T in Ts]
    import pickle 
    with open("/tmp/new_result.pickle", "wb") as f:
        pickle.dump((Is, tms),f)
    
    import _t_tnu_vein
    Ts_ = []
    Tnus_ = []
    tu = []
    for (T_, Xn_) in _t_tnu_vein.t_tnu.items():
        tu.append((T_, Xn_))

    tu = sorted(tu, reverse=True)
    tu
    for (T_, Xn_) in tu:
        Ts_.append(T_)
        Tnus_.append(Xn_)
        
    plt.xscale('log')
    # plt.ylim([-2, 4])
    # plt.yscale('log')
    plt.plot(Ts_, Tnus_)
    plt.plot(constants.to_norm_tempreture(Ts, units="K"), [T/TnuFromT(T) for T in Ts], 
            linewidth=2.0, label=r'X_n')
    plt.show()
    tu

    plt.cla()
    plt.clf()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e8, 0.99*1e10])
    plt.xlim([1, 1e4])
    plt.xlabel(r'\textbf{time} (s)')
    plt.ylabel(r'\textbf{tempreture} (K)')
    print(constants.to_norm_tempreture(Ts, units="K"))
    plt.plot(constants.to_norm_time(ts), constants.to_norm_tempreture(Ts, units="K"),
        'r', linewidth=2.0, label=r't(T) modeling result')
    plt.plot(constants.to_norm_time(ts), constants.to_norm_tempreture(np.array([Tfromt(t) for t in ts]), units="K"),
        'b', linewidth=2.0, label=r't(T) modeling result new')

    plt.plot([0.994*(constants.to_norm_tempreture(T, units="K")/10**10)**(-2)-0.994*(10)**(-2) for T in Ts], constants.to_norm_tempreture(Ts, units="K"), 
        linewidth=1.0, label=r'$t \to 0$')
    plt.plot([1.78*(constants.to_norm_tempreture(T, units="K")/10**10)**(-2)-1.78*(10)**(-2) for T in Ts], constants.to_norm_tempreture(Ts, units="K"),
        linewidth=1.0, label=r'$t \to \infty$')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)
    plt.show()

    # plt.cla()
    # plt.clf()
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # # plt.xscale('log')
    # # plt.yscale('log')
    # # plt.ylim([1e8, 0.99*1e10])
    # # plt.xlim([1, 1e4])
    # plt.xlabel(r'\textbf{time} (s)')
    # plt.ylabel(r'\textbf{tempreture} (K)')
    # print(constants.to_norm_tempreture(Ts, units="K"))
    # plt.plot(Ts, [derriviate_T_from_t(T) for T in Ts],
    #     'r', linewidth=2.0, label=r'derriv')
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #        ncol=3, mode="expand", borderaxespad=0.)
    # plt.show()
    # print([1/derriviate_T_from_t(T) for T in Ts])


    # plt.cla()
    # plt.clf()
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel(r'\textbf{time} (s)')
    # plt.ylabel(r'\textbf{tempreture} (K)')
    # for t in Ts:
    #     print("T = {:.2E}".format( constants.to_norm_tempreture(t, units="K")))
    # plt.plot(-derriviate_T_from_t(Ts), Ts,
    #     'r--', label=r'dT/dt modeling result')
    # f1 = lambda T: 0.994*(constants.to_norm_tempreture(T, units="K")/10**10)**(-2)-0.994*(10)**(-2)
    # plt.plot([-1/derivative(f1, T) for T in Ts], constants.to_norm_tempreture(Ts, units="K"), 
        # linewidth=2.0, label=r'$t \to 0$')
    # f2 = lambda T: 1.78*(constants.to_norm_tempreture(T, units="K")/10**10)**(-2)-1.78*(10)**(-2)
    # plt.plot([-1/derivative(f2, T) for T in Ts], constants.to_norm_tempreture(Ts, units="K"),
        # linewidth=2.0, label=r'$t \to \infty$')


    # plt.legend()
    # plt.show()

    with open("tempreture_data.dat", "w") as f:
        for i in range(len(Ts)):
            f.write("{:.1E}: {:.4E}\n".format(Ts[i], ts[i]))