r"""
Высчитывание времени с момента БВ по температуре.
Формулы взяты из книги "Космология" Стивена Вейнберга,
Переводное издание 2013 года.
Глава 3.1
"""

import math
import numpy as np
import scipy
from constants import *
from scipy import integrate
from scipy.misc import derivative
import Cacher

cacher = Cacher.Cacher()


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
    # print (res)
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
    integ = t_e*scipy.integrate.quad(__under_int__, 0, m_e/(T))[0]
    return integ

__t0__ = __tfromT__(k_b*10**11)

@cacher.sql_tempreture_cache
def tfromT(T, *, units="eV"):
    r"""
    зависимость времени от температуры. По умолчанию в эВ, 
       доступные единицы
       eV, K
    нормировка на 10**11 K
    """
    if units=="K":
        T = k_b*T
    return __tfromT__(T) - __t0__


@cacher.sql_tempreture_cache
def TnuFromT(T, *, units="eV"):
    r"""
    температура нейтрино. При T>>m_e равна температуре фотонов, 
    при T<<m_e T/T_{\nu}=(11/4)^(1/3)
    в общем случае
    T_\nu=(4/11)^{1/3}TS^{1/3}(m_e/T)
    """
    Tu = T
    if units=="K":
        Tu = k_b*T
    return (4*__S_beaut__(m_e/T)/11.0)**(1/3)*Tu


if __name__ == '__main__':
    import sys
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    if len(sys.argv)-1:
        T = sys.argv[1]
        T = eval(T)
        print("{:.1E}: {:.4E}\n".format(T, tfromT(T, units="K")))
        exit()
    Ts = np.logspace(math.log10(10**8), math.log10(10**11), num=100)
    ts = [tfromT(T, units="K") for T in Ts]
    print("DONE")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'\textbf{time} (s)')
    plt.ylabel(r'\textbf{tempreture} (K)')
    plt.plot([0.994*(T/10**10)**(-2)-0.994*(10)**(-2) for T in Ts], Ts, 
        linewidth=2.0, label=r'$t \to 0$')
    plt.plot([1.78*(T/10**10)**(-2)-1.78*(10)**(-2) for T in Ts], Ts,
        linewidth=2.0, label=r'$t \to \infty$')
    plt.plot(ts, Ts,
        'r--', label=r't(T) modeling result')
    plt.legend()
    # time.sleep(2)
    plt.show()
    time.sleep(2)
    with open("tempreture_data.dat", "w") as f:
        for i in range(len(Ts)):
            f.write("{:.1E}: {:.4E}\n".format(Ts[i], ts[i]))