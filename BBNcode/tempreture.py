r"""
Высчитывание времени с момента БВ по температуре.
Формулы взяты из книги "Космология" Стивена Вейнберга,
Переводное издание 2013 года.
Глава 3.1
"""

import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import sys
from constants import *
from scipy import integrate
from scipy.misc import derivative

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

def tfromT(T, *, units="eV"):
    r"""
    зависимость времени от температуры. По умолчанию в эВ, доступные единицы
    eV, K
    """
    if units=="K":
        T = k_b*T
    return __tfromT__(T) - __t0__