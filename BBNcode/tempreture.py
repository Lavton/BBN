import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import sys
from scipy import integrate
from scipy.misc import derivative
from tempreture import *
from constants import *

def __s_beaut_integrand_f__(x, y):
    try:
        res = (y**2)*(math.sqrt(y**2+x**2)+(y**2)/(3*math.sqrt(y**2+x**2)))*1/(math.exp(math.sqrt(y**2+x**2))+1)
    except OverflowError:
        res = 0.0
    # print (res)
    return res

def __S_beaut__(x, max_y=np.inf):
    """3.1.18"""
    integ = scipy.integrate.quad(lambda y: __s_beaut_integrand_f__(x, y), 0, max_y)
    return 1 + (45/(2*math.pi**4))*integ[0]


def __eps_beaut_integrand_f__(x, y):
    try:
        res = y**2*math.sqrt(y**2+x**2)/(math.exp(math.sqrt(y**2+x**2))+1)
    except OverflowError:
        res = 0.0
    return res

def __epsilon__(x, max_y=np.inf):
    """3.1.23"""
    integ = scipy.integrate.quad(lambda y: __eps_beaut_integrand_f__(x, y), 0, max_y)
    return 1+(21/8.)*(4.*__S_beaut__(x)/11)**(4./3)+(30/(math.pi**4))*integ[0]

def __under_int__(x):
    return (3-x*derivative(__S_beaut__, x)/__S_beaut__(x))*(__epsilon__(x)**(-1/2.))*x

def __tfromT__(T):
    """3.1.24"""
    t_e = 4.3694
    integ = t_e*scipy.integrate.quad(__under_int__, 0, m_e/(k_b*T))[0]
    return integ

__t0__ = __tfromT__(10**11)

def tfromT(T):
    return __tfromT__(T) - __t0__