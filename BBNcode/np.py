from constants import *
from scipy import integrate
from math import pi

__lambda0__ = 1.63616
__a0__ = (4*g**2)/(pi**3*(c*h*2*pi)**6*(2*pi*h))

def __x__(T):
    return m_e/(k_b*T)

def __x_nu__(T_nu):
    pass

def __q__():
    return Q/m_e

def __eps__(E_e):
    return E_e/m_e

def __Kp_under_int__(T, T_nu, eps):

def __Kp__(T, T_nu, max_y=np.inf):
    return integrate.quad(lambda eps: __Kp_under_int__(T, T_nu, eps), 1, max_y)
